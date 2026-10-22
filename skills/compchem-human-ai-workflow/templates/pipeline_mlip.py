"""MLIP 版の計算パイプライン（compchem-human-ai-workflow テンプレート）。

    python pipeline_mlip.py plan <試行>   計画レビュー用の表と初期構造を出す（計算しない）
    python pipeline_mlip.py run  <試行>   参照と候補を順に緩和し、results/results.db に保存して結果の表を作る

1 プロセスで全構造を順に計算するので、ase.db に直接書いてよい（多数のジョブから同時に書かない）。
必要なパッケージは ase と pyyaml（mace-mp を使うなら mace-torch、画像を出すなら matplotlib）。
"""

from __future__ import annotations

import argparse
import inspect
import os
import socket
import subprocess
import warnings
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import ase
import ase.db
import ase.optimize
import numpy as np
import yaml
from ase.calculators.singlepoint import SinglePointCalculator
from ase.geometry import find_mic
from ase.io import write

import recipes

warnings.filterwarnings("ignore", message="Skipping unhashable information adsorbate_info")  # 吸着サイトの情報は extxyz に書けないだけで無害

HERE = Path(__file__).resolve().parent


# ======== 全体の流れ ========


def pipeline(trial: Path, plan: dict) -> None:
    """モデルの準備 → 構造作成 → 緩和 → 確認 → 保存 → 結果の表。"""
    calc = make_calculator(plan["model"])  # 1. モデルの準備
    structures, candidates = build_structures(trial, plan)  # 2. 構造作成（参照と候補）
    (trial / "results").mkdir(exist_ok=True)
    db = ase.db.connect(trial / "results" / "results.db")
    slab_size = len(structures[plan["candidates"]["slab"]])  # 候補では、これ以降の原子が吸着種
    for label, atoms in structures.items():
        if any(row.converged for row in db.select(label=label)):  # 保証 2：完了済みは計算し直さない
            print("SKIP（完了済み）", label)
            continue
        relaxed, result = relax(atoms, calc, plan["relax"])  # 3. 構造緩和
        observed, warnings = check(plan, atoms, relaxed, slab_size if label in candidates else None)  # 4. 確認
        db.write(  # 5. 保存（保証 1・3・4：収束、条件、実行記録を一緒に）
            relaxed,
            label=label,
            converged=result["converged"],
            energy_eV=result["energy_eV"],
            fmax_eV_A=result["fmax_eV_A"],
            steps=result["steps"],
            data=dict(observed=observed, warnings=warnings, model=plan["model"], relax=plan["relax"], provenance=provenance(plan)),
        )
        print("DONE" if result["converged"] else "NOT CONVERGED", label, f"{result['energy_eV']:.6f} eV")
    report(trial, plan, db, structures, candidates)  # 6. 結果の表


def make_calculator(model: dict):
    """計画の model から ASE calculator を作る。"""
    if model["name"] == "emt":
        from ase.calculators.emt import EMT

        return EMT()  # 動作確認用。物理的な結論には使わない
    if model["name"] == "mace-mp":
        from mace.calculators import mace_mp

        return mace_mp(model=model.get("size", "medium"), device=model.get("device", "cpu"), default_dtype=model.get("dtype", "float64"))
    raise ValueError(f"未対応のモデル {model['name']!r}。make_calculator に書き足してください")


def build_structures(trial: Path, plan: dict) -> tuple[dict, list[str]]:
    """参照と候補の初期構造を作る。候補は名前付きの生成手順（candidates.recipe）で作る。"""
    ctx = recipes.Context(trial, {name: ref["build"] for name, ref in plan["references"].items()})
    structures = {name: ctx.build(name) for name in plan["references"]}
    spec = {k: v for k, v in plan["candidates"].items() if k not in ("recipe", "max_lateral_shift")}
    generated = recipes.GENERATORS[plan["candidates"]["recipe"]](ctx, **spec)
    structures.update(generated)
    return structures, list(generated)


def relax(atoms, calc, cfg: dict):
    """ASE の最適化手法で緩和し、最終の残留力で収束を判定する（固定原子の力は除く）。"""
    atoms = atoms.copy()
    atoms.calc = calc
    optimizer = getattr(ase.optimize, cfg["optimizer"])(atoms, logfile=None)
    optimizer.run(fmax=cfg["fmax"], steps=cfg["steps"])
    energy = float(atoms.get_potential_energy())
    forces = atoms.get_forces()
    fmax = float(np.linalg.norm(forces, axis=1).max())
    atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=forces)  # モデル本体を db に持ち込まない
    return atoms, dict(energy_eV=energy, fmax_eV_A=fmax, steps=optimizer.nsteps, converged=fmax <= cfg["fmax"])


def check(plan: dict, start, final, adsorbate_from: int | None) -> tuple[dict, list[str]]:
    """最大変位と、候補なら吸着種（index adsorbate_from 以降の原子）の面内の移動を調べる。"""
    vectors, _ = find_mic(final.positions - start.positions, final.cell, final.pbc)
    observed = dict(max_displacement_A=float(np.linalg.norm(vectors, axis=1).max()))
    warnings = []
    limit = plan["candidates"].get("max_lateral_shift")
    if adsorbate_from is not None and limit is not None:
        shift = float(np.linalg.norm(vectors[adsorbate_from:, :2], axis=1).max())
        observed["lateral_shift_A"] = shift
        if shift > limit:
            warnings.append(f"吸着種が面内に {shift:.2f} Å 動いた（初期サイトから移動した可能性）")
    return observed, warnings


# ======== 計画と結果 ========


def show_plan(trial: Path, plan: dict) -> None:
    """計画レビューの土台：モデル、緩和条件、参照と候補、使う手順、規模を表示し、初期構造を書き出す。"""
    structures, candidates = build_structures(trial, plan)
    folder = trial / "results" / "initial"
    folder.mkdir(parents=True, exist_ok=True)
    out = [f"# 計画の実効条件：{trial.name}", "", f"目的：{plan.get('purpose', '').strip()}", ""]
    out += [f"- モデル：{fmt(plan['model'])}", f"- 緩和：{fmt(plan['relax'])}", ""]
    out += ["|構造|種類|原子数|役割|", "|---|---|---:|---|"]
    for label, atoms in structures.items():
        kind = "候補" if label in candidates else "参照"
        role = plan["references"].get(label, {}).get("role", "")
        out.append(f"|{label}|{kind}|{len(atoms)}|{role}|")
        write(folder / f"{safe(label)}.extxyz", atoms)
        write_images(atoms, folder / safe(label))
    used = {ref["build"]["recipe"] for ref in plan["references"].values()}
    used |= {plan["candidates"]["recipe"]} | {spec["recipe"] for spec in plan.get("results", {}).values()}
    docs = recipes.recipe_docs()
    out += ["", "## 使う手順", ""] + [f"- `{n}`：{docs.get(n, '（説明なし）')}" for n in sorted(used)]
    out += ["", "## 既知の近似と未検証の事項", ""] + [f"- {c}" for c in plan.get("caveats", [])]
    out += ["", f"構造の数：{len(structures)}（参照 {len(structures) - len(candidates)}、候補 {len(candidates)}）", f"初期構造：{folder}"]
    print("\n".join(out))


def report(trial: Path, plan: dict, db, structures: dict, candidates: list[str]) -> None:
    """最新の行から結果の表と確認用の一覧を作る。未収束の構造は値の計算に使わない。"""
    latest = {}
    for row in db.select():
        if row.label in structures and (row.label not in latest or row.id > latest[row.label].id):
            latest[row.label] = row
    energies = {label: row.energy_eV for label, row in latest.items() if row.converged}
    text = [f"# 結果：{trial.name}", "", f"更新 UTC：{now()}、モデル：{fmt(plan['model'])}", ""]
    for name, spec in plan.get("results", {}).items():
        func = recipes.ANALYSES[spec["recipe"]]
        args = {k: v for k, v in spec.items() if k not in ("recipe", "meaning", "compare")}
        if "targets" in inspect.signature(func).parameters:
            args.setdefault("targets", candidates)
        values = func(energies, **args)
        text += [f"## {name}", "", spec.get("meaning", ""), "", "|構造|値 eV|状態|警告|", "|---|---:|---|---|"]
        for label, value in sorted(values.items(), key=lambda item: item[1]):
            warnings = latest[label].data["warnings"]
            text.append(f"|{label}|{value:.4f}|{'要確認' if warnings else '確認待ち'}|{'; '.join(warnings)}|")
        missing = [label for label in candidates if label not in values]
        if missing:
            text.append(f"\n未完了または未収束：{', '.join(missing)}")
        text.append("")
    text += ["## 構造ごとの記録", "", "|構造|収束|最大残留力 eV/Å|ステップ|観測値|", "|---|---|---:|---:|---|"]
    for label, row in latest.items():
        text.append(f"|{label}|{'✓' if row.converged else '✗'}|{row.fmax_eV_A:.3f}|{row.steps}|{fmt(row.data['observed'])}|")
        write_images(row.toatoms(), trial / "results" / "structures" / safe(label))
    (trial / "results" / "RESULTS.md").write_text("\n".join(text) + "\n", encoding="utf-8")
    print("WROTE", trial / "results" / "RESULTS.md")


# ======== 記録と小物 ========


def provenance(plan: dict) -> dict:
    """保証 4：何を、どこで、どの版で実行したか。"""

    def git(*args):
        try:
            return subprocess.run(["git", *args], cwd=HERE, capture_output=True, text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return None

    package = {"mace-mp": "mace-torch"}.get(plan["model"]["name"])
    commit = git("rev-parse", "HEAD")
    return dict(
        utc=now(),
        host=socket.gethostname(),
        slurm_job_id=os.environ.get("SLURM_JOB_ID"),
        ase_version=ase.__version__,
        model_package_version=metadata.version(package) if package else None,
        script=str(Path(__file__).resolve()),
        git_commit=commit,
        git_dirty=bool(git("status", "--porcelain")) if commit else None,
    )


def write_images(atoms, stem: Path) -> None:
    """上面図と側面図の PNG を書く。matplotlib がなければ何もしない。"""
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return
    stem.parent.mkdir(parents=True, exist_ok=True)
    write(f"{stem}_top.png", atoms)
    write(f"{stem}_side.png", atoms, rotation="-90x")


def safe(label: str) -> str:
    return label.replace("+", "_").replace("@", "_")


def fmt(value) -> str:
    if isinstance(value, dict):
        return ", ".join(f"{k}={fmt(v)}" for k, v in value.items())
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["plan", "run"])
    parser.add_argument("trial", type=Path, help="plan.yaml のある試行フォルダ")
    args = parser.parse_args()
    trial = args.trial.resolve()
    plan = yaml.safe_load((trial / "plan.yaml").read_text(encoding="utf-8"))
    if args.command == "plan":
        show_plan(trial, plan)
    else:
        pipeline(trial, plan)


if __name__ == "__main__":
    main()
