"""VASP 版の計算パイプライン（compchem-human-ai-workflow テンプレート）。

    python pipeline_vasp.py plan    <試行>            計画レビュー用の表と初期構造を出す（計算しない）
    python pipeline_vasp.py submit  <試行> [系 ...]   系ごとのジョブと集計ジョブを投入する（人間の OK の後）
    python pipeline_vasp.py run     <試行> <系>       ジョブの中で 1 系を計算する
    python pipeline_vasp.py collect <試行>            結果の表・確認カード・results.db を作る

試行フォルダに plan.yaml を置く。必要なパッケージは ase と pyyaml（画像を出すなら matplotlib）。
VASP の起動コマンドは環境変数 ASE_VASP_COMMAND、PAW の場所は VASP_PP_PATH で与える（job.slurm）。
"""

from __future__ import annotations

import argparse
import json
import shlex
import signal
import socket
import subprocess
import os
import warnings
from datetime import datetime, timezone
from pathlib import Path

import ase
import ase.db
import numpy as np
import yaml
from ase.calculators.vasp import Vasp
from ase.geometry import find_mic
from ase.io import read, write

import recipes

warnings.filterwarnings("ignore", message="Skipping unhashable information adsorbate_info")  # 吸着サイトの情報は extxyz に書けないだけで無害

HERE = Path(__file__).resolve().parent
STEPS = ("relax", "static")  # このテンプレートが扱う段階。増やすなら pipeline() に書き足す
ACTIVE_FOLDER: Path | None = None  # 終了予告のとき STOPCAR を書く先


# ======== 1 系ぶんの流れ ========


def pipeline(trial: Path, plan: dict, name: str) -> None:
    """構造作成 → 構造緩和 → 1 点計算。各段階で収束を確かめてから次へ進む。"""
    steps = plan["systems"][name]["steps"]
    atoms = build_structure(trial, plan, name)  # 1. 構造作成
    if "relax" in steps:
        atoms = calculate(trial, plan, name, "relax", atoms)  # 2. 構造緩和（VASP 内部の最適化）
    if "static" in steps:
        calculate(trial, plan, name, "static", atoms)  # 3. 1 点計算（緩和後の構造で）


def build_structure(trial: Path, plan: dict, name: str):
    """計画の build（名前付きの手順）から初期構造を作り、初期磁気モーメントを設定する。"""
    ctx = recipes.Context(trial, plan["_build_specs"])
    system = plan["systems"][name]
    atoms = recipes.build(system["build"], ctx)
    if "magmom" in system:
        atoms.set_initial_magnetic_moments(system["magmom"])  # SCF の出発点。結果ではない
    return atoms


def calculate(trial: Path, plan: dict, name: str, step: str, atoms):
    """1 段階を VASP で計算し、収束を確かめて保存する。完了済みならその最終構造を返す。"""
    global ACTIVE_FOLDER
    folder = trial / "runs" / name / step
    done = load_json(folder / "result.json")
    if done and done["converged"]:  # 保証 2：完了済みは計算し直さない
        print("SKIP（完了済み）", folder)
        return restore(read(folder / "final.extxyz"), atoms)
    previous = archive_attempt(folder)  # 途中で終わった出力は消さずに退避する
    if previous and step == "relax" and (previous / "final.extxyz").exists():
        atoms = restore(read(previous / "final.extxyz"), atoms)  # 前回の途中の構造から緩和を続ける
    folder.mkdir(parents=True, exist_ok=True)

    start = atoms.copy()
    params = effective_params(plan, name, step)
    atoms = atoms.copy()
    atoms.calc = Vasp(directory=str(folder), txt="vasp.out", **params)
    ACTIVE_FOLDER = folder
    try:
        energy = atoms.get_potential_energy()  # ここで VASP を 1 回起動する。緩和もこの中で VASP が行う
    finally:
        ACTIVE_FOLDER = None
    calc = atoms.calc
    forces = atoms.get_forces()  # 固定原子の力は 0 として扱われる
    result = dict(
        system=name,
        step=step,
        energy_eV=float(energy),  # σ→0 外挿のエネルギー。比較にはこちらを使う
        free_energy_eV=float(calc.results["free_energy"]),  # TOTEN（電子エントロピーを含む）
        fmax_eV_A=float(np.linalg.norm(forces, axis=1).max()),
        magmom_muB=float(calc.results.get("magmom", 0.0)),
        restarted_from=str(previous) if previous else None,
        params=params,
        provenance=provenance(),
    )
    result["converged"], result["convergence_notes"] = check_convergence(folder, params, result["fmax_eV_A"])
    result["observed"], result["warnings"] = check_expectations(plan["systems"][name], start, atoms, result)
    final = atoms.copy()
    write(folder / "final.extxyz", final)
    save_json(folder / "result.json", result)  # 保証 3・4：条件と実行記録を結果と一緒に保存
    print("DONE" if result["converged"] else "NOT CONVERGED", folder, f"{energy:.6f} eV")
    if not result["converged"]:  # 保証 1：未収束なら次の段階へ進まない
        raise RuntimeError(f"未収束です: {folder} {result['convergence_notes']}")
    return final


# ======== 計算条件 ========


def effective_params(plan: dict, name: str, step: str) -> dict:
    """実効パラメータ = 共通条件 < 段階の条件 < 系の条件 < 系の段階ごとの条件 < 収束試験の値。"""
    system = plan["systems"][name]
    params: dict = {}
    for layer in (
        plan["method"],
        plan["steps"][step],
        system.get("params", {}),
        system.get("step_params", {}).get(step, {}),
        system.get("scan_override", {}),
    ):
        params.update(layer)
    for key, value in params.items():
        if isinstance(value, str) and looks_like_number(value):
            raise ValueError(f"{name}/{step}: {key}={value!r} が文字列として読まれました。YAML では 1.0e-7 のように書いてください")
    return params


def looks_like_number(text: str) -> bool:
    try:
        float(text)
        return True
    except ValueError:
        return False


# ======== 確認 ========


def check_convergence(folder: Path, params: dict, fmax: float) -> tuple[bool, list[str]]:
    """全 SCF ループ（VASP 6 の OUTCAR）と、緩和なら最終の残留力を確かめる。問題がなければ収束。"""
    problems = []
    loops = [line for line in (folder / "OUTCAR").read_text(errors="replace").splitlines() if "aborting loop" in line]
    if not loops:
        problems.append("OUTCAR に SCF の終了理由がありません（VASP 6 以外か、途中で終了）")
    failed = [i for i, line in enumerate(loops, 1) if "EDIFF is reached" not in line]
    if failed:
        problems.append(f"SCF が収束しなかったループ: {failed}")
    if params.get("nsw", 0) > 0 and params.get("ibrion", -1) in (1, 2, 3):
        if params.get("ediffg", 0) >= 0:
            problems.append("EDIFFG を負（力の基準、eV/Å）にしてください。エネルギー基準は判定できません")
        elif fmax > abs(params["ediffg"]):
            problems.append(f"最大残留力 {fmax:.3f} eV/Å > {abs(params['ediffg'])}")
    return not problems, problems


def check_expectations(system: dict, start, final, result: dict) -> tuple[dict, list[str]]:
    """計画の expect と照らし合わせ、観測値と警告を返す。警告は確認カードに載る。"""
    expect = system.get("expect", {})
    observed: dict = {}
    warnings = []
    vectors, _ = find_mic(final.positions - start.positions, final.cell, final.pbc)
    observed["max_displacement_A"] = float(np.linalg.norm(vectors, axis=1).max())
    if "magmom" in expect and abs(result["magmom_muB"] - expect["magmom"]) > 0.1:
        warnings.append(f"磁化 {result['magmom_muB']:.2f} μB（期待 {expect['magmom']}）")
    if "distance" in expect:
        i, j = expect["distance"]["atoms"]
        low, high = expect["distance"]["range"]
        observed["distance_A"] = float(final.get_distance(i, j, mic=True))
        if not low <= observed["distance_A"] <= high:
            warnings.append(f"原子 {i}-{j} の距離 {observed['distance_A']:.3f} Å（期待 {low}〜{high}）")
    if "max_displacement" in expect and observed["max_displacement_A"] > expect["max_displacement"]:
        warnings.append(f"最大変位 {observed['max_displacement_A']:.2f} Å（サイトの移動や再構成の可能性）")
    return observed, warnings


# ======== 計画 ========


def load_plan(trial: Path) -> dict:
    """plan.yaml を読み、収束試験の展開と、段階・パラメータの検証をする。"""
    plan = yaml.safe_load((trial / "plan.yaml").read_text(encoding="utf-8"))
    plan["_build_specs"] = {name: system["build"] for name, system in plan["systems"].items()}
    if "scan" in plan:
        expand_scan(plan)
    for name, system in plan["systems"].items():
        unknown = set(system["steps"]) - set(STEPS)
        if unknown:
            raise ValueError(f"{name}: 未対応の段階 {unknown}。pipeline() に処理を書き足してください")
        for step in system["steps"]:
            effective_params(plan, name, step)
    return plan


def expand_scan(plan: dict) -> None:
    """収束試験：scan.systems の各系を、振る値ごとの系に展開する（名前は <系>__<キー>-<値>）。"""
    scan = plan["scan"]
    for base in scan["systems"]:
        system = plan["systems"].pop(base)
        for value in scan["values"]:
            copy = json.loads(json.dumps(system))
            copy.update(scan_base=base, scan_value=value, scan_override={scan["key"]: value})
            plan["systems"][scan_name(base, scan["key"], value)] = copy


def scan_name(base: str, key: str, value) -> str:
    tag = "x".join(map(str, value)) if isinstance(value, list) else str(value)
    return f"{base}__{key}-{tag}"


def show_plan(trial: Path, plan: dict) -> None:
    """計画レビューの土台：共通条件、段階の条件、系ごとの差分、使う手順、規模を表示し、初期構造を書き出す。"""
    out = [f"# 計画の実効条件：{trial.name}", "", f"目的：{plan.get('purpose', '').strip()}", ""]
    out += ["## 共通条件", "", "|キー|値|", "|---|---|"] + [f"|{k}|{v}|" for k, v in plan["method"].items()]
    out += ["", "## 段階の条件", ""] + [f"- {step}: {fmt(params)}" for step, params in plan["steps"].items()]
    out += ["", "## 計算対象", "", "|系|役割|構造の作り方|段階|系ごとの差分|期待値|", "|---|---|---|---|---|---|"]
    folder = trial / "results" / "initial"
    folder.mkdir(parents=True, exist_ok=True)
    used = set()
    for name, system in plan["systems"].items():
        diff = dict(system.get("params", {}), **system.get("scan_override", {}))
        extra = "; ".join(f"{s}: {fmt(p)}" for s, p in system.get("step_params", {}).items())
        build = system["build"]
        used.add(build["recipe"])
        out.append(
            f"|{name}|{system.get('role', '')}|{build['recipe']}|{' → '.join(system['steps'])}"
            f"|{fmt(diff)}{'; ' + extra if extra else ''}|{fmt(system.get('expect', {}))}|"
        )
        atoms = build_structure(trial, plan, name)  # 手順が動くことの確認も兼ねる
        write(folder / f"{name}.extxyz", atoms)
        write_images(atoms, folder / name)
    used |= {spec["recipe"] for spec in plan.get("results", {}).values()}
    docs = recipes.recipe_docs()
    out += ["", "## 使う手順", ""] + [f"- `{n}`：{docs.get(n, '（説明なし）')}" for n in sorted(used)]
    out += ["", "## 既知の近似と未検証の事項", ""] + [f"- {c}" for c in plan.get("caveats", [])]
    runs = sum(len(s["steps"]) for s in plan["systems"].values())
    out += ["", "## 規模", "", f"- 系の数（= 計算ジョブ数）：{len(plan['systems'])}、VASP の実行回数：{runs}"]
    if "resources" in plan:
        out.append(f"- 資源：{fmt(plan['resources'])}")
    out += ["", f"初期構造：{folder}"]
    print("\n".join(out))


# ======== 投入と集計 ========


def submit(trial: Path, plan: dict, names: list[str], dry_run: bool) -> None:
    """系ごとのジョブと、全ジョブ終了後（afterany）の集計ジョブを投入する。"""
    res = plan["resources"]
    common = [f"--{k}={res[k]}" for k in ("account", "partition", "qos") if k in res]
    logs = trial / "logs"
    logs.mkdir(exist_ok=True)
    job, me = str(HERE / "job.slurm"), str(Path(__file__).resolve())
    ids = []
    for name in names or list(plan["systems"]):
        ids.append(
            sbatch(
                [*common, f"--ntasks={res['ntasks']}", f"--time={res['time']}", f"--job-name={name}",
                 f"--output={logs}/{name}-%j.out", job, me, "run", str(trial), name],
                logs, dry_run,
            )
        )
    sbatch(
        [*common, "--ntasks=1", "--time=00:15:00", "--job-name=collect", f"--output={logs}/collect-%j.out",
         "--dependency=afterany:" + ":".join(ids), job, me, "collect", str(trial)],
        logs, dry_run,
    )


def sbatch(args: list[str], logs: Path, dry_run: bool) -> str:
    """sbatch を 1 回実行し、JOBID を logs/submissions.jsonl に記録する。"""
    command = ["sbatch", "--parsable", *args]
    print(shlex.join(command))
    if dry_run:
        return "DRYRUN"
    job_id = subprocess.run(command, check=True, capture_output=True, text=True).stdout.strip().split(";")[0]
    with (logs / "submissions.jsonl").open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(dict(utc=now(), job_id=job_id, argv=command), ensure_ascii=False) + "\n")
    return job_id


def collect(trial: Path, plan: dict) -> None:
    """結果の表・確認カード・構造の画像・results.db を作る。DFT は起動しない。"""
    out = trial / "results"
    out.mkdir(exist_ok=True)
    db_path = out / "results.db"
    db_path.unlink(missing_ok=True)  # collect のたびに作り直す派生物
    finals, cards = {}, []
    with ase.db.connect(db_path) as db:
        for name, system in plan["systems"].items():
            for step in system["steps"]:
                folder = trial / "runs" / name / step
                result = load_json(folder / "result.json")
                cards.append(card(name, step, system, result))
                finals[name] = result
                if result is None:
                    continue
                atoms = read(folder / "final.extxyz")
                db.write(atoms, system=name, step=step, converged=result["converged"], energy_eV=result["energy_eV"],
                         data=dict(params=result["params"], warnings=result["warnings"], provenance=result["provenance"]))
                write_images(atoms, out / "structures" / f"{name}_{step}")
    energies = {n: r["energy_eV"] for n, r in finals.items() if r and r["converged"]}
    text = [f"# 結果：{trial.name}", "", f"更新 UTC：{now()}", "", "電子エネルギー差（eV）。ZPE・熱補正なし。",
            "状態の `確認待ち` は、人間が確認カードを見て `妥当` か `除外` に変える。", ""]
    text += results_table(plan, energies, finals)
    text += ["", "## 確認カード", ""] + cards
    (out / "RESULTS.md").write_text("\n".join(text) + "\n", encoding="utf-8")
    print("WROTE", out / "RESULTS.md")


def results_table(plan: dict, energies: dict, finals: dict) -> list[str]:
    """計画の results を、名前付きの解析で計算して表にする。未収束の項を含む値は出さない。"""
    lines = ["|量|値 eV|比較対象|状態|意味|", "|---|---:|---|---|---|"]
    for label, spec in plan.get("results", {}).items():
        args = {k: v for k, v in spec.items() if k not in ("recipe", "meaning", "compare")}
        func = recipes.ANALYSES[spec["recipe"]]
        variants = [("", energies, finals)]
        if "scan" in plan:  # 収束試験：振った値ごとに同じ量を計算する
            variants = [(f" [{plan['scan']['key']}={v}]", *scan_slice(plan, v, energies, finals)) for v in plan["scan"]["values"]]
        rows = []
        for suffix, e, f in variants:
            names = list(args.get("terms", {})) + list(args.get("reference", {})) + list(args.get("targets", []))
            status = term_status(names, f)
            rows.append((suffix, func(e, **args) if status in ("確認待ち", "要確認") else None, status))
        last = rows[-1][1] if "scan" in plan and isinstance(rows[-1][1], float) else None  # 収束試験では最後の値を基準にする
        for suffix, value, status in rows:
            shown = f"{value:.6f}" if isinstance(value, float) else ("—" if value is None else fmt(value))
            if last is not None and isinstance(value, float):
                shown += f"（最後の値との差 {value - last:+.4f}）"
            lines.append(f"|{label}{suffix}|{shown}|{spec.get('compare', '')}|{status}|{spec.get('meaning', '')}|")
    return lines


def scan_slice(plan: dict, value, energies: dict, finals: dict) -> tuple[dict, dict]:
    """収束試験の 1 つの値について、展開した名前を元の系の名前に戻す。"""
    scan = plan["scan"]
    rename = {scan_name(b, scan["key"], value): b for b in scan["systems"]}
    return ({rename.get(n, n): e for n, e in energies.items()}, {rename.get(n, n): r for n, r in finals.items()})


def term_status(names: list[str], finals: dict) -> str:
    if any(finals.get(n) is None for n in names):
        return "未完了"
    if any(not finals[n]["converged"] for n in names):
        return "未収束"
    return "要確認" if any(finals[n]["warnings"] for n in names) else "確認待ち"


def card(name: str, step: str, system: dict, result: dict | None) -> str:
    """確認カード 1 枚（review-points.md の型）。"""
    if result is None:
        return f"```text\n[{name}/{step}]  未完了\n```"
    conv = "✓" if result["converged"] else "✗ " + "; ".join(result["convergence_notes"])
    return "\n".join([
        "```text",
        f"[{name}/{step}]  役割: {system.get('role', '')}",
        f"  収束      : {conv}（最大残留力 {result['fmax_eV_A']:.3f} eV/Å）",
        f"  磁化      : {result['magmom_muB']:.2f} μB" + (f"（期待 {system['expect']['magmom']}）" if "magmom" in system.get("expect", {}) else ""),
        f"  観測値    : {fmt(result['observed'])}",
        f"  エネルギー: {result['energy_eV']:.6f} eV",
        f"  警告      : {'; '.join(result['warnings']) or 'なし'}",
        f"  構造      : results/structures/{name}_{step}_top.png, _side.png",
        "```",
    ])


# ======== 記録と小物 ========


def provenance() -> dict:
    """保証 4：何を、どこで、どの版で実行したか。"""

    def git(*args):
        try:
            return subprocess.run(["git", *args], cwd=HERE, capture_output=True, text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return None

    commit = git("rev-parse", "HEAD")
    return dict(
        utc=now(),
        host=socket.gethostname(),
        slurm_job_id=os.environ.get("SLURM_JOB_ID"),
        ase_version=ase.__version__,
        script=str(Path(__file__).resolve()),
        git_commit=commit,
        git_dirty=bool(git("status", "--porcelain")) if commit else None,
    )


def archive_attempt(folder: Path) -> Path | None:
    """前回の出力が残っていれば <段階>.attempt-<時刻> に移して返す。"""
    if not folder.exists() or not any(folder.iterdir()):
        return None
    target = folder.with_name(f"{folder.name}.attempt-{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}")
    folder.rename(target)
    print("ARCHIVED", target)
    return target


def restore(atoms, template):
    """ファイルから読んだ構造に、元の拘束と初期磁気モーメントを戻す（原子の順番は同じ前提）。"""
    atoms.set_constraint(template.constraints)
    atoms.set_initial_magnetic_moments(template.get_initial_magnetic_moments())
    return atoms


def write_images(atoms, stem: Path) -> None:
    """上面図と側面図の PNG を書く。matplotlib がなければ何もしない。"""
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return
    stem.parent.mkdir(parents=True, exist_ok=True)
    write(f"{stem}_top.png", atoms)
    write(f"{stem}_side.png", atoms, rotation="-90x")


def request_stop(signum, frame) -> None:
    """終了予告（SIGTERM）を受けたら、実行中の VASP に STOPCAR で停止を頼む。"""
    if ACTIVE_FOLDER is not None:
        (ACTIVE_FOLDER / "STOPCAR").write_text("LSTOP = .TRUE.\n")
        print("STOPCAR を書きました", ACTIVE_FOLDER, flush=True)


def fmt(value) -> str:
    if isinstance(value, dict):
        return ", ".join(f"{k}={fmt(v)}" for k, v in value.items())
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_json(path: Path) -> dict | None:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def save_json(path: Path, data: dict) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(path)  # 書きかけのファイルを残さない


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["plan", "submit", "run", "collect"])
    parser.add_argument("trial", type=Path, help="plan.yaml のある試行フォルダ")
    parser.add_argument("names", nargs="*", help="run では 1 系、submit では投入する系（省略で全系）")
    parser.add_argument("--dry-run", action="store_true", help="submit で sbatch を表示するだけにする")
    args = parser.parse_args()
    trial = args.trial.resolve()
    plan = load_plan(trial)
    if args.command == "plan":
        show_plan(trial, plan)
    elif args.command == "submit":
        submit(trial, plan, args.names, args.dry_run)
    elif args.command == "run":
        if len(args.names) != 1:
            parser.error("run には系の名前を 1 つ指定してください")
        signal.signal(signal.SIGTERM, request_stop)
        pipeline(trial, plan, args.names[0])
    else:
        collect(trial, plan)


if __name__ == "__main__":
    main()
