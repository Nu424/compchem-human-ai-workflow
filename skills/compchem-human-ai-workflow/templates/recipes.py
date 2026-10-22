"""構造の作り方と解析の「名前付きの手順」（compchem-human-ai-workflow テンプレート）。

計画ファイルからは `recipe: <関数名>` と引数だけで呼ぶ。新しい手順は関数を足して
デコレータで登録し、計画レビューの「新しく足した手順」で知らせる。
長さの単位は Å、エネルギーは eV。
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
from ase import Atoms
from ase.build import add_adsorbate, fcc111, molecule
from ase.constraints import FixAtoms
from ase.io import read

BUILDERS: dict[str, Callable] = {}
GENERATORS: dict[str, Callable] = {}
ANALYSES: dict[str, Callable] = {}


def builder(func):
    BUILDERS[func.__name__] = func
    return func


def generator(func):
    GENERATORS[func.__name__] = func
    return func


def analysis(func):
    ANALYSES[func.__name__] = func
    return func


class Context:
    """手順が、試行フォルダや他の系の構造を参照するための窓口。"""

    def __init__(self, trial: Path, specs: dict[str, dict]):
        self.trial = trial
        self.specs = specs  # 系の名前 → build の指定

    def build(self, name: str) -> Atoms:
        return build(self.specs[name], self)


def build(spec: dict, ctx: Context) -> Atoms:
    """build の指定（recipe と引数）から初期構造を作る。"""
    args = dict(spec)
    name = args.pop("recipe")
    if name not in BUILDERS:
        raise KeyError(f"構造の手順 {name!r} がありません。候補: {sorted(BUILDERS)}")
    atoms = BUILDERS[name](ctx, **args)
    atoms.pbc = True  # 平面波 DFT・多くの MLIP は分子も周期セルで扱う
    return atoms


def recipe_docs() -> dict[str, str]:
    """計画レビュー用に、手順の名前と docstring の 1 行目を返す。"""
    table = {**BUILDERS, **GENERATORS, **ANALYSES}
    return {name: (func.__doc__ or "").strip().splitlines()[0] for name, func in table.items()}


# ---- 構造の作り方 ----


@builder
def molecule_in_box(ctx, name, box):
    """ASE の分子データベースの分子を、箱（Å。数値か [a, b, c]）の中央に置く。"""
    atoms = molecule(name)
    atoms.set_cell([box] * 3 if np.isscalar(box) else box)
    atoms.center()
    return atoms


@builder
def atom_in_box(ctx, symbol, cell):
    """孤立原子を箱の中央に置く。立方体の箱は軌道の縮退を残すので、非立方体を推奨。"""
    atoms = Atoms(symbol, cell=cell)
    atoms.center()
    return atoms


@builder
def fcc111_slab(ctx, symbol, size, a, vacuum, fix_layers=0):
    """fcc(111) スラブ。size=[nx, ny, 層数]、a は格子定数、vacuum は上下それぞれの真空（Å）。下から fix_layers 層を固定。"""
    slab = fcc111(symbol, size=tuple(size), a=a, vacuum=vacuum)
    if fix_layers:
        layers = size[2]
        slab.set_constraint(FixAtoms(mask=slab.get_tags() > layers - fix_layers))  # tag は上の層から 1, 2, ...
    return slab


@builder
def adsorbate_on(ctx, slab, adsorbate, site, height):
    """別の系（slab）の初期構造に吸着種を置く。site は ontop / bridge / fcc / hcp、height は表面からの高さ（Å）。"""
    atoms = ctx.build(slab)
    add_adsorbate(atoms, adsorbate, height, position=site)
    return atoms


@builder
def graphite_ab(ctx, a, c):
    """AB 積層 graphite の単位胞（C4）。a は面内、c は積層方向の格子定数（Å）。"""
    return Atoms(
        "C4",
        scaled_positions=[[0, 0, 0], [1 / 3, 2 / 3, 0], [0, 0, 0.5], [2 / 3, 1 / 3, 0.5]],
        cell=[[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]],
    )


@builder
def from_file(ctx, path):
    """試行フォルダからの相対パスで構造ファイルを読む（.traj なら拘束と磁気モーメントも保持）。"""
    return read(ctx.trial / path)


# ---- 候補構造の生成（複数の構造を返す） ----


@generator
def adsorbate_sites(ctx, slab, adsorbate, sites):
    """1 つのスラブの各サイトに吸着種を置いた候補を作る。sites は {サイト名: 高さ Å}。"""
    candidates = {}
    for site, height in sites.items():
        atoms = ctx.build(slab)
        add_adsorbate(atoms, adsorbate, height, position=site)
        atoms.pbc = True
        candidates[f"{slab}+{adsorbate}@{site}"] = atoms
    return candidates


# ---- 解析 ----


@analysis
def linear_combination(energies, terms):
    """係数付きの線形結合 Σ c_i E_i（eV）。足りない項があれば KeyError。"""
    return sum(coefficient * energies[name] for name, coefficient in terms.items())


@analysis
def relative_to_references(energies, targets, reference):
    """各 target について E(target) + Σ c_i E(reference_i) を返す（吸着エネルギーなど）。"""
    offset = sum(coefficient * energies[name] for name, coefficient in reference.items())
    return {name: energies[name] + offset for name in targets if name in energies}
