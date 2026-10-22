# 計画ファイルの書き方

目的は「汎用的な設定体系」を作ることではなく、**その計画が人間にも機械にも分かる形で書けていること**。すべての計算を1つの型で表そうとしない。ただし、同じ種類の計算は同じ型にする（毎回形が変わると読む側が困る）。

## 形式

- 人間が書く計画は **YAML**（`plan.yaml`）。コメントを書けるので、理由や注記を添えられる。
- 値の単位は、キー名かコメントで必ず分かるようにする（例：`box: 16  # Å`）。
- YAML（PyYAML）の読み違いに注意する。`1e-7` は文字列になる（`1.0e-7` と書く）。`on` `off` `yes` `no` は、キーでも値でも真偽値になる（キー名に使わない。文字列なら引用符で囲む）。
- 実効パラメータ（共通条件にすべての上書きを適用した結果）はコードが計算し、`plan` コマンドで表示し、結果ファイルにも保存する。人間が手で展開しない。

## 共通の骨格

どの型の計画も、次の見出しを持つ。

```yaml
type: <計画の型の名前>        # reaction-energy / convergence / mlip-screening / 新しい型
purpose: <何を知りたいか、1〜2 行>
caveats:                      # 既知の近似と未検証の事項
  - <...>
resources: {...}              # Slurm の資源（使わない計算では省略）
method: {...}                 # 全対象に共通の計算条件
steps: {...}                  # 処理の段階ごとの条件（relax / static など）
targets / systems: {...}      # 計算対象。中身の形は型ごとに決める
results: {...}                # 何を求めるか（反応式、比較、図など）
```

各対象には、できるだけ次を書く。

| キー | 中身 |
|---|---|
| `role` | この対象が何のためにあるか（「O₂ 解離の参照」など） |
| `build` | 構造の作り方。`recipe`（手順の名前）と引数 |
| `steps` | 処理の段階の並び（例：`[relax, static]`） |
| `params` / `step_params` | 共通条件からの差分、段階ごとの差分 |
| `expect` | 確認カードで自動確認する値（磁化、結合長の範囲など） |
| `notes` | この対象に固有の判断の解説（learning.md の書き方で、1〜3 行） |

## 名前付きの手順

構造の作り方（`build.recipe`）、候補の生成（`candidates.recipe`）、解析（`results` の `recipe`）は、計画からは**名前と引数だけ**で呼ぶ。中身は `recipes.py` の関数として自由に書き、デコレータで登録する。

```yaml
build: {recipe: adsorbate_on, slab: pt111, adsorbate: H, site: fcc, height: 1.0}
```

```python
@builder
def adsorbate_on(ctx, slab, adsorbate, site, height):
    """別の系（slab）の初期構造に吸着種を置く。site は ontop / bridge / fcc / hcp、height は表面からの高さ（Å）。"""
    atoms = ctx.build(slab)
    add_adsorbate(atoms, adsorbate, height, position=site)
    return atoms
```

docstring の 1 行目は `plan` コマンドの「使う手順」に表示されるので、何をするかと単位をそこに書く。

- 足し引きで表せない解析（状態方程式のフィット、最小値の選択、平均など）も、同じように名前付きの手順にする。
- 計画で表しにくい特殊な処理を、YAML の条件分岐や繰り返しで無理に表現しない。Python の手順にして、名前と docstring で何をするか分かるようにする。
- 新しい手順を書いたら、計画レビューの「新しく足した手順」に載せる。

## 型の見本

| 型 | テンプレート | 対象の形 |
|---|---|---|
| reaction-energy | `templates/plan_reaction-energy.yaml` | 系の一覧と、反応式（係数付きの線形結合） |
| convergence | `templates/plan_convergence.yaml` | 固定する構造と、`scan`（振るキー・値の列・対象の系）。`pipeline_vasp.py` が値ごとの系に展開する |
| mlip-screening | `templates/plan_mlip-screening.yaml` | 参照（`references`）と候補の生成規則（`candidates`）、モデル・緩和条件 |

reaction-energy と convergence は `pipeline_vasp.py`、mlip-screening は `pipeline_mlip.py` で動く。

## 新しい型を作るとき

1. 共通の骨格から始め、対象（`systems` / `targets`）と `results` の形だけをその計算に合わせて決める。
2. 型の名前を `type` に書き、試行の README に「新しい型」と明記する。
3. 同じ種類の計算を 2 回目に行うときに、その計画を `templates/plan_<型>.yaml` に昇格させる。
4. 既存の型を少し広げれば済むなら、新しい型を作らずに広げる。
