---
name: compchem-human-ai-workflow
description: >-
  計算化学シミュレーション（ASE・VASP・MLIP・Slurm）を人間と AI エージェントが協働して
  計画・実装・実行・確認・記録するための理念、進め方、計画ファイルとコードのテンプレート。
  DFT や MLIP の計算コードを書くとき、計算を計画・投入するとき、結果を確認・まとめるとき、
  試行フォルダや解析結果を整理するときに使う。Use for computational chemistry workflows:
  ASE, VASP, INCAR, MLIP, ase.db, Slurm/sbatch jobs, reaction/adsorption energies,
  convergence tests, structure relaxation.
---

# 計算化学の人間–AI 協働

## 理念

1. **計算するのは AI、問いを立てて結論に責任を持つのは人間。** 何を計算し、どの近似を使い、結果をどう解釈するかを、人間が自分の言葉で説明できる状態にする。
2. **コードではなく、計画と証拠をレビューしてもらう。** 人間が読むのは計画レビューと確認カード。配管のコードを読ませない作りにする。
3. **計算者の理解と洞察が広がるように進める。** 説明は、その計算に固有の判断に添える短い助言にする。一般論の講義、予想欄、形式的な問いかけはしない。
4. **近似と未検証の事項を隠さない。** 「動いた」と「正しい」を区別する。仮説を事実として書かない。
5. **すべての数値と図は、元の計算までたどれるようにする。** 棄却したものも理由付きで残す。
6. **シンプルに保つ。** 仕組みは必要になってから足す。作法は規模に応じて提案し、強制しない。

## 役割分担

| 人間 | AI |
|---|---|
| 問い・計算対象・近似の選択を決める | 計画を起こし、選択肢と根拠を示す |
| 投入前に計画を確認して OK を出す | 実装・投入・監視・記録・集計 |
| 確認カードと結果を見て判断する | 確認材料を作り、異常を指摘し、原因の候補を示す |
| 結論を確定する | まとめを書き、主張と証拠を対応させる |

## 進め方

```text
- [ ] 0. 規模を見て、計画を立てるか提案する
- [ ] 1. 計画ファイルを書き、計画レビューを提示する → 人間の OK を待つ（必須）
- [ ] 2. テンプレートの型で実装する
- [ ] 3. （提案）代表的な 1〜2 系で試し計算 → 確認カードを見せる
- [ ] 4. 本計算 → 確認カードと結果の表を作る
- [ ] 5. まとめを書く →（提案）結論前の確認
- [ ] 6. 試行の README とテーマの索引を更新する。出てきた概念を学習メモに残す
```

**0. 計画を立てるかの目安。** 次のどれかに当てはまれば「計画を立ててから進めたほうがよさそうです」と提案する。決めるのは人間。

- 初めて扱う種類の系・手法・計算エンジンがある
- ジョブ数が多い、または合計の計算時間が長い（目安：数百コア時間以上）
- 結果を結論や図に使う予定がある
- 作業フォルダが変わるくらい、計算の内容が大きく変わる

当てはまらない小さな試し計算は、試行の README に「探索。結論には使わない」と書いて軽く進めてよい。

**1. 計画レビュー（必須）。** `sbatch` や長時間の計算を始める前に、[review-points.md](references/review-points.md) の型で計画レビューを提示し、人間の OK を待つ。10 分で読める分量にする。計画ファイルの書き方は [plan-format.md](references/plan-format.md)。

**3〜5. 確認ポイント。** 試し計算の後は「残りを同じ条件で流してよいか」、結論の前は「主張と証拠が対応しているか」を確認してもらう。見せるものは [review-points.md](references/review-points.md)。

**6. 記録と学び。** フォルダと README の書き方は [project-layout.md](references/project-layout.md)、助言と学習メモの書き方は [learning.md](references/learning.md)。

## 5つの保証

壊れると誤った結論に直結するものだけを、すべての計算コードで守る。これ以外の仕組み（入力ハッシュ、ロック機構、承認ファイルなど）は、実際に問題が起きるまで作らない。

1. **未収束の結果を使わない。** 電子（全 SCF ループ）とイオン（残留力）の収束を判定して結果に保存し、集計で除外する。終了コード 0 を収束とみなさない。
2. **既存の結果を上書きしない。** 完了済みはスキップする。やり直しは、古い出力を退避してから明示的に行う。
3. **結果に計算条件を添えて保存する。** エネルギーと一緒に、実効パラメータ・最終構造・磁化・残留力を保存する。
4. **何を実行したかを記録する。** git のコミット、Slurm の JOBID、実行時刻、ASE のバージョンを結果に書く。
5. **投入前に人間の OK をもらう。** これはスクリプトではなく、エージェントの振る舞いとして守る。

## 実装の型

- **処理の流れを一直線に書く。** 1 系ぶんの流れを、上から読めば「構造作成 → 緩和 → 1 点計算 → 確認 → 保存」と分かる関数にする。全体は `plan`（計画表示）→ `submit` → `run`（ジョブ内）→ `collect` のコマンドに分ける。
- **系ごとの違いは計画ファイルに書く。** コードの中で系の名前による分岐（`if name == "pt111"` など）をしない。
- **構造の作り方と解析は「名前付きの手順」にする。** 計画ファイルから名前で呼び、手順の中身は Python で自由に書く。新しい手順を足したら計画レビューで知らせる。
- **部品は使う、枠組みには縛られない。** ASE・pymatgen・packmol などの定評ある部品は積極的に使う。重いワークフロー管理の枠組みは、数百件規模になるまで導入しない。
- **保存の単位。** VASP のようにファイル入出力が必要なエンジンは `runs/<系>/<段階>/` に 1 計算 1 フォルダ。Python で完結する計算（MLIP など）は `ase.db`。集約は `collect` が 1 プロセスで行う（ネットワーク上のファイルシステムで、多数のジョブから同じ SQLite に同時に書かない）。
- **コメントには、意味・単位・理由・既知の限界を書く。** 構文の説明は書かない。物理的な判断の解説は、コードではなく計画ファイルの `notes` に書く。
- **共通スクリプトは引数で試行を受け取る。** 試行固有の値をスクリプトに書き込まない。

## テンプレート

新しい試行では、`templates/` から必要なものをコピーして使う（コピーした時点の版が試行に残る）。

| ファイル | 用途 |
|---|---|
| [plan_reaction-energy.yaml](templates/plan_reaction-energy.yaml) | 反応・生成・吸着エネルギー型の計画（VASP） |
| [plan_convergence.yaml](templates/plan_convergence.yaml) | ENCUT・k 点などの収束試験型の計画 |
| [plan_mlip-screening.yaml](templates/plan_mlip-screening.yaml) | MLIP による構造緩和・スクリーニング型の計画 |
| [pipeline_vasp.py](templates/pipeline_vasp.py) | VASP 版の処理（`plan` / `submit` / `run` / `collect`）。人間が流れを読むファイル |
| [pipeline_mlip.py](templates/pipeline_mlip.py) | MLIP 版の処理（`plan` / `run`。`ase.db` に保存） |
| [recipes.py](templates/recipes.py) | 構造の作り方・候補の生成・解析の名前付きの手順。AI が手順を足すファイル |
| [job.slurm](templates/job.slurm) | Slurm ジョブの雛形（終了予告で VASP を止める処理付き） |
| [trial_README.md](templates/trial_README.md) | 試行の README |
| [theme_README.md](templates/theme_README.md) | 研究テーマの README（試行の索引） |

計画の型がどれにも合わないときは、[plan-format.md](references/plan-format.md) の共通骨格から、その計算に適した型を作る。

## やってはいけないこと

- 人間の OK の前に `sbatch` する、または長時間の計算を始める
- 未収束の値を結果の表や図に載せる（載せる場合は「未収束」と明記して区別する）
- 既存の計算結果を削除・上書きする
- 推測した原因を、検証済みの事実として書く
- 計画レビューに一般論の解説を長々と書く
- 既知の落とし穴（[pitfalls.md](references/pitfalls.md)）を確認せずに条件を決める

## 参照

- [references/review-points.md](references/review-points.md)：計画レビュー、確認カード、結果の表、結論前の確認の型
- [references/plan-format.md](references/plan-format.md)：計画ファイルの共通骨格と型の作り方
- [references/project-layout.md](references/project-layout.md)：フォルダ構成、README、状態、スクリプトの使い回し
- [references/learning.md](references/learning.md)：助言と学習メモの書き方
- [references/pitfalls.md](references/pitfalls.md)：DFT・VASP・ASE・Slurm・MLIP の落とし穴
