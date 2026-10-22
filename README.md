# compchem-human-ai-workflow

計算化学シミュレーション（ASE・VASP・MLIP・Slurm）を、人間と AI エージェントが協働して計画・実装・実行・確認・記録するための Agent Skill です。

**計算するのは AI、問いを立てて結論に責任を持つのは人間。** 人間はコードではなく計画と証拠をレビューし、AI は実装・投入・監視・集計・記録を担います。

## 特徴

- **投入前の計画レビュー**：`sbatch` や長時間の計算の前に、AI が 10 分で読める計画レビューを提示し、人間の OK を待ちます
- **確認カード**：試し計算の後や結論の前に、判断に必要な材料だけをまとめて見せます
- **5 つの保証**：未収束の結果を使わない、既存の結果を上書きしない、計算条件を添えて保存する、何を実行したかを記録する、投入前に人間の OK をもらう
- **テンプレート**：反応・吸着エネルギー、収束試験、MLIP スクリーニングの計画ファイルと、VASP 版・MLIP 版のパイプライン、Slurm ジョブの雛形
- **落とし穴集**：DFT・VASP・ASE・Slurm・MLIP でよくある失敗を、条件を決める前に確認します

## インストール

[skills](https://github.com/vercel-labs/skills) CLI でインストールできます。

```bash
npx skills add Nu424/compchem-human-ai-workflow
```

使っているエージェントを選ぶ画面が出ます。直接指定する場合は `-a` を、すべてのプロジェクトで使う場合は `-g` を付けます。

```bash
npx skills add Nu424/compchem-human-ai-workflow -a claude-code
npx skills add Nu424/compchem-human-ai-workflow -a cursor -g
```

手動で入れる場合は、`skills/compchem-human-ai-workflow/` をエージェントのスキル用フォルダ（Claude Code なら `.claude/skills/`）にコピーしてください。

## 使い方

インストール後は、計算化学の作業を頼むとエージェントが自動でこのスキルを参照します。たとえば次のように話しかけます。

- 「Pt(111) への H 吸着エネルギーを VASP で計算したい。計画を立てて」
- 「ENCUT と k 点の収束試験をしたい」
- 「MACE で候補構造を緩和してスクリーニングしたい」
- 「終わった計算の結果を集計して、確認カードを見せて」

進め方はおおむね次のとおりです。

1. AI が規模を見て、計画を立てるかを提案する
2. AI が計画ファイルを書き、計画レビューを提示する → **人間が OK を出す**
3. AI がテンプレートの型で実装し、代表的な系で試し計算をする → 確認カードを見る
4. 本計算 → 確認カードと結果の表を見る
5. AI がまとめを書き、人間が結論を確定する
6. 試行の README とテーマの索引、学習メモを更新する

詳しくは [SKILL.md](skills/compchem-human-ai-workflow/SKILL.md) を参照してください。

## 構成

```text
skills/compchem-human-ai-workflow/
  SKILL.md        理念、役割分担、進め方、5 つの保証、実装の型
  CHANGELOG.md    変更履歴
  references/     確認ポイント、計画ファイルの書き方、フォルダ構成、学びの助言、落とし穴
  templates/      計画ファイル、パイプライン、手順集、Slurm ジョブ、README の雛形
```

## 動作確認の状況

テンプレートは ASE の EMT 計算機による代用でのみ動作を確認しています。実機の VASP・MACE では未確認です。使う前に、お使いの計算環境に合わせて確認・調整してください。

## ライセンス

[MIT License](LICENSE)
