# <YYMMDD>_<試行名>

**状態：** 計画中 <!-- 計画中 / 実行中 / 完了 / 採用 / 探索 / 置き換え済み（→ 試行名） / 棄却（理由） -->

## 目的

<何を知りたいか。1〜3 行>

## 計画の要約

- 計画ファイル：`plan.yaml`（型：<reaction-energy など>）
- 対象：<系の一覧を 1 行で>
- 主な条件：<汎関数、カットオフ、k 点など、議論に効くものだけ>
- 処理の流れ：<構造作成 → 緩和 → 1 点計算 → 集計>
- 既知の近似と未検証：<caveats の要約>
- 計画レビューの OK：<日付、確認した人>

## 実行

- スクリプト：`<scripts/pipeline_vasp.py>`（git <コミット>）
- 生データの場所：<クラスタ上のパス。ローカルにない場合>

```bash
python scripts/pipeline_vasp.py plan    trials/<試行>
python scripts/pipeline_vasp.py submit  trials/<試行>
python scripts/pipeline_vasp.py collect trials/<試行>   # 途中経過の確認にも使える
```

## 結果

<results/RESULTS.md の要点。値、比較対象、状態>

| 量 | 値 | 状態 |
|---|---:|---|

## 結論

<結論前の確認を済ませた主張だけ。仮説は「仮説」と明記する>

## 残った疑問・次の試行

- <...>

## ファイル

| パス | 中身 |
|---|---|
| `plan.yaml` | 計画 |
| `runs/<系>/<段階>/` | 計算の入出力、`result.json`（結果・条件・実行記録）、`final.extxyz`（最終構造） |
| `results/RESULTS.md` | 結果の表と確認カード |
| `results/results.db` | 全結果の ase.db（collect が作り直す） |
| `results/structures/` | 最終構造の画像 |
| `logs/` | Slurm の出力、`submissions.jsonl` |
