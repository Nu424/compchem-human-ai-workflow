#!/bin/bash -l
# 計算ジョブの雛形（compchem-human-ai-workflow）。資源（account, ntasks, time など）は pipeline の submit が
# sbatch のコマンドラインで渡す。ここには、どのジョブでも同じ指定と、クラスタ固有の環境だけを書く。
#   使い方: sbatch [資源] job.slurm <pipeline の .py> <コマンド> <試行> [系]
#SBATCH --nodes=1
#SBATCH --signal=B:TERM@120
set -euo pipefail

# ---- クラスタ固有の環境（新しいクラスタ・新しい VASP では必ず確認する） ----
module unload openmpi || true                       # 既定の MPI を外す
module load <compiler>/<version> <mpi>/<version>    # VASP をビルドしたコンパイラと MPI に揃える
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1      # MPI の各プロセスがスレッドを立ててコアを奪い合わないように
export PYTHONNOUSERSITE=1                            # ~/.local のパッケージを無視し、venv だけを使う
PYTHON=/path/to/.venv/bin/python
export VASP_PP_PATH=/path/to/potentials              # ASE が POTCAR を組み立てる PAW の場所
export ASE_VASP_COMMAND="mpiexec -np ${SLURM_NTASKS:-1} /path/to/vasp_std"   # Python は 1 つ、VASP だけ MPI

# ---- 終了予告で VASP を止める ----
# Slurm は --signal=B:TERM@120 で、制限時間の 120 秒前にこの bash へ SIGTERM を送る。
# bash はフォアグラウンドのコマンドの実行中は trap を処理しないので、Python を & で起動して wait で待つ。
# 転送された SIGTERM を受けた Python が、実行中の計算フォルダに STOPCAR（LSTOP）を書く。
"$PYTHON" "$@" &
worker=$!
trap 'kill -TERM "$worker" 2>/dev/null || true' TERM
# 失敗・中断なら終了コードを保存し、Python の後始末を待ってから、その終了コードで終わる。
wait "$worker" || { code=$?; wait "$worker" || true; exit "$code"; }
