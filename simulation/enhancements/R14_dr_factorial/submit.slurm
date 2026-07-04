#!/bin/bash
# =====================================================================
# submit.slurm  —  R14_dr_factorial (double-robustness factorial, spec A6)
# NON-DESTRUCTIVE enhancement run. Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch codes/arc_runs/R14_dr_factorial/submit.slurm
#
# One array task = one (q_spec x g_spec x sampling x rung x rep-chunk) cell.
#   #tasks = 16 cells x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL : SIM_N_REPS=1000, SIM_CHUNK=100 -> 16*10 = 160 -> --array=1-160
#   TEST : --array=1-160:10  (chunk 1 of each of the 16 cells; grid is
#          CELL-MAJOR, so tasks 1-10 = cell 1, 11-20 = cell 2, ...)
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
#
# Cell order (expand.grid: q_spec fastest, then g_spec, sampling, rung):
#   cells  1- 4 : (C,C) (W,C) (C,W) (W,W)  x info    x L1_param
#   cells  5- 8 : the same four            x noninfo x L1_param
#   cells  9-16 : the same eight                     x L2_smooth
# All cells are CHEAP (GLM / smooth library; no deep RF) -> --time=02:00:00.
# Re-submitting the FULL array is idempotent (arc_skip_if_done per chunk).
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R14_dr_factorial
#SBATCH --array=1-160
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=sim_output/logs/R14_dr_factorial_%A_%a.out
#SBATCH --error=sim_output/logs/R14_dr_factorial_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; mkdir -p sim_output/logs
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: everything organized UNDER the submit dir (scratch); no clobber ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export R14_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R14_dr_factorial"
export R14_OUT="$DATA_ROOT/arc_runs/R14_dr_factorial"

# --- run size (edit with --array above) ---
export SIM_N_REPS=1000
export SIM_CHUNK=100

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R14_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
