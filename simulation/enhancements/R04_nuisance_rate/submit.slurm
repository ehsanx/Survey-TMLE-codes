#!/bin/bash
# =====================================================================
# submit.slurm  —  R04_nuisance_rate (nuisance-rate diagnostic ladder)
# NON-DESTRUCTIVE enhancement run. Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch simulation/enhancements/R04_nuisance_rate/submit.slurm
#
# One array task = one (scenario x rung x rep-chunk) cell  (mirrors run_sim.R).
#   #tasks = 2 scenarios x 4 rungs x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL : SIM_N_REPS=200, SIM_CHUNK=50  -> 2*4*4 = 32  -> --array=1-32
#   (L4 deep-RF is the heaviest rung; 200 reps is ample for stable mean L2 errors
#    -- this is a diagnostic, not a coverage table, so fewer reps than the main sim.)
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
#
# Grid task order (expand.grid scenario fastest, then rung, then chunk):
#   tasks  1- 8 : standard/R1 x L1_param  (chunks 1-4 interleaved by scenario)
#   ... up to ...
#   tasks 25-32 : standard/R1 x L4_aggressive
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R04_nuisance_rate
#SBATCH --array=1-32
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/R04_nuisance_rate_%A_%a.out
#SBATCH --error=sim_output/logs/R04_nuisance_rate_%A_%a.err
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
export R04_DIR="$SLURM_SUBMIT_DIR/simulation/enhancements/R04_nuisance_rate"
export R04_OUT="$DATA_ROOT/arc_runs/R04_nuisance_rate"

# --- run size (edit with --array above) ---
export SIM_N_REPS=200
export SIM_CHUNK=50

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R04_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
