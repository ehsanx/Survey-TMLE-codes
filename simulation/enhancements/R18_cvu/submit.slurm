#!/bin/bash
# =====================================================================
# submit.slurm  —  R18_cvu (unweighted-nuisance internal-CV sim arm, spec A4)
# NON-DESTRUCTIVE enhancement run. Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch codes/arc_runs/R18_cvu/submit.slurm
#
# One array task = one (scenario x rung x rep-chunk) cell; ONE cheap arm per rep.
#   #tasks = 2 scenarios x 2 rungs (L2_smooth, L3_adaptive) x ceil(SIM_N_REPS/SIM_CHUNK)
#   FULL : SIM_N_REPS=1000, SIM_CHUNK=100 -> 2*2*10 = 40  -> --array=1-40
#   TEST : --array=1,11,21,31   (one chunk per cell; verify before the full run)
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
#
# Grid task order (CELL-MAJOR; see run.R):
#   tasks  1-10 : standard x L2_smooth    (chunks 1-10)
#   tasks 11-20 : R1       x L2_smooth
#   tasks 21-30 : standard x L3_adaptive
#   tasks 31-40 : R1       x L3_adaptive
# Per-chunk checkpointing: re-submitting the FULL array is idempotent; finished
# chunks skip in seconds (arc_skip_if_done).
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R18_cvu
#SBATCH --array=1-40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=sim_output/logs/R18_cvu_%A_%a.out
#SBATCH --error=sim_output/logs/R18_cvu_%A_%a.err
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
export R18_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R18_cvu"
export R18_OUT="$DATA_ROOT/arc_runs/R18_cvu"

# --- run size (edit with --array above) ---
export SIM_N_REPS=1000
export SIM_CHUNK=100

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R18_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
