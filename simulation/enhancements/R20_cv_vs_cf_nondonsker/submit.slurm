#!/bin/bash
# =====================================================================
# submit.slurm  —  R20_cv_vs_cf_nondonsker (CV vs CF at a non-Donsker
# multi-learner library L5 = glm + earth + deep RF; reviewer follow-up).
# NON-DESTRUCTIVE. Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch codes/arc_runs/R20_cv_vs_cf_nondonsker/submit.slurm
#
# One array task = one (scenario x rep-chunk) cell (run.R grid).
#   #tasks = 2 scenarios x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL : SIM_N_REPS=1000, SIM_CHUNK=100 -> 2*10 = 20 -> --array=1-20
#   TEST : SIM_N_REPS=100,  SIM_CHUNK=100 -> 2*1  =  2 -> --array=1-2
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
#   tasks 1-10 standard L5 ; tasks 11-20 R1 L5
#
# Walltime: 06:00:00 — L5 includes the deep RF, so per-rep cost matches the
# heaviest headline rung (L4) PLUS the extra internal-CV arm; 6 h is a
# conservative cap. Re-submitting the FULL array is idempotent (per-chunk
# checkpoint). If a chunk times out, resubmit ONLY that task id at the same
# SIM_CHUNK with a longer --time (never change SIM_CHUNK against existing files).
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R20_cv_vs_cf_nondonsker
#SBATCH --array=1-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=sim_output/logs/R20_cv_vs_cf_nondonsker_%A_%a.out
#SBATCH --error=sim_output/logs/R20_cv_vs_cf_nondonsker_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; mkdir -p sim_output/logs

module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export R20_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R20_cv_vs_cf_nondonsker"
export R20_OUT="$DATA_ROOT/arc_runs/R20_cv_vs_cf_nondonsker"

export SIM_N_REPS=1000                        # TEST=100 ; FULL=1000
export SIM_CHUNK=100                          # TEST=100 ; FULL=100

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R20_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
