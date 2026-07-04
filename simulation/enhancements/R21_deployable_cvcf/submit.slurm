#!/bin/bash
# =====================================================================
# submit.slurm  —  R21_deployable_cvcf (deployable non-Donsker library L6 +
# the new AIPW cluster-CV arm; 6 arms/cell: TMLE-{SF,CV,CF} + AIPW-{SF,CV,CF}).
# NON-DESTRUCTIVE. Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch codes/arc_runs/R21_deployable_cvcf/submit.slurm
#
# One array task = one (scenario x rep-chunk) cell (run.R grid).
#   #tasks = 2 scenarios x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL : SIM_N_REPS=1000, SIM_CHUNK=100 -> 2*10 = 20 -> --array=1-20
#   TEST : SIM_N_REPS=100,  SIM_CHUNK=100 -> 2*1  =  2 -> --array=1-2
#   tasks 1-10 standard L6 ; tasks 11-20 R1 L6
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
#
# OPTIONAL m-sweep: set R21_MSWEEP=1 -> grid doubles to 4 cells (base_m {6,20})
#   -> --array=1-40 (gate it OFF for the core run; ship the 20-task core first).
# OPTIONAL HAL: set R21_WITH_HAL=1 to append SL.hal9001 (slow / may fail to
#   install on some ARC images; never make the headline depend on it).
#
# Walltime: 12:00:00 — L6 is the HEAVIEST run yet: a 5-member CV-SuperLearner
# (tuned xgboost + deep RF carry the cost) fit THREE ways on the TMLE side
# (single-fit, internal-CV, 5-fold CF) PLUS a full independent AIPW estimator
# (SF + 5-fold CF + internal-CV), all per rep x 100 reps/chunk. That is
# materially more than R20's 6 h (L5, TMLE-only). 12 h is a conservative cap;
# PROBE first (sbatch --array=1 --time=12:00:00 ...; read sacct Elapsed, then
# run the rest at ~3x). Re-submitting the FULL array is idempotent (per-chunk
# checkpoint). If a chunk times out, resubmit ONLY that task id at the same
# SIM_CHUNK with a longer --time (never change SIM_CHUNK against existing files).
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R21_deployable_cvcf
#SBATCH --array=1-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=sim_output/logs/R21_deployable_cvcf_%A_%a.out
#SBATCH --error=sim_output/logs/R21_deployable_cvcf_%A_%a.err
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
export R21_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R21_deployable_cvcf"
export R15_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R15_aipw_benchmark"   # aipw_helpers.R
export R21_OUT="$DATA_ROOT/arc_runs/R21_deployable_cvcf"

export SIM_N_REPS=1000                        # TEST=100 ; FULL=1000
export SIM_CHUNK=100                          # TEST=100 ; FULL=100
# export R21_MSWEEP=1                          # uncomment + set --array=1-40 for the base_m sweep
# export R21_WITH_HAL=1                        # uncomment to add SL.hal9001 (slow; optional)
# export R21_XGB_TUNE=1                        # uncomment to CV-tune xgb depth {3,6,8} (default: fixed 6)

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R21_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
