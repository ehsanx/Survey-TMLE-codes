#!/bin/bash
# =====================================================================
# submit.slurm  —  R19_rate_sweep (sqrt(m)-scaled product-error diagnostic
#                  across the F.5 m-sweep; crosses R04 instrumentation x R02 m knob)
# NON-DESTRUCTIVE enhancement run. Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch codes/arc_runs/R19_rate_sweep/submit.slurm
#
# One array task = one (base_m x rung x rep-chunk) cell at FIXED scenario
# 'standard' (mirrors run_sim.R/R04 task shape).
#   #tasks = 4 base_m x 4 rungs x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL : SIM_N_REPS=200, SIM_CHUNK=50 -> 16 cells x 4 chunks = 64 -> --array=1-64
#   TEST : --array=1 (cheapest cell: L1_param x base_m=6, chunk 1) or
#          --array=61 (heaviest: L4_aggressive x base_m=30, n=7590, chunk 1)
#          to time the binding L4 corner. (Task 49 is L4 x base_m=6 -- the
#          LIGHTEST L4 cell, n=1518; the walltime-binding cells are 61-64.)
#   (200 reps suffices: this is a mean-L2-error diagnostic, not a coverage table
#    -- same rationale and size as R04.)
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
#
# Grid task order (cells row-major with chunk fastest; expand.grid puts base_m
# fastest within rung):
#   tasks  1-16 : L1_param      (m6 c1-4, m12 c1-4, m20 c1-4, m30 c1-4)
#   tasks 17-32 : L2_smooth      (same m/chunk layout)
#   tasks 33-48 : L3_adaptive
#   tasks 49-64 : L4_aggressive  (deep RF; heaviest, esp. m30 = tasks 61-64)
# Re-submitting the FULL array is idempotent (arc_skip_if_done per chunk).
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R19_rate_sweep
#SBATCH --array=1-64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/R19_rate_sweep_%A_%a.out
#SBATCH --error=sim_output/logs/R19_rate_sweep_%A_%a.err
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
export R19_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R19_rate_sweep"
export R19_OUT="$DATA_ROOT/arc_runs/R19_rate_sweep"
# R04's helpers are SOURCED read-only from R04's own folder:
export R04_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R04_nuisance_rate"

# --- run size (edit with --array above) ---
export SIM_N_REPS=200
export SIM_CHUNK=50

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R19_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
