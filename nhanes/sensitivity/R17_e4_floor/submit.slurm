#!/bin/bash
# =====================================================================
# R17_e4_floor/submit.slurm  —  E4 floor sensitivity + share-at-floor for N6
#   (spec item A9b) as a SLURM ARRAY job.
# Mirrors Nhanes/arc_runs/R05_harmonized_floor/submit.slurm resources/env;
# changes ONLY array/driver/R17 paths. Submit from the repo root
# (= SLURM_SUBMIT_DIR = scratch checkout):
#   sbatch Nhanes/arc_runs/R17_e4_floor/submit.slurm
#
# CELL GRID (one task = one (example, floor) cell):
#   task 1: E2, floor 0.05    (N6 share row)
#   task 2: E3, floor 0.05    (N6 share row)
#   task 3: E4, floor 0.05    (N6 share row + sensitivity anchor)
#   task 4: E4, floor 0.025   (sensitivity)
#   task 5: E4, floor 0.01    (sensitivity)
#   FULL : --array=1-5   (the committed value below)
#   TEST : --array=3     (single E4/f050 cell; cheapest end-to-end ARC check)
# Re-submitting the FULL array is idempotent: completed cells skip in seconds
# via arc_skip_if_done().
#
# Walltime: each cell = B=20 splits of ONE arm (primary Fully-Aware-CF,
# LIB = SL.glm+SL.earth+SL.glmnet) -> roughly 1/5 of the locked 5-arm ladder
# cell (~1-2 h there). E2 (n_domain~29k) is the heaviest cell; E3/E4 (~11k)
# are light. 4 h is a safe limit; after a test run `seff <jobid>` and tighten.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R17_e4_floor
#SBATCH --array=1-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/R17_e4_floor_%A_%a.out
#SBATCH --error=sim_output/logs/R17_e4_floor_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p sim_output/logs

# --- proven ARC R environment (same as codes/sim.slurm, Nhanes/nhanes.slurm) ---
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: engine in codes/, data in Nhanes/analytic, R17 outputs isolated ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export R17_DIR="$SLURM_SUBMIT_DIR/Nhanes/arc_runs/R17_e4_floor"
export NH_ANA="$SLURM_SUBMIT_DIR/Nhanes/analytic"
export R17_OUT="$SLURM_SUBMIT_DIR/Nhanes/nhanes_output/arc_runs/R17"        # NON-destructive
export R17_MANIFEST="$SLURM_SUBMIT_DIR/Nhanes/nhanes_output/arc_runs/R17/manifest"

# --- analysis size (primary library is hardcoded in run.R; no override here) ---
export NHANES_B=20                      # cross-fit-split repeats (uses the cores)

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R17_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
