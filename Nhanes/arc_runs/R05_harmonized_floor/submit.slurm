#!/bin/bash
# =====================================================================
# R05_harmonized_floor/submit.slurm  —  NHANES harmonized-floor + weighted-OOF-CF
#   ladder on E3/E4 (de-confound the real-data ladder) as a SLURM ARRAY job.
# Mirrors Nhanes/nhanes.slurm resources/env; changes ONLY array/driver/NH_OUT.
# Submit from the repo root (= SLURM_SUBMIT_DIR = scratch checkout):
#   sbatch Nhanes/arc_runs/R05_harmonized_floor/submit.slurm
#
# GRID: {E3,E4} x {L1_param,L2_smooth,L3_adaptive,L4_aggressive} = 8 cells.
#   task t: example = c("E3","E4")[((t-1)%%2)+1]; rung = RUNGS[((t-1)%/%2)+1]
#   -> tasks 1,2 = L1 (E3,E4); 3,4 = L2; 5,6 = L3; 7,8 = L4.
# To run only the data-adaptive rungs L3,L4 (where divergence is largest):
#   change --array=1-8 to --array=5-8.
#
# Walltime: each cell runs B random splits of the FOUR arms (FA, FA-h05, FA-CF,
# FA-CF-wOOF) -> roughly the locked nhanes_arc cost + one extra cross-fit pass
# (the weighted-OOF arm). L4 deep-RF on the ~11k-row E3/E4 domains is heaviest.
# 4 h is a safe limit; after a test `seff <jobid>` and tighten, or lower NHANES_B.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R05_harmonized_floor
#SBATCH --array=1-8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/R05_harmonized_floor_%A_%a.out
#SBATCH --error=sim_output/logs/R05_harmonized_floor_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p sim_output/logs

# --- proven ARC R environment (same as codes/sim.slurm, Nhanes/nhanes.slurm) ---
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: engine in codes/, data in Nhanes/analytic, R05 outputs isolated ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export R05_DIR="$SLURM_SUBMIT_DIR/Nhanes/arc_runs/R05_harmonized_floor"
export NH_ANA="$SLURM_SUBMIT_DIR/Nhanes/analytic"
export NH_OUT="$SLURM_SUBMIT_DIR/Nhanes/nhanes_output/arc_runs/R05"          # NON-destructive
export NH_MANIFEST="$SLURM_SUBMIT_DIR/Nhanes/nhanes_output/arc_runs/R05/manifest"

# --- analysis size + harmonized floors ---
export NHANES_B=20                      # cross-fit-split repeats (uses the cores)
export R05_FA_GBOUND=0.05               # harmonized Fully-Aware EIF floor (= OOF floor)
export R05_OOF_BOUND=0.05               # OOF propensity floor (matches the locked CF arm)

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R05_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
