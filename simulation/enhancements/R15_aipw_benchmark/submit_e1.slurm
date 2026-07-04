#!/bin/bash
# =====================================================================
# submit_e1.slurm  —  R15_aipw_benchmark NHANES E1 driver (spec A7, single task)
# Survey-weighted AIPW (SF + CF arms) on the E1 domain with JKn replicate +
# Eq-8 linearized SEs, B=NHANES_B repeated cross-fit splits.
# Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch codes/arc_runs/R15_aipw_benchmark/submit_e1.slurm
#
# SINGLE task (no --array). The B splits parallelize across the 32 cores.
# E1 is the largest domain (~31k rows); with the primary library
# (SL.glm+SL.earth+SL.glmnet) one split is roughly comparable to one locked
# run_estimators() call, and all B=20 splits run concurrently -> well inside
# the 2 h walltime. If it times out, lower NHANES_B (the seeds are
# 20260607+1:B, so a lower B is a prefix of the same splits).
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R15_aipw_E1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=sim_output/logs/R15_aipw_E1_%j.out
#SBATCH --error=sim_output/logs/R15_aipw_E1_%j.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; mkdir -p sim_output/logs
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: engine in codes/, data in Nhanes/analytic, R15 outputs isolated ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export R15_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R15_aipw_benchmark"
export NH_ANA="$SLURM_SUBMIT_DIR/Nhanes/analytic"
export R15_E1_OUT="$SLURM_SUBMIT_DIR/Nhanes/nhanes_output/arc_runs/R15_e1"   # NON-destructive

# --- analysis size ---
export NHANES_B=20                       # cross-fit-split repeats (uses the cores)

echo "Job $SLURM_JOB_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R15_DIR/nhanes_e1.R"
echo "Job $SLURM_JOB_ID finished $(date)"
