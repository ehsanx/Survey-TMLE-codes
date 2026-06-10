#!/bin/bash
# =====================================================================
# Nhanes/arc_runs/R06_mi/submit.slurm  —  R06_mi enhancement run on UBC ARC
#
# Multiple imputation (m=40) headline NHANES Table-2 intervals for the CERTIFIED
# PRIMARY estimator (3-learner near-Donsker LIB = SL.glm/SL.earth/SL.glmnet).
# m=40 follows White, Royston & Wood (2011): m >= the percentage of incomplete
# cases (max 29.2% for E4). The estimated FMI never exceeds 0.14, so this is
# conservative -- the headline numbers are insensitive to m beyond ~15.
#
# !! BEFORE RE-SUBMITTING at a new m: delete the cached per-example results, or
#    the per-example checkpoint in run.R will reuse the old (m=10) files:
#       rm -f Nhanes/nhanes_output/arc_runs/R06/mi_E*.rds
#
# One array task = one example. Each task imputes m=40 completed datasets and
# runs the 5-arm suite (LIB only, NOT the L1-L4 ladder) on each, then pools by
# Rubin's rules. Submit from the repo root (SLURM submit dir / scratch):
#     sbatch Nhanes/arc_runs/R06_mi/submit.slurm
#
# Resources mirror nhanes.slurm exactly (YOUR_ALLOCATION, gcc/9.4.0 + r/4.4.0,
# 32 cpus, 64G); only array/time/output paths differ. The 32 cores parallelize
# the m=40 imputations (capped to min(cores, m)=32 workers per task, ~2 waves).
# 5 h is a safe limit; E1/E2 (~30k-row domains) are the heaviest. After a test run
# use `seff <jobid>` and tighten.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R06_mi
#SBATCH --array=1-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=05:00:00
#SBATCH --output=sim_output/logs/R06_mi_%A_%a.out
#SBATCH --error=sim_output/logs/R06_mi_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p sim_output/logs

# --- proven ARC R environment (same as nhanes.slurm / sim.slurm) ---
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: engine in codes/, NHANES data + this run's outputs under Nhanes/ ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export NH_ANA="$SLURM_SUBMIT_DIR/Nhanes/analytic"
export R06_OUT="$SLURM_SUBMIT_DIR/Nhanes/nhanes_output/arc_runs/R06"
export R06_RESULTS="$SLURM_SUBMIT_DIR/results/arc"
export R06_MANIFEST="$SLURM_SUBMIT_DIR/Nhanes/nhanes_output/arc_runs/R06/manifest"

# --- run-specific knobs ---
export MI_M=40                          # number of imputations (>= % incomplete cases; White et al. 2011)
export MI_SEED=20260607                 # base mice seed (== 02_impute.R)

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$SLURM_SUBMIT_DIR/Nhanes/arc_runs/R06_mi/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
