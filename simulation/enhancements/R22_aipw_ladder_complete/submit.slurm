#!/bin/bash
# =====================================================================
# submit.slurm  —  R22_aipw_ladder_complete (the alpha completion: AIPW arms
# at the ladder/ensemble rungs L2/L3/L5 that R15/R20 skipped, so the unified
# L1-L7 x {TMLE,AIPW} appendix figure is gap-free except the principled
# L1/L4 internal-CV cells). 6 fold-protocol arms/cell: TMLE-{SF,CV,CF} +
# AIPW-{SF,CV,CF}. NON-DESTRUCTIVE (sources engine + R15 + R21 helpers
# read-only; no xgboost/HAL).
#
# Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch codes/arc_runs/R22_aipw_ladder_complete/submit.slurm
#
# One array task = one (scenario x rung x rep-chunk) cell (run.R grid).
#   #tasks = 2 scenarios x 3 rungs x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL : SIM_N_REPS=1000, SIM_CHUNK=100 -> 2*3*10 = 60 -> --array=1-60
#   TEST : SIM_N_REPS=100,  SIM_CHUNK=100 -> 2*3*1  =  6 -> --array=1-6
#   cell order (scenario fastest within rung): tasks
#     1-10  std-L2 | 11-20 R1-L2 | 21-30 std-L3 | 31-40 R1-L3 |
#     41-50 std-L5 | 51-60 R1-L5      <-- L5 is the expensive block (deep RF)
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
#
# Walltime 12:00:00 (conservative). L2/L3 are cheap (no deep RF); L5 ~ R20's L5
# cost x ~2 (TMLE 3-way + a full AIPW estimator), still well under 12 h for a
# 100-rep chunk. PROBE the expensive block first:
#   sbatch --array=41 --time=12:00:00 codes/arc_runs/R22_aipw_ladder_complete/submit.slurm
#   sacct -j <id> --format=Elapsed   # then run the rest at ~3x observed
# Re-submitting the FULL array is idempotent (per-chunk checkpoint). If a chunk
# times out, resubmit ONLY that task id at the SAME SIM_CHUNK with a longer
# --time (never change SIM_CHUNK against existing files).
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R22_aipw_ladder_complete
#SBATCH --array=1-60
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=sim_output/logs/R22_aipw_ladder_complete_%A_%a.out
#SBATCH --error=sim_output/logs/R22_aipw_ladder_complete_%A_%a.err
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
export R22_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R22_aipw_ladder_complete"
export R21_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R21_deployable_cvcf"     # r21_helpers.R (aipw_cv_arm/aipw_arms_all)
export R15_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R15_aipw_benchmark"      # aipw_helpers.R (aipw_arms, .aipw_jkn)
export R22_OUT="$DATA_ROOT/arc_runs/R22_aipw_ladder_complete"

export SIM_N_REPS=1000                        # TEST=100 ; FULL=1000
export SIM_CHUNK=100                          # TEST=100 ; FULL=100

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R22_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
