#!/bin/bash
# =====================================================================
# nhanes_finish.slurm  —  re-run ONLY the four E1/E2 deep-learner cells that
# the main run (job 11622292) cancelled at the 2 h wall-clock limit.
#
# Background: the 16-task array is grid = expand.grid(example, rung) with
# example varying fastest, so:
#     task  9 = E1 x L3_adaptive      task 13 = E1 x L4_aggressive
#     task 10 = E2 x L3_adaptive      task 14 = E2 x L4_aggressive
# These four TIMED OUT (E1/E2 domains are ~30k rows, ~3x E3/E4); the deep
# random-forest rungs on 30k rows did not finish in 2 h. Nothing was wrong
# with the code or data — E3/E4 (the same rungs) completed in 40-57 min.
#
# This job is NON-DESTRUCTIVE: it writes only the four missing intermediates
#   nh_E1_L3_adaptive.rds, nh_E2_L3_adaptive.rds,
#   nh_E1_L4_aggressive.rds, nh_E2_L4_aggressive.rds
# alongside the 12 existing files. The completed E3/E4 + L1/L2 cells are untouched.
#
# Two changes vs nhanes.slurm:  (1) --array restricted to 9,10,13,14;
# (2) --time raised to 8 h and NHANES_B lowered 20 -> 5. With B=5 only 5
# split-workers run concurrently (min(cores,B)), so each deep-RF fit gets the
# node's memory without contention and runs at full speed; 8 h is ample headroom.
# (If you want E1/E2 deep rungs to match E3/E4's B=20 exactly, set NHANES_B=20
# and --time=16:00:00 — but B=5 is plenty for the point-estimate ladder figure.)
#
# Submit from /scratch/YOUR_ALLOCATION/GitHub/survey-tmle2/ :
#     sbatch Nhanes/nhanes_finish.slurm
# IMPORTANT: ensure the repo on ARC is current first (the very first attempt,
# job 11622228, failed because the nodes had a stale codes/estimators.R lacking
# the W_cols/nest/inpop args). `git -C /scratch/.../survey-tmle2 pull` before submit.
#
# AFTER it finishes, regenerate the aggregates + figures (now a full 4x4 grid):
#     Rscript Nhanes/R/aggregate_nhanes.R     # -> nhanes_output/results/*.csv + nhanes_combined.rds
#     Rscript Nhanes/R/figures_arc.R          # -> Nhanes/figures/ladder_FA_vs_CF.png, overlap_L4.png
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=nhtmle_fin
#SBATCH --array=9,10,13,14
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=nhanes_output/logs/nhtmle_%A_%a.out
#SBATCH --error=nhanes_output/logs/nhtmle_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p nhanes_output/logs

# --- proven ARC R environment (arc/job8.sh) ---
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: engine in codes/, data + outputs under Nhanes/ ---
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export NH_ANA="$SLURM_SUBMIT_DIR/Nhanes/analytic"
export NH_OUT="$SLURM_SUBMIT_DIR/Nhanes/nhanes_output/intermediate"
export NH_MANIFEST="$SLURM_SUBMIT_DIR/Nhanes/nhanes_output/manifest"

# --- analysis size: fewer split-repeats so the heavy E1/E2 deep rungs finish ---
export NHANES_B=5                       # was 20; 5 concurrent workers -> no memory contention

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$SIM_CODE/../Nhanes/R/nhanes_arc.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
