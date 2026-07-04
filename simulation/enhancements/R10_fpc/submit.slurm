#!/bin/bash
# =====================================================================
# R10_fpc/submit.slurm — FPC / first-stage-fraction sensitivity (ARC array job)
# Submit from the repo root (== SLURM_SUBMIT_DIR == scratch):
#     sbatch simulation/enhancements/R10_fpc/submit.slurm
#
# Task grid (run.R builds it; 5 tasks):
#   task 1      PART A : R1 design, rungs {L1_param, L2_smooth}, fpc vs no-fpc
#   task 2..5   PART B : first-stage-fraction sweep J in {8,16,32,64}, L1_param CF
#
# Resources mirror codes/sim.slurm exactly (32 cpus, 64G). Time bumped to 02:00:00
# because PART A runs the engine TWICE per rep (L1+L2) plus a local CF re-fit, and
# the sweep grows N at J=64 (N = 50*64*300 = 960k). Each task is single-cell and
# embarrassingly parallel over its reps, so 2 h is a safe backfill-friendly cap.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R10_fpc
#SBATCH --array=1-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=sim_output/logs/R10_fpc_%A_%a.out
#SBATCH --error=sim_output/logs/R10_fpc_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; mkdir -p sim_output/logs
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/
export REPO_ROOT="$SLURM_SUBMIT_DIR"; export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"; export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
# run-specific env vars
export FPC_RUN_DIR="$SLURM_SUBMIT_DIR/simulation/enhancements/R10_fpc"
export FPC_N_REPS=1000

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$FPC_RUN_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
