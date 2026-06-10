#!/bin/bash
# =====================================================================
# submit.slurm  —  R07_informative (informative-sampling-beyond-C stress test)
# Mirrors codes/sim.slurm resources EXACTLY (YOUR_ALLOCATION, gcc/9.4.0 + r/4.4.0,
# 32 cpus, 64G). Only the array/time and the driver path differ.
# Submit from the repo root (= SLURM submit dir / scratch):
#     sbatch codes/arc_runs/R07_informative/submit.slurm
#
# One array task = one (rho x rung x rep-chunk) cell.
#   #tasks = |rho| x |rung| x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL: rho in {0,0.3,0.6,0.9} (4) x rung in {L2_smooth,L3_adaptive} (2)
#         x ceil(1000/100)=10  =  4*2*10 = 80  -> --array=1-80
#   (rho/rung/rep order is the cell grid built inside run.R; keep --array,
#    SIM_N_REPS, SIM_CHUNK consistent -- run.R stops on a bad index.)
#
# Walltime: each rep runs TWO cross-fitted arms (unwt + wt OOF). The heavier
# rung is L3_adaptive (default RF in the SL library). Budget conservatively at
# ~2x a single canonical CF arm; 4 h is a safe, backfill-friendly ceiling for a
# 100-rep chunk on 32 cores. Tighten after the first task with `seff <jobid>`.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R07_informative
#SBATCH --array=1-80
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/R07_informative_%A_%a.out
#SBATCH --error=sim_output/logs/R07_informative_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p sim_output/logs

# --- proven ARC R environment (codes/sim.slurm) ---
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: everything organized UNDER the submit dir (scratch) ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export R07_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R07_informative"
export R07_OUT="$DATA_ROOT/arc_runs/R07_informative"

# --- run size (edit with --array above) ---
export SIM_N_REPS=1000                        # FULL=1000
export SIM_CHUNK=100                          # FULL=100

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R07_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
