#!/bin/bash
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R02_largem_sweep
#SBATCH --array=1-80
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=sim_output/logs/R02_largem_sweep_%A_%a.out
#SBATCH --error=sim_output/logs/R02_largem_sweep_%A_%a.err
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
# ---- run-specific env vars ----
# Array = 4 m-levels {6,12,20,30} x 20 rep-chunks (1000 reps / 50 per chunk) = 80 tasks.
# ALL 4 rungs (L4_aggressive, RF_shallow, L1_param, L3_adaptive) run inside each task.
# Chunk halved 100->50 (and --array doubled 40->80, --time 6h->8h) so the deep-RF L4
# cell at m=300 fits the walltime; the per-chunk checkpoint makes re-submits idempotent.
export SIM_N_REPS=1000
export SIM_CHUNK=50
Rscript "$SLURM_SUBMIT_DIR/simulation/enhancements/R02_largem_sweep/run.R"
