#!/bin/bash
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R03_isolation_2x2
#SBATCH --array=1-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/R03_isolation_2x2_%A_%a.out
#SBATCH --error=sim_output/logs/R03_isolation_2x2_%A_%a.err
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

# ---- run-specific env vars --------------------------------------------------
# Grid built in run.R: scenario=standard x rung{L4_aggressive,L1_param} x 10 chunks
#   = 20 array tasks. Tasks 1-10 = L4_aggressive chunks 1-10; 11-20 = L1_param.
export SIM_N_REPS=1000
export SIM_CHUNK=100
export R03_GFLOOR=0.05
export R03_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R03_isolation_2x2"
export R03_OUT="$SLURM_SUBMIT_DIR/sim_output/arc_runs/R03_isolation_2x2"

Rscript "$SLURM_SUBMIT_DIR/codes/arc_runs/R03_isolation_2x2/run.R"
