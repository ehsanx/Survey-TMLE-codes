#!/bin/bash
# J5_vfolds: number-of-folds (V) sensitivity for the CF arm
# FULL: --array=1-60 ; TEST: --array=1-6 with SIM_N_REPS=50 SIM_CHUNK=50
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=J5_vfolds
#SBATCH --array=1-60
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/J5_vfolds_%A_%a.out
#SBATCH --error=sim_output/logs/J5_vfolds_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=FAIL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p sim_output/logs sim_output/enhancements/J5_vfolds

# load your cluster's R module here, e.g.:
# module load gcc/9.4.0 r/4.4.0
# export R_LIBS=/path/to/your/R/library

export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/R"
export J_OUT="$SLURM_SUBMIT_DIR/sim_output/enhancements/J5_vfolds"
export SIM_N_REPS=${SIM_N_REPS:-1000}
export SIM_CHUNK=${SIM_CHUNK:-100}

Rscript simulation/enhancements/J5_vfolds/run.R
