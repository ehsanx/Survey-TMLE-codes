#!/bin/bash
# =====================================================================
# Nhanes/arc_runs/R09_sensitivity/submit.slurm  —  OPTIONAL convenience submit.
#
# R09_sensitivity is a LOCAL run: it is a covariate-set swap loop + one design
# tabulation on the primary library, and finishes on a LAPTOP in a few
# minutes (see NOTES.md). You do NOT need SLURM. This script is provided only for
# convenience if you would rather run it on ARC (e.g. to keep all enhancement
# runs uniform). It is a SINGLE task (no array): the three sensitivity cells +
# the design check run serially inside one Rscript call.
#
# Submit from the repo root (SLURM submit dir / scratch):
#   sbatch Nhanes/arc_runs/R09_sensitivity/submit.slurm
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R09_sensitivity
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/R09_sensitivity_%A_%a.out
#SBATCH --error=sim_output/logs/R09_sensitivity_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p sim_output/logs

# --- proven ARC R environment (same as codes/sim.slurm, Nhanes/nhanes.slurm) ---
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: engine in codes/, imputed data under Nhanes/analytic ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export NH_ANA="$SLURM_SUBMIT_DIR/Nhanes/analytic"
export R09_OUT="$SLURM_SUBMIT_DIR/Nhanes/nhanes_output/arc_runs/R09"
export R09_RESULTS="$SLURM_SUBMIT_DIR/results/arc"
export R09_MANIFEST="$R09_OUT/manifest"

echo "Job $SLURM_JOB_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$SLURM_SUBMIT_DIR/Nhanes/arc_runs/R09_sensitivity/run.R"
echo "finished $(date)"
