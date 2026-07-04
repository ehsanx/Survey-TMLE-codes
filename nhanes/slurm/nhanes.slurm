#!/bin/bash
# =====================================================================
# nhanes.slurm  —  NHANES application as a SLURM ARRAY job on UBC ARC (Sockeye)
# Same proven environment as codes/sim.slurm (YOUR_ALLOCATION, gcc/9.4.0 + r/4.4.0).
# Submit from /scratch/YOUR_ALLOCATION/GitHub/survey-tmle2/ :  sbatch Nhanes/nhanes.slurm
#
# One array task = one (example x SL-library-rung) cell.
#   #tasks = 4 examples x 4 rungs = 16  ->  --array=1-16
# Each task runs the full 5-arm suite at that rung and repeats it over NHANES_B
# random cross-fitting / SuperLearner splits across the 32 cores (split-stability
# diagnostic). The example x rung grid is the APPLIED echo of simulation Figure 1.
#
# Walltime: most cells finish in minutes; the deep-RF rung (L4) on the ~30k-row
# examples (E1/E2) is the heaviest. 2 h is a safe limit; after a test run `seff
# <jobid>` and tighten, or lower NHANES_B / split the heavy cells out.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=nhtmle
#SBATCH --array=1-16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
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

# --- analysis size ---
export NHANES_B=20                      # cross-fit-split repeats (uses the cores; lower if memory-tight)

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$SIM_CODE/../Nhanes/R/nhanes_arc.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
