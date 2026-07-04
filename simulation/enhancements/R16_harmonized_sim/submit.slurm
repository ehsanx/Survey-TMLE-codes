#!/bin/bash
# =====================================================================
# submit.slurm  —  R16_harmonized_sim (harmonized-truncation headline grid, A9a)
# NON-DESTRUCTIVE enhancement run. Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch codes/arc_runs/R16_harmonized_sim/submit.slurm
#
# One array task = one (scenario x rung x rep-chunk) cell  (mirrors run_sim.R).
#   #tasks = 2 scenarios x 4 rungs x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL : SIM_N_REPS=1000, SIM_CHUNK=100 -> 2*4*10 = 80  -> --array=1-80
#   TEST : --array=1,11,21,31,41,51,61,71   (chunk 1 of each of the 8 cells)
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
# Re-submitting the FULL array is idempotent (arc_skip_if_done per-chunk
# checkpoint): completed chunks skip in seconds.
#
# Grid task order (cell-major; cells = expand.grid scenario-fastest, then rung):
#   tasks  1-10 : standard x L1_param        chunks 1-10
#   tasks 11-20 : R1       x L1_param        chunks 1-10
#   tasks 21-30 : standard x L2_smooth       chunks 1-10
#   tasks 31-40 : R1       x L2_smooth       chunks 1-10
#   tasks 41-50 : standard x L3_adaptive     chunks 1-10
#   tasks 51-60 : R1       x L3_adaptive     chunks 1-10
#   tasks 61-70 : standard x L4_aggressive   chunks 1-10
#   tasks 71-80 : R1       x L4_aggressive   chunks 1-10
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R16_harmonized_sim
#SBATCH --array=1-80
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/R16_harmonized_sim_%A_%a.out
#SBATCH --error=sim_output/logs/R16_harmonized_sim_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; mkdir -p sim_output/logs
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: everything organized UNDER the submit dir (scratch); no clobber ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export R16_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R16_harmonized_sim"
export R16_OUT="$DATA_ROOT/arc_runs/R16_harmonized_sim"

# --- run size (edit with --array above) ---
export SIM_N_REPS=1000
export SIM_CHUNK=100
# --- harmonized knobs (defaults match run.R; uncomment to override) ---
# NOTE: changing the floor/q knobs requires a FRESH R16_OUT (or deleting the old
# chunks) -- the per-chunk checkpoint cannot distinguish knob values, so a
# re-submit into the same dir would skip old-knob chunks and silently pool a
# mixed-floor run (aggregate.R refuses such pools).
# export R16_GFLOOR=0.05
# export R16_QLO=0.001

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R16_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
