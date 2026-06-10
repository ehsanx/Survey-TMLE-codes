#!/bin/bash
# =====================================================================
# sim.slurm  —  survey-TMLE simulation as a SLURM ARRAY job on UBC ARC (Sockeye)
# Matches the proven arc/job8.sh environment (YOUR_ALLOCATION, gcc/9.4.0 + r/4.4.0).
# Submit from /scratch/YOUR_ALLOCATION/GitHub/survey-tmle2/ :  sbatch codes/sim.slurm
#
# One array task = one (scenario x library-rung x rep-chunk) cell.
#   #tasks = 2 scenarios x 4 rungs x ceil(SIM_N_REPS / SIM_CHUNK)
#   TEST: SIM_N_REPS=50,   SIM_CHUNK=50   -> 2*4*1  =  8  -> --array=1-8
#   FULL: SIM_N_REPS=1000, SIM_CHUNK=100  -> 2*4*10 = 80  -> --array=1-80
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run_sim.R stops on a bad index).
#
# Walltime: the heaviest cell (L4 deep-RF / R1) is ~440 s/rep single-core; at
# CHUNK=100 over 32 cores that's ~23 min/task (allow ~2x for ARC CPUs => ~45 min).
# 1 h is a safe, backfill-friendly limit -- do NOT inflate it (over-requesting
# delays scheduling). After the TEST, run `seff <jobid>` and tighten if desired.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=svytmle
#SBATCH --array=1-80                         # TEST=1-8 ; FULL=1-80
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=sim_output/logs/svytmle_%A_%a.out
#SBATCH --error=sim_output/logs/svytmle_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p sim_output/logs

# --- proven ARC R environment (arc/job8.sh) ---
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: everything organized UNDER the submit dir (scratch) ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export SIM_OUT="$DATA_ROOT/intermediate"

# --- run size (edit with --array above) ---
export SIM_N_REPS=1000                        # TEST=50 ; FULL=1000
export SIM_CHUNK=100                         # TEST=50 ; FULL=100

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$SIM_CODE/run_sim.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
