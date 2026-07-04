#!/bin/bash
# =====================================================================
# R01_simple_control / submit.slurm
#   Pipeline-correctness control: L1 (+L2) under model_type='simple' for both
#   scenarios, plus L1 under 'complex' side-by-side. All five arms, 1000 reps.
#   Mirrors codes/sim.slurm resources EXACTLY; only array/time/mem-comment and
#   the driver path differ. Submit from the repo root (= SLURM_SUBMIT_DIR):
#       sbatch simulation/enhancements/R01_simple_control/submit.slurm
#
# GRID (run.R): cells = (2 scenarios) x (3 model_type/rung pairs) = 6 cells.
#   pairs: simple/L1_param, simple/L2_smooth, complex/L1_param.
#   #tasks = 6 * ceil(SIM_N_REPS / SIM_CHUNK)
#     FULL: SIM_N_REPS=1000, SIM_CHUNK=100 -> 6*10 = 60 -> --array=1-60
#     TEST: SIM_N_REPS=100,  SIM_CHUNK=100 -> 6*1  =  6 -> --array=1-6
#   Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
#
# Walltime: all cells here are L1 (GLM) or L2 (glm+gam+earth) -- far cheaper than
# the headline L4 deep-RF. 1 h is generous; tighten via `seff <jobid>` if wanted.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R01_simple_control
#SBATCH --array=1-60                          # TEST=1-6 ; FULL=1-60
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=sim_output/logs/R01_simple_control_%A_%a.out
#SBATCH --error=sim_output/logs/R01_simple_control_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p sim_output/logs

# --- proven ARC R environment (matches codes/sim.slurm) ---
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: everything organized UNDER the submit dir (scratch) ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
# private per-task output dir; NEVER the locked sim_output/intermediate
export ARC_OUT="$DATA_ROOT/arc_runs/R01_simple_control"

# --- run size (edit with --array above) ---
export SIM_N_REPS=1000                        # TEST=100 ; FULL=1000
export SIM_CHUNK=100                          # TEST=100 ; FULL=100

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$SIM_CODE/arc_runs/R01_simple_control/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
