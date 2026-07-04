#!/bin/bash
# =====================================================================
# submit.slurm  —  R13_null_typeI (null / type-I + small-effect power; spec A5)
# NON-DESTRUCTIVE enhancement run. Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch codes/arc_runs/R13_null_typeI/submit.slurm
#
# One array task = one (te x scenario x rung x rep-chunk) cell (run.R grid).
#   #tasks = 2 te x 2 scenarios x 4 rungs x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL : SIM_N_REPS=1000, SIM_CHUNK=100 -> 16*10 = 160 -> --array=1-160
#   TEST : SIM_N_REPS=100,  SIM_CHUNK=100 -> 16*1  =  16 -> --array=1-16
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
#
# Task order (run.R): cells are te-major (tasks 1-80 te=0, 81-160 te=0.3);
# within te, expand.grid(scenario, rung) scenario-fastest; chunks fastest:
#   tasks   1- 10 te=0.0 standard L1_param ... tasks  71- 80 te=0.0 R1 L4_aggressive
#   tasks  81- 90 te=0.3 standard L1_param ... tasks 151-160 te=0.3 R1 L4_aggressive
#
# Walltime: 04:00:00 — the L4 deep-RF cells are binding (headline L4 chunk was
# ~23-45 min at 32 cores; 4 h is a conservative cap for ARC CPU variability).
# Cheaper cells exit early. Re-submitting the FULL array is idempotent: the
# per-chunk checkpoint (arc_skip_if_done) skips completed chunks in seconds.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R13_null_typeI
#SBATCH --array=1-160
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/R13_null_typeI_%A_%a.out
#SBATCH --error=sim_output/logs/R13_null_typeI_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"; mkdir -p sim_output/logs

# --- proven ARC R environment (matches codes/sim.slurm) ---
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: everything organized UNDER the submit dir (scratch); no clobber ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export R13_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R13_null_typeI"
export R13_OUT="$DATA_ROOT/arc_runs/R13_null_typeI"

# --- run size (edit with --array above) ---
export SIM_N_REPS=1000                        # TEST=100 ; FULL=1000
export SIM_CHUNK=100                          # TEST=100 ; FULL=100

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R13_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
