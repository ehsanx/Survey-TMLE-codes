#!/bin/bash
# =====================================================================
# R12_highdeff / submit.slurm  —  High-DEFF "Design C" stress run on UBC ARC
#
# One array task = one (Design C rung x rep-chunk) cell.
#   #tasks = 3 rungs x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL : SIM_N_REPS=1000, SIM_CHUNK=100 -> 3*10 = 30 -> --array=1-30
#   (For 500 reps: SIM_N_REPS=500, SIM_CHUNK=100 -> 3*5 = 15 -> --array=1-15)
#   Rung order in run.R: 1=L1_param, 2=L3_adaptive, 3=L4_aggressive (chunks within).
#
# Walltime: L4 deep-RF is the heavy rung AND Design C uses a larger
# within-PSU n (base_n0 45 vs 25) -> bigger samples per rep. Allow 5 h.
# mem 96G (larger samples + SL fits per worker).
#
# Non-destructive: writes ONLY to sim_output/arc_runs/R12_highdeff and
# results/arc; the locked headline outputs are untouched.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R12_highdeff
#SBATCH --array=1-30
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=05:00:00
#SBATCH --output=sim_output/logs/R12_highdeff_%A_%a.out
#SBATCH --error=sim_output/logs/R12_highdeff_%A_%a.err
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
export SIM_OUT="$DATA_ROOT/arc_runs/R12_highdeff"

# --- run size (keep consistent with --array above) ---
export SIM_N_REPS=1000                        # 500 also fine (set --array=1-15)
export SIM_CHUNK=100

# --- Design C knobs (defaults baked into design_c.R; override here to RETUNE) ---
# If the smoke gate showed realized deff_clust < ~2.5, RAISE these and resubmit:
# export DC_SIGMA2_Y=4.0
# export DC_SIGMA2_A=2.0
# export DC_BASE_N0=50
# (default Design C: DC_SIGMA2_Y=3.5, DC_SIGMA2_A=1.7, DC_BASE_N0=45)

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$SIM_CODE/arc_runs/R12_highdeff/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
