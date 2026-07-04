#!/bin/bash
# =====================================================================
# submit.slurm  —  R15_aipw_benchmark (survey-weighted AIPW benchmark, spec A7)
# NON-DESTRUCTIVE enhancement run. Submit from the repo root / SLURM_SUBMIT_DIR:
#   sbatch codes/arc_runs/R15_aipw_benchmark/submit.slurm
#
# One array task = one (scenario x rung x rep-chunk) cell  (mirrors run_sim.R).
#   #tasks = 2 scenarios x 4 rungs x ceil(SIM_N_REPS / SIM_CHUNK)
#   FULL : SIM_N_REPS=1000, SIM_CHUNK=100 -> 2*4*10 = 80  -> --array=1-80
#   TEST : --array=1,61   (one cheap standard/L1 chunk + one heavy standard/L4
#          chunk; seff both, then submit the full array -- completed chunks
#          re-skip instantly via the checkpoint)
# Keep --array, SIM_N_REPS, SIM_CHUNK consistent (run.R stops on a bad index).
# WALLTIME RECOVERY (e.g. L4 chunks hit the 4 h limit): resubmit ONLY the
# failed task ids at the UNCHANGED SIM_CHUNK=100 with a longer --time:
#   sbatch --array=<failed ids> --time=08:00:00 codes/arc_runs/R15_aipw_benchmark/submit.slurm
# That is truly idempotent: checkpoint filenames are keyed on the chunk NUMBER,
# whose rep range is derived from SIM_CHUNK at runtime. DO NOT change SIM_CHUNK
# on a partially complete run -- old CHUNK=100 files alias the new chunk
# numbers (different rep ranges), so missing reps are never computed and
# surviving reps are double-counted. If you MUST change SIM_CHUNK (e.g. halve
# to 50, --array=1-160), FIRST delete/move ALL r15_<scenario>_<rung>_chunk*.rds
# (and their manifests) for the affected cells. See NOTES.md.
#
# Grid task order (expand.grid scenario fastest; chunks innermost):
#   tasks  1-10  standard x L1_param      11-20  R1 x L1_param
#   tasks 21-30  standard x L2_smooth     31-40  R1 x L2_smooth
#   tasks 41-50  standard x L3_adaptive   51-60  R1 x L3_adaptive
#   tasks 61-70  standard x L4_aggressive 71-80  R1 x L4_aggressive
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R15_aipw_benchmark
#SBATCH --array=1-80
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=sim_output/logs/R15_aipw_benchmark_%A_%a.out
#SBATCH --error=sim_output/logs/R15_aipw_benchmark_%A_%a.err
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
export R15_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R15_aipw_benchmark"
export R15_OUT="$DATA_ROOT/arc_runs/R15_aipw_benchmark"

# --- run size (edit with --array above) ---
export SIM_N_REPS=1000
export SIM_CHUNK=100

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R15_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"
