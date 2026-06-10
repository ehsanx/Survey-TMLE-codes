#!/bin/bash
# =====================================================================
# submit.slurm  —  ARC run R11_resampling_eff (Web Appendix D)
#   (a) resampling (survey JKn jackknife) SE check vs design linearization
#   (b) IPW/svyglm efficiency benchmark vs cross-fitted TMLE
#
# GRID  = rung x scenario x rep-chunk, with the driver ordering rungs
#         (L1_param, L3_adaptive) so L1 cells come FIRST:
#   R11_N_REPS=500, R11_CHUNK=100 -> 5 chunks/cell
#   2 rungs x 2 scenarios x 5 chunks = 20 tasks -> --array=1-20
#     tasks  1-10  = L1_param  (standard chunks 1-5, then R1 chunks 1-5)
#     tasks 11-20  = L3_adaptive
#
# SMOKE-GATE: run the SMOKE locally first; if the L1 jackknife/linearization
# ratio is ~1 (see NOTES.md), submit all 20. If you want to gate on the
# cluster too, submit only --array=1-10 (L1), check, then --array=11-20 (L3).
#
# Resources mirror codes/sim.slurm exactly (32 cpus, 64G). Time bumped to
# 03:00:00 because each L3 rep does the engine arms PLUS a CF-EIF refit PLUS
# two survey JKn replicate-variance passes (delete-one-PSU; many replicates
# in the `standard` scenario) PLUS the svyglm baseline.
# =====================================================================
#SBATCH --account=YOUR_ALLOCATION
#SBATCH --job-name=R11_resampling_eff
#SBATCH --array=1-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=sim_output/logs/R11_resampling_eff_%A_%a.out
#SBATCH --error=sim_output/logs/R11_resampling_eff_%A_%a.err
#SBATCH --mail-user=YOUR_EMAIL@EXAMPLE.COM
#SBATCH --mail-type=ALL
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p sim_output/logs

module load gcc/9.4.0
module load r/4.4.0
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/

# --- paths: everything under the submit dir (scratch), outputs isolated ---
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export DATA_ROOT="$SLURM_SUBMIT_DIR/sim_output"
export SIM_CODE="$SLURM_SUBMIT_DIR/codes"
export R11_DIR="$SLURM_SUBMIT_DIR/codes/arc_runs/R11_resampling_eff"
export R11_OUT="$SLURM_SUBMIT_DIR/sim_output/arc_runs/R11_resampling_eff"

# --- run size (keep consistent with --array above) ---
export R11_N_REPS=500
export R11_CHUNK=100

echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID started $(date) on $(hostname) [$SLURM_CPUS_PER_TASK cores]"
Rscript "$R11_DIR/run.R"
echo "task $SLURM_ARRAY_TASK_ID finished $(date)"

# After ALL tasks complete, aggregate (run once, e.g. via a dependency job):
#   Rscript "$R11_DIR/aggregate.R"
