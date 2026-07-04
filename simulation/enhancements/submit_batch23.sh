#!/bin/bash
# =====================================================================
# submit_batch23.sh — submit the Phase-4 batch-2+3 enhancement runs (R13-R20)
# from the repo root on the cluster scratch (= SLURM_SUBMIT_DIR).
#
#   bash simulation/enhancements/submit_batch23.sh            # all 8 (+ R15 E1)
#   bash simulation/enhancements/submit_batch23.sh R17 R20    # any subset
#
# Order = cheap science first (no run depends on another). Each run's
# NOTES.md documents its smoke command, TEST-vs-FULL array, aggregate
# command, and decision gates.
#
# WALLTIME RECOVERY (all runs): resubmit ONLY failed task ids at the
# UNCHANGED SIM_CHUNK with a longer --time. Never change SIM_CHUNK
# against existing chunk files (checkpoint aliasing).
# =====================================================================
set -e
mkdir -p sim_output/logs

LIST=("$@"); [ ${#LIST[@]} -eq 0 ] && LIST=(R17 R19 R18 R14 R13 R16 R20 R15)

for WHAT in "${LIST[@]}"; do
  case "$WHAT" in
    R13) sbatch simulation/enhancements/R13_null_typeI/submit.slurm ;;
    R14) sbatch simulation/enhancements/R14_dr_factorial/submit.slurm ;;
    R15) sbatch simulation/enhancements/R15_aipw_benchmark/submit.slurm
         sbatch simulation/enhancements/R15_aipw_benchmark/submit_e1.slurm ;;
    R16) sbatch simulation/enhancements/R16_harmonized_sim/submit.slurm ;;
    R17) sbatch Nhanes/sensitivity/R17_e4_floor/submit.slurm ;;
    R18) sbatch simulation/enhancements/R18_cvu/submit.slurm ;;
    R19) sbatch simulation/enhancements/R19_rate_sweep/submit.slurm ;;
    R20) sbatch simulation/enhancements/R20_cv_vs_cf_nondonsker/submit.slurm ;;
    *) echo "unknown run: $WHAT (use R13..R20)"; exit 2 ;;
  esac
done
echo "submitted: ${LIST[*]} — watch with squeue -u \$USER"
