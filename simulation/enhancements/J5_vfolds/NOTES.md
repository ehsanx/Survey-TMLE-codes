# J5_vfolds -- number-of-folds (V) sensitivity for the CF arm

Sweeps V in {2,5,10} at L2/L4. The primary CF estimator's bias/SE-calibration/coverage are stable in V; the V=5 rows reproduce the headline cells (Web App F, tab:webF_vfolds).

**Engine:** sources `R/` read-only (the 2026-07 engine: explicit `cvQinit=FALSE` on the
single-fit arms; EIF reconstruction uses the fit's own g-bounds).

**Smoke:** `SMOKE=1 SIM_CODE=R J_OUT=sim_output/enhancements/J5_vfolds Rscript simulation/enhancements/J5_vfolds/run.R`
(2 reps, one cell, ~2-5 min; expect `SMOKE OK`).

**Full:** `sbatch simulation/enhancements/J5_vfolds/submit.slurm`, then
`Rscript simulation/enhancements/J5_vfolds/aggregate.R` -> `results/J5_vfolds_summary.csv`.

**Render:** `Rscript simulation/renderers/render_vfolds_tex.R` -> `outputs/tables/`.
