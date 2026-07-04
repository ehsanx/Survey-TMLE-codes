# J3_calib -- calibrated (raked) vs design weights

Analyzes each Design-A replicate twice -- true design (HT) weights vs those weights raked to two known population margins (stratum counts x L1 tertiles) -- across the 5-arm suite at L1/L2. Shows the primary CF estimator's coverage is unchanged by raking (Web App F, tab:webF_calib).

**Engine:** sources `R/` read-only (the 2026-07 engine: explicit `cvQinit=FALSE` on the
single-fit arms; EIF reconstruction uses the fit's own g-bounds).

**Smoke:** `SMOKE=1 SIM_CODE=R J_OUT=sim_output/enhancements/J3_calib Rscript simulation/enhancements/J3_calib/run.R`
(2 reps, one cell, ~2-5 min; expect `SMOKE OK`).

**Full:** `sbatch simulation/enhancements/J3_calib/submit.slurm`, then
`Rscript simulation/enhancements/J3_calib/aggregate.R` -> `results/J3_calib_summary.csv`.

**Render:** `Rscript simulation/renderers/render_calib_tex.R` -> `outputs/tables/`.
