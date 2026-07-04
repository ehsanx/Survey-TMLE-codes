# J4_overlap -- known-truth limited-overlap (positivity) stress

Raises the confounding odds ratio e^alpha_g from 1.3 to 2.5/5.0 with treatment prevalence held fixed (truth Psi invariant by construction), at L2 (3 levels) and L4 (2 levels). Shows the headline contrast reproduces at base overlap and that BOTH arms degrade under genuine near-positivity violation -- the floor acts as real truncation; cross-fitting dominates but does not rescue absent overlap (Web App F, tab:webF_overlapsim).

**Engine:** sources `R/` read-only (the 2026-07 engine: explicit `cvQinit=FALSE` on the
single-fit arms; EIF reconstruction uses the fit's own g-bounds).

**Smoke:** `SMOKE=1 SIM_CODE=R J_OUT=sim_output/enhancements/J4_overlap Rscript simulation/enhancements/J4_overlap/run.R`
(2 reps, one cell, ~2-5 min; expect `SMOKE OK`).

**Full:** `sbatch simulation/enhancements/J4_overlap/submit.slurm`, then
`Rscript simulation/enhancements/J4_overlap/aggregate.R` -> `results/J4_overlap_summary.csv`.

**Render:** `Rscript simulation/renderers/render_overlapsim_tex.R` -> `outputs/tables/`.
