# ARC re-runs — round 2 (2026-06-08): m=40 MI + CV-arm fix

Two re-runs decided this session. **Both need the engine change below on ARC first.**

---

## ⚠ ENGINE CHANGED — sync `R/estimators.R` to ARC before either re-run

`R/estimators.R` was edited (the Fully-Aware-CV arm now fits its nuisances
**unweighted on the real-data domain path** — `inpop` supplied — identical to the
CF arm; the simulation path is unchanged). This fixes the deterministic
E2/L2 `Fully-Aware-CV = -0.065` sign-flip (weighted full-sample SuperLearner
separation).

The old `COPY_MANIFEST.txt` says "engine UNCHANGED — do not copy". **That is now
stale.** Get the new engine onto ARC by EITHER:
- `git pull` on ARC (if `R/estimators.R` is committed + pushed), **or**
- copying the single file `R/estimators.R` to
  `/scratch/YOUR_ALLOCATION/GitHub/survey-tmle2/R/estimators.R`.

Verify on ARC: `grep -n "w_cv <- if (is.null(inpop))" R/estimators.R` (should hit).

---

## Re-run 1 — R06 multiple imputation at m=40 (headline Table 2)

`m=40` follows White, Royston & Wood (2011): m ≥ % incomplete cases
(E1 18.2%, E2 8.6%, E3 13.3%, **E4 29.2%**). FMI ≤ 0.14, so it is conservative.
`submit.slurm` now sets `MI_M=40` and `--time=05:00:00`.

```bash
cd /scratch/YOUR_ALLOCATION/GitHub/survey-tmle2
rm -f nhanes/nhanes_output/arc_runs/R06/mi_E*.rds   # else the checkpoint reuses m=10
sbatch nhanes/arc_runs/R06_mi/submit.slurm
```
Outputs: `results/R06_mi_summary.csv` + `R06_mi_primary.csv` (overwritten),
`nhanes/nhanes_output/arc_runs/R06/mi_E*.rds`. rsync those back, then wire the
Fully-Aware-CF MI intervals into Table 2; single imputation -> appendix sensitivity.

## Re-run 2 — CV-arm fix on the NHANES ladder (L2 + L3 only)

The CV arm exists only at the multi-learner rungs L2_smooth and L3_adaptive.
In `nhanes_arc.R` the grid is `expand.grid(example=E1..E4, rung=L1..L4)` with
example varying fastest, so:
  L1 = tasks 1–4, **L2 = tasks 5–8, L3 = tasks 9–12**, L4 = tasks 13–16.

Re-run just the 8 CV cells (runner has no checkpoint -> it overwrites cleanly;
back up the 8 locked files first if you want a safety copy):

```bash
cd /scratch/YOUR_ALLOCATION/GitHub/survey-tmle2
# (optional backup) cp nhanes/nhanes_output/intermediate/nh_E?_L{2_smooth,3_adaptive}.rds /tmp/
sbatch --array=5-12 nhanes/nhanes.slurm
```
Then rsync back the 8 `nh_E?_L{2_smooth,3_adaptive}.rds`, re-aggregate
(`nhanes/aggregate_nhanes.R`), regenerate figures, and re-run the number audit
(`verify_numbers.R`). Only the CV column should move; CF/FA/PA/NA
shift at most by split-Monte-Carlo noise (b_split_sd ~ 0.001–0.006). The primary
results and the ladder story are unchanged.
