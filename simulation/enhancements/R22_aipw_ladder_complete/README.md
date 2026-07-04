# R22_aipw_ladder_complete — RUN README (copy / run / download)

**What this run does (one line).** Fills the only fillable gaps in the unified
L1–L7 × {TMLE, AIPW} × {single-fit, cluster-CV, CF} appendix figure by computing
the AIPW cells nobody has run yet: **AIPW-{SF,CV,CF} @ L5** and **AIPW-CV @ L2, L3**.
No TMLE, no duplicate cells, no xgboost/HAL. **SMOKE-gated locally: PASS.**

> **Standing rule:** this is **>10 min/cell → ARC; the author submits.** Do NOT
> run the full array locally. Only `SMOKE=1` (8 reps, glm+earth) runs locally.

---

## 0. TL;DR (three commands on ARC, via the bundle)

```bash
# from the bundle root on ARC scratch (output/ must be empty for R22):
sbatch --array=41 --time=12:00:00 job.slurm R22     # PROBE one L5 task first (read its Elapsed)
bash submit_all.sh R22                               # then the full --array=1-60
bash aggregate_all.sh R22                            # after all tasks finish -> summary CSV
# download the whole output/ folder back to your laptop (see §4)
```

---

## 1. What R22 computes (and why so little)

| rung | library | AIPW arms computed | why |
|---|---|---|---|
| `L2_smooth` | glm + gam + earth | **AIPW-CV** only | R15 already has AIPW-SF/CF @ L2 |
| `L3_adaptive` | glm + earth + ranger.t1 | **AIPW-CV** only | R15 already has AIPW-SF/CF @ L3 |
| `L5_nondonsker` | glm + earth + ranger.deep | **AIPW-SF, AIPW-CV, AIPW-CF** | R20 ran TMLE only at L5 |

Everything else in the figure is already locked (sim_full L1–L4 TMLE; R15 L1–L4
AIPW-SF/CF; R20 L5 TMLE; R21 L6/L7 both). `L1/L4 internal-CV` is **N/A by
construction** (a lone learner has nothing to weight). So after R22 the figure is
gap-free except those two principled N/A cells. **No locked number changes.**

Trust in the new cells = **code identity**: the AIPW machinery is the SAME
(`R15/aipw_helpers.R::aipw_arms`, `R21/r21_helpers.R::aipw_cv_arm`), already
validated at L1–L4 (R15) and L6/L7 (R21); R22 only applies it to the L2/L3/L5
libraries on the SAME seeds.

---

## 2. Prerequisites (already satisfied in the bundle)

The bundle's `engine/codes/arc_runs/` must contain **R15_aipw_benchmark/** (for
`aipw_helpers.R`) and **R21_deployable_cvcf/** (for `r21_helpers.R`) — both are
already there. R22 sources them read-only; under `REPO_ROOT=engine` the default
paths resolve with no extra wiring. Packages needed on ARC: SuperLearner, tmle,
survey, surveyCV, earth, gam, ranger (all already used by R15/R20). **No xgboost,
no hal9001** — so none of R21's HAL install fragility.

---

## 3. EXACT steps

### Path A — via the ARC bundle `<ARC-bundle-dir>` (recommended)

R22 is **already wired into the bundle** (this session): the package is at
`engine/codes/arc_runs/R22_aipw_ladder_complete/`, and `job.slurm`,
`submit_all.sh`, `aggregate_all.sh` each have an `R22` case. So:

1. **Make sure `output/` is clean for R22.** R22 writes only under
   `output/sim_output/arc_runs/R22_aipw_ladder_complete/` (a brand-new dir), so
   leftover R21 results do NOT contaminate it. If you want a fully lean upload,
   empty `output/` first (standard practice) — but it is not required for R22.
2. **Copy the whole bundle folder to ARC scratch**, e.g.
   `/scratch/YOUR_ALLOCATION/<ARC-bundle-dir>/`.
   (sim-only run — you do NOT need `data/analytic/`.)
3. **Probe the expensive block first** (L5 = tasks 41–60):
   ```bash
   cd /scratch/YOUR_ALLOCATION/<ARC-bundle-dir>
   sbatch --array=41 --time=12:00:00 job.slurm R22
   sacct -j <jobid> --format=JobID,Elapsed,State    # read Elapsed when it ends
   ```
4. **Submit the full array** (60 tasks: 1–10 std-L2, 11–20 R1-L2, 21–30 std-L3,
   31–40 R1-L3, 41–50 std-L5, 51–60 R1-L5):
   ```bash
   bash submit_all.sh R22
   ```
   (R22 is OPT-IN — it is NOT in the default `submit_all.sh` list; you must name it.)
5. **Monitor:** `squeue -u $USER` and/or `bash progress.sh` (DONE/EXP chunks).
6. **Aggregate after all tasks finish:**
   ```bash
   bash aggregate_all.sh R22
   ```
   This writes `output/results_arc/arc/R22_aipw_ladder_complete_summary.csv` and
   prints the deliverable cells + locked-neighbour plausibility context.

### Path B — via a plain repo checkout (alternative, no bundle)

```bash
# from the repo root on ARC (engine = the repo itself):
sbatch --array=41 --time=12:00:00 codes/arc_runs/R22_aipw_ladder_complete/submit.slurm   # probe L5
sbatch codes/arc_runs/R22_aipw_ladder_complete/submit.slurm                              # full 1-60
Rscript codes/arc_runs/R22_aipw_ladder_complete/aggregate.R
```
(`submit.slurm` already exports `R22_DIR`, `R21_DIR`, `R15_DIR`, `R22_OUT`.)

---

## 4. What to download, and where it goes

After `aggregate_all.sh R22`, **download the whole `output/` folder** from ARC
back to `<ARC-bundle-dir>\output\`. The pieces that matter:

| ARC path | what | local destination |
|---|---|---|
| `output/results_arc/arc/R22_aipw_ladder_complete_summary.csv` | the aggregated summary | → `<ARC-bundle-dir>/output/results_arc/arc/` |
| `output/sim_output/arc_runs/R22_aipw_ladder_complete/*.rds` | the per-task raw RDS (+ `manifest/`) | → same relative path |
| `output/logs/R22_*` | SLURM logs (keep for provenance) | → same |

You do **not** need to copy anything else back for R22.

---

## 5. After the download (next session, in the working repo)

1. **Re-derive the summary locally** from the raw RDS (cheap; best provenance):
   ```powershell
   $env:R22_OUT="<ARC-bundle-dir>\output\sim_output\arc_runs\R22_aipw_ladder_complete"
   & "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" codes/arc_runs/R22_aipw_ladder_complete/aggregate.R
   ```
   Confirm `n_reps==1000` in every deliverable cell and that AIPW-CV @ L2/L3 sits
   near nominal (≈ the locked TMLE-CV there) and AIPW @ L5 sits between the L4
   collapse and the L6 deployable values (the plausibility context printout).
2. **Fill the appendix figure** — the production render script pulls
   AIPW-CV@L2/L3 + AIPW@L5 from `R22_..._summary.csv` (everything else from the
   locked runs); the figure auto-completes on recompile.
3. **Extend `verify_inline_numbers.R`** for any newly-quoted R22 number; recompile
   **both** PDFs; confirm audit green (0 undefined refs).
4. **Freeze** R22 raw RDS + summary into the internal `survey-tmle-reproduce`
   archive. **Do NOT push — author signs off.**

---

## 6. Cost / walltime / recovery

- L2/L3 are cheap (AIPW-CV only). **L5 (tasks 41–60) is the expensive block**
  (deep RF + a full AIPW estimator), but lighter than R21 (no xgboost/HAL, no
  TMLE 3-way). The 12 h cap is very conservative — probe task 41, then run the
  rest at ≈3× observed Elapsed.
- **Idempotent:** re-submitting the full array skips finished chunks
  (`arc_skip_if_done` keyed on chunk number). If a chunk times out, resubmit
  ONLY that task id at the SAME `SIM_CHUNK` with a longer `--time`. NEVER change
  `SIM_CHUNK` against existing chunk files.

## 7. Local validation already done

`SMOKE=1 Rscript codes/arc_runs/R22_aipw_ladder_complete/run.R` →
`[SMOKE-GATE] PASS` (standard, glm+earth, 8 reps; AIPW-SF/CV/CF all present,
finite b/se). Plumbing is validated; the deep-RF science path is the same code
already validated by R15/R20/R21.
