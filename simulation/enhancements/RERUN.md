# RERUN.md — re-running the ARC enhancement jobs (updated 2026-06-08)

Post-mortem of the **first** ARC submission, the fixes now baked into every run package,
and exactly how to copy the code, submit with checkpointing, and resume after a timeout.
Companion to `README.md` (which has the per-run *scientific* detail).

---

## 1. Why the first runs failed (diagnosed from the 760 SLURM logs)

| Run | squeue status | Root cause (actual log text) | Fix applied |
|---|---|---|---|
| **R04** | FAILED (32/32) | `object 'pop' not found` — the population was built on the master but never sent to the 32 PSOCK workers | `clusterExport(cl, c("pop","learners","rung","scenario"))` in `run.R` |
| **R06** | FAILED (4/4) | `there is no package called 'mice'` | **install `mice` on ARC** (see §2) — no code bug |
| **R01** | Mixed (1/60) | `no applicable method for 'predict' applied to ... "NULL"` — a degenerate learner fit returned NULL | per-rep `tryCatch` → drop that rep |
| **R07** | Mixed (3/80) | same NULL-learner crash | per-rep `tryCatch` → drop that rep |
| **R02** | Mixed (10/40) | deep-RF L4 at m=300, 100-rep chunks > 6 h walltime | chunk 100→50, `--array=1-80`, 8 h |
| **R09** | TIMEOUT | single job > 1 h — but R09 is meant to run **locally** (a few minutes) | run it **locally** (§5); 4 h ARC fallback |
| **R10** | Mixed (1/5) | the J=64 cell (N = 960k) > 2 h | walltime 2 h → 8 h |
| *(all)* | — | `sh: line 1: git: command not found` — **benign** (manifest git-SHA on a node with no git) | `arc_git_sha()` records NA silently |

---

## 2. One-time ARC setup

Install **`mice`** (R06's only blocker). On ARC:

```bash
module load gcc/9.4.0 r/4.4.0
R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4/ \
  Rscript -e 'install.packages("mice", repos="https://cloud.r-project.org")'
```

Everything else (SuperLearner, tmle, survey, surveyCV, earth, gam, glmnet, ranger) is already installed.

---

## 3. Copy the updated code to ARC (manual, no git)

**Only the run packages changed.** The canonical engine `R/*.R` (config, dgp, estimators,
diagnostics, learners) is **UNCHANGED** — do **not** delete or replace it.

On ARC, delete the two old run-package folders, then paste the local ones in their place:

```bash
# on ARC
cd /scratch/YOUR_ALLOCATION/GitHub/survey-tmle2
rm -rf simulation/enhancements nhanes/arc_runs
```

| Copy this local folder | → to this ARC path |
|---|---|
| `<repo-root>\codes\arc_runs\`  | `/scratch/YOUR_ALLOCATION/GitHub/survey-tmle2/simulation/enhancements/` |
| `<repo-root>\Nhanes\arc_runs\` | `/scratch/YOUR_ALLOCATION/GitHub/survey-tmle2/nhanes/arc_runs/` |

(A ready-to-drag snapshot of exactly these two folders is staged at `..\<ARC-bundle-dir>\`.)

---

## 4. What changed in every run package

- **Checkpointing** — `_checkpoint.R::arc_skip_if_done(fn, task)` at the top of each `run.R`:
  if a task's output `.rds` already exists and is a readable file, the task prints `SKIP checkpoint`
  and exits 0. **Re-submitting the whole `--array` is now idempotent** — finished chunks skip in
  seconds, only missing / failed / timed-out chunks actually run. A truncated file from a killed
  task is detected (unreadable RDS) and recomputed.
- **Robust git** — `arc_git_sha()` replaces the raw `system("git rev-parse HEAD")`; no more
  `git: command not found` on the compute nodes.
- **Per-rep error guard** (sim runs R01, R03, R04, R07, R11, R12) — a degenerate learner fit now
  drops that single replicate (the aggregators already use `na.rm` / count reps) instead of
  crashing the whole task.
- **Walltime / chunk** — R02 (`SIM_CHUNK=50`, `--array=1-80`, 8 h), R10 (8 h), R09 (4 h fallback).

---

## 5. Submit / resume

From the repo root on scratch (`= $SLURM_SUBMIT_DIR`):

```bash
# the runs that failed — checkpointing makes a full re-submit safe
sbatch simulation/enhancements/R04_nuisance_rate/submit.slurm
sbatch simulation/enhancements/R01_simple_control/submit.slurm
sbatch simulation/enhancements/R07_informative/submit.slurm
sbatch simulation/enhancements/R02_largem_sweep/submit.slurm
sbatch simulation/enhancements/R10_fpc/submit.slurm
sbatch nhanes/arc_runs/R06_mi/submit.slurm          # only after installing mice (§2)
```

- **R09 — run LOCALLY** (its own header says so; it is minutes on a laptop and timed out oddly on
  ARC). See `nhanes/arc_runs/R09_sensitivity/NOTES.md` for the env vars, then
  `Rscript nhanes/arc_runs/R09_sensitivity/run.R`. (A 4 h ARC fallback walltime is set if you
  prefer the cluster.)
- **Resume after a timeout** — just `sbatch` the same `submit.slurm` again. Finished chunks skip via
  the checkpoint; only the unfinished ones re-run.
- **Re-run only specific failed task ids** (optional) — `sbatch --array=<ids> simulation/enhancements/<RID>/submit.slurm`.

---

## 6. Aggregate (after each array finishes) — SET THE ENV VARS FIRST

`aggregate.R` run **by hand** (not via the slurm script) has no `REPO_ROOT`/`DATA_ROOT`
set, so it falls back to the dev-machine Windows paths and reports `no *.rds found in
<repo-root>/...`. Export the paths once, from the repo root on scratch:

```bash
cd /scratch/YOUR_ALLOCATION/GitHub/survey-tmle2
export REPO_ROOT="$PWD" DATA_ROOT="$PWD/sim_output" SIM_CODE="$PWD/codes"

# already-completed runs that were never aggregated:
Rscript simulation/enhancements/R03_isolation_2x2/aggregate.R
Rscript simulation/enhancements/R11_resampling_eff/aggregate.R
Rscript nhanes/arc_runs/R05_harmonized_floor/aggregate_R05.R
# each re-run, after it completes -> results/<RID>_summary.csv:
Rscript simulation/enhancements/<RID>/aggregate.R          # R01, R02, R04, R07, R10, R06
```

Then rsync `results/`, `sim_output/arc_runs/`, and `nhanes/nhanes_output/` back to your machine.

---

## 7. Walltimes (post-fix; unchanged unless noted)

| Run | array | mem | walltime | note |
|---|---|---|---|---|
| R01 | 1-60 | 64G | 1 h | light |
| **R02** | **1-80** | 64G | **8 h** | chunk 50; deep-RF L4 @ m=300 is the binding cell |
| R03 | 1-20 | 64G | 4 h | |
| R04 | 1-32 | 64G | 4 h | |
| R05 | 1-8 | 64G | 4 h | nhanes |
| R06 | 1-4 | 64G | 3 h | needs `mice` |
| R07 | 1-80 | 64G | (orig) | |
| **R09** | — | — | **local** | run on a laptop; 4 h ARC fallback |
| **R10** | 1-5 | 64G | **8 h** | J=64 cell, N=960k |
| R11 | 1-20 | 64G | (orig) | |
| R12 | 1-30 | 64G | (orig) | |
