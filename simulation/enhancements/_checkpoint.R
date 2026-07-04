# =====================================================================
# _checkpoint.R  —  shared helpers for the ARC enhancement runs
#
# Sourced (read-only) by every <RID>/run.R. Two utilities:
#
#   arc_skip_if_done(fn, task)
#     Per-chunk checkpoint. Call it EARLY in run.R (right after the task's
#     scenario/rung/chunk -> output filename `fn` is known, BEFORE any heavy
#     compute). If `fn` already exists as a valid, non-trivial, READABLE .rds,
#     it prints a SKIP line and exits the task with status 0. Effect:
#     re-submitting the FULL --array is idempotent -- completed chunks skip in
#     seconds, only missing / failed / timed-out chunks actually run. A
#     truncated file from a killed task (unreadable RDS) is treated as NOT done
#     and recomputed.
#
#   arc_git_sha()
#     Records the git commit for the manifest WITHOUT the
#     `sh: line 1: git: command not found` stderr noise seen in every ARC log
#     (compute nodes have no git on PATH). Returns NA cleanly when git is absent.
# =====================================================================

arc_skip_if_done <- function(fn, task = NA, min_bytes = 200L) {
  if (file.exists(fn) && isTRUE(file.info(fn)$size > min_bytes)) {
    ok <- tryCatch({ invisible(readRDS(fn)); TRUE }, error = function(e) FALSE)
    if (ok) {
      cat(sprintf("[task %s] SKIP checkpoint: %s already complete\n",
                  as.character(task), basename(fn)))
      quit(save = "no", status = 0)
    }
    cat(sprintf("[task %s] checkpoint %s is corrupt/truncated -> recomputing\n",
                as.character(task), basename(fn)))
  }
  invisible(NULL)
}

arc_git_sha <- function() tryCatch({
  x <- suppressWarnings(system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE))
  if (length(x) >= 1L && nzchar(x[1])) x[1] else NA_character_
}, error = function(e) NA_character_)
