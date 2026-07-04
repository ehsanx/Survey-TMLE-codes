# =====================================================================
# nhanes/00_download_cache.R  —  reproducible raw NHANES download + cache
#
# Downloads (via nhanesA) every component table needed across the four examples
# (nhanes/examples.R) for cycles 2007-2008 .. 2017-2018, caches each raw table
# to nhanes/raw/<TABLE>.rds (download-if-missing, so re-runs are offline), and
# writes a provenance manifest. Raw INTEGER codes are kept (translated = FALSE):
# they are stable across cycles; decoding happens in a controlled step in 01.
#
# Key fix vs the legacy pipeline: BP table is BPXO_J (oscillometric) in 2017-18,
# not BPX_J -> handled by examples.R::tbl_name(). Tables absent in a cycle are
# recorded as "missing" and skipped, not fatal.
#
# Usage:  Rscript nhanes/00_download_cache.R        (full union; idempotent)
#         SIM_NH_TABLES="DEMO_E,BMX_E,BPX_E,BPXO_J" Rscript nhanes/00_download_cache.R   (subset, for validation)
# =====================================================================

suppressMessages({ library(nhanesA) })

REPO_ROOT <- Sys.getenv("REPO_ROOT", ".")
SRC <- file.path(REPO_ROOT, "nhanes", "examples.R")
if (!file.exists(SRC)) stop("cannot find examples.R at ", SRC)
source(SRC)

raw_dir <- file.path(REPO_ROOT, "nhanes", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

## ---- which tables to fetch -------------------------------------------------
tbl_override <- Sys.getenv("SIM_NH_TABLES", "")
if (nzchar(tbl_override)) {
  tables <- trimws(strsplit(tbl_override, ",")[[1]])
} else {
  tables <- unlist(lapply(COMPONENTS_ALL, function(comp)
    vapply(CYCLE_SUFFIXES, function(suf) tbl_name(comp, suf), character(1))))
  tables <- unique(tables)
}
cat("Will ensure", length(tables), "raw NHANES tables are cached under", raw_dir, "\n")

## ---- robust single-table downloader (version-tolerant) --------------------
fetch_raw <- function(tbl) {
  # try the modern signature first, then fall back for older nhanesA
  out <- tryCatch(nhanesA::nhanes(tbl, translated = FALSE, cleanse_numeric = TRUE),
                  error = function(e) tryCatch(nhanesA::nhanes(tbl, translated = FALSE),
                                               error = function(e2) NULL))
  out
}

nhanesA_ver <- as.character(utils::packageVersion("nhanesA"))
has_digest  <- requireNamespace("digest", quietly = TRUE)
manifest <- list()

for (tbl in tables) {
  f <- file.path(raw_dir, paste0(tbl, ".rds"))
  if (file.exists(f)) {
    x <- readRDS(f); status <- "cached"
  } else {
    x <- fetch_raw(tbl)
    if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) {
      status <- "missing"; x <- NULL
    } else {
      saveRDS(x, f); status <- "downloaded"
    }
  }
  manifest[[tbl]] <- data.frame(
    table = tbl, status = status,
    nrow = if (is.null(x)) NA_integer_ else nrow(x),
    ncol = if (is.null(x)) NA_integer_ else ncol(x),
    sha1 = if (!is.null(x) && has_digest) digest::digest(x, algo = "sha1") else NA_character_,
    nhanesA = nhanesA_ver,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE)
  cat(sprintf("  %-10s %-12s  n=%s\n", tbl, status,
              ifelse(is.null(x), "-", nrow(x))))
}

man <- do.call(rbind, manifest); rownames(man) <- NULL
write.csv(man, file.path(raw_dir, "manifest.csv"), row.names = FALSE)

cat("\nDONE. ",
    sum(man$status %in% c("downloaded", "cached")), " tables present, ",
    sum(man$status == "missing"), " missing (absent that cycle).\n",
    "Manifest: ", file.path(raw_dir, "manifest.csv"), "\n", sep = "")
