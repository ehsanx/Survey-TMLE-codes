# =====================================================================
# render_table1_tex.R  —  turn the saved svy-weighted Table 1 objects
# (nhanes/results/<id>/table1_descriptives.rds $df) into clean booktabs
# LaTeX for Web Appendix E.  Deterministic: numbers come straight from $df
# (survey-weighted % / mean(SD)); unweighted per-arm n from the imputed data.
# Output: outputs/tables/webE_table1s.tex  (four \begin{table} blocks).
# =====================================================================
options(stringsAsFactors = FALSE)

EXPOSED_LAB <- c(E1 = "Short sleep", E2 = "Food insecure",
                 E3 = "Ever-vaped", E4 = "GDM history")
UNEXP_LAB   <- c(E1 = "7--9 h sleep", E2 = "Food secure",
                 E3 = "Never vaped", E4 = "No GDM")
OUTCOME_LAB <- c(E1 = "Obesity (BMI $\\ge$ 30), \\%",
                 E2 = "Depression (PHQ-9 $\\ge$ 10), \\%",
                 E3 = "Hypertension, \\%", E4 = "Hypertension, \\%")
CAPTION <- c(
  E1 = "Survey-weighted characteristics of the E1 analytic domain (adults aged $\\ge$ 20, NHANES 2007--2018), by short-sleep exposure. Cells are survey-weighted percentages or means (standard deviations); SMD is the standardized mean difference between exposure groups.",
  E2 = "Survey-weighted characteristics of the E2 analytic domain (adults aged $\\ge$ 20, NHANES 2007--2018), by household food-insecurity exposure. Cells are survey-weighted percentages or means (SD); SMD is the standardized mean difference.",
  E3 = "Survey-weighted characteristics of the E3 analytic domain (adults aged $\\ge$ 20, NHANES 2015--2018), by ever-use of electronic cigarettes. Cells are survey-weighted percentages or means (SD); SMD is the standardized mean difference.",
  E4 = "Survey-weighted characteristics of the E4 analytic domain (ever-pregnant, non-diabetic women aged $\\ge$ 20, NHANES 2007--2018), by history of gestational diabetes. Cells are survey-weighted percentages or means (SD); SMD is the standardized mean difference.")

# readable label for a variable / level string ------------------------------
clean_label <- function(v, id) {
  v <- trimws(v)
  # continuous "x (mean (SD))"
  if (grepl("\\(mean \\(SD\\)\\)", v)) {
    nm <- trimws(sub("\\(mean \\(SD\\)\\)", "", v))
    map <- c(age = "Age (years)", pir = "Poverty--income ratio",
             phq = "PHQ-9 score", hhsize = "Household size",
             bmi = "Body-mass index (kg/m$^2$)", cigday = "Cigarettes per day",
             npreg = "No. of pregnancies", parity = "Parity",
             agefirst = "Age at first birth (years)")
    lab <- if (nm %in% names(map)) map[[nm]] else nm
    return(list(type = "cont", label = paste0(lab, ", mean (SD)")))
  }
  # binary "x = Level (%)"
  if (grepl(" = .*\\(%\\)$", v)) {
    stub <- trimws(sub("\\(%\\)$", "", v))           # "sex = Female"
    nm   <- trimws(sub(" =.*$", "", stub))
    lvl  <- trimws(sub("^.*= ", "", stub))
    bmap <- c(sex = "Female, \\%", married = "Married or partnered, \\%",
              pa = "Physically active, \\%", insured = "Health insurance, \\%",
              diab = "Diabetes, \\%", menop = "Premenopausal, \\%",
              famdm = "Family history of diabetes, \\%",
              cycle = "2017--2018 cycle, \\%")
    if (nm == "Yf") return(list(type = "bin", label = OUTCOME_LAB[[id]]))
    lab <- if (nm %in% names(bmap)) bmap[[nm]] else paste0(nm, " = ", lvl, ", \\%")
    return(list(type = "bin", label = lab))
  }
  # categorical header "x (%)"
  if (grepl("\\(%\\)$", v)) {
    nm <- trimws(sub("\\(%\\)$", "", v))
    hmap <- c(race = "Race/ethnicity, \\%", educ = "Education, \\%",
              cycle = "NHANES cycle, \\%", smk = "Smoking status, \\%")
    lab <- if (nm %in% names(hmap)) hmap[[nm]] else paste0(nm, ", \\%")
    return(list(type = "head", label = lab))
  }
  # n row
  if (v == "n") return(list(type = "n", label = "n"))
  # otherwise a category level
  lmap <- c(NHWhite = "\\quad Non-Hispanic White", NHBlack = "\\quad Non-Hispanic Black",
            MexAm = "\\quad Mexican American", OtherHisp = "\\quad Other Hispanic",
            Other = "\\quad Other",
            "<HS" = "\\quad $<$ High school", HS = "\\quad High school",
            "SomeColl+" = "\\quad Some college or more",
            Never = "\\quad Never", Former = "\\quad Former", Current = "\\quad Current",
            "5" = "\\quad 2007--2008", "6" = "\\quad 2009--2010", "7" = "\\quad 2011--2012",
            "8" = "\\quad 2013--2014", "9" = "\\quad 2015--2016", "10" = "\\quad 2017--2018")
  lab <- if (v %in% names(lmap)) lmap[[v]] else paste0("\\quad ", v)
  list(type = "lvl", label = lab)
}

# extract the trailing parenthesized number (the weighted %) from "count (pct)"
pct_of <- function(cell) {
  cell <- trimws(cell)
  if (cell == "" || is.na(cell)) return("")
  m <- regmatches(cell, regexpr("\\(([0-9.]+)\\)\\s*$", cell))
  if (length(m) == 0) return(cell)
  sub("[()]", "", gsub("[()]", "", m))
}

fmt_smd <- function(s) { s <- trimws(s); if (s == "" || is.na(s)) "" else s }

render_one <- function(id) {
  obj <- readRDS(file.path("nhanes","results",id,"table1_descriptives.rds"))
  d   <- obj$df
  imp <- readRDS(file.path("nhanes","analytic", paste0(id, "_imputed.rds")))
  dom <- imp[isTRUE2(imp$inpop) & !is.na(imp$A), , drop = FALSE]
  n_un <- sum(dom$A == 0); n_ex <- sum(dom$A == 1); n_to <- nrow(dom)

  out <- c()
  out <- c(out, "\\begin{table}[H]\\centering\\small")
  out <- c(out, sprintf("\\caption{%s}", CAPTION[[id]]))
  out <- c(out, sprintf("\\label{tab:webE_%s}", id))
  out <- c(out, "\\begin{tabular}{@{}l rrr c@{}}")
  out <- c(out, "\\toprule")
  out <- c(out, sprintf("Characteristic & Overall & %s & %s & SMD \\\\",
                        UNEXP_LAB[[id]], EXPOSED_LAB[[id]]))
  out <- c(out, "\\midrule")
  out <- c(out, sprintf("Unweighted $n$ & %s & %s & %s & \\\\",
                        format(n_to, big.mark=","), format(n_un, big.mark=","),
                        format(n_ex, big.mark=",")))
  for (i in seq_len(nrow(d))) {
    v <- d$variable[i]
    info <- clean_label(v, id)
    ov <- d$Overall[i]; un <- d$Unexposed[i]; ex <- d$Exposed[i]; smd <- fmt_smd(d$SMD[i])
    if (info$type == "n") next                              # replaced above
    if (info$type == "cont") {
      out <- c(out, sprintf("%s & %s & %s & %s & %s \\\\",
                            info$label, trimws(ov), trimws(un), trimws(ex), smd))
    } else if (info$type == "bin") {
      out <- c(out, sprintf("%s & %s & %s & %s & %s \\\\",
                            info$label, pct_of(ov), pct_of(un), pct_of(ex), smd))
    } else if (info$type == "head") {
      out <- c(out, sprintf("%s & & & & %s \\\\", info$label, smd))
    } else {  # level
      out <- c(out, sprintf("%s & %s & %s & %s & \\\\",
                            info$label, pct_of(ov), pct_of(un), pct_of(ex)))
    }
  }
  out <- c(out, "\\bottomrule")
  out <- c(out, "\\end{tabular}")
  out <- c(out, "\\end{table}")
  paste(out, collapse = "\n")
}

isTRUE2 <- function(x) !is.na(x) & x == TRUE

blocks <- vapply(c("E1","E2","E3","E4"), render_one, character(1))
tex <- paste(c("% Auto-generated by nhanes/render_table1_tex.R -- do not edit by hand.",
               "% Survey-weighted Table 1s for Web Appendix E (numbers from",
               "% nhanes/results/<id>/table1_descriptives.rds; unweighted n from analytic/<id>_imputed.rds).",
               blocks), collapse = "\n\n")
writeLines(tex, "outputs/tables/webE_table1s.tex")
cat(tex)
cat("\n\n--- wrote outputs/tables/webE_table1s.tex ---\n")
