# make_tables4.R -- render battery4 outputs as markdown tables for paper.qmd.
# Output: output/battery4_tables.md (all numbers are read from battery4_table.csv
# and battery4_selection.csv; nothing is typed by hand).
tab <- read.csv("output/battery4_table.csv", comment.char = "#", stringsAsFactors = FALSE)
sel <- read.csv("output/battery4_selection.csv", comment.char = "#", stringsAsFactors = FALSE)
w <- reshape(tab[, c("cell", "arm", "measure", "value")], idvar = c("cell", "arm"), timevar = "measure", direction = "wide")
names(w) <- sub("^value\\.", "", names(w))
m <- reshape(tab[, c("cell", "arm", "measure", "mcse")], idvar = c("cell", "arm"), timevar = "measure", direction = "wide")
names(m) <- sub("^mcse\\.", "", names(m))
w$bias_mcse <- m$bias[match(paste(w$cell, w$arm), paste(m$cell, m$arm))]
w$cov_stk_mcse <- m$cov_stk[match(paste(w$cell, w$arm), paste(m$cell, m$arm))]
w$cov_emp_mcse <- m$cov_emp[match(paste(w$cell, w$arm), paste(m$cell, m$arm))]
cells_main <- c("I", "II", "III", "M", "P1", "P2", "C", "Bd", "H", "IV6")
cell_lab <- c(I = "I ($n{=}500, p{=}20$)", II = "II ($n{=}200, p{=}50$)", III = "III (local boundary)",
              M = "M (misspecified OM)", P1 = "P1 ($n{=}500, p{=}100$)", P2 = "P2 ($n{=}500, p{=}200$)",
              C = "C (AR(1) $\\rho{=}.5$)", Bd = "Bd (bounded)", H = "H (heteroscedastic)", IV6 = "IV6 (six instruments)")
arm_lab <- c(oracle = "oracle set + AIPW", full = "full set + AIPW", out_thr = "outcome-only thr.", out_bic = "outcome-only BIC",
             oal_fix = "OAL arm alone", oal_wamd = "naive OAL sel. + AIPW", oads_thr = "ODS-thr", oads_bic = "ODS-bic",
             trt_bic = "treatment-only lasso", bch = "BCH-DS", oracle_ols = "oracle OLS", ols_oadsthr = "naive OLS on ODS-thr set",
             oal_ipw = "naive OAL-IPW", goal_ipw = "GOAL-IPW", ry2 = "RY-ODS v2", ctmle = "C-TMLE")
pc <- function(x, d = 1) ifelse(is.na(x), "–", formatC(100 * x, format = "f", digits = d))
f3 <- function(x, d = 3) ifelse(is.na(x), "–", formatC(x, format = "f", digits = d))
get <- function(cell, arm, meas) { v <- w[w$cell == cell & w$arm == arm, meas]; if (length(v)) v else NA }
sink("output/battery4_tables.md")

cat("## Main-text Table 1: coverage (stacked sandwich), percent\n\n")
arms1 <- c("oal_ipw", "goal_ipw", "out_thr", "oads_thr", "oads_bic", "bch", "full")
cat("| Cell |", paste(arm_lab[arms1], collapse = " | "), "|\n|---|", paste(rep("---", length(arms1)), collapse = "|"), "|\n", sep = "")
for (cl in cells_main) {
  vals <- sapply(arms1, function(a) { v <- get(cl, a, if (a %in% c("oal_ipw", "goal_ipw")) "cov_emp" else "cov_stk"); pc(v) })
  vals["oads_thr"] <- paste0("**", vals["oads_thr"], "**")
  cat("| ", cell_lab[cl], " | ", paste(vals, collapse = " | "), " |\n", sep = "")
}
cat("\nn_reps per cell: ", paste(sapply(cells_main, function(cl) paste0(cl, "=", get(cl, "oads_thr", "n_reps"))), collapse = ", "), "\n\n")

cat("## Width relative to BCH (mean stacked SE ratio) and to oracle\n\n")
cat("| Cell | ODS-thr / BCH | ODS-bic / BCH | out_thr / BCH | ODS-thr / oracle OLS | P(S_oal \\\\ S_out non-empty), thr | same, bic |\n|---|---|---|---|---|---|---|\n")
for (cl in cells_main) {
  r <- function(a, b) f3(get(cl, a, "mean_se_stk") / get(cl, b, "mean_se_stk"))
  s1 <- sel[sel$cell == cl & sel$arm == "oads_thr", "p_oal_extra"]; s2 <- sel[sel$cell == cl & sel$arm == "oads_bic", "p_oal_extra"]
  cat("| ", cell_lab[cl], " | ", r("oads_thr", "bch"), " | ", r("oads_bic", "bch"), " | ", r("out_thr", "bch"), " | ",
      f3(get(cl, "oads_thr", "mean_se_stk") / get(cl, "oracle_ols", "mean_se_emp")), " | ", pc(s1), " | ", pc(s2), " |\n", sep = "")
}

cat("\n## Web Appendix C: full tables per cell\n\n")
arms_all <- c("oracle_ols", "oracle", "full", "out_thr", "out_bic", "oal_fix", "oal_wamd", "oads_thr", "oads_bic", "trt_bic", "bch",
              "ols_oadsthr", "oal_ipw", "goal_ipw", "ry2", "ctmle")
for (cl in unique(w$cell)) {
  cat("### Cell ", cl, ": ", unique(tab$design[tab$cell == cl]), " (n_reps = ", get(cl, "oads_thr", "n_reps"), ")\n\n", sep = "")
  cat("| Arm | bias (MCSE) | SD | RMSE | mean SE emp | mean SE stk | SE/SD emp | SE/SD stk | cov emp | cov stk (MCSE) | width stk | fail |\n")
  cat("|---|---|---|---|---|---|---|---|---|---|---|---|\n")
  for (a in arms_all) {
    if (!any(w$cell == cl & w$arm == a)) next
    g <- function(mm) get(cl, a, mm)
    cat("| ", arm_lab[a], " | ", f3(g("bias")), " (", f3(g("bias_mcse")), ") | ", f3(g("emp_sd")), " | ", f3(g("rmse")), " | ",
        f3(g("mean_se_emp")), " | ", f3(g("mean_se_stk")), " | ", f3(g("se_ratio_emp"), 2), " | ", f3(g("se_ratio_stk"), 2), " | ",
        pc(g("cov_emp")), " | ", pc(g("cov_stk")), " (", pc(g("cov_stk_mcse")), ") | ", f3(g("width_stk"), 2), " | ",
        ifelse(is.na(g("n_fail")), "–", g("n_fail")), " |\n", sep = "")
  }
  cat("\nSelection:\n\n| Arm | mean size | P(exact oracle) | P(superset) | P(IV excluded) | P(boundary retained) | mean noise | P(S_oal \\\\ S_out non-empty) |\n|---|---|---|---|---|---|---|---|\n")
  ss <- sel[sel$cell == cl, ]
  for (a in arms_all) { r <- ss[ss$arm == a, ]; if (!nrow(r)) next
    cat("| ", arm_lab[a], " | ", f3(r$mean_size, 1), " | ", pc(r$p_exact_oracle), " | ", pc(r$p_superset_oracle), " | ", pc(r$p_iv_excluded), " | ",
        pc(r$p_bdry_retained), " | ", f3(r$mean_noise, 2), " | ", pc(r$p_oal_extra), " |\n", sep = "") }
  cat("\n")
}

cat("\n## Boundary sequence\n\n")
bd <- grep("^B[0-9]", unique(w$cell), value = TRUE)
cat("| n | kappa | b_n | b_n/kappa_n(tr) | P(bdry retained) out_thr | oads_thr | out_bic | bias out_thr | bias/SD out_thr | cov stk: oracle | out_thr | oads_thr | oads_bic | BCH |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (cl in bd) {
  n <- unique(tab$n[tab$cell == cl]); kp <- as.numeric(sub(".*_k", "", cl)); bn <- log(n)^kp / sqrt(n); ntr <- 0.8 * n; kn <- log(ntr) / sqrt(ntr)
  ss <- sel[sel$cell == cl, ]; pb <- function(a) pc(ss[ss$arm == a, "p_bdry_retained"])
  cat("| ", n, " | ", kp, " | ", f3(bn), " | ", f3(bn / kn, 2), " | ", pb("out_thr"), " | ", pb("oads_thr"), " | ", pb("out_bic"), " | ",
      f3(get(cl, "out_thr", "bias")), " | ", f3(get(cl, "out_thr", "bias") / get(cl, "out_thr", "emp_sd"), 2), " | ",
      pc(get(cl, "oracle", "cov_stk")), " | ", pc(get(cl, "out_thr", "cov_stk")), " | ", pc(get(cl, "oads_thr", "cov_stk")), " | ",
      pc(get(cl, "oads_bic", "cov_stk")), " | ", pc(get(cl, "bch", "cov_stk")), " |\n", sep = "")
}
cat("\n## Key numbers\n\n")
for (cl in cells_main) {
  cat(cl, ": oads_thr cov_stk ", pc(get(cl, "oads_thr", "cov_stk")), " cov_emp ", pc(get(cl, "oads_thr", "cov_emp")),
      "; out_thr cov_stk ", pc(get(cl, "out_thr", "cov_stk")), " SD ", f3(get(cl, "out_thr", "emp_sd")), " vs oads_thr SD ", f3(get(cl, "oads_thr", "emp_sd")),
      "; oads_bic cov_stk ", pc(get(cl, "oads_bic", "cov_stk")), "; bch ", pc(get(cl, "bch", "cov_stk")),
      "; oal_ipw ", pc(get(cl, "oal_ipw", "cov_emp")), "; goal ", pc(get(cl, "goal_ipw", "cov_emp")),
      "; ry2 stk ", pc(get(cl, "ry2", "cov_stk")), " width/oracleOLS ", f3(get(cl, "ry2", "mean_se_stk") / get(cl, "oracle_ols", "mean_se_emp"), 2),
      "; widened B1 oads_thr ", pc(get(cl, "oads_thr", "cov_w1_stk")), " width ratio ", f3(get(cl, "oads_thr", "width_w1_stk") / get(cl, "oracle_ols", "width_emp"), 2),
      "; se_ratio_emp/stk oads_thr ", f3(get(cl, "oads_thr", "se_ratio_emp"), 2), "/", f3(get(cl, "oads_thr", "se_ratio_stk"), 2),
      " full ", f3(get(cl, "full", "se_ratio_emp"), 2), "/", f3(get(cl, "full", "se_ratio_stk"), 2),
      " bch ", f3(get(cl, "bch", "se_ratio_emp"), 2), "/", f3(get(cl, "bch", "se_ratio_stk"), 2),
      "; trim_frac full ", f3(get(cl, "full", "trim_frac")), " oads_thr ", f3(get(cl, "oads_thr", "trim_frac")), "\n", sep = "")
}
sink()
cat(readLines("output/battery4_tables.md"), sep = "\n")
