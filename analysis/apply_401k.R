# apply_401k.R -- 401(k) eligibility illustration on a documented public input (added 2026-09-16,
# review round 2, comment R2-M7).
#
# Earlier versions of the paper quoted a 401(k) analysis produced by apply_realdata.R from a
# historical file (pension.rda) that is not part of the replication archive; those numbers were
# withdrawn from the manuscript. This script reruns the CURRENT pipeline of apply_realdata2.R
# (K = 5 cross-fitted AIPW, unpenalised refits, clipping at 0.05, both variance estimators, the
# same arms) on the public `pension` data of the CRAN package hdm (n = 9,915; Chernozhukov and
# Hansen 2004; the data set used in the 401(k) example of Chernozhukov et al. 2018), with the
# same 15 candidate controls as before: age, inc, educ, fsize, marr, twoearn, db, pira, hown,
# their squares (age, inc, educ, fsize) and the products inc x age and marr x fsize, all
# standardised. Treatment D = e401 (eligibility), outcome Y = net_tfa (net financial assets).
#
# Nothing in the paper depends on this output; it is provided so that a reader can reproduce
# the 401(k) comparison from a public, versioned input. Requires the hdm package (CRAN; not
# installable in the author's container, run on the Mac):
#   Rscript -e 'install.packages("hdm", repos = "https://cloud.r-project.org")'
#   Rscript apply_401k.R          # from analysis/; about one minute
# Outputs: output/applications_401k.csv (unrounded), output/applications_401k.txt.
suppressMessages(library(glmnet))
if (!requireNamespace("hdm", quietly = TRUE)) stop("package 'hdm' is required: install.packages('hdm')")

# ---- pipeline functions of apply_realdata2.R (everything before its data section) ---------
src <- readLines("apply_realdata2.R")
cut <- grep("^# -+ data$", src)
stopifnot(length(cut) == 1)
eval(parse(text = src[seq_len(cut - 1)]))          # pilot(), sel_*(), refit_fold(), run_pipelines(), K, TRIM

# ---- data ---------------------------------------------------------------------------------
data("pension", package = "hdm"); d <- pension
covs <- c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown")
stopifnot(all(c(covs, "e401", "net_tfa") %in% names(d)))
X <- d[, covs]
X$age2 <- d$age^2; X$inc2 <- d$inc^2; X$educ2 <- d$educ^2; X$fsize2 <- d$fsize^2
X$inc_age <- d$inc * d$age; X$marr_fsize <- d$marr * d$fsize
Z <- scale(as.matrix(X)); attr(Z, "scaled:center") <- NULL; attr(Z, "scaled:scale") <- NULL
A <- as.integer(d$e401); Y <- as.numeric(d$net_tfa)
cat("401(k): n =", nrow(Z), " p =", ncol(Z), " (hdm::pension, package version", as.character(packageVersion("hdm")), ")\n")

# ---- analysis ------------------------------------------------------------------------------
t0 <- Sys.time()
r <- run_pipelines(Z, A, Y, "gaussian", seed = 20260903)
tb <- r$table; tb$outcome <- "401(k): net financial assets (dollars), D = eligibility"; tb$n <- nrow(Z); tb$p <- ncol(Z)
tb$lo_stk <- tb$estimate - 1.96 * tb$se_stk; tb$hi_stk <- tb$estimate + 1.96 * tb$se_stk
tb <- tb[, c("outcome", "n", "p", "arm", "estimate", "se_emp", "se_stk", "lo_stk", "hi_stk", "mean_size", "trim_frac")]
dir.create("output", showWarnings = FALSE)
write.csv(tb, "output/applications_401k.csv", row.names = FALSE)
sink("output/applications_401k.txt")
cat("apply_401k.R --", format(Sys.time(), "%Y-%m-%d %H:%M"), "; input hdm::pension version", as.character(packageVersion("hdm")),
    "; n =", nrow(Z), "; p =", ncol(Z), "; time", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("Arms as in apply_realdata2.R (K = 5 cross-fitted AIPW, clipping 0.05); se_emp = empirical variance, se_stk = stacked sandwich (B.1).\n\n")
print(format(tb, digits = 4), row.names = FALSE)
cat("\nselection frequency over folds (columns retained by each arm):\n"); print(round(r$selmat, 2))
sink()
cat(readLines("output/applications_401k.txt"), sep = "\n")
