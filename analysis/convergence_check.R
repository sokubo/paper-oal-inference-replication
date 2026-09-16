# convergence_check.R -- arm-specific numerical diagnostics for the cross-fitted AIPW arms
# of sim_battery4.R (added 2026-09-16, review round 2, comment R2-M4).
#
# The battery counts a replication as a failure only when an arm's estimate or standard
# error is nonfinite ("n_fail" in battery4_table.csv).  Warnings from the logistic refits
# are suppressed there, aliased (rank-deficient) coefficients are set to zero, and a
# singular sandwich block is solved with a generalized inverse.  This script re-creates,
# for a fixed subset of replications of every cell, exactly the data, folds, selected
# sets and refits of the battery (same per-replication seeds, same functions, sourced
# from sim_battery4.R), and records for every propensity and outcome refit of every
# AIPW arm:
#   ps_nonconv   the logistic refit stopped at glm.fit's iteration limit (converged = FALSE)
#   ps_boundary  fitted probabilities numerically at 0 or 1 (glm.fit boundary flag)
#   ps_rankdef   aliased propensity coefficients (set to zero by the battery)
#   out_rankdef  aliased outcome-refit coefficients (set to zero by the battery)
#   ginv         a sandwich block (J_gamma, J_1, J_0) was singular and MASS::ginv was used
#   ntrim        number of test-fold propensities clipped to [0.05, 0.95]
#   dpsi_maxit   max |psi - psi'| and max |corr - corr'| when the propensity refit is
#                repeated with glm.control(maxit = 200) (only meaningful for ps_nonconv fits)
#   d_est, d_se_*, cov_flip   per replication and arm, the change in the arm's estimate, standard errors and
#                95 percent (stacked) coverage indicator when every non-converged propensity refit is replaced
#                by its maxit = 200 refit (zero when no refit of the arm was non-converged)
#   est_check    per replication, the largest difference between the arm estimates rebuilt here
#                and those returned by the battery's aipw_from_sets() on the same folds and sets
#                (an identity check; should be 0)
# and summarises them per cell and arm (output/convergence_check.csv, .txt).  The diagnostic
# uses the first n_reps_main (n_reps_bdry) replications of each cell, so it is a fixed subset,
# not an estimate of the overall frequencies; its purpose is to make the numerical events
# that the failure count does not record visible, arm by arm.
#
# Usage:  Rscript convergence_check.R [n_reps_main] [n_reps_bdry] [cores] [cells]
#   (2026-09-16 container run: 25 10 2)

suppressMessages({ library(glmnet); library(parallel) })
cc_args <- commandArgs(trailingOnly = TRUE)
RD_MAIN  <- if (length(cc_args) >= 1) as.integer(cc_args[1]) else 25L
RD_BDRY  <- if (length(cc_args) >= 2) as.integer(cc_args[2]) else 10L
RD_CORES <- if (length(cc_args) >= 3) as.integer(cc_args[3]) else max(1L, detectCores())
RD_ONLY  <- if (length(cc_args) >= 4) strsplit(cc_args[4], ",")[[1]] else NULL

# ---- source the battery's definitions (everything before its run section) --------------
src <- readLines("sim_battery4.R")
cut <- grep("^# -+ run$", src)
stopifnot(length(cut) == 1)
commandArgs <- function(...) character(0)          # neutralise the battery's argument parsing
eval(parse(text = src[seq_len(cut - 1)]))          # defines pilot(), sel_*(), refit_fold(), CELL_DEFS, ...
rm(commandArgs)
R_MAIN <- RD_MAIN; R_BDRY <- RD_BDRY; CORES <- RD_CORES   # the battery's own N_MAIN/CORES are not used here
CELL_DEFS_ALL <- CELL_DEFS
if (!is.null(RD_ONLY)) CELL_DEFS_ALL <- Filter(function(cd) cd$id %in% RD_ONLY, CELL_DEFS_ALL)

# ---- instrumented copy of refit_fold (identical numerics; diagnostics added) ------------
refit_fold_diag <- function(Z, A, Y, tr, te, S, trim = TRIM) {
  xtr <- cbind(1, Z[tr, S, drop = FALSE]); xte <- cbind(1, Z[te, S, drop = FALSE])
  Atr <- A[tr]; Ytr <- Y[tr]; Ate <- A[te]; Yte <- Y[te]
  ntr <- nrow(xtr); nte <- nrow(xte)
  ps <- suppressWarnings(glm.fit(xtr, Atr, family = binomial()))
  g <- ps$coefficients; ps_rankdef <- sum(is.na(g)); g[is.na(g)] <- 0
  b1 <- qr.coef(qr(xtr[Atr == 1, , drop = FALSE]), Ytr[Atr == 1]); out_rankdef <- sum(is.na(b1)); b1[is.na(b1)] <- 0
  b0 <- qr.coef(qr(xtr[Atr == 0, , drop = FALSE]), Ytr[Atr == 0]); out_rankdef <- out_rankdef + sum(is.na(b0)); b0[is.na(b0)] <- 0
  eraw <- as.vector(plogis(xte %*% g)); trimmed <- eraw < trim | eraw > 1 - trim
  e <- pmin(pmax(eraw, trim), 1 - trim)
  m1 <- as.vector(xte %*% b1); m0 <- as.vector(xte %*% b0)
  psi <- m1 - m0 + Ate * (Yte - m1) / e - (1 - Ate) * (Yte - m0) / (1 - e)
  dg <- -(Ate * (Yte - m1) * (1 - e) / e + (1 - Ate) * (Yte - m0) * e / (1 - e)); dg[trimmed] <- 0
  Dg <- colMeans(xte * dg); D1 <- colMeans(xte * (1 - Ate / e)); D0 <- colMeans(-xte * (1 - (1 - Ate) / (1 - e)))
  etr <- as.vector(plogis(xtr %*% g))
  Jg <- crossprod(xtr * sqrt(etr * (1 - etr))) / ntr
  J1 <- crossprod(xtr[Atr == 1, , drop = FALSE]) / ntr
  J0 <- crossprod(xtr[Atr == 0, , drop = FALSE]) / ntr
  ginv_used <- 0L
  solve_diag <- function(M, v) { out <- tryCatch(solve(M, v), error = function(e) NULL)
    if (is.null(out)) { ginv_used <<- ginv_used + 1L; out <- MASS::ginv(M) %*% v }; out }
  vg <- solve_diag(Jg, Dg); v1 <- solve_diag(J1, D1); v0 <- solve_diag(J0, D0)
  Ug <- xtr * (Atr - etr); U1 <- xtr * (Atr * (Ytr - as.vector(xtr %*% b1))); U0 <- xtr * ((1 - Atr) * (Ytr - as.vector(xtr %*% b0)))
  corr <- as.vector(Ug %*% vg + U1 %*% v1 + U0 %*% v0) * (nte / ntr)
  # repeat the propensity refit with a larger iteration cap (only the propensity part changes)
  dpsi <- NA_real_; dcorr <- NA_real_; psi2 <- psi; corr2 <- corr
  if (!isTRUE(ps$converged)) {
    ps2 <- suppressWarnings(glm.fit(xtr, Atr, family = binomial(), control = glm.control(maxit = 200)))
    g2 <- ps2$coefficients; g2[is.na(g2)] <- 0
    eraw2 <- as.vector(plogis(xte %*% g2)); trimmed2 <- eraw2 < trim | eraw2 > 1 - trim
    e2 <- pmin(pmax(eraw2, trim), 1 - trim)
    psi2 <- m1 - m0 + Ate * (Yte - m1) / e2 - (1 - Ate) * (Yte - m0) / (1 - e2)
    dg2 <- -(Ate * (Yte - m1) * (1 - e2) / e2 + (1 - Ate) * (Yte - m0) * e2 / (1 - e2)); dg2[trimmed2] <- 0
    Dg2 <- colMeans(xte * dg2); D12 <- colMeans(xte * (1 - Ate / e2)); D02 <- colMeans(-xte * (1 - (1 - Ate) / (1 - e2)))
    etr2 <- as.vector(plogis(xtr %*% g2)); Jg2 <- crossprod(xtr * sqrt(etr2 * (1 - etr2))) / ntr
    vg2 <- safe_solve(Jg2, Dg2); v12 <- safe_solve(J1, D12); v02 <- safe_solve(J0, D02)
    Ug2 <- xtr * (Atr - etr2)
    corr2 <- as.vector(Ug2 %*% vg2 + U1 %*% v12 + U0 %*% v02) * (nte / ntr)
    dpsi <- max(abs(psi2 - psi)); dcorr <- max(abs(corr2 - corr))
  }
  list(psi = psi, corr = corr, psi2 = psi2, corr2 = corr2,
       diag = c(ps_nonconv = as.integer(!isTRUE(ps$converged)), ps_boundary = as.integer(isTRUE(ps$boundary)),
                ps_rankdef = ps_rankdef, out_rankdef = out_rankdef, ginv = ginv_used, ntrim = sum(trimmed),
                d = ncol(xtr), dpsi_maxit = dpsi, dcorr_maxit = dcorr))
}

one_rep_diag <- function(seed, dgp, arms_aipw) {
  set.seed(seed)
  d <- dgp(); Z <- d$Z; A <- d$A; Y <- d$Y; n <- nrow(Z); p <- ncol(Z)
  folds <- sample(rep(1:K, length.out = n))       # identical to one_rep(): data, then folds
  rows <- list(); setlist <- vector("list", K)
  psi_diag <- matrix(NA_real_, n, length(arms_aipw), dimnames = list(NULL, arms_aipw))
  psi2_diag <- psi_diag
  corr_diag <- matrix(0, n, length(arms_aipw), dimnames = list(NULL, arms_aipw)); corr2_diag <- corr_diag
  for (k in 1:K) {
    tr <- folds != k; te <- !tr; Ztr <- Z[tr, , drop = FALSE]; Atr <- A[tr]; Ytr <- Y[tr]; ntr <- sum(tr)
    pl <- pilot(Ztr, Atr, Ytr)
    s_thr <- sel_thr(pl$b, ntr)
    s_bic <- sel_bic(Ztr, Ytr, "gaussian", unpen = cbind(Atr))
    s_oal <- sel_oal_fixed(Ztr, Atr, pl$b)
    s_trt <- if (any(c("trt_bic", "bch") %in% arms_aipw)) sel_bic(Ztr, Atr, "binomial") else integer(0)
    s_wamd <- if ("oal_wamd" %in% arms_aipw) sel_oal_wamd(Ztr, Atr, pl$b)$sel else integer(0)
    all_sets <- list(oracle = d$Sstar, full = seq_len(p), out_thr = s_thr, out_bic = s_bic,
                     oal_fix = s_oal, oal_wamd = s_wamd, oads_thr = sort(union(s_thr, s_oal)),
                     oads_bic = sort(union(s_bic, s_oal)), trt_bic = s_trt,
                     bch = sort(union(s_bic, s_trt)))
    setlist[[k]] <- all_sets[arms_aipw]
    cache <- list()
    for (a in arms_aipw) {
      S <- all_sets[[a]]; key <- paste0("S", paste(S, collapse = ","))
      if (is.null(cache[[key]])) cache[[key]] <- refit_fold_diag(Z, A, Y, tr, te, S)
      psi_diag[te, a] <- cache[[key]]$psi; psi2_diag[te, a] <- cache[[key]]$psi2
      corr_diag[tr, a] <- corr_diag[tr, a] + cache[[key]]$corr; corr2_diag[tr, a] <- corr2_diag[tr, a] + cache[[key]]$corr2
      rows[[length(rows) + 1]] <- data.frame(seed = seed, fold = k, arm = a, pilot = pl$type,
                                             t(cache[[key]]$diag), row.names = NULL)
    }
  }
  # identity check against the battery's own estimator function on the same sets and folds
  fit <- aipw_from_sets(Z, A, Y, folds, setlist)
  out <- do.call(rbind, rows)
  out$est_check <- max(abs(colMeans(psi_diag) - fit$est))
  # per-arm effect of the larger iteration cap on the replication's estimate, SEs and coverage indicator
  est1 <- colMeans(psi_diag); est2 <- colMeans(psi2_diag)
  se_stk1 <- sqrt(colMeans(sweep(psi_diag + corr_diag, 2, est1)^2) / n)
  se_stk2 <- sqrt(colMeans(sweep(psi2_diag + corr2_diag, 2, est2)^2) / n)
  se_emp1 <- apply(psi_diag, 2, sd) / sqrt(n); se_emp2 <- apply(psi2_diag, 2, sd) / sqrt(n)
  cov1 <- abs(est1 - TAU) <= 1.96 * se_stk1; cov2 <- abs(est2 - TAU) <= 1.96 * se_stk2
  eff <- data.frame(arm = arms_aipw, d_est = est2 - est1, d_se_emp = se_emp2 - se_emp1, d_se_stk = se_stk2 - se_stk1,
                    cov_flip = as.integer(cov1 != cov2), row.names = NULL)
  out <- merge(out, eff, by = "arm", sort = FALSE)
  out
}

dir.create("output", showWarnings = FALSE)
t0 <- Sys.time(); all <- list()
for (cd in CELL_DEFS_ALL) {
  R <- switch(cd$R, main = R_MAIN, big = R_MAIN, bdry = R_BDRY)
  arms <- if (is.null(cd$arms_aipw)) ARMS_AIPW else cd$arms_aipw
  seeds <- cd$seed * 1000L + seq_len(R)
  dgp <- cd$dgp()
  res <- mclapply(seeds, function(s) tryCatch(one_rep_diag(s, dgp, arms), error = function(e) NULL), mc.cores = CORES)
  ok <- !vapply(res, is.null, TRUE)
  df <- do.call(rbind, res[ok]); df$cell <- cd$id; df$design <- cd$design; df$n <- cd$n; df$p <- cd$p
  all[[cd$id]] <- df
  cat(sprintf("[%s] %d replications (%d failed to run)  %.1f min\n", cd$id, sum(ok), sum(!ok),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
fits <- do.call(rbind, all)
write.csv(fits, "output/convergence_check_fits.csv", row.names = FALSE)

# ---- per-cell, per-arm summary ------------------------------------------------------------
summ <- do.call(rbind, lapply(split(fits, list(fits$cell, fits$arm), drop = TRUE), function(x) {
  data.frame(cell = x$cell[1], design = x$design[1], n = x$n[1], p = x$p[1], arm = x$arm[1],
             arm_label = unname(ARM_LABELS[x$arm[1]]), pilot = paste(unique(x$pilot), collapse = "/"),
             reps = length(unique(x$seed)), fits = nrow(x), mean_dim = mean(x$d),
             ps_nonconv = sum(x$ps_nonconv), ps_boundary = sum(x$ps_boundary),
             ps_rankdef_fits = sum(x$ps_rankdef > 0), out_rankdef_fits = sum(x$out_rankdef > 0),
             ginv_fits = sum(x$ginv > 0), trim_frac = sum(x$ntrim) / (length(unique(x$seed)) * x$n[1]),
             max_dpsi_maxit = if (any(!is.na(x$dpsi_maxit))) max(x$dpsi_maxit, na.rm = TRUE) else NA_real_,
             max_dcorr_maxit = if (any(!is.na(x$dcorr_maxit))) max(x$dcorr_maxit, na.rm = TRUE) else NA_real_,
             est_check = max(x$est_check),
             max_d_est_maxit = max(abs(x$d_est)), max_d_se_stk_maxit = max(abs(x$d_se_stk)),
             cov_flips_maxit = sum(x$cov_flip[!duplicated(x$seed)]),
             row.names = NULL)
}))
ord <- order(match(summ$cell, vapply(CELL_DEFS_ALL, `[[`, "", "id")), match(summ$arm, ARMS_AIPW))
summ <- summ[ord, ]
hdr <- c(sprintf("# convergence_check.R  run on %s; replications per cell: main %d, boundary %d; cores %d; %.1f min",
                 format(Sys.time(), "%Y-%m-%d %H:%M"), R_MAIN, R_BDRY, CORES, as.numeric(difftime(Sys.time(), t0, units = "mins"))),
         "# fits = propensity/outcome refits examined (folds x replications); ps_nonconv = logistic refit at glm.fit's iteration limit;",
         "# ps_boundary = fitted probabilities numerically at 0/1; *_rankdef_fits = refits with aliased coefficients (set to zero);",
         "# ginv_fits = refits in which a sandwich block was solved with a generalized inverse; trim_frac = share of test-fold",
         "# propensities clipped; max_dpsi_maxit / max_dcorr_maxit = largest change in a pseudo-observation / correction term when",
         "# a non-converged propensity refit is repeated with maxit = 200 (NA when no refit was non-converged);",
         "# max_d_est_maxit / max_d_se_stk_maxit / cov_flips_maxit = largest change in a replication's arm estimate / stacked SE,",
         "# and number of replications whose 95 percent stacked-interval coverage indicator changes, under the maxit = 200 refits.")
writeLines(hdr, "output/convergence_check.csv")
suppressWarnings(write.table(summ, "output/convergence_check.csv", sep = ",", row.names = FALSE, append = TRUE))
sink("output/convergence_check.txt")
cat(paste(hdr, collapse = "\n"), "\n\n")
print(format(summ[, c("cell", "arm", "reps", "fits", "mean_dim", "ps_nonconv", "ps_boundary", "ps_rankdef_fits",
                      "out_rankdef_fits", "ginv_fits", "trim_frac", "max_dpsi_maxit", "max_d_est_maxit", "max_d_se_stk_maxit",
                      "cov_flips_maxit", "est_check")], digits = 3), row.names = FALSE)
cat("\nest_check = largest |estimate from the instrumented refits - estimate from aipw_from_sets()| over the replications (identity check).\n")
cat("\nTotals over all cells and arms:\n")
print(colSums(summ[, c("fits", "ps_nonconv", "ps_boundary", "ps_rankdef_fits", "out_rankdef_fits", "ginv_fits")]))
sink()
cat(readLines("output/convergence_check.txt"), sep = "\n")
