# sim_battery4.R -- Paper 3, revised simulation battery (Phase B).
#
# ADEMP summary (details in the paper's Web Appendix C).
#  Aims: (1) compare the thresholded ODS estimator analysed in Theorem 2
#        (ODS-thr) with the BIC-lasso implementation variant (ODS-bic);
#        (2) isolate the contribution of the selector by giving every arm
#        the SAME downstream estimator (K = 5 cross-fitted AIPW with
#        unpenalised refits, propensity trimming at 0.05);
#        (3) report both variance estimators used in the paper: the
#        empirical variance of the AIPW pseudo-observations ("emp", the one
#        that produced Table 1 of the previous version, sim_battery3.R) and
#        the stacked sandwich (B.1) of Lemma A4 ("stk");
#        (4) measure how often the OAL arm adds anything to the outcome arm
#        (S_oal \ S_out non-empty) -- the outcome-only baseline;
#        (5) trace the boundary sequence b_n = (log n)^kappa / sqrt(n),
#        kappa in {1.5, 1.1, 0.75}: kappa > 1 satisfies Definition 1's rate
#        sqrt(n) b_n / log n -> infinity, kappa = 0.75 violates it.
#  Data-generating mechanisms: Y = tau A + Z'beta + eps, A ~ Bern(expit(Z'alpha)),
#        tau = 1, alpha = (.8,.8,a3,a3,0,0,1,1,0,...), beta = (.6,.6,b3,b3,.6,.6,0,...):
#        Z1,Z2 strong confounders; Z3,Z4 "boundary" confounders (a3, b3 vary);
#        Z5,Z6 outcome-only predictors; Z7,Z8 instruments; the rest noise.
#        Designs: Gaussian independent (cells I, II, III, P1, P2, boundary
#        sequence B*), Hermite-cubic outcome misspecification (M), AR(1)
#        correlated Gaussian (C, rho = .5), bounded Rademacher covariates
#        (Bd, coefficients scaled so that e(Z) in [.057,.943]), heteroscedastic
#        errors (H, sd depends on Z1+Z2), instrument-heavy (IV6, six
#        instruments).
#  Estimand: ATE tau = 1.
#  Methods: see ARMS below.  Selection is done inside each training fold;
#        the pilot is OLS when p <= n_tr/3 and ridge (lambda = .1) otherwise;
#        kappa_n = log(n_tr)/sqrt(n_tr) and the OAL lambda_n = n_tr^{1/4}
#        (sum-log-likelihood scale) use the training-fold size n_tr.
#  Performance: bias, empirical SD, RMSE, mean SE (emp, stk), SE/SD ratio,
#        coverage of nominal 95 percent Wald intervals (emp, stk), mean
#        interval width, fit failures, Monte Carlo SEs for bias and coverage;
#        selection frequencies (exact oracle set, superset of the oracle
#        set, instruments excluded, boundary confounders retained, mean |S|,
#        P(S_oal \ S_out non-empty)).
#
# Usage:  Rscript sim_battery4.R [n_main] [n_big] [n_bdry] [cores] [cells]
#   defaults here (2-core container, <= 60 min): see the header of the output.
#   Full size on the author's machine: see RUN_ON_MAC.md.
# Outputs: output/battery4_table.csv (long: design, n, arm, measure, value, mcse)
#          output/battery4_selection.csv, output/battery4_summary.txt

suppressMessages({ library(glmnet); library(parallel) })
TAU <- 1; K <- 5; TRIM <- 0.05
HAVE_CTMLE <- requireNamespace("ctmle", quietly = TRUE)

args   <- commandArgs(trailingOnly = TRUE)
N_MAIN <- if (length(args) >= 1) as.integer(args[1]) else 500
N_BIG  <- if (length(args) >= 2) as.integer(args[2]) else 300
N_BDRY <- if (length(args) >= 3) as.integer(args[3]) else 300
CORES  <- if (length(args) >= 4) as.integer(args[4]) else max(1L, detectCores())
CELLS  <- if (length(args) >= 5) strsplit(args[5], ",")[[1]] else NULL

# ---------------------------------------------------------------- selectors
pilot <- function(Z, A, Y) {
  n <- nrow(Z); p <- ncol(Z)
  if (p <= n / 3) {
    b <- coef(lm.fit(cbind(1, A, Z), Y))[-(1:2)]; b[is.na(b)] <- 0
    list(b = b, type = "ols")
  } else {
    f <- suppressWarnings(glmnet(cbind(A, Z), Y, alpha = 0, lambda = 0.1,
                                 penalty.factor = c(0, rep(1, p))))
    list(b = as.vector(coef(f))[-(1:2)], type = "ridge")
  }
}
sel_thr <- function(b, n) which(abs(b) >= log(n) / sqrt(n))
sel_bic <- function(x, y, family, unpen = NULL) {
  if (is.null(unpen)) { X <- x; pf <- rep(1, ncol(x)); off <- 0 }
  else { X <- cbind(unpen, x); pf <- c(rep(0, ncol(unpen)), rep(1, ncol(x))); off <- ncol(unpen) }
  fit <- suppressWarnings(glmnet(X, y, family = family, nlambda = 40, penalty.factor = pf))
  bic <- deviance(fit) + fit$df * log(length(y))
  cf <- as.vector(coef(fit, s = fit$lambda[which.min(bic)]))[-1]
  which(cf[(off + 1):length(cf)] != 0)
}
sel_oal_fixed <- function(Z, A, b) {
  # adaptive logistic lasso, weights |b_j|^{-2}, lambda_sum = n^{1/4}
  # (glmnet mean-log-likelihood scale: n^{1/4}/n = n^{-3/4})
  n <- nrow(Z)
  Zs <- sweep(Z, 2, pmax(abs(b)^2, 1e-8), "*"); lam <- n^(-0.75)
  fit <- suppressWarnings(glmnet(Zs, A, family = "binomial", standardize = FALSE,
                                 lambda = c(2 * lam, lam)))
  which(as.vector(coef(fit, s = lam))[-1] != 0)
}
sel_oal_wamd <- function(Z, A, b, alpha = 1) {
  # Shortreed-Ertefaie wAMD tuning over their lambda grid (practitioner default)
  n <- nrow(Z); p <- ncol(Z)
  Zs <- sweep(Z, 2, pmax(abs(b)^2, 1e-8), "*")
  lg <- sort(n^c(-10, -5, -1, -.75, -.5, -.25, .25, .49) / n, decreasing = TRUE)
  fit <- suppressWarnings(glmnet(Zs, A, family = "binomial", standardize = FALSE,
                                 lambda = lg, alpha = alpha))
  crit <- rep(Inf, length(lg)); es <- vector("list", length(lg))
  sdz <- apply(Z, 2, sd)
  for (k in seq_along(lg)) {
    e <- tryCatch(as.vector(predict(fit, Zs, s = lg[k], type = "response")),
                  error = function(x) NULL)
    if (is.null(e)) next
    e <- pmin(pmax(e, .01), .99); es[[k]] <- e; ipw <- A / e + (1 - A) / (1 - e)
    w1 <- ipw * A; w0 <- ipw * (1 - A)
    smd <- abs(colSums(Z * w1) / sum(w1) - colSums(Z * w0) / sum(w0)) / sdz
    crit[k] <- sum(abs(b) * smd)
  }
  k <- which.min(crit)
  list(sel = which(as.vector(coef(fit, s = lg[k]))[-1] != 0), e = es[[k]])
}
hajek <- function(A, Y, e) {
  ipw <- A / e + (1 - A) / (1 - e); X <- cbind(1, A)
  br <- solve(crossprod(X * sqrt(ipw))); bh <- br %*% crossprod(X * ipw, Y)
  r <- Y - X %*% bh; V <- br %*% crossprod(X * (ipw * as.vector(r))) %*% br
  c(bh[2], sqrt(V[2, 2]))
}

# ---------------------------------------------- cross-fitted AIPW, both SEs
safe_solve <- function(M, v) {
  out <- tryCatch(solve(M, v), error = function(e) NULL)
  if (is.null(out)) out <- MASS::ginv(M) %*% v
  out
}
refit_fold <- function(Z, A, Y, tr, te, S, trim = TRIM) {
  # Refit on training rows tr with set S, evaluate on test rows te.
  # Returns psi (test), and the stacked-sandwich correction for training rows.
  xtr <- cbind(1, Z[tr, S, drop = FALSE]); xte <- cbind(1, Z[te, S, drop = FALSE])
  Atr <- A[tr]; Ytr <- Y[tr]; Ate <- A[te]; Yte <- Y[te]
  d <- ncol(xtr); ntr <- nrow(xtr); nte <- nrow(xte)
  ps <- suppressWarnings(glm.fit(xtr, Atr, family = binomial()))
  g <- ps$coefficients; g[is.na(g)] <- 0
  b1 <- qr.coef(qr(xtr[Atr == 1, , drop = FALSE]), Ytr[Atr == 1]); b1[is.na(b1)] <- 0
  b0 <- qr.coef(qr(xtr[Atr == 0, , drop = FALSE]), Ytr[Atr == 0]); b0[is.na(b0)] <- 0
  # test-fold pseudo-observations
  eraw <- as.vector(plogis(xte %*% g)); trimmed <- eraw < trim | eraw > 1 - trim
  e <- pmin(pmax(eraw, trim), 1 - trim)
  m1 <- as.vector(xte %*% b1); m0 <- as.vector(xte %*% b0)
  psi <- m1 - m0 + Ate * (Yte - m1) / e - (1 - Ate) * (Yte - m0) / (1 - e)
  # D_k = test-fold mean of d psi / d theta, theta = (gamma, b1, b0)
  dg <- -(Ate * (Yte - m1) * (1 - e) / e + (1 - Ate) * (Yte - m0) * e / (1 - e))
  dg[trimmed] <- 0
  Dg <- colMeans(xte * dg)
  D1 <- colMeans(xte * (1 - Ate / e))
  D0 <- colMeans(-xte * (1 - (1 - Ate) / (1 - e)))
  # J_k = training mean of -dU/dtheta' (block diagonal), U_i training scores
  etr <- as.vector(plogis(xtr %*% g))
  Jg <- crossprod(xtr * sqrt(etr * (1 - etr))) / ntr
  J1 <- crossprod(xtr[Atr == 1, , drop = FALSE]) / ntr
  J0 <- crossprod(xtr[Atr == 0, , drop = FALSE]) / ntr
  vg <- safe_solve(Jg, Dg); v1 <- safe_solve(J1, D1); v0 <- safe_solve(J0, D0)
  Ug <- xtr * (Atr - etr)
  U1 <- xtr * (Atr * (Ytr - as.vector(xtr %*% b1)))
  U0 <- xtr * ((1 - Atr) * (Ytr - as.vector(xtr %*% b0)))
  corr <- as.vector(Ug %*% vg + U1 %*% v1 + U0 %*% v0) * (nte / ntr)
  list(psi = psi, corr = corr, ntrim = sum(trimmed))
}
aipw_from_sets <- function(Z, A, Y, folds, setlist) {
  # setlist[[k]] : named list of index sets (one per arm) for fold k
  n <- nrow(Z); arms <- names(setlist[[1]])
  psi <- matrix(NA_real_, n, length(arms), dimnames = list(NULL, arms))
  corr <- matrix(0, n, length(arms), dimnames = list(NULL, arms))
  ntrim <- setNames(numeric(length(arms)), arms)
  for (k in 1:K) {
    tr <- folds != k; te <- !tr; cache <- list()
    for (a in arms) {
      S <- setlist[[k]][[a]]; key <- paste0("S", paste(S, collapse = ","))
      if (is.null(cache[[key]])) cache[[key]] <- refit_fold(Z, A, Y, tr, te, S)
      r <- cache[[key]]
      psi[te, a] <- r$psi; corr[tr, a] <- corr[tr, a] + r$corr; ntrim[a] <- ntrim[a] + r$ntrim
    }
  }
  est <- colMeans(psi)
  se_emp <- apply(psi, 2, sd) / sqrt(n)
  phi <- psi + corr
  se_stk <- sqrt(colMeans(sweep(phi, 2, est)^2) / n)
  list(est = est, se_emp = se_emp, se_stk = se_stk, ntrim = ntrim / n)
}

# ------------------------------------------------------------------- C-TMLE
ctmle_arm <- function(Z, A, Y) {
  # Reference implementation (Ju, Gruber, van der Laan; CRAN package ctmle), C-TMLE1 as in the
  # package vignette: the lambda sequence for the propensity lasso starts at the cv.glmnet
  # minimiser and continues towards less penalisation; Q is the BIC-lasso outcome fit (the same
  # initial estimator as the other arms). Y and Q are rescaled to [0, 1] for the logistic
  # fluctuation and the estimate is rescaled back. Skipped where the package is unavailable
  # (RUN_ON_MAC.md); the first error messages are written to output/ctmle_errors.txt.
  if (!HAVE_CTMLE) return(c(NA_real_, NA_real_))
  tryCatch({
    W <- as.data.frame(Z); names(W) <- paste0("W", seq_len(ncol(Z)))
    SQ <- sel_bic(Z, Y, "gaussian", unpen = cbind(A))
    X1 <- cbind(1, A, Z[, SQ, drop = FALSE])
    qc <- qr.coef(qr(X1), Y); qc[is.na(qc)] <- 0
    Q <- cbind(as.vector(cbind(1, 0, Z[, SQ, drop = FALSE]) %*% qc),
               as.vector(cbind(1, 1, Z[, SQ, drop = FALSE]) %*% qc))
    lo <- min(Y, Q); hi <- max(Y, Q); rg <- hi - lo
    Ys <- (Y - lo) / rg; Qs <- (Q - lo) / rg
    cvg <- glmnet::cv.glmnet(x = Z, y = A, family = "binomial", nlambda = 20)
    lambdas <- cvg$lambda[which(cvg$lambda == cvg$lambda.min):length(cvg$lambda)]
    fit <- ctmle::ctmleGlmnet(Y = Ys, A = A, W = W, Q = Qs, lambdas = lambdas, ctmletype = 1,
                              family = "gaussian", gbound = TRIM, V = 5)
    c(fit$est * rg, sqrt(fit$var.psi) * rg)
  }, error = function(e) {
    f <- "output/ctmle_errors.txt"
    if (!file.exists(f) || length(readLines(f)) < 20) cat(conditionMessage(e), "\n", file = f, append = TRUE)
    c(NA_real_, NA_real_)
  })
}

# --------------------------------------------------------------------- arms
ARMS_AIPW <- c("oracle", "full", "out_thr", "out_bic", "oal_fix", "oal_wamd",
               "oads_thr", "oads_bic", "trt_bic", "bch")
ARMS_OTHER <- c("oracle_ols", "ols_oadsthr", "oal_ipw", "goal_ipw", "ry2", "ctmle")
ARM_LABELS <- c(
  oracle = "oracle set + AIPW", full = "full set + AIPW",
  out_thr = "outcome-only thresholded + AIPW", out_bic = "outcome-only BIC-lasso + AIPW",
  oal_fix = "OAL arm alone (fixed lambda) + AIPW", oal_wamd = "naive OAL selection (wAMD) + AIPW",
  oads_thr = "ODS-thr (theory object) + AIPW", oads_bic = "ODS-bic (software variant) + AIPW",
  trt_bic = "treatment-only BIC-lasso + AIPW", bch = "BCH double selection + AIPW",
  oracle_ols = "oracle OLS (textbook SE)", ols_oadsthr = "naive OLS on ODS-thr set (textbook SE)",
  oal_ipw = "naive OAL + Hajek IPW (practitioner pipeline)",
  goal_ipw = "GOAL (elastic net) + Hajek IPW", ry2 = "RY-ODS v2 (exploratory randomisation variant)",
  ctmle = "C-TMLE (ctmle package)")

ry2_arm <- function(Z, A, Y, folds) {
  # Web Appendix D: selection on U = Y + W, estimation on V = Y - W, W ~ N(0, sigma_hat^2),
  # sigma_hat^2 = df-corrected residual variance of the full OLS of Y on (1, A, Z);
  # the ODS-bic union is selected once on (Z, A, U) and held fixed; the AIPW refits
  # use V and the same K folds (both variance estimators computed on V).
  n <- nrow(Z); p <- ncol(Z)
  ff <- if (p < n - 2) lm.fit(cbind(1, A, Z), Y) else lm.fit(cbind(1, A, Z[, seq_len(floor(n / 2)), drop = FALSE]), Y)
  sig <- sqrt(sum(ff$residuals^2) / max(ff$df.residual, 1))
  W <- rnorm(n, sd = sig); U <- Y + W; V <- Y - W
  pl <- pilot(Z, A, U)
  S <- sort(union(sel_bic(Z, U, "gaussian", unpen = cbind(A)), sel_oal_fixed(Z, A, pl$b)))
  f <- aipw_from_sets(Z, A, V, folds, lapply(1:K, function(k) list(ry2 = S)))
  c(f$est, f$se_emp, f$se_stk, length(S))
}

one_rep <- function(seed, dgp, arms_aipw, arms_other) {
  set.seed(seed)
  d <- dgp(); Z <- d$Z; A <- d$A; Y <- d$Y; n <- nrow(Z); p <- ncol(Z)
  folds <- sample(rep(1:K, length.out = n))
  setlist <- vector("list", K); selstat <- vector("list", K)
  for (k in 1:K) {
    tr <- folds != k; Ztr <- Z[tr, , drop = FALSE]; Atr <- A[tr]; Ytr <- Y[tr]; ntr <- sum(tr)
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
    selstat[[k]] <- t(sapply(arms_aipw, function(a) {
      S <- all_sets[[a]]
      c(size = length(S), exact = as.numeric(setequal(S, d$Sstar)),
        superset = as.numeric(all(d$Sstar %in% S)),
        iv_excluded = as.numeric(!any(d$iv %in% S)),
        bdry_retained = mean(d$bdry %in% S),
        noise_count = sum(!(S %in% c(d$Sstar, d$iv))),
        oal_extra = switch(a, oads_thr = as.numeric(length(setdiff(s_oal, s_thr)) > 0),
                           oads_bic = as.numeric(length(setdiff(s_oal, s_bic)) > 0), NA_real_),
        oal_extra_n = switch(a, oads_thr = length(setdiff(s_oal, s_thr)),
                             oads_bic = length(setdiff(s_oal, s_bic)), NA_real_))
    }))
  }
  fit <- aipw_from_sets(Z, A, Y, folds, setlist)
  sel <- Reduce("+", selstat) / K
  out <- data.frame(arm = arms_aipw, est = fit$est, se_emp = fit$se_emp, se_stk = fit$se_stk,
                    trim_frac = fit$ntrim, sel, row.names = NULL)
  # full-sample (non-cross-fitted) arms
  oth <- lapply(arms_other, function(a) {
    v <- tryCatch(switch(a,
      oracle_ols = d$oracle_ols(d),
      ols_oadsthr = {
        pl <- pilot(Z, A, Y); S <- sort(union(sel_thr(pl$b, n), sel_oal_fixed(Z, A, pl$b)))
        m <- if (length(S)) lm(Y ~ A + Z[, S, drop = FALSE]) else lm(Y ~ A)
        c(summary(m)$coefficients["A", 1:2], length(S))
      },
      oal_ipw = { pl <- pilot(Z, A, Y); w <- sel_oal_wamd(Z, A, pl$b); c(hajek(A, Y, w$e), length(w$sel)) },
      goal_ipw = { pl <- pilot(Z, A, Y); w <- sel_oal_wamd(Z, A, pl$b, alpha = .5); c(hajek(A, Y, w$e), length(w$sel)) },
      ry2 = { r <- ry2_arm(Z, A, Y, folds); c(r[1], r[2], r[4], r[3]) },
      ctmle = c(ctmle_arm(Z, A, Y), NA)), error = function(e) c(NA, NA, NA))
    data.frame(arm = a, est = v[1], se_emp = v[2], se_stk = if (a == "ry2") v[4] else NA_real_,
               trim_frac = NA_real_,
               size = v[3], exact = NA, superset = NA, iv_excluded = NA, bdry_retained = NA,
               noise_count = NA, oal_extra = NA, oal_extra_n = NA)
  })
  rbind(out, do.call(rbind, oth))
}

summarise_cell <- function(res, design, n, p, R) {
  # res: list of per-rep data.frames (NULL for failed reps)
  ok <- !vapply(res, is.null, TRUE)
  df <- do.call(rbind, res[ok]); arms <- unique(df$arm)
  rows <- list(); sel <- list()
  for (a in arms) {
    x <- df[df$arm == a, ]
    est <- x$est; fin <- is.finite(est) & is.finite(x$se_emp)
    nf <- sum(!fin) + sum(!ok); est <- est[fin]; se1 <- x$se_emp[fin]; se2 <- x$se_stk[fin]
    m <- length(est)
    if (m < 2) {
      ms <- c("bias", "emp_sd", "rmse", "mean_se_emp", "mean_se_stk", "se_ratio_emp", "se_ratio_stk",
              "cov_emp", "cov_stk", "width_emp", "width_stk", "trim_frac", "n_fail", "n_reps")
      rows[[a]] <- data.frame(design = design, n = n, p = p, arm = a, measure = ms,
                              value = ifelse(ms == "n_fail", nf, ifelse(ms == "n_reps", m, NA)), mcse = NA)
      next }
    bias <- mean(est) - TAU; sd_e <- sd(est)
    cov1 <- mean(abs(est - TAU) <= 1.96 * se1)
    cov2 <- if (all(is.na(se2))) NA else mean(abs(est - TAU) <= 1.96 * se2)
    meas <- c(bias = bias, emp_sd = sd_e, rmse = sqrt(mean((est - TAU)^2)),
              mean_se_emp = mean(se1), mean_se_stk = if (all(is.na(se2))) NA else mean(se2),
              se_ratio_emp = mean(se1) / sd_e, se_ratio_stk = if (all(is.na(se2))) NA else mean(se2) / sd_e,
              cov_emp = cov1, cov_stk = cov2,
              width_emp = mean(2 * 1.96 * se1), width_stk = if (all(is.na(se2))) NA else mean(2 * 1.96 * se2),
              cov_w1_stk = if (all(is.na(se2))) NA else mean(abs(est - TAU) <= (1.96 + 1) * se2),
              cov_w2_stk = if (all(is.na(se2))) NA else mean(abs(est - TAU) <= (1.96 + 2) * se2),
              width_w1_stk = if (all(is.na(se2))) NA else mean(2 * (1.96 + 1) * se2),
              width_w2_stk = if (all(is.na(se2))) NA else mean(2 * (1.96 + 2) * se2),
              trim_frac = mean(x$trim_frac[fin]), n_fail = nf, n_reps = m)
    mcse <- c(bias = sd_e / sqrt(m), emp_sd = sd_e / sqrt(2 * (m - 1)),
              rmse = NA, mean_se_emp = sd(se1) / sqrt(m), mean_se_stk = if (all(is.na(se2))) NA else sd(se2) / sqrt(m),
              se_ratio_emp = NA, se_ratio_stk = NA,
              cov_emp = sqrt(cov1 * (1 - cov1) / m), cov_stk = if (is.na(cov2)) NA else sqrt(cov2 * (1 - cov2) / m),
              width_emp = NA, width_stk = NA, cov_w1_stk = NA, cov_w2_stk = NA, width_w1_stk = NA,
              width_w2_stk = NA, trim_frac = NA, n_fail = NA, n_reps = NA)
    rows[[a]] <- data.frame(design = design, n = n, p = p, arm = a, measure = names(meas),
                            value = unname(meas), mcse = unname(mcse[names(meas)]))
    sel[[a]] <- data.frame(design = design, n = n, p = p, arm = a,
                           mean_size = mean(x$size, na.rm = TRUE),
                           p_exact_oracle = mean(x$exact), p_superset_oracle = mean(x$superset),
                           p_iv_excluded = mean(x$iv_excluded), p_bdry_retained = mean(x$bdry_retained),
                           mean_noise = mean(x$noise_count),
                           p_oal_extra = mean(x$oal_extra), mean_oal_extra = mean(x$oal_extra_n),
                           n_reps = m)
  }
  list(table = do.call(rbind, rows), selection = do.call(rbind, sel))
}

run_cell <- function(design, dgp, n, p, R, seed, arms_aipw = ARMS_AIPW, arms_other = ARMS_OTHER) {
  t0 <- Sys.time()
  seeds <- seed * 1000L + seq_len(R)
  res <- mclapply(seeds, function(s) tryCatch(one_rep(s, dgp, arms_aipw, arms_other),
                                              error = function(e) NULL), mc.cores = CORES)
  out <- summarise_cell(res, design, n, p, R)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("[%s] n=%d p=%d R=%d  %.1f min  (failed reps: %d)\n", design, n, p, R,
              el, sum(vapply(res, is.null, TRUE))))
  out$table$minutes <- el
  out
}

# --------------------------------------------------------------------- DGPs
lin_dgp <- function(n, p, a3, b3, sd_y = 2, design = "gauss", rho = .5,
                    alpha = NULL, beta = NULL, hetero = FALSE, nonlin = FALSE) {
  if (is.null(alpha)) alpha <- c(.8, .8, a3, a3, 0, 0, 1, 1, rep(0, p - 8))
  if (is.null(beta))  beta  <- c(.6, .6, b3, b3, .6, .6, 0, 0, rep(0, p - 8))
  Sstar <- which(beta != 0); iv <- which(beta == 0 & alpha != 0); bdry <- c(3, 4)
  chol_S <- if (design == "ar1") chol(rho^abs(outer(1:p, 1:p, "-"))) else NULL
  force(sd_y)
  function() {
    Z <- switch(design,
      gauss = matrix(rnorm(n * p), n, p),
      ar1   = matrix(rnorm(n * p), n, p) %*% chol_S,
      rademacher = matrix(sample(c(-1, 1), n * p, replace = TRUE), n, p))
    A <- rbinom(n, 1, plogis(Z %*% alpha))
    nl <- if (nonlin) 0.6 * (Z[, 1]^3 - 3 * Z[, 1] + Z[, 2]^3 - 3 * Z[, 2]) / sqrt(6) else 0
    sdv <- if (hetero) sd_y * exp(0.3 * (Z[, 1] + Z[, 2]) / sqrt(2) - 0.09) else sd_y
    Y <- TAU * A + as.vector(Z %*% beta) + nl + rnorm(n, sd = sdv)
    list(Z = Z, A = A, Y = Y, Sstar = Sstar, iv = iv, bdry = bdry, beta = beta,
         oracle_ols = if (nonlin) function(d) {
           H <- (d$Z[, 1]^3 - 3 * d$Z[, 1] + d$Z[, 2]^3 - 3 * d$Z[, 2]) / sqrt(6)
           c(summary(lm(d$Y ~ d$A + d$Z[, d$Sstar, drop = FALSE] + H))$coefficients[2, 1:2], length(d$Sstar))
         } else function(d)
           c(summary(lm(d$Y ~ d$A + d$Z[, d$Sstar, drop = FALSE]))$coefficients[2, 1:2], length(d$Sstar)))
  }
}
bn <- function(n, kappa) log(n)^kappa / sqrt(n)

CELL_DEFS <- list(
  list(id = "I",   design = "I (n=500,p=20)",          n = 500,  p = 20,  R = "main", seed = 801,
       dgp = function() lin_dgp(500, 20, .8, .6)),
  list(id = "II",  design = "II (n=200,p=50)",         n = 200,  p = 50,  R = "main", seed = 802,
       dgp = function() lin_dgp(200, 50, .7, .25)),
  list(id = "III", design = "III (n=2000,p=20, local boundary)", n = 2000, p = 20, R = "main", seed = 803,
       dgp = function() lin_dgp(2000, 20, 1, 3 / sqrt(2000))),
  list(id = "M",   design = "M (misspecified OM)",     n = 500,  p = 20,  R = "main", seed = 804,
       dgp = function() lin_dgp(500, 20, .8, .6, nonlin = TRUE)),
  list(id = "P1",  design = "P1 (n=500,p=100)",        n = 500,  p = 100, R = "main", seed = 805,
       dgp = function() lin_dgp(500, 100, .7, .25)),
  list(id = "P2",  design = "P2 (n=500,p=200)",        n = 500,  p = 200, R = "big",  seed = 806,
       dgp = function() lin_dgp(500, 200, .7, .25),
       arms_aipw = setdiff(ARMS_AIPW, "full")),   # a 201-parameter logistic refit on 400 obs. is degenerate
  list(id = "C",   design = "C (AR(1) rho=.5, n=500,p=20)", n = 500, p = 20, R = "main", seed = 807,
       dgp = function() lin_dgp(500, 20, .8, .6, design = "ar1")),
  list(id = "Bd",  design = "Bd (bounded Rademacher, n=500,p=20)", n = 500, p = 20, R = "main", seed = 808,
       dgp = function() lin_dgp(500, 20, NA, NA, design = "rademacher",
                                alpha = c(.4, .4, .4, .4, 0, 0, .6, .6, rep(0, 12)),
                                beta = c(.6, .6, .6, .6, .6, .6, 0, 0, rep(0, 12)))),
  list(id = "H",   design = "H (heteroscedastic, n=500,p=20)", n = 500, p = 20, R = "main", seed = 809,
       dgp = function() lin_dgp(500, 20, .8, .6, hetero = TRUE)),
  list(id = "IV6", design = "IV6 (six instruments, n=500,p=20)", n = 500, p = 20, R = "main", seed = 810,
       dgp = function() lin_dgp(500, 20, NA, NA,
                                alpha = c(.8, .8, .8, .8, 0, 0, rep(1, 6), rep(0, 8)),
                                beta = c(.6, .6, .6, .6, .6, .6, rep(0, 14)))))
BDRY_ARMS <- c("oracle", "full", "out_thr", "out_bic", "oal_fix", "oads_thr", "oads_bic", "bch")
for (nn in c(500, 2000, 8000)) for (kp in c(1.5, 1.1, 0.75)) {
  CELL_DEFS[[length(CELL_DEFS) + 1]] <- local({
    nn <- nn; kp <- kp
    list(id = sprintf("B%d_k%.2f", nn, kp),
         design = sprintf("B (boundary seq., kappa=%.2f, b_n=%.3f)", kp, bn(nn, kp)),
         n = nn, p = 20, R = "bdry", seed = 900 + round(kp * 100) + nn %/% 1000,
         dgp = function() lin_dgp(nn, 20, 1, bn(nn, kp)),
         arms_aipw = BDRY_ARMS, arms_other = c("oracle_ols", "ols_oadsthr"))
  })
}

# ---------------------------------------------------------------------- run
dir.create("output", showWarnings = FALSE)
if (!is.null(CELLS)) CELL_DEFS <- Filter(function(cd) cd$id %in% CELLS, CELL_DEFS)
t0 <- Sys.time(); tabs <- list(); sels <- list()
# Per-cell checkpoints: a finished cell is stored in output/battery4_cells/<id>_R<reps>.rds and
# reused on a rerun with the same replication count, so an interrupted run resumes where it stopped.
dir.create("output/battery4_cells", showWarnings = FALSE)
for (cd in CELL_DEFS) {
  R <- switch(cd$R, main = N_MAIN, big = N_BIG, bdry = N_BDRY)
  ck <- sprintf("output/battery4_cells/%s_R%d.rds", cd$id, R)
  if (file.exists(ck)) {
    out <- readRDS(ck); cat(sprintf("[%s] reused checkpoint %s\n", cd$design, ck))
  } else {
    out <- run_cell(cd$design, cd$dgp(), cd$n, cd$p, R, cd$seed,
                    arms_aipw = if (is.null(cd$arms_aipw)) ARMS_AIPW else cd$arms_aipw,
                    arms_other = if (is.null(cd$arms_other)) ARMS_OTHER else cd$arms_other)
    saveRDS(out, ck)
  }
  out$table$cell <- cd$id; out$selection$cell <- cd$id
  tabs[[cd$id]] <- out$table; sels[[cd$id]] <- out$selection
}
tab <- do.call(rbind, tabs); sel <- do.call(rbind, sels)
tab <- tab[, c("cell", "design", "n", "p", "arm", "measure", "value", "mcse", "minutes")]
tab$arm_label <- ARM_LABELS[tab$arm]
total_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
hdr <- c(sprintf("# sim_battery4.R  run on %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
         sprintf("# n_reps: main cells %d, growing-p cells %d, boundary-sequence cells %d; cores %d; total %.1f min",
                 N_MAIN, N_BIG, N_BDRY, CORES, total_min),
         sprintf("# ctmle package available: %s%s", HAVE_CTMLE,
                 if (HAVE_CTMLE) "" else " -- C-TMLE arm written as NA (not run in this environment; see RUN_ON_MAC.md)"),
         "# cell P2 omits the full-set arm (201-parameter logistic refit on 400 training observations is degenerate).",
         "# Boundary-sequence cells run the reduced arm set oracle/full/out_thr/out_bic/oal_fix/oads_thr/oads_bic/bch (+ oracle_ols, ols_oadsthr).",
         "# se_emp = empirical variance of AIPW pseudo-observations (the estimator used for Table 1 of the previous version, sim_battery3.R);",
         "# se_stk = stacked sandwich (B.1) of Lemma A4. Coverage of nominal 95 percent Wald intervals. mcse = Monte Carlo SE.")
writeLines(hdr, "output/battery4_table.csv")
suppressWarnings(write.table(tab, "output/battery4_table.csv", sep = ",", row.names = FALSE, append = TRUE))
writeLines(hdr, "output/battery4_selection.csv")
suppressWarnings(write.table(sel, "output/battery4_selection.csv", sep = ",", row.names = FALSE, append = TRUE))

# wide summary for reading
wide <- reshape(tab[, c("cell", "arm", "measure", "value")], idvar = c("cell", "arm"),
                timevar = "measure", direction = "wide")
names(wide) <- sub("^value\\.", "", names(wide))
sink("output/battery4_summary.txt")
cat(paste(hdr, collapse = "\n"), "\n\n")
for (cid in unique(wide$cell)) {
  cat("==== cell", cid, ":", unique(tab$design[tab$cell == cid]), "\n")
  w <- wide[wide$cell == cid, c("arm", "bias", "emp_sd", "rmse", "mean_se_emp", "mean_se_stk",
                                 "se_ratio_emp", "se_ratio_stk", "cov_emp", "cov_stk", "width_emp", "n_fail")]
  print(format(w, digits = 3), row.names = FALSE)
  s <- sel[sel$cell == cid, c("arm", "mean_size", "p_exact_oracle", "p_superset_oracle", "p_iv_excluded",
                              "p_bdry_retained", "mean_noise", "p_oal_extra", "mean_oal_extra")]
  print(format(s, digits = 3), row.names = FALSE); cat("\n")
}
sink()
cat(sprintf("battery4 done: %.1f min total\n", total_min))
