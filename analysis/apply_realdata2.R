# apply_realdata2.R -- Paper 3, revised RHC application (Phase B).
#
# Data: SUPPORT right-heart-catheterisation study (Connors et al. 1996),
#   rhc_full.rda, n = 5,735 ICU admissions; treatment = RHC within the first
#   24 hours of ICU admission (swang1); baseline covariates measured on day 1.
# Outcomes: (a) 30-day survival, Y = 1{dth30 == "No"} (risk difference);
#           (b) hospital length of stay in days, dschdte - sadmdte (n = 5,734
#           with a recorded discharge date; descriptive estimand that ignores
#           the competing event death).
# Candidate pool: the SAME dummy-coded baseline pool for both outcomes (see
#   output/applications2_pool.txt for the variable list with SUPPORT labels).
# Arms (all with the same downstream estimator: K = 5 cross-fitted AIPW,
#   unpenalised refits, propensity trimming at 0.05, both variance estimators):
#   full, out_thr, out_bic, oal_fix, oal_wamd, oads_thr, oads_bic, trt_bic, bch;
#   plus the practitioner pipeline naive OAL (wAMD) + Hajek IPW.
#   For the binary outcome the pilot, the outcome refits and the stacked
#   sandwich are logistic and the threshold kappa_n is applied on the
#   log-odds scale; for LOS everything is linear.  kappa_n and lambda_n use
#   the training-fold size.
# Outputs: output/applications2_table.csv (unrounded), applications2_table.txt
#   (4 significant digits), applications2_selection.csv (per-variable selection
#   frequencies across folds), applications2_pool.txt, applications2_splits.csv
#   (stability over 5 fold splits), applications2_overlap.txt.
suppressMessages(library(glmnet))
K <- 5; TRIM <- 0.05; N_SPLITS <- 5
DATA_PATH <- local({                       # first existing copy of the public RHC file
  cand <- c("data/rhc_full.rda", "../data/rhc_full.rda",
            "../../../9br-update/data/rhc_full.rda", "../../../9br/data/rhc_full.rda")
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("rhc_full.rda not found; looked in: ", paste(cand, collapse = ", "))
  hit[1]
})
dir.create("output", showWarnings = FALSE)

# ------------------------------------------------------------ selectors
pilot <- function(Z, A, Y, family) {
  if (family == "gaussian") { b <- coef(lm.fit(cbind(1, A, Z), Y))[-(1:2)] }
  else b <- suppressWarnings(coef(glm.fit(cbind(1, A, Z), Y, family = binomial()))[-(1:2)])
  b[is.na(b)] <- 0; b
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
  n <- nrow(Z); Zs <- sweep(Z, 2, pmax(abs(b)^2, 1e-8), "*"); lam <- n^(-0.75)
  fit <- suppressWarnings(glmnet(Zs, A, family = "binomial", standardize = FALSE, lambda = c(2 * lam, lam)))
  which(as.vector(coef(fit, s = lam))[-1] != 0)
}
sel_oal_wamd <- function(Z, A, b) {
  n <- nrow(Z); Zs <- sweep(Z, 2, pmax(abs(b)^2, 1e-8), "*")
  lg <- sort(n^c(-10, -5, -1, -.75, -.5, -.25, .25, .49) / n, decreasing = TRUE)
  fit <- suppressWarnings(glmnet(Zs, A, family = "binomial", standardize = FALSE, lambda = lg))
  crit <- rep(Inf, length(lg)); es <- vector("list", length(lg)); sdz <- apply(Z, 2, sd) + 1e-12
  for (k in seq_along(lg)) {
    e <- tryCatch(as.vector(predict(fit, Zs, s = lg[k], type = "response")), error = function(x) NULL)
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

# -------------------------------------------- cross-fitted AIPW, both SEs
safe_solve <- function(M, v) { out <- tryCatch(solve(M, v), error = function(e) NULL)
  if (is.null(out)) out <- MASS::ginv(M) %*% v; out }
refit_fold <- function(Z, A, Y, tr, te, S, family, trim = TRIM) {
  xtr <- cbind(1, Z[tr, S, drop = FALSE]); xte <- cbind(1, Z[te, S, drop = FALSE])
  Atr <- A[tr]; Ytr <- Y[tr]; Ate <- A[te]; Yte <- Y[te]; ntr <- nrow(xtr); nte <- nrow(xte)
  ps <- suppressWarnings(glm.fit(xtr, Atr, family = binomial())); g <- ps$coefficients; g[is.na(g)] <- 0
  fit_out <- function(rows) {
    if (family == "gaussian") { b <- qr.coef(qr(xtr[rows, , drop = FALSE]), Ytr[rows]) }
    else b <- suppressWarnings(glm.fit(xtr[rows, , drop = FALSE], Ytr[rows], family = binomial()))$coefficients
    b[is.na(b)] <- 0; b
  }
  b1 <- fit_out(Atr == 1); b0 <- fit_out(Atr == 0)
  link <- if (family == "gaussian") function(u) u else plogis
  dlink <- if (family == "gaussian") function(u) rep(1, length(u)) else function(u) plogis(u) * (1 - plogis(u))
  eraw <- as.vector(plogis(xte %*% g)); trimmed <- eraw < trim | eraw > 1 - trim
  e <- pmin(pmax(eraw, trim), 1 - trim)
  l1 <- as.vector(xte %*% b1); l0 <- as.vector(xte %*% b0); m1 <- link(l1); m0 <- link(l0)
  psi <- m1 - m0 + Ate * (Yte - m1) / e - (1 - Ate) * (Yte - m0) / (1 - e)
  dg <- -(Ate * (Yte - m1) * (1 - e) / e + (1 - Ate) * (Yte - m0) * e / (1 - e)); dg[trimmed] <- 0
  Dg <- colMeans(xte * dg)
  D1 <- colMeans(xte * (dlink(l1) * (1 - Ate / e)))
  D0 <- colMeans(-xte * (dlink(l0) * (1 - (1 - Ate) / (1 - e))))
  etr <- as.vector(plogis(xtr %*% g)); lt1 <- as.vector(xtr %*% b1); lt0 <- as.vector(xtr %*% b0)
  Jg <- crossprod(xtr * sqrt(etr * (1 - etr))) / ntr
  J1 <- crossprod(xtr[Atr == 1, , drop = FALSE] * sqrt(dlink(lt1)[Atr == 1])) / ntr
  J0 <- crossprod(xtr[Atr == 0, , drop = FALSE] * sqrt(dlink(lt0)[Atr == 0])) / ntr
  vg <- safe_solve(Jg, Dg); v1 <- safe_solve(J1, D1); v0 <- safe_solve(J0, D0)
  Ug <- xtr * (Atr - etr); U1 <- xtr * (Atr * (Ytr - link(lt1))); U0 <- xtr * ((1 - Atr) * (Ytr - link(lt0)))
  corr <- as.vector(Ug %*% vg + U1 %*% v1 + U0 %*% v0) * (nte / ntr)
  list(psi = psi, corr = corr, ntrim = sum(trimmed), e = eraw)
}
run_pipelines <- function(Z, A, Y, family, seed, trim = TRIM) {
  set.seed(seed); n <- nrow(Z); p <- ncol(Z); folds <- sample(rep(1:K, length.out = n))
  arms <- c("full", "out_thr", "out_bic", "oal_fix", "oal_wamd", "oads_thr", "oads_bic", "trt_bic", "bch")
  psi <- matrix(NA_real_, n, length(arms), dimnames = list(NULL, arms)); corr <- matrix(0, n, length(arms), dimnames = list(NULL, arms))
  ntrim <- setNames(numeric(length(arms)), arms)
  selmat <- matrix(0, p, length(arms), dimnames = list(colnames(Z), arms)); extra <- list(); iv_cand <- list()
  for (k in 1:K) {
    tr <- folds != k; te <- !tr; ntr <- sum(tr)
    b <- pilot(Z[tr, ], A[tr], Y[tr], family)
    s_thr <- sel_thr(b, ntr); s_bic <- sel_bic(Z[tr, ], Y[tr], family, unpen = cbind(A[tr]))
    s_oal <- sel_oal_fixed(Z[tr, ], A[tr], b); s_w <- sel_oal_wamd(Z[tr, ], A[tr], b)$sel
    s_trt <- sel_bic(Z[tr, ], A[tr], "binomial")
    sets <- list(full = seq_len(p), out_thr = s_thr, out_bic = s_bic, oal_fix = s_oal, oal_wamd = s_w,
                 oads_thr = sort(union(s_thr, s_oal)), oads_bic = sort(union(s_bic, s_oal)),
                 trt_bic = s_trt, bch = sort(union(s_bic, s_trt)))
    extra[[k]] <- list(thr = colnames(Z)[setdiff(s_oal, s_thr)], bic = colnames(Z)[setdiff(s_oal, s_bic)])
    iv_cand[[k]] <- colnames(Z)[setdiff(s_trt, sets$oads_thr)]
    cache <- list()
    for (a in arms) {
      S <- sets[[a]]; key <- paste0("S", paste(S, collapse = ","))
      if (is.null(cache[[key]])) cache[[key]] <- refit_fold(Z, A, Y, tr, te, S, family, trim = trim)
      r <- cache[[key]]; psi[te, a] <- r$psi; corr[tr, a] <- corr[tr, a] + r$corr; ntrim[a] <- ntrim[a] + r$ntrim
      selmat[S, a] <- selmat[S, a] + 1 / K
    }
  }
  est <- colMeans(psi); se_emp <- apply(psi, 2, sd) / sqrt(n)
  se_stk <- sqrt(colMeans(sweep(psi + corr, 2, est)^2) / n)
  tab <- data.frame(arm = arms, estimate = est, se_emp = se_emp, se_stk = se_stk,
                    mean_size = colSums(selmat), trim_frac = ntrim / n, row.names = NULL)
  # practitioner pipeline: naive OAL (wAMD) + Hajek IPW, fixed-weight sandwich SE
  b <- pilot(Z, A, Y, family); w <- sel_oal_wamd(Z, A, b); e <- pmin(pmax(w$e, trim), 1 - trim)
  h <- hajek(A, Y, e)
  tab <- rbind(tab, data.frame(arm = "oal_ipw", estimate = h[1], se_emp = h[2], se_stk = NA, mean_size = length(w$sel),
                               trim_frac = mean(w$e < trim | w$e > 1 - trim)))
  list(table = tab, selmat = selmat, extra = extra, iv_cand = iv_cand, folds = folds)
}

# --------------------------------------------------------------- data
load(DATA_PATH); f <- rhc
labs <- sapply(names(f), function(v) { l <- attr(f[[v]], "label"); if (is.null(l)) NA_character_ else l })
D <- as.integer(f$swang1 == "RHC")
Y_surv <- as.integer(f$dth30 == "No")
LOS <- as.numeric(f$dschdte - f$sadmdte)
drop <- c("ptid", "sadmdte", "dschdte", "dthdte", "lstctdte", "death", "dth30", "t3d30", "swang1", "adld3p")
X <- f[, setdiff(names(f), drop)]
X$cat2 <- factor(ifelse(is.na(as.character(X$cat2)), "None", as.character(X$cat2)))
X$urin1_miss <- as.integer(is.na(X$urin1)); X$urin1[is.na(X$urin1)] <- median(X$urin1, na.rm = TRUE)
X$wtkilo1_zero <- as.integer(X$wtkilo1 == 0)
for (v in names(X)) if (is.numeric(X[[v]]) && any(is.na(X[[v]]))) X[[v]][is.na(X[[v]])] <- median(X[[v]], na.rm = TRUE)
X <- droplevels(X)
Xm <- model.matrix(~ ., data = X)[, -1]
Xm <- Xm[, apply(Xm, 2, sd) > 0]
Zs <- scale(Xm); attr(Zs, "scaled:center") <- NULL; attr(Zs, "scaled:scale") <- NULL
cat("RHC pool: n =", nrow(Zs), " p =", ncol(Zs), "\n")
# pool listing with SUPPORT labels
src_var <- sapply(colnames(Xm), function(cn) { m <- names(X)[startsWith(cn, names(X))]; m[which.max(nchar(m))] })
pool <- data.frame(column = colnames(Xm), source = src_var, label = labs[src_var], row.names = NULL)
pool$label[is.na(pool$label)] <- c(cat1 = "primary disease category", cat2 = "secondary disease category (None if absent)",
  ca = "cancer (none / yes / metastatic)", sex = "sex", dnr1 = "DNR status on day 1", ninsclas = "insurance class",
  resp = "respiratory diagnosis", card = "cardiovascular diagnosis", neuro = "neurological diagnosis",
  gastr = "gastrointestinal diagnosis", renal = "renal diagnosis", meta = "metabolic diagnosis",
  hema = "hematologic diagnosis", seps = "sepsis diagnosis", trauma = "trauma diagnosis", ortho = "orthopedic diagnosis",
  race = "race", income = "income", urin1_miss = "urine output day 1 missing (indicator)",
  wtkilo1_zero = "weight recorded as 0 (missing weight indicator)")[pool$source[is.na(pool$label)]]
writeLines(c(sprintf("Candidate pool for both RHC outcomes: %d dummy-coded baseline columns (standardised), n = %d.", ncol(Zs), nrow(Zs)),
             "Excluded: identifiers and dates, death/dth30/t3d30 (outcomes), swang1 (treatment), adld3p (day-3 ADL, post-baseline).",
             "urin1 median-imputed with a missingness indicator; wtkilo1 = 0 (missing weight) flagged; cat2 NA coded as 'None'.", "",
             sprintf("%-28s %-12s %s", "column", "source", "SUPPORT label"),
             sprintf("%-28s %-12s %s", pool$column, pool$source, pool$label)), "output/applications2_pool.txt")

# ---------------------------------------------------------- analyses
run_outcome <- function(Y, Z, family, label, seed_main) {
  keep <- !is.na(Y); Z <- Z[keep, , drop = FALSE]; A <- D[keep]; Y <- Y[keep]
  main <- run_pipelines(Z, A, Y, family, seed_main)
  tb <- main$table; tb$outcome <- label; tb$n <- length(Y); tb$p <- ncol(Z)
  tb$lo_emp <- tb$estimate - 1.96 * tb$se_emp; tb$hi_emp <- tb$estimate + 1.96 * tb$se_emp
  tb$lo_stk <- tb$estimate - 1.96 * tb$se_stk; tb$hi_stk <- tb$estimate + 1.96 * tb$se_stk
  tb$lo_w1 <- tb$estimate - (1.96 + 1) * tb$se_stk; tb$hi_w1 <- tb$estimate + (1.96 + 1) * tb$se_stk
  # stability over fold splits
  sp <- do.call(rbind, lapply(seq_len(N_SPLITS), function(s) {
    r <- if (s == 1) main else run_pipelines(Z, A, Y, family, seed_main + s)
    data.frame(outcome = label, split = s, r$table[, c("arm", "estimate", "se_emp", "se_stk", "mean_size")]) }))
  # sensitivity: trimming at 0.01 (the level used by apply_realdata.R in the previous version), main split
  sens <- run_pipelines(Z, A, Y, family, seed_main, trim = 0.01)$table; sens$outcome <- label; sens$trim <- 0.01
  # selection frequencies and the S_oal \ S_out sets
  sel <- data.frame(outcome = label, column = rownames(main$selmat), main$selmat, row.names = NULL)
  ex_thr <- table(unlist(lapply(main$extra, `[[`, "thr"))) / K
  ex_bic <- table(unlist(lapply(main$extra, `[[`, "bic"))) / K
  ivc <- table(unlist(main$iv_cand)) / K
  # overlap diagnostic (full-set logistic on the whole sample)
  ps <- suppressWarnings(glm.fit(cbind(1, Z), A, family = binomial()))$fitted.values
  ov <- c(quantile(ps, c(0, .01, .05, .25, .5, .75, .95, .99, 1)), below05 = mean(ps < .05), above95 = mean(ps > .95))
  list(table = tb, splits = sp, sel = sel, ex_thr = ex_thr, ex_bic = ex_bic, ivc = ivc, overlap = ov, sens = sens)
}
t0 <- Sys.time()
r_surv <- run_outcome(Y_surv, Zs, "binomial", "30-day survival (risk difference)", 20260901)
r_los  <- run_outcome(LOS, Zs, "gaussian", "length of stay (days)", 20260902)
cat("time:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

tab <- rbind(r_surv$table, r_los$table)
tab <- tab[, c("outcome", "n", "p", "arm", "estimate", "se_emp", "se_stk", "lo_emp", "hi_emp", "lo_stk", "hi_stk",
               "lo_w1", "hi_w1", "mean_size", "trim_frac")]
write.csv(tab, "output/applications2_table.csv", row.names = FALSE)
write.csv(rbind(r_surv$splits, r_los$splits), "output/applications2_splits.csv", row.names = FALSE)
write.csv(rbind(r_surv$sens, r_los$sens), "output/applications2_trim01.csv", row.names = FALSE)
write.csv(rbind(r_surv$sel, r_los$sel), "output/applications2_selection.csv", row.names = FALSE)
sink("output/applications2_table.txt")
cat("apply_realdata2.R --", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("se_emp: empirical variance of AIPW pseudo-observations (Hajek fixed-weight sandwich for oal_ipw); se_stk: stacked sandwich (B.1).\n")
cat("lo/hi_w1: Theorem-3 widened interval with B = one stacked SE. mean_size: mean selected-set size over the 5 folds.\n\n")
print(format(tab, digits = 4), row.names = FALSE)
for (r in list(r_surv, r_los)) {
  cat("\n---", r$table$outcome[1], "---\n")
  cat("S_oal \\ S_out(thr), frequency over folds:\n"); print(r$ex_thr)
  cat("S_oal \\ S_out(bic), frequency over folds:\n"); print(r$ex_bic)
  cat("selected by the treatment lasso but not by ODS-thr (candidate treatment-only predictors), frequency over folds:\n"); print(sort(r$ivc, decreasing = TRUE))
  cat("full-set propensity score quantiles and tail mass:\n"); print(round(r$overlap, 4))
  cat("sensitivity: trimming at 0.01 (main split):\n")
  print(format(r$sens[, c("arm", "estimate", "se_emp", "se_stk", "mean_size", "trim_frac")], digits = 4), row.names = FALSE)
  cat("split stability (estimate range over", N_SPLITS, "fold splits):\n")
  sp <- r$splits; arms_u <- unique(sp$arm)
  print(format(data.frame(arm = arms_u,
                          est_min = tapply(sp$estimate, sp$arm, min)[arms_u], est_max = tapply(sp$estimate, sp$arm, max)[arms_u],
                          se_stk_min = tapply(sp$se_stk, sp$arm, min)[arms_u], se_stk_max = tapply(sp$se_stk, sp$arm, max)[arms_u]),
               digits = 4), row.names = FALSE)
}
sink()
writeLines(capture.output({ for (r in list(r_surv, r_los)) { cat(r$table$outcome[1], "\n"); print(round(r$overlap, 4)) } }),
           "output/applications2_overlap.txt")
cat(readLines("output/applications2_table.txt"), sep = "\n")
