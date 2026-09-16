# proposed_estimator.R — Paper 3: prototype of Outcome-Adaptive Double
# Selection (OADS) and comparators, evaluated in the two demo regimes.
#
# OADS: K-fold cross-fitting. In each training set:
#   S_out = outcome lasso (gaussian, BIC over glmnet path)
#   S_oal = OAL (adaptive logistic lasso, outcome-model weights, wAMD tuning)
#   S     = S_out UNION S_oal   (instruments excluded: neither arm admits them)
#   Refit unpenalized: logistic PS on S; OLS outcome models on S (by arm).
#   Held-out fold: AIPW pseudo-outcomes. SE = sd(EIF)/sqrt(n).
# BCH-DS: same but S = S_out UNION S_trt (treatment lasso, BIC) --
#   uniformly valid, but re-admits instruments (variance/positivity cost).
# Reported alongside the naive pipelines from undercoverage_demo.R.

suppressMessages(library(glmnet))
dir.create("output", showWarnings = FALSE)
TAU <- 1; R <- 500; K <- 5

bic_select <- function(x, y, family, unpen = NULL) {
  # unpen: matrix of always-included, unpenalized columns (e.g., treatment A
  # in the outcome arm). Selection reported for x columns only.
  if (is.null(unpen)) { X <- x; pf <- rep(1, ncol(x)); off <- 0 }
  else { X <- cbind(unpen, x); pf <- c(rep(0, ncol(unpen)), rep(1, ncol(x)))
         off <- ncol(unpen) }
  fit <- suppressWarnings(glmnet(X, y, family = family, nlambda = 40,
                                 penalty.factor = pf))
  n <- length(y)
  bic <- deviance(fit) + fit$df * log(n)
  cf <- as.vector(coef(fit, s = fit$lambda[which.min(bic)]))[-1]
  which(cf[(off + 1):length(cf)] != 0)
}

oal_select <- function(Z, A, Y, tune = "wamd") {
  n <- nrow(Z); p <- ncol(Z)
  b_out <- coef(lm(Y ~ A + Z))[-(1:2)]
  scl <- pmax(abs(b_out)^2, 1e-8)          # gamma = 2; x*_j = x_j * |b_j|^2
  Zs <- sweep(Z, 2, scl, "*")
  lam_grid <- sort(n^c(-10,-5,-1,-.75,-.5,-.25,.25,.49), decreasing = TRUE)
  fit <- suppressWarnings(glmnet(Zs, A, family = "binomial",
                                 standardize = FALSE, lambda = lam_grid))
  crit <- rep(Inf, length(lam_grid))
  for (k in seq_along(lam_grid)) {
    e <- tryCatch(as.vector(predict(fit, Zs, s = lam_grid[k], type = "response")),
                  error = function(x) NULL)
    if (is.null(e)) next
    e <- pmin(pmax(e, .01), .99)
    if (tune == "bic") {
      df <- sum(as.vector(coef(fit, s = lam_grid[k]))[-1] != 0)
      crit[k] <- -2*sum(A*log(e) + (1-A)*log(1-e)) + df*log(n)
    } else {
      ipw <- A/e + (1-A)/(1-e)
      smd <- vapply(1:p, function(j)
        abs(weighted.mean(Z[A==1,j], ipw[A==1]) -
            weighted.mean(Z[A==0,j], ipw[A==0])) / sd(Z[,j]), 0)
      crit[k] <- sum(abs(b_out) * smd)
    }
  }
  which(as.vector(coef(fit, s = lam_grid[which.min(crit)]))[-1] != 0)
}

aipw_crossfit <- function(Z, A, Y, select_fun) {
  n <- nrow(Z); folds <- sample(rep(1:K, length.out = n))
  psi <- numeric(n)
  for (k in 1:K) {
    tr <- folds != k; te <- !tr
    S <- select_fun(Z[tr,,drop=FALSE], A[tr], Y[tr])
    Ztr <- Z[tr, S, drop = FALSE]; Zte <- Z[te, S, drop = FALSE]
    # PS (unpenalized refit)
    ps_fit <- if (length(S)) glm.fit(cbind(1, Ztr), A[tr],
                 family = binomial()) else glm.fit(matrix(1, sum(tr)), A[tr],
                 family = binomial())
    e <- plogis(cbind(1, Zte) %*% ps_fit$coefficients)
    if (!length(S)) e <- rep(plogis(ps_fit$coefficients[1]), sum(te))
    e <- pmin(pmax(as.vector(e), .05), .95)
    # outcome models by arm
    for (a in 0:1) {
      ia <- tr & A == a
      cf <- qr.coef(qr(cbind(1, Z[ia, S, drop = FALSE])), Y[ia])
      cf[is.na(cf)] <- 0
      m <- as.vector(cbind(1, Zte) %*% cf)
      if (a == 1) m1 <- m else m0 <- m
    }
    psi[te] <- m1 - m0 + A[te]*(Y[te]-m1)/e - (1-A[te])*(Y[te]-m0)/(1-e)
  }
  c(mean(psi), sd(psi)/sqrt(n))
}

oal_select_fixed <- function(Z, A, Y) {
  # S-E oracle-condition tuning for gamma = 2: lambda_sum = n^{1/4}
  # (lambda_sum / sqrt(n) -> 0 and lambda_sum * n^{gamma/2 - 1} -> Inf).
  # glmnet uses the mean-log-likelihood scale: lambda = n^{1/4} / n = n^{-3/4}.
  n <- nrow(Z)
  b_out <- coef(lm(Y ~ A + Z))[-(1:2)]
  scl <- pmax(abs(b_out)^2, 1e-8)
  Zs <- sweep(Z, 2, scl, "*")
  lam <- n^(-0.75)
  fit <- suppressWarnings(glmnet(Zs, A, family = "binomial",
                                 standardize = FALSE,
                                 lambda = c(2*lam, lam)))
  which(as.vector(coef(fit, s = lam))[-1] != 0)
}

sel_oads <- function(Z, A, Y) union(bic_select(Z, Y, "gaussian", unpen = cbind(A)), oal_select_fixed(Z, A, Y))
sel_bch  <- function(Z, A, Y) union(bic_select(Z, Y, "gaussian", unpen = cbind(A)),
                                    bic_select(Z, A, "binomial"))

run_regime <- function(n, p, a3, b3, sd_y, seed, label) {
  set.seed(seed)
  alpha <- c(.8,.8,a3,a3,0,0,1,1, rep(0, p-8))
  beta  <- c(.6,.6,b3,b3,.6,.6,0,0, rep(0, p-8))
  res <- matrix(NA_real_, R, 6)
  for (r in 1:R) {
    Z <- matrix(rnorm(n*p), n, p)
    A <- rbinom(n, 1, plogis(Z %*% alpha))
    Y <- TAU*A + as.vector(Z %*% beta) + rnorm(n, sd = sd_y)
    res[r, 1:2] <- aipw_crossfit(Z, A, Y, sel_oads)
    res[r, 3:4] <- aipw_crossfit(Z, A, Y, sel_bch)
    # instrument retention check (full-data selection, for reporting)
    if (r <= 100) {
      so <- sel_oads(Z, A, Y); sb <- sel_bch(Z, A, Y)
      res[r, 5] <- mean(c(7, 8) %in% so)   # IVs are Z7, Z8
      res[r, 6] <- mean(c(7, 8) %in% sb)
      if (r == 1) extra <<- NULL
      extra <<- rbind(extra, c(mean(c(3, 4) %in% so), mean(c(3, 4) %in% sb)))
    }
  }
  meth <- c("OADS + AIPW (proposed)", "BCH double selection + AIPW")
  out <- do.call(rbind, lapply(1:2, function(i) {
    est <- res[, (i-1)*2 + 1]; se <- res[, (i-1)*2 + 2]
    data.frame(regime = label, method = meth[i],
               bias = mean(est) - TAU, emp_sd = sd(est),
               se_ratio = mean(se)/sd(est),
               coverage95 = mean(abs(est - TAU) <= 1.96*se))
  }))
  out$iv_retention <- c(mean(res[1:100, 5], na.rm = TRUE),
                        mean(res[1:100, 6], na.rm = TRUE))
  out$boundary_retention <- colMeans(extra)
  out
}

t0 <- Sys.time()
tabs <- rbind(
  run_regime(500, 20, .8, .6,  2, 301, "I: benign (n=500, p=20)"),
  run_regime(200, 50, .7, .25, 2, 302, "II: many controls (n=200, p=50)"))
cat("time:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
print(tabs, digits = 3, row.names = FALSE)
write.csv(tabs, "output/proposed_table.csv", row.names = FALSE)
capture.output(print(tabs, digits = 3), file = "output/proposed_summary.txt")
