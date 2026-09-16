# sim_battery.R — Paper 3, Section 6: full simulation battery.
# Regimes:
#  I   benign            n=500,  p=20, boundary vars strong (a=.8, b=.6)
#  II  many controls     n=200,  p=50, boundary vars (a=.7, b=.25)
#  III local boundary    n=2000, p=20, boundary vars on the sqrt(n) scale:
#      a=1.0, b = 3/sqrt(n)  -> separation FAILS by construction (Thm 1 zone)
# Estimators: oracle OLS; naive OAL+IPW; naive OAL+refit OLS; OADS+AIPW;
#  BCH-DS+AIPW; RY-OADS (Rasines-Young U/V randomization, gamma=1);
#  OADS widened intervals (Theorem 3, B in {1,2} x sigma-hat units).
# R = 500 reps. Outputs: output/battery_table.csv/.txt

suppressMessages(library(glmnet))
set.seed(20260819); R <- 500; K <- 5; TAU <- 1

bic_select <- function(x, y, family, unpen = NULL) {
  if (is.null(unpen)) { X <- x; pf <- rep(1, ncol(x)); off <- 0 }
  else { X <- cbind(unpen, x); pf <- c(rep(0, ncol(unpen)), rep(1, ncol(x)))
         off <- ncol(unpen) }
  fit <- suppressWarnings(glmnet(X, y, family = family, nlambda = 40,
                                 penalty.factor = pf))
  bic <- deviance(fit) + fit$df * log(length(y))
  cf <- as.vector(coef(fit, s = fit$lambda[which.min(bic)]))[-1]
  which(cf[(off + 1):length(cf)] != 0)
}
oal_fixed <- function(Z, A, Y) {
  n <- nrow(Z)
  b_out <- coef(lm(Y ~ A + Z))[-(1:2)]; b_out[is.na(b_out)] <- 0
  scl <- pmax(abs(b_out)^2, 1e-8)
  Zs <- sweep(Z, 2, scl, "*"); lam <- n^(-0.75)
  fit <- suppressWarnings(glmnet(Zs, A, family = "binomial",
                                 standardize = FALSE, lambda = c(2*lam, lam)))
  which(as.vector(coef(fit, s = lam))[-1] != 0)
}
oal_wamd <- function(Z, A, Y) {
  n <- nrow(Z); p <- ncol(Z)
  b_out <- coef(lm(Y ~ A + Z))[-(1:2)]; b_out[is.na(b_out)] <- 0
  scl <- pmax(abs(b_out)^2, 1e-8); Zs <- sweep(Z, 2, scl, "*")
  lam_grid <- sort(n^c(-10,-5,-1,-.75,-.5,-.25,.25,.49)/n, decreasing = TRUE)
  fit <- suppressWarnings(glmnet(Zs, A, family="binomial",
                                 standardize = FALSE, lambda = lam_grid))
  crit <- rep(Inf, length(lam_grid)); es <- vector("list", length(lam_grid))
  for (k in seq_along(lam_grid)) {
    e <- tryCatch(as.vector(predict(fit, Zs, s=lam_grid[k], type="response")),
                  error = function(x) NULL)
    if (is.null(e)) next
    e <- pmin(pmax(e, .01), .99); es[[k]] <- e
    ipw <- A/e + (1-A)/(1-e)
    smd <- vapply(1:p, function(j)
      abs(weighted.mean(Z[A==1,j], ipw[A==1]) -
          weighted.mean(Z[A==0,j], ipw[A==0]))/sd(Z[,j]), 0)
    crit[k] <- sum(abs(b_out) * smd)
  }
  k <- which.min(crit)
  list(sel = which(as.vector(coef(fit, s = lam_grid[k]))[-1] != 0),
       e = es[[k]])
}
hajek_ipw <- function(A, Y, e) {
  ipw <- A/e + (1-A)/(1-e); X <- cbind(1, A)
  br <- solve(crossprod(X * sqrt(ipw))); bh <- br %*% crossprod(X * ipw, Y)
  r <- Y - X %*% bh
  V <- br %*% crossprod(X * (ipw * as.vector(r))) %*% br
  c(bh[2], sqrt(V[2,2]))
}
aipw_cf <- function(Z, A, Y, select_fun, trim = .05) {
  n <- nrow(Z); folds <- sample(rep(1:K, length.out = n)); psi <- numeric(n)
  for (k in 1:K) {
    tr <- folds != k; te <- !tr
    S <- select_fun(Z[tr,,drop=FALSE], A[tr], Y[tr])
    Ztr <- Z[tr, S, drop=FALSE]; Zte <- Z[te, S, drop=FALSE]
    ps <- suppressWarnings(glm.fit(cbind(1, Ztr), A[tr], family = binomial()))
    e <- pmin(pmax(as.vector(plogis(cbind(1, Zte) %*% ps$coefficients)),
                   trim), 1 - trim)
    m1c <- qr.coef(qr(cbind(1, Z[tr & A==1, S, drop=FALSE])), Y[tr & A==1]); m1c[is.na(m1c)] <- 0
    m0c <- qr.coef(qr(cbind(1, Z[tr & A==0, S, drop=FALSE])), Y[tr & A==0]); m0c[is.na(m0c)] <- 0
    m1 <- as.vector(cbind(1, Zte) %*% m1c); m0 <- as.vector(cbind(1, Zte) %*% m0c)
    psi[te] <- m1 - m0 + A[te]*(Y[te]-m1)/e - (1-A[te])*(Y[te]-m0)/(1-e)
  }
  c(mean(psi), sd(psi)/sqrt(n))
}
aipw_full <- function(Z, A, Y, S, Ysel = Y, trim = .05) {
  # single-sample refit AIPW using outcome Ysel for estimation (RY: V-part)
  ps <- suppressWarnings(glm.fit(cbind(1, Z[,S,drop=FALSE]), A, family=binomial()))
  e <- pmin(pmax(as.vector(plogis(cbind(1, Z[,S,drop=FALSE]) %*% ps$coefficients)),
                 trim), 1-trim)
  m1c <- qr.coef(qr(cbind(1, Z[A==1, S, drop=FALSE])), Ysel[A==1]); m1c[is.na(m1c)] <- 0
  m0c <- qr.coef(qr(cbind(1, Z[A==0, S, drop=FALSE])), Ysel[A==0]); m0c[is.na(m0c)] <- 0
  m1 <- as.vector(cbind(1, Z[,S,drop=FALSE]) %*% m1c)
  m0 <- as.vector(cbind(1, Z[,S,drop=FALSE]) %*% m0c)
  psi <- m1 - m0 + A*(Ysel-m1)/e - (1-A)*(Ysel-m0)/(1-e)
  c(mean(psi), sd(psi)/sqrt(length(Ysel)))
}
sel_oads <- function(z,a,y) union(bic_select(z, y, "gaussian", unpen = cbind(a)),
                                  oal_fixed(z, a, y))
sel_bch  <- function(z,a,y) union(bic_select(z, y, "gaussian", unpen = cbind(a)),
                                  bic_select(z, a, "binomial"))

run_regime <- function(n, p, a3, b3fun, sd_y, label, seed) {
  set.seed(seed)
  b3 <- b3fun(n)
  alpha <- c(.8,.8,a3,a3,0,0,1,1, rep(0, p-8))
  beta  <- c(.6,.6,b3,b3,.6,.6,0,0, rep(0, p-8))
  meth <- c("oracle","oal_ipw","oal_ols","oads","bch","ry_oads")
  res <- array(NA_real_, c(R, length(meth), 2))
  for (r in 1:R) {
    Z <- matrix(rnorm(n*p), n, p)
    A <- rbinom(n, 1, plogis(Z %*% alpha))
    Y <- TAU*A + as.vector(Z %*% beta) + rnorm(n, sd = sd_y)
    # oracle
    mo <- lm(Y ~ A + Z[, which(beta != 0), drop = FALSE])
    res[r,1,] <- summary(mo)$coefficients["A", 1:2]
    # naive OAL
    ow <- oal_wamd(Z, A, Y)
    res[r,2,] <- hajek_ipw(A, Y, ow$e)
    m2 <- if (length(ow$sel)) lm(Y ~ A + Z[, ow$sel, drop=FALSE]) else lm(Y ~ A)
    res[r,3,] <- summary(m2)$coefficients["A", 1:2]
    # OADS / BCH (cross-fitted)
    res[r,4,] <- aipw_cf(Z, A, Y, sel_oads)
    res[r,5,] <- aipw_cf(Z, A, Y, sel_bch)
    # RY-OADS: U/V randomization (gamma = 1), selection on U, estimation on V
    sig <- sd(resid(lm(Y ~ A + Z)))
    W <- rnorm(n, sd = sig)
    U <- Y + W; V <- Y - W
    S <- sel_oads(Z, A, U)
    res[r,6,] <- aipw_full(Z, A, Y, S, Ysel = V)
  }
  out <- do.call(rbind, lapply(seq_along(meth), function(i) {
    est <- res[,i,1]; se <- res[,i,2]
    data.frame(regime = label, method = meth[i],
               bias = mean(est) - TAU, emp_sd = sd(est),
               se_ratio = mean(se)/sd(est),
               cov95 = mean(abs(est - TAU) <= 1.96*se),
               width_rel_oracle = mean(se)/mean(res[,1,2]))
  }))
  # Theorem-3 widened intervals from the OADS run
  for (b in c(1, 2)) {
    est <- res[,4,1]; se <- res[,4,2]
    out <- rbind(out, data.frame(regime = label,
      method = sprintf("oads_widened_B%dsd", b),
      bias = mean(est) - TAU, emp_sd = sd(est), se_ratio = NA,
      cov95 = mean(abs(est - TAU) <= (1.96 + b)*se),
      width_rel_oracle = mean((1.96 + b)*se)/(1.96*mean(res[,1,2]))))
  }
  out
}

t0 <- Sys.time()
tabs <- rbind(
  run_regime(500, 20, .8, function(n) .6,        2, "I: benign (n=500,p=20)", 501),
  run_regime(200, 50, .7, function(n) .25,       2, "II: many controls (n=200,p=50)", 502),
  run_regime(2000,20, 1.0, function(n) 3/sqrt(n),2, "III: local boundary (n=2000,p=20)", 503))
cat("battery time:", round(difftime(Sys.time(), t0, units="mins"),1), "min\n")
print(tabs, digits = 3, row.names = FALSE)
write.csv(tabs, "output/battery_table.csv", row.names = FALSE)
capture.output(print(tabs, digits = 3), file = "output/battery_summary.txt")
cat("battery done\n")
