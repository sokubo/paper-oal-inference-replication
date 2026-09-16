# sim_battery3.R — Paper 3: competitor arms + misspecification regime.
# Arms: oracle / naive OAL-IPW / OAL-refit-OLS / OADS / BCH / RY-OADS v2 /
#   GOAL-IPW (generalized OAL: elastic-net alpha=.5, ridge pilot, wAMD
#   tuning; Hajek IPW + sandwich, the field's native practice) /
#   C-TMLE (lasso-path collaborative TMLE, 2-fold CV loss, EIF-based SE).
# Cells: Regimes I-III as before; Regime M = misspecified outcome model
#   (Hermite cubic in the two strong confounders; PS form stays correct,
#   so doubly robust arms survive and singly robust refits break);
#   growing-p P1/P2.
suppressMessages(library(glmnet)); TAU <- 1; K <- 5

bic_select <- function(x, y, family, unpen = NULL) {
  if (is.null(unpen)) { X <- x; pf <- rep(1, ncol(x)); off <- 0 }
  else { X <- cbind(unpen, x); pf <- c(rep(0, ncol(unpen)), rep(1, ncol(x))); off <- ncol(unpen) }
  fit <- suppressWarnings(glmnet(X, y, family = family, nlambda = 40, penalty.factor = pf))
  bic <- deviance(fit) + fit$df * log(length(y))
  cf <- as.vector(coef(fit, s = fit$lambda[which.min(bic)]))[-1]
  which(cf[(off+1):length(cf)] != 0)
}
pilot <- function(Z, A, Y) {
  n <- nrow(Z); p <- ncol(Z)
  if (p <= n/3) { b <- coef(lm(Y ~ A + Z))[-(1:2)]; b[is.na(b)] <- 0; b }
  else { f <- suppressWarnings(glmnet(cbind(A, Z), Y, alpha = 0, lambda = 0.1,
           penalty.factor = c(0, rep(1, p)))); as.vector(coef(f))[-(1:2)] }
}
oal_fixed <- function(Z, A, Y) {
  n <- nrow(Z); b <- pilot(Z, A, Y)
  Zs <- sweep(Z, 2, pmax(abs(b)^2, 1e-8), "*"); lam <- n^(-0.75)
  fit <- suppressWarnings(glmnet(Zs, A, family="binomial", standardize=FALSE, lambda=c(2*lam,lam)))
  which(as.vector(coef(fit, s = lam))[-1] != 0)
}
oal_wamd <- function(Z, A, Y, alpha = 1) {
  n <- nrow(Z); p <- ncol(Z); b <- pilot(Z, A, Y)
  Zs <- sweep(Z, 2, pmax(abs(b)^2, 1e-8), "*")
  lg <- sort(n^c(-10,-5,-1,-.75,-.5,-.25,.25,.49)/n, decreasing = TRUE)
  fit <- suppressWarnings(glmnet(Zs, A, family="binomial", standardize=FALSE,
                                 lambda=lg, alpha = alpha))
  crit <- rep(Inf, length(lg)); es <- vector("list", length(lg))
  for (k in seq_along(lg)) {
    e <- tryCatch(as.vector(predict(fit, Zs, s=lg[k], type="response")), error=function(x) NULL)
    if (is.null(e)) next
    e <- pmin(pmax(e,.01),.99); es[[k]] <- e; ipw <- A/e + (1-A)/(1-e)
    smd <- vapply(1:p, function(j) abs(weighted.mean(Z[A==1,j],ipw[A==1]) -
             weighted.mean(Z[A==0,j],ipw[A==0]))/sd(Z[,j]), 0)
    crit[k] <- sum(abs(b)*smd)
  }
  k <- which.min(crit)
  list(sel = which(as.vector(coef(fit, s=lg[k]))[-1] != 0), e = es[[k]])
}
hajek <- function(A, Y, e) {
  ipw <- A/e + (1-A)/(1-e); X <- cbind(1, A)
  br <- solve(crossprod(X*sqrt(ipw))); bh <- br %*% crossprod(X*ipw, Y)
  r <- Y - X %*% bh; V <- br %*% crossprod(X*(ipw*as.vector(r))) %*% br
  c(bh[2], sqrt(V[2,2]))
}
aipw_cf <- function(Z, A, Y, sel_fun, Yest = NULL, trim = .05) {
  if (is.null(Yest)) Yest <- Y
  n <- nrow(Z); folds <- sample(rep(1:K, length.out = n)); psi <- numeric(n)
  for (k in 1:K) {
    tr <- folds != k; te <- !tr
    S <- sel_fun(Z[tr,,drop=FALSE], A[tr], Y[tr])
    ps <- suppressWarnings(glm.fit(cbind(1, Z[tr,S,drop=FALSE]), A[tr], family=binomial()))
    pc <- ps$coefficients; pc[is.na(pc)] <- 0
    e <- pmin(pmax(as.vector(plogis(cbind(1, Z[te,S,drop=FALSE]) %*% pc)), trim), 1-trim)
    c1 <- qr.coef(qr(cbind(1, Z[tr & A==1, S, drop=FALSE])), Yest[tr & A==1]); c1[is.na(c1)] <- 0
    c0 <- qr.coef(qr(cbind(1, Z[tr & A==0, S, drop=FALSE])), Yest[tr & A==0]); c0[is.na(c0)] <- 0
    m1 <- as.vector(cbind(1, Z[te,S,drop=FALSE]) %*% c1)
    m0 <- as.vector(cbind(1, Z[te,S,drop=FALSE]) %*% c0)
    psi[te] <- m1 - m0 + A[te]*(Yest[te]-m1)/e - (1-A[te])*(Yest[te]-m0)/(1-e)
  }
  c(mean(psi), sd(psi)/sqrt(n))
}
sel_oads <- function(z,a,y) union(bic_select(z,y,"gaussian",unpen=cbind(a)), oal_fixed(z,a,y))
sel_bch  <- function(z,a,y) union(bic_select(z,y,"gaussian",unpen=cbind(a)), bic_select(z,a,"binomial"))
ry2 <- function(Z, A, Y) {
  n <- nrow(Z)
  fullfit <- if (ncol(Z) < n - 2) lm(Y ~ A + Z) else
    lm(Y ~ A + Z[, seq_len(floor(n/2)), drop=FALSE])
  sig <- sqrt(sum(resid(fullfit)^2) / max(fullfit$df.residual, 1))
  W <- rnorm(n, sd = sig); U <- Y + W; V <- Y - W
  S <- sel_oads(Z, A, U)
  aipw_cf(Z, A, V, function(z,a,y) S, Yest = V)
}
ctmle <- function(Z, A, Y, nlam = 8, trim = .05) {
  n <- nrow(Z)
  SQ <- bic_select(Z, Y, "gaussian", unpen = cbind(A))
  qc <- qr.coef(qr(cbind(1, A, Z[, SQ, drop=FALSE])), Y); qc[is.na(qc)] <- 0
  QA <- as.vector(cbind(1, A, Z[,SQ,drop=FALSE]) %*% qc)
  Q1 <- as.vector(cbind(1, 1, Z[,SQ,drop=FALSE]) %*% qc)
  Q0 <- as.vector(cbind(1, 0, Z[,SQ,drop=FALSE]) %*% qc)
  ps_path <- suppressWarnings(glmnet(Z, A, family = "binomial"))
  lam_idx <- unique(round(seq(1, length(ps_path$lambda), length.out = nlam)))
  lams <- ps_path$lambda[lam_idx]
  half <- sample(rep(1:2, length.out = n))
  cvloss <- rep(Inf, length(lams))
  for (k in seq_along(lams)) {
    g <- pmin(pmax(as.vector(predict(ps_path, Z, s = lams[k],
                                     type = "response")), trim), 1-trim)
    H <- A/g - (1-A)/(1-g)
    loss <- 0
    for (h in 1:2) {
      inh <- half == h; outh <- !inh
      eps <- tryCatch(coef(lm.fit(cbind(H[inh]), Y[inh] - QA[inh]))[1],
                      error = function(e) 0)
      if (!is.finite(eps)) eps <- 0
      loss <- loss + sum((Y[outh] - (QA[outh] + eps * H[outh]))^2)
    }
    cvloss[k] <- loss
  }
  k <- which.min(cvloss)
  g <- pmin(pmax(as.vector(predict(ps_path, Z, s = lams[k],
                                   type = "response")), trim), 1-trim)
  H <- A/g - (1-A)/(1-g)
  eps <- coef(lm.fit(cbind(H), Y - QA))[1]; if (!is.finite(eps)) eps <- 0
  Q1s <- Q1 + eps / g; Q0s <- Q0 - eps / (1-g); QAs <- QA + eps * H
  psi <- mean(Q1s - Q0s)
  IC <- H * (Y - QAs) + Q1s - Q0s - psi
  c(psi, sd(IC)/sqrt(n))
}

METH <- c("oracle","oal_ipw","oal_ols","oads","bch","ry2","goal_ipw","ctmle")
run_cell <- function(n, p, label, seed, R, dgp, oracle_fit) {
  set.seed(seed)
  res <- array(NA_real_, c(R, length(METH), 2))
  for (r in 1:R) {
    d <- dgp(n, p); Z <- d$Z; A <- d$A; Y <- d$Y
    for (mi in seq_along(METH)) {
      res[r,mi,] <- switch(METH[mi],
        oracle  = oracle_fit(d),
        oal_ipw = { ow <- oal_wamd(Z,A,Y); hajek(A, Y, ow$e) },
        oal_ols = { ow <- oal_wamd(Z,A,Y)
                    m2 <- if (length(ow$sel)) lm(Y ~ A + Z[,ow$sel,drop=FALSE]) else lm(Y ~ A)
                    summary(m2)$coefficients["A",1:2] },
        oads    = aipw_cf(Z, A, Y, sel_oads),
        bch     = aipw_cf(Z, A, Y, sel_bch),
        ry2     = ry2(Z, A, Y),
        goal_ipw= { gw <- oal_wamd(Z,A,Y, alpha=.5); hajek(A, Y, gw$e) },
        ctmle   = ctmle(Z, A, Y))
    }
  }
  do.call(rbind, lapply(seq_along(METH), function(i) {
    est <- res[,i,1]; se <- res[,i,2]
    data.frame(cell = label, method = METH[i], bias = mean(est)-TAU,
               emp_sd = sd(est), se_ratio = mean(se)/sd(est),
               cov95 = mean(abs(est-TAU) <= 1.96*se),
               width_rel_oracle = mean(se)/mean(res[,1,2]))
  }))
}

lin_dgp <- function(a3, b3, sd_y) function(n, p) {
  alpha <- c(.8,.8,a3,a3,0,0,1,1, rep(0, p-8))
  beta  <- c(.6,.6,b3,b3,.6,.6,0,0, rep(0, p-8))
  Z <- matrix(rnorm(n*p), n, p)
  A <- rbinom(n, 1, plogis(Z %*% alpha))
  list(Z=Z, A=A, Y = TAU*A + as.vector(Z %*% beta) + rnorm(n, sd=sd_y),
       beta=beta)
}
lin_oracle <- function(d)
  summary(lm(d$Y ~ d$A + d$Z[, which(d$beta!=0), drop=FALSE]))$coefficients[2,1:2]

mis_dgp <- function(n, p) {
  alpha <- c(.8,.8,.8,.8,0,0,1,1, rep(0, p-8))
  beta  <- c(.6,.6,.6,.6,.6,.6,0,0, rep(0, p-8))
  Z <- matrix(rnorm(n*p), n, p)
  A <- rbinom(n, 1, plogis(Z %*% alpha))
  nl <- 0.6*(Z[,1]^3 - 3*Z[,1] + Z[,2]^3 - 3*Z[,2])/sqrt(6)
  list(Z=Z, A=A, Y = TAU*A + as.vector(Z %*% beta) + nl + rnorm(n, sd=2),
       beta=beta, nl=nl)
}
mis_oracle <- function(d) {   # knows the true functional form
  H <- (d$Z[,1]^3 - 3*d$Z[,1] + d$Z[,2]^3 - 3*d$Z[,2])/sqrt(6)
  summary(lm(d$Y ~ d$A + d$Z[, which(d$beta!=0), drop=FALSE] + H))$coefficients[2,1:2]
}

args <- commandArgs(trailingOnly = TRUE)
Rmain <- if (length(args)) as.integer(args[1]) else 500
Rbig  <- if (length(args) > 1) as.integer(args[2]) else 300
t0 <- Sys.time()
tabs <- rbind(
  run_cell(500, 20, "I (n=500,p=20)",    701, Rmain, lin_dgp(.8,.6,2),  lin_oracle),
  run_cell(200, 50, "II (n=200,p=50)",   702, Rmain, lin_dgp(.7,.25,2), lin_oracle),
  run_cell(2000,20, "III (n=2000,p=20)", 703, Rmain, lin_dgp(1,3/sqrt(2000),2), lin_oracle),
  run_cell(500, 20, "M (misspec OM)",    704, Rmain, mis_dgp,           mis_oracle),
  run_cell(500, 100,"P1 (n=500,p=100)",  705, Rbig,  lin_dgp(.7,.25,2), lin_oracle),
  run_cell(500, 200,"P2 (n=500,p=200)",  706, Rbig,  lin_dgp(.7,.25,2), lin_oracle))
cat("time:", round(difftime(Sys.time(), t0, units="mins"),1), "min\n")
print(tabs, digits = 3, row.names = FALSE)
write.csv(tabs, "output/battery3_table.csv", row.names = FALSE)
capture.output(print(tabs, digits=3), file = "output/battery3_summary.txt")
cat("battery3 done\n")
