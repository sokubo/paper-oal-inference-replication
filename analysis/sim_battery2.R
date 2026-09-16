# sim_battery2.R — RY variant v2 (df-corrected sigma + V-side cross-fitted
# refits) in Regimes I-III, plus growing-p cells (n=500, p in {100, 200})
# for naive OAL / OADS / BCH / RY2. Ridge pilot when p > n/3.
suppressMessages(library(glmnet)); set.seed(20260820)
R1 <- 500; R2 <- 300; K <- 5; TAU <- 1

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
oal_wamd <- function(Z, A, Y) {
  n <- nrow(Z); p <- ncol(Z); b <- pilot(Z, A, Y)
  Zs <- sweep(Z, 2, pmax(abs(b)^2, 1e-8), "*")
  lg <- sort(n^c(-10,-5,-1,-.75,-.5,-.25,.25,.49)/n, decreasing = TRUE)
  fit <- suppressWarnings(glmnet(Zs, A, family="binomial", standardize=FALSE, lambda=lg))
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
  # selection uses Y; refits+evaluation use Yest (defaults to Y)
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
  n <- nrow(Z); p <- ncol(Z)
  fullfit <- lm(Y ~ A + Z)
  sig <- sqrt(sum(resid(fullfit)^2) / max(fullfit$df.residual, 1))  # df-corrected
  W <- rnorm(n, sd = sig); U <- Y + W; V <- Y - W
  S <- sel_oads(Z, A, U)
  # V-side cross-fitted refits with FIXED selection S
  aipw_cf(Z, A, V, function(z,a,y) S, Yest = V)
}

run_cell <- function(n, p, a3, b3, sd_y, label, seed, R, methods) {
  set.seed(seed)
  alpha <- c(.8,.8,a3,a3,0,0,1,1, rep(0, p-8))
  beta  <- c(.6,.6,b3,b3,.6,.6,0,0, rep(0, p-8))
  res <- array(NA_real_, c(R, length(methods), 2))
  for (r in 1:R) {
    Z <- matrix(rnorm(n*p), n, p)
    A <- rbinom(n, 1, plogis(Z %*% alpha))
    Y <- TAU*A + as.vector(Z %*% beta) + rnorm(n, sd = sd_y)
    for (mi in seq_along(methods)) {
      m <- methods[mi]
      res[r,mi,] <- switch(m,
        oracle  = summary(lm(Y ~ A + Z[, which(beta!=0), drop=FALSE]))$coefficients["A",1:2],
        oal_ipw = { ow <- oal_wamd(Z,A,Y); hajek(A, Y, ow$e) },
        oal_ols = { ow <- oal_wamd(Z,A,Y)
                    m2 <- if (length(ow$sel)) lm(Y ~ A + Z[,ow$sel,drop=FALSE]) else lm(Y ~ A)
                    summary(m2)$coefficients["A",1:2] },
        oads    = aipw_cf(Z, A, Y, sel_oads),
        bch     = aipw_cf(Z, A, Y, sel_bch),
        ry2     = ry2(Z, A, Y))
    }
  }
  do.call(rbind, lapply(seq_along(methods), function(i) {
    est <- res[,i,1]; se <- res[,i,2]
    data.frame(cell = label, method = methods[i], bias = mean(est)-TAU,
               emp_sd = sd(est), se_ratio = mean(se)/sd(est),
               cov95 = mean(abs(est-TAU) <= 1.96*se))
  }))
}

t0 <- Sys.time()
m6 <- c("oracle","oal_ipw","oal_ols","oads","bch","ry2")
tabs <- rbind(
  run_cell(500, 20, .8, .6,        2, "I (n=500,p=20)",   601, R1, c("oracle","ry2")),
  run_cell(200, 50, .7, .25,       2, "II (n=200,p=50)",  602, R1, c("oracle","ry2")),
  run_cell(2000,20, 1.0, 3/sqrt(2000),2,"III (n=2000,p=20)",603, R1, c("oracle","ry2")),
  run_cell(500, 100, .7, .25,      2, "P1 (n=500,p=100)", 604, R2, m6),
  run_cell(500, 200, .7, .25,      2, "P2 (n=500,p=200)", 605, R2, m6))
cat("time:", round(difftime(Sys.time(), t0, units="mins"),1), "min\n")
print(tabs, digits = 3, row.names = FALSE)
write.csv(tabs, "output/battery2_table.csv", row.names = FALSE)
capture.output(print(tabs, digits=3), file = "output/battery2_summary.txt")
cat("battery2 done\n")
