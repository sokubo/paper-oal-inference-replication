# ry_target_check.R -- Web Appendix D: coverage of the exploratory RY-ODS v2
# variant for (a) the ATE and (b) its own projection target
# tau(e*_S, m*_S) for the selected set S (the reduced pseudo-true AIPW
# functional of Lemma A3), in cells I, II and III of sim_battery4.R.
# The projection target is computed per replication from a population
# sample of N_POP draws of the same DGP (reduced logistic fit and per-arm
# least squares on x_S, then the population AIPW average).
# Usage: Rscript ry_target_check.R [n_reps] [cores]
# Output: output/ry_target_check.csv / .txt
suppressMessages({ library(glmnet); library(parallel) })
src <- readLines("sim_battery4.R"); cut <- grep("^# -+ run$", src); eval(parse(text = src[1:(cut - 1)]))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 300
CORES <- if (length(args) >= 2) as.integer(args[2]) else max(1L, detectCores())
N_POP <- 200000

proj_target <- function(pop, S) {
  x <- cbind(1, pop$Z[, S, drop = FALSE]); A <- pop$A; Y <- pop$Y
  g <- suppressWarnings(glm.fit(x, A, family = binomial()))$coefficients; g[is.na(g)] <- 0
  e <- as.vector(plogis(x %*% g))
  b1 <- qr.coef(qr(x[A == 1, , drop = FALSE]), Y[A == 1]); b1[is.na(b1)] <- 0
  b0 <- qr.coef(qr(x[A == 0, , drop = FALSE]), Y[A == 0]); b0[is.na(b0)] <- 0
  m1 <- as.vector(x %*% b1); m0 <- as.vector(x %*% b0)
  mean(m1 - m0 + A * (Y - m1) / e - (1 - A) * (Y - m0) / (1 - e))   # untrimmed population functional
}
ry_rep <- function(seed, dgp, pop) {
  set.seed(seed); d <- dgp(); Z <- d$Z; A <- d$A; Y <- d$Y; n <- nrow(Z); p <- ncol(Z)
  folds <- sample(rep(1:K, length.out = n))
  ff <- lm.fit(cbind(1, A, Z), Y); sig <- sqrt(sum(ff$residuals^2) / max(ff$df.residual, 1))
  W <- rnorm(n, sd = sig); U <- Y + W; V <- Y - W
  pl <- pilot(Z, A, U)
  S <- sort(union(sel_bic(Z, U, "gaussian", unpen = cbind(A)), sel_oal_fixed(Z, A, pl$b)))
  f <- aipw_from_sets(Z, A, V, folds, lapply(1:K, function(k) list(ry2 = S)))
  c(est = unname(f$est), se_emp = unname(f$se_emp), se_stk = unname(f$se_stk),
    target = proj_target(pop, S), size = length(S), superset = as.numeric(all(d$Sstar %in% S)))
}
cells <- list(I = list(dgp = lin_dgp(500, 20, .8, .6), seed = 801, popdgp = lin_dgp(N_POP, 20, .8, .6)),
              II = list(dgp = lin_dgp(200, 50, .7, .25), seed = 802, popdgp = lin_dgp(N_POP, 50, .7, .25)),
              III = list(dgp = lin_dgp(2000, 20, 1, 3 / sqrt(2000)), seed = 803, popdgp = lin_dgp(N_POP, 20, 1, 3 / sqrt(2000))))
out <- list()
for (cn in names(cells)) {
  cl <- cells[[cn]]; set.seed(cl$seed); pop <- cl$popdgp()
  res <- do.call(rbind, mclapply(cl$seed * 1000L + seq_len(R), function(s) ry_rep(s, cl$dgp, pop), mc.cores = CORES))
  cov <- function(center, se) mean(abs(res[, "est"] - center) <= 1.96 * se)
  out[[cn]] <- data.frame(cell = cn, n_reps = R,
    bias_ate = mean(res[, "est"]) - TAU, mean_target_minus_ate = mean(res[, "target"]) - TAU,
    sd_target = sd(res[, "target"]), emp_sd = sd(res[, "est"]),
    mean_se_emp = mean(res[, "se_emp"]), mean_se_stk = mean(res[, "se_stk"]),
    cov_ate_emp = cov(TAU, res[, "se_emp"]), cov_ate_stk = cov(TAU, res[, "se_stk"]),
    cov_target_emp = cov(res[, "target"], res[, "se_emp"]), cov_target_stk = cov(res[, "target"], res[, "se_stk"]),
    mcse_cov = sqrt(.95 * .05 / R), mean_size = mean(res[, "size"]), p_superset = mean(res[, "superset"]))
  cat(cn, "done\n")
}
tab <- do.call(rbind, out)
write.csv(tab, "output/ry_target_check.csv", row.names = FALSE)
capture.output(print(format(tab, digits = 3), row.names = FALSE), file = "output/ry_target_check.txt")
print(format(tab, digits = 3), row.names = FALSE)
