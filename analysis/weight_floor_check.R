# weight_floor_check.R -- Paper 3, review round 4 (R4-m3): the adaptive-weight floor.
#
# The OAL arm rescales covariate j by max(b_j^2, eps_n) with eps_n = min(1e-8, 1/n)
# (n = training size), i.e. it uses the stabilized adaptive weight {max(b_j^2, eps_n)}^{-1}.
# Earlier code used the fixed floor 1e-8. This script records three facts:
#   (1) eps_n == 1e-8 exactly (as a double) at every training size used in the paper,
#       so the shipped fits are unchanged by the definition;
#   (2) the rescaling and the OAL selection are identical() under the two floors on
#       simulated data from the battery's DGP, and the rescaling on an RHC pilot (numeric
#       columns); this is implied by (1) and is recorded as a direct check, not as evidence
#       about the shipped fits beyond (1);
#   (3) the rate bookkeeping of Lemma A2(ii): on the event |b_j| <= log(n)/sqrt(n), the
#       stabilized denominator is <= log(n)^2/n (because eps_n <= 1/n <= log(n)^2/n for
#       n >= 3), so lambda_n * w_j >= n^{5/4}/log(n)^2; and the contrast with a floor fixed
#       at 1e-8 for all n, whose penalty n^{1/4} * 1e8 falls below n once n^{3/4} > 1e8.
# Output: output/weight_floor_check.txt. Run from analysis/: Rscript weight_floor_check.R
suppressPackageStartupMessages(library(glmnet))
src <- readLines("sim_battery4.R"); cut <- grep("^# -+ run$", src); eval(parse(text = src[1:(cut - 1)]))
old_scale <- function(b, n) pmax(abs(b)^2, 1e-8)     # the pre-round-4 rescaling
sel_old <- function(Z, A, b) {                        # the pre-round-4 sel_oal_fixed, verbatim except for the floor
  n <- nrow(Z); Zs <- sweep(Z, 2, old_scale(b, n), "*"); lam <- n^(-0.75)
  fit <- suppressWarnings(glmnet(Zs, A, family = "binomial", standardize = FALSE, lambda = c(2 * lam, lam)))
  which(as.vector(coef(fit, s = lam))[-1] != 0)
}
out <- character(0); say <- function(...) out <<- c(out, sprintf(...))

# (1) eps_n at the training sizes used: battery n in {500, 2000, 8000} with K = 5 (n_tr = 0.8 n),
#     RHC n = 5735 (n_tr = 4588), Section 2 prototypes (n up to 8000).
n_used <- c(0.8 * c(500, 2000, 8000), 4588, 500, 2000, 8000, 5735)
say("(1) eps_n = min(1e-8, 1/n) at the training sizes used:")
for (n in n_used) say("    n = %5d: eps_n = %.17g; identical(eps_n, 1e-8) = %s", n, oal_floor(n), identical(oal_floor(n), 1e-8))
stopifnot(all(sapply(n_used, function(n) identical(oal_floor(n), 1e-8))))

# (2) identical rescaling and selection under the two floors
set.seed(20260917)
same_sel <- TRUE; same_scl <- TRUE
for (cell in list(list(n = 500, p = 20, a3 = .5, b3 = .25), list(n = 2000, p = 20, a3 = .5, b3 = .25),
                  list(n = 8000, p = 20, a3 = .5, b3 = .25), list(n = 500, p = 50, a3 = .5, b3 = .25))) {
  for (r in 1:20) {
    d <- lin_dgp(cell$n, cell$p, cell$a3, cell$b3)()
    ntr <- floor(0.8 * cell$n); Z <- d$Z[1:ntr, ]; A <- d$A[1:ntr]; Y <- d$Y[1:ntr]
    b <- coef(lm(Y ~ A + Z))[-(1:2)]
    same_scl <- same_scl && identical(oal_scale(b, ntr), old_scale(b, ntr))
    Zs_new <- sweep(Z, 2, oal_scale(b, ntr), "*"); Zs_old <- sweep(Z, 2, old_scale(b, ntr), "*")
    same_sel <- same_sel && identical(Zs_new, Zs_old) && identical(sel_oal_fixed(Z, A, b), sel_old(Z, A, b))
  }
}
# an exactly zero pilot coefficient (the case the floor exists for)
b0 <- c(0, 1e-5, 0.3); say("    zero coefficient: oal_scale(c(0, 1e-5, 0.3), 400) = %s (old: %s)",
                            paste(format(oal_scale(b0, 400), digits = 3), collapse = " "),
                            paste(format(old_scale(b0, 400), digits = 3), collapse = " "))
same_scl <- same_scl && identical(oal_scale(b0, 400), old_scale(b0, 400))
say("(2) rescaled design matrices and OAL selections identical under both floors (80 simulated training sets): %s; rescaling vectors identical: %s", same_sel, same_scl)
stopifnot(same_sel, same_scl)

# RHC pilot coefficients, both outcomes, on the whole sample (the fold pilots use subsets of it)
rhc_path <- Filter(file.exists, c("../data/rhc_full.rda", "../../../replication/paper-oal-inference-replication/data/rhc_full.rda"))
if (length(rhc_path)) {
  load(rhc_path[1]); f <- rhc
  Y_los <- as.numeric(f$dschdte - f$sadmdte); A <- as.integer(f$swang1 == "RHC")
  num <- sapply(f, is.numeric); Z <- scale(as.matrix(f[, num & !(names(f) %in% c("ptid", "sadmdte", "dschdte", "dthdte", "lstctdte", "t3d30", "adld3p"))]))
  Z[is.na(Z)] <- 0; ok <- !is.na(Y_los)
  b <- coef(lm(Y_los[ok] ~ A[ok] + Z[ok, ]))[-(1:2)]; b[is.na(b)] <- 0
  say("    RHC (numeric columns, LOS pilot, n = %d): min |b_j| = %.3e; rescaling identical under both floors: %s",
      sum(ok), min(abs(b)), identical(oal_scale(b, sum(ok)), old_scale(b, sum(ok))))
}

# (3) rate bookkeeping
say("(3) Lemma A2(ii): max(b^2, eps_n) <= log(n)^2/n on the event |b| <= log(n)/sqrt(n), since eps_n <= 1/n <= log(n)^2/n for n >= 3:")
for (n in c(3, 400, 1600, 6400, 1e6, 1e8, 1e10, 1e12)) {
  eps <- oal_floor(n); bound <- log(n)^2 / n
  say("    n = %.0e: eps_n = %.1e <= log(n)^2/n = %.3e: %s; shrinking-floor penalty lower bound n^{5/4}/log(n)^2 = %.3e; fixed-floor (1e-8) penalty at an exact zero, n^{1/4}*1e8 = %.3e, versus the gradient scale n = %.0e -> %s",
      n, eps, bound, eps <= bound, n^1.25 / log(n)^2, n^0.25 * 1e8, n, if (n^0.25 * 1e8 > n) "exclusion argument holds" else "fixed floor too weak (n^{3/4} > 1e8)")
}
say("    the fixed floor 1e-8 stops dominating the gradient at n = 1e8^(4/3) = %.2e", 1e8^(4 / 3))
writeLines(out, "output/weight_floor_check.txt"); cat(out, sep = "\n")
