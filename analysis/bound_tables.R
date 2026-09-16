# bound_tables.R -- Paper 3, Web Appendix B: population (quadrature) checks.
#
# Table 1 (bias): the true bias of the reduced pseudo-true AIPW target,
#   tau(e*_S, m*_S) - tau, when a weak confounder Z2 is omitted (S = {1},
#   O = {2}), against (a) the leading Stein term of Proposition A1'
#   (Gaussian independent design only) and (b) the Holder bound of Lemma A3,
#   for five two-covariate designs: Gaussian independent, correlated Gaussian
#   (rho = .5), bounded Rademacher, bounded uniform, and a nonlinear outcome
#   model (Hermite cubic in Z1).  Heteroscedasticity does not change the
#   pseudo-true projections (E[eps | A, Z] = 0), so it enters Table 2 only.
# Table 2 (r_mis): the misspecification term of Proposition 1(b),
#   r_mis = Var(phi*) - V_eff(S*), when a pure instrument Z2 is omitted
#   (beta_2 = 0, S* = {1}), together with the uncorrected asymptotic variance
#   Var(psi(e*_S, m)) (the "empirical" variance estimand), the corrected
#   variance Var(phi*) (Lemma A4), and V_eff(S* u I) for the same designs plus a
#   heteroscedastic one (sd depends on Z1).
# Design: e(z) = expit(a0 + a1 z1 + a2 z2), m_a(z) = tau a + b1 z1 + b2 z2 [+ c h(z1)],
#   eps | A, Z ~ (0, sigma^2(z)).  Quadrature: Gauss-Hermite 80 x 80 (Gaussian),
#   Gauss-Legendre 80 x 80 (uniform), exact 2 x 2 sums (Rademacher).
# Outputs: output/bound_tables_bias.csv, output/bound_tables_rmis.csv, output/bound_tables.txt
suppressMessages(library(statmod))
TAU <- 1; NQ <- 80

grid2 <- function(design, rho = .5) {
  if (design %in% c("gauss", "ar1", "hetero", "nonlin")) {
    q <- gauss.quad.prob(NQ, "normal"); g <- expand.grid(u1 = q$nodes, u2 = q$nodes)
    w <- as.vector(outer(q$weights, q$weights))
    if (design == "ar1") { z1 <- g$u1; z2 <- rho * g$u1 + sqrt(1 - rho^2) * g$u2 } else { z1 <- g$u1; z2 <- g$u2 }
  } else if (design == "unif") {
    q <- gauss.quad.prob(NQ, "uniform", l = -sqrt(3), u = sqrt(3)); g <- expand.grid(u1 = q$nodes, u2 = q$nodes)
    w <- as.vector(outer(q$weights, q$weights)); z1 <- g$u1; z2 <- g$u2
  } else if (design == "rademacher") {
    g <- expand.grid(u1 = c(-1, 1), u2 = c(-1, 1)); w <- rep(.25, 4); z1 <- g$u1; z2 <- g$u2
  }
  list(z1 = z1, z2 = z2, w = w / sum(w))
}
population <- function(design, a0, a1, a2, b1, b2, cnl = 0, het = 0, sig2 = 1) {
  G <- grid2(design); z1 <- G$z1; z2 <- G$z2; w <- G$w
  eta <- a0 + a1 * z1 + a2 * z2; e <- plogis(eta); e0 <- plogis(-eta)   # e0 = 1 - e without cancellation
  h <- if (design == "nonlin") cnl * (z1^3 - 3 * z1) / sqrt(6) else 0
  m1 <- TAU + b1 * z1 + b2 * z2 + h; m0 <- b1 * z1 + b2 * z2 + h
  s2 <- sig2 * exp(het * z1 - het^2 / 2)      # E[s2] = sig2 under N(0,1) z1
  list(z1 = z1, z2 = z2, w = w, e = e, e0 = e0, m1 = m1, m0 = m0, s2 = s2, x = cbind(1, z1))
}
E <- function(f, w) sum(f * w)
reduced_logit <- function(P) {                    # population logistic fit of A on (1, z1)
  g <- c(qlogis(E(P$e, P$w)), 0)
  for (it in 1:50) {
    pi <- as.vector(plogis(P$x %*% g)); sc <- crossprod(P$x, P$w * (P$e - pi))
    H <- crossprod(P$x * sqrt(P$w * pi * (1 - pi))); step <- as.vector(solve(H, sc)); g <- g + step
    if (max(abs(step)) < 1e-13) break
  }
  eta <- as.vector(P$x %*% g); list(es = plogis(eta), es0 = plogis(-eta))
}
arm_projection <- function(P, a) {                # population LS of Y on (1, z1) under P(. | A = a)
  wa <- if (a == 1) P$e else P$e0; ma <- if (a == 1) P$m1 else P$m0
  M <- crossprod(P$x * sqrt(P$w * wa)); v <- crossprod(P$x, P$w * wa * ma)
  as.vector(P$x %*% solve(M, v))
}
analyse_bias <- function(design, a0, a1, a2, b1, b2, cnl = 0) {
  P <- population(design, a0, a1, a2, b1, b2, cnl = cnl); w <- P$w; e <- P$e; e0 <- P$e0
  rl <- reduced_logit(P); es <- rl$es; es0 <- rl$es0
  ms1 <- arm_projection(P, 1); ms0 <- arm_projection(P, 0)
  bias <- E(ms1 - ms0 + e * (P$m1 - ms1) / es - e0 * (P$m0 - ms0) / es0, w) - TAU
  # Lemma A3 identity (check) and Holder bound with exponents (2, 4, 4)
  ident <- E((e - es) * ((P$m1 - ms1) / es + (P$m0 - ms0) / es0), w)
  n2 <- function(f) sqrt(E(f^2, w)); n4 <- function(f) E(f^4, w)^(1 / 4)
  holder <- n2(e - es) * (n4(P$m1 - ms1) * n4(1 / es) + n4(P$m0 - ms0) * n4(1 / es0))
  # Proposition A1' decomposition (Stein leading term is exact only for Gaussian independent z)
  hS <- 1 / es + 1 / es0
  T1 <- E((e - es) * hS * b2 * P$z2, w)
  cPS <- E(hS * e * e0, w)
  q1 <- (P$m1 - ms1) - b2 * P$z2; q0 <- (P$m0 - ms0) - b2 * P$z2     # -q_a(z_S) (+ nonlinear residual)
  T2 <- -E((e - es) * (q1 / es + q0 / es0), w)
  data.frame(design = design, a0 = a0, a1 = a1, a2 = a2, b1 = b1, b2 = b2, c_nl = cnl,
             bias = bias, identity_check = ident - bias, a2b2 = a2 * b2, bias_over_a2b2 = bias / (a2 * b2),
             stein_T1 = T1, stein_c = cPS, a2b2_c = a2 * b2 * cPS, T2 = T2,
             holder_bound = holder, bias_over_holder = abs(bias) / holder,
             max_e_minus_estar = max(abs(e - es)))
}
analyse_rmis <- function(design, a0, a1, a2, b1, het = 0, sig2 = 1) {
  P <- population(design, a0, a1, a2, b1, 0, het = het, sig2 = sig2); w <- P$w; e <- P$e; e0 <- P$e0; x <- P$x; s2 <- P$s2
  rl <- reduced_logit(P); es <- rl$es; es0 <- rl$es0
  # coarsened propensity e_S = E[e | z1] and arm-specific conditional variances given z1
  key <- round(P$z1, 10); eS <- ave(e * w, key, FUN = sum) / ave(w, key, FUN = sum)
  s2_1 <- ave(e * s2 * w, key, FUN = sum) / ave(e * w, key, FUN = sum)
  s2_0 <- ave(e0 * s2 * w, key, FUN = sum) / ave(e0 * w, key, FUN = sum)
  eS0 <- ave(e0 * w, key, FUN = sum) / ave(w, key, FUN = sum)
  Veff_S <- E(s2_1 / eS + s2_0 / eS0, w)                    # Var(m1 - m0 | z_S) = 0 (constant effect)
  Veff_SI <- E(s2 / e + s2 / e0, w)
  Vnaive <- E(e * s2 / es^2 + e0 * s2 / es0^2, w)            # Var psi(e*_S, m)
  D1 <- crossprod(x, w * (1 - e / es)); D0 <- -crossprod(x, w * (1 - e0 / es0))
  J1 <- crossprod(x * sqrt(w * e)); J0 <- crossprod(x * sqrt(w * e0))
  v1 <- as.vector(x %*% solve(J1, D1)); v0 <- as.vector(x %*% solve(J0, D0))
  Vcorr <- E(e * s2 * (1 / es + v1)^2 + e0 * s2 * (-1 / es0 + v0)^2, w)
  data.frame(design = design, a0 = a0, a1 = a1, a2 = a2, b1 = b1, het = het,
             Veff_S = Veff_S, Veff_SI = Veff_SI, ratio_SI_over_S = Veff_SI / Veff_S,
             Var_psi_naive = Vnaive, Var_phi_corrected = Vcorr,
             r_mis = Vcorr - Veff_S, r_mis_pct = 100 * (Vcorr - Veff_S) / Veff_S,
             sd_ratio_naive_over_corrected = sqrt(Vnaive / Vcorr),
             D_b1 = paste(signif(D1, 3), collapse = ";"), max_eS_minus_estar = max(abs(eS - es)))
}

dir.create("output", showWarnings = FALSE)
coefs <- rbind(c(.5, 1, 1, 1, .3), c(0, 1, 2, 1, .3), c(.5, 1, .5, 1, .3), c(1, 1.5, 1.5, 1, .3))
colnames(coefs) <- c("a0", "a1", "a2", "b1", "b2")
tab1 <- do.call(rbind, lapply(c("gauss", "ar1", "rademacher", "unif", "nonlin"), function(d)
  do.call(rbind, lapply(seq_len(nrow(coefs)), function(i) {
    cf <- coefs[i, ]; analyse_bias(d, cf["a0"], cf["a1"], cf["a2"], cf["b1"], cf["b2"],
                                   cnl = if (d == "nonlin") 0.6 else 0) }))))
# nonlinear outcome model with Z2 a pure instrument (b2 = 0): both reduced models misspecified
tab1 <- rbind(tab1, do.call(rbind, lapply(seq_len(nrow(coefs)), function(i) {
  cf <- coefs[i, ]; analyse_bias("nonlin", cf["a0"], cf["a1"], cf["a2"], cf["b1"], 0, cnl = 0.6) })))
tab1$bias_over_a2b2[tab1$a2b2 == 0] <- NA
tab2 <- rbind(
  analyse_rmis("gauss", .5, 1, 1, 1), analyse_rmis("gauss", 0, 1, 2, 1), analyse_rmis("gauss", .5, 1, .5, 1),
  analyse_rmis("gauss", 1, 1.5, 1.5, 1),
  analyse_rmis("ar1", .5, 1, 1, 1), analyse_rmis("ar1", 0, 1, 2, 1),
  analyse_rmis("rademacher", .5, 1, 1, 1), analyse_rmis("rademacher", 0, 1, 2, 1),
  analyse_rmis("unif", .5, 1, 1, 1), analyse_rmis("unif", 0, 1, 2, 1),
  analyse_rmis("hetero", .5, 1, 1, 1, het = .5), analyse_rmis("hetero", 0, 1, 2, 1, het = .5),
  analyse_rmis("hetero", 0, 1, 2, 1, het = 1))
write.csv(rbind(cbind(table = "bias", tab1[, setdiff(names(tab1), "")], stringsAsFactors = FALSE)),
          "output/bound_tables_bias.csv", row.names = FALSE)
write.csv(tab2, "output/bound_tables_rmis.csv", row.names = FALSE)
sink("output/bound_tables.txt")
cat("bound_tables.R --", format(Sys.time(), "%Y-%m-%d %H:%M"), " quadrature nodes per dimension:", NQ, "\n\n")
cat("Table 1: bias of the reduced pseudo-true AIPW target when Z2 (b2 != 0) is omitted, S = {1}\n")
cat("  bias_over_a2b2: bias / (a2 b2); stein_T1: E[(e - e*) h z2] b2 (equals a2 b2 c_PS under Gaussian independence);\n")
cat("  T2 = E[(e - e*){q1/e* + q0/(1 - e*)}] (arm-projection remainder; bias = T1 - T2 exactly, all designs);\n")
cat("  holder_bound: ||e - e*||_2 sum_a ||m_a - m*_a||_4 ||1/e*_a||_4 (Lemma A3).\n\n")
print(format(tab1[, c("design", "a0", "a1", "a2", "b1", "b2", "c_nl", "bias", "a2b2", "bias_over_a2b2", "stein_T1", "T2",
                      "holder_bound", "bias_over_holder", "identity_check")], digits = 4), row.names = FALSE)
cat("\nTable 2: r_mis = Var(phi*) - V_eff(S*) when a pure instrument Z2 is omitted (b2 = 0, S* = {1})\n")
cat("  Var_psi_naive = Var psi(e*_S, m) (estimand of the empirical variance); Var_phi_corrected = Var(phi*) (Lemma A4).\n\n")
print(format(tab2[, c("design", "a0", "a1", "a2", "b1", "het", "Veff_S", "Veff_SI", "ratio_SI_over_S", "Var_psi_naive",
                      "Var_phi_corrected", "r_mis", "r_mis_pct", "sd_ratio_naive_over_corrected", "max_eS_minus_estar")],
             digits = 4), row.names = FALSE)
sink()
cat(readLines("output/bound_tables.txt"), sep = "\n")
