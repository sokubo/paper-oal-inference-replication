# undercoverage_demo.R — Paper 3, Section 2: how naive inference after the
# outcome-adaptive lasso (Shortreed & Ertefaie 2017, Biometrics) behaves.
#
# Two regimes, same pipelines:
#  Regime I  ("benign", S-E-like): n=500, p=20, strong signals, wAMD tuning.
#  Regime II ("many candidate controls"): n=200, p=50, two boundary
#            confounders (alpha=.7, beta=.25) whose selection is a coin
#            flip under BIC tuning; wAMD variant also reported.
# Estimators: oracle OLS; full OLS; OAL->Hajek IPW with fixed-weight
# sandwich SE; OAL-selected -> refit OLS with textbook SE.
# R = 500 replications per regime. Message: calibration is erratic across
# regimes/tuning in ways the analyst cannot see -- valid post-selection
# inference, not luck, is required.

suppressMessages(library(glmnet))
dir.create("output", showWarnings = FALSE)
TAU <- 1; R <- 500

run_regime <- function(n, p, a3, b3, sd_y, tune, seed, label) {
  set.seed(seed)
  alpha <- c(.8,.8,a3,a3,0,0,1,1, rep(0, p-8))
  beta  <- c(.6,.6,b3,b3,.6,.6,0,0, rep(0, p-8))
  lam_grid <- sort(n^c(-10,-5,-1,-.75,-.5,-.25,.25,.49), decreasing = TRUE)
  ncols <- 8; res <- matrix(NA_real_, R, 8 + 2)
  for (r in 1:R) {
    Z <- matrix(rnorm(n*p), n, p)
    A <- rbinom(n, 1, plogis(Z %*% alpha))
    Y <- TAU*A + as.vector(Z %*% beta) + rnorm(n, sd = sd_y)
    # oracle / full
    mo <- lm(Y ~ A + Z[, which(beta != 0)])
    mf <- lm(Y ~ A + Z)
    co <- summary(mo)$coefficients["A", 1:2]
    cf <- summary(mf)$coefficients["A", 1:2]
    # OAL
    b_out <- coef(mf)[-(1:2)]
    w <- pmin(abs(b_out)^(-2), 1e4)
    fit <- suppressWarnings(glmnet(Z, A, family = "binomial",
                                   penalty.factor = w, lambda = lam_grid))
    crit <- rep(Inf, length(lam_grid))
    ee <- vector("list", length(lam_grid))
    for (k in seq_along(lam_grid)) {
      e <- tryCatch(as.vector(predict(fit, Z, s = lam_grid[k],
                                      type = "response")),
                    error = function(x) NULL)
      if (is.null(e)) next
      e <- pmin(pmax(e, .01), .99); ee[[k]] <- e
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
    k <- which.min(crit); e <- ee[[k]]
    sel <- which(as.vector(coef(fit, s = lam_grid[k]))[-1] != 0)
    # OAL -> Hajek IPW, fixed-weight sandwich
    ipw <- A/e + (1-A)/(1-e); X <- cbind(1, A)
    br <- solve(crossprod(X * sqrt(ipw)))
    bh <- br %*% crossprod(X * ipw, Y)
    rr <- Y - X %*% bh
    V  <- br %*% crossprod(X * (ipw * as.vector(rr))) %*% br
    # OAL -> refit OLS
    m2 <- if (length(sel)) lm(Y ~ A + Z[, sel, drop = FALSE]) else lm(Y ~ A)
    c2 <- summary(m2)$coefficients["A", 1:2]
    res[r, ] <- c(co, cf, bh[2], sqrt(V[2,2]), c2,
                  3 %in% sel, 4 %in% sel)
  }
  meth <- c("oracle","full","oal_ipw","oal_ols")
  tab <- do.call(rbind, lapply(1:4, function(i) {
    est <- res[, (i-1)*2 + 1]; se <- res[, (i-1)*2 + 2]
    data.frame(regime = label, tune = tune, method = meth[i],
               bias = mean(est) - TAU, emp_sd = sd(est),
               se_ratio = mean(se)/sd(est),
               coverage95 = mean(abs(est - TAU) <= 1.96*se))
  }))
  tab$z3_rate <- mean(res[, 9]); tab$z4_rate <- mean(res[, 10])
  tab
}

t0 <- Sys.time()
tabs <- rbind(
  run_regime(500, 20, .8, .6, 2, "wamd", 101, "I: benign (n=500, p=20)"),
  run_regime(200, 50, .7, .25, 2, "bic",  202, "II: many controls (n=200, p=50)"),
  run_regime(200, 50, .7, .25, 2, "wamd", 202, "II: many controls (n=200, p=50)"))
cat("time:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
print(tabs, digits = 3, row.names = FALSE)
write.csv(tabs, "output/undercoverage_table.csv", row.names = FALSE)
capture.output(print(tabs, digits = 3), file = "output/undercoverage_summary.txt")

# ---- Figure: coverage by regime ------------------------------------------
suppressMessages(library(ggplot2))
sub <- tabs[(tabs$regime == "I: benign (n=500, p=20)" & tabs$tune == "wamd") |
            (tabs$regime == "II: many controls (n=200, p=50)" & tabs$tune == "bic"), ]
lab <- c(oracle = "Oracle OLS", full = "Full OLS",
         oal_ipw = "OAL + IPW, naive SE", oal_ols = "OAL + refit OLS, naive SE")
sub$mlab <- factor(lab[sub$method], levels = rev(lab))
fig <- ggplot(sub, aes(x = coverage95, y = mlab)) +
  geom_vline(xintercept = .95, color = "#0b0b0b", linetype = "dashed",
             linewidth = .35) +
  geom_col(width = .55, fill = "#2a78d6") +
  geom_text(aes(label = sprintf("%.0f%%", 100*coverage95)), hjust = -.15,
            color = "#0b0b0b", size = 3.2, fontface = "bold") +
  facet_wrap(~regime) +
  scale_x_continuous(limits = c(0, 1.09), expand = c(0, 0),
                     labels = function(v) sprintf("%.0f%%", 100*v)) +
  labs(title = "The same pipeline, two calibrations",
       subtitle = paste("95% CI coverage of the ATE over 500 replications; dashed line = nominal.",
                        "Regime II: two boundary confounders selected in only 49% / 47% of replications (BIC tuning).",
                        sep = "\n"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "#e1e0d9", linewidth = .3),
        plot.background = element_rect(fill = "#fcfcfb", color = NA),
        text = element_text(color = "#0b0b0b"),
        strip.text = element_text(color = "#52514e", face = "bold"),
        plot.subtitle = element_text(color = "#52514e", size = 8.5))
ggsave("output/fig1_undercoverage.png", fig, width = 7.5, height = 3.4,
       dpi = 200, bg = "#fcfcfb")
cat("Demo complete.\n")
