# convergence_table.R -- Table 11 (numerical diagnostics) read from the diagnostic output.
#
# Reads output/convergence_check.csv and output/convergence_check_fits.csv (written by
# convergence_check.R 25 10 <cores>) and writes output/convergence_table.md: the rows of
# Table 11 and the five prose figures of the "Numerical diagnostics" paragraph of Web
# Appendix C, each tied to its source column. The manuscript's table is copied from this
# file; running this script after any rerun of convergence_check.R shows at once whether
# the manuscript entries are stale. No simulation is run here.
#
# Usage: Rscript convergence_table.R      (from analysis/)

summ <- read.csv("output/convergence_check.csv", comment.char = "#", stringsAsFactors = FALSE)
fits <- read.csv("output/convergence_check_fits.csv", stringsAsFactors = FALSE)

stopifnot(nrow(fits) == sum(summ$fits))
ev <- summ[summ$ps_nonconv > 0 | summ$ginv_fits > 0, ]
ev <- ev[order(ev$cell), ]
other <- summ[!(paste(summ$cell, summ$arm) %in% paste(ev$cell, ev$arm)), ]
stopifnot(all(other$ps_nonconv == 0), all(other$ginv_fits == 0),
          all(summ$ps_boundary == 0), all(summ$ps_rankdef_fits == 0), all(summ$out_rankdef_fits == 0))

# replications affected (from the per-fit records: distinct seeds with a non-converged refit)
reps_affected <- sapply(seq_len(nrow(ev)), function(i) {
  f <- fits[fits$cell == ev$cell[i] & fits$arm == ev$arm[i] & fits$ps_nonconv > 0, ]
  length(unique(f$seed))
})
ev$reps_affected <- reps_affected

# arm labels as printed in the manuscript's Table 11 (the CSV labels carry "+ AIPW" and "(wAMD)")
paper_label <- function(l) sub(" \\(wAMD\\)", "", sub("^full set \\+ AIPW$", "full set", l))
fmt <- function(x, d) formatC(x, format = "f", digits = d)
row_md <- function(r) {
  sprintf("| %s | %s | %d | %d | %s | %d | %d | %s | %s | %s | %d |",
          r$cell, paper_label(r$arm_label), r$reps, r$fits,
          if (abs(r$mean_dim - round(r$mean_dim)) < 1e-9) as.character(round(r$mean_dim)) else fmt(r$mean_dim, 1),
          r$ps_nonconv, r$ginv_fits, fmt(100 * r$trim_frac, 1),
          fmt(r$max_d_est_maxit, 3), fmt(r$max_d_se_stk_maxit, 3), r$cov_flips_maxit)
}
hdr <- c("| Cell | Arm | Reps | Refits | Dim. | Iter. limit | Gen. inverse | Clipped (pct) | Max Δ estimate | Max Δ stacked SE | Coverage changes |",
         "|---|---|---|---|---|---|---|---|---|---|---|")
rows <- vapply(seq_len(nrow(ev)), function(i) row_md(ev[i, ]), character(1))
other_row <- sprintf("| other | all other %d cell–arm pairs | %s | %s | %s–%s | 0 | 0 | — | — | — | — |",
                     nrow(other), paste(sort(unique(other$reps), decreasing = TRUE), collapse = " or "),
                     format(sum(other$fits), big.mark = ","),
                     fmt(min(other$mean_dim), 1), as.character(round(max(other$mean_dim))))

prose <- c(
  sprintf("- total refits examined: %s (column `fits`, summed)", format(sum(summ$fits), big.mark = ",")),
  sprintf("- refits that stopped at the iteration limit: %d (column `ps_nonconv`, summed)", sum(summ$ps_nonconv)),
  sprintf("- generalized-inverse fallbacks: %d (column `ginv_fits`, summed)", sum(summ$ginv_fits)),
  sprintf("- %s / %s: %d non-converged refits in %d of %d replications, %d generalized inverse(s), %s percent clipped; max |Δ estimate| %s, max |Δ stacked SE| %s, coverage changes %d",
          ev$cell, paper_label(ev$arm_label), ev$ps_nonconv, ev$reps_affected, ev$reps, ev$ginv_fits, fmt(100 * ev$trim_frac, 0),
          fmt(ev$max_d_est_maxit, 3), fmt(ev$max_d_se_stk_maxit, 3), ev$cov_flips_maxit))
hdr_line <- readLines("output/convergence_check.csv", n = 1)
out <- c("# Table 11 entries generated from output/convergence_check.csv", "",
         paste0("Source run: ", sub("^# ", "", hdr_line)), "", hdr, rows, other_row, "",
         "Prose figures (Web Appendix C, numerical diagnostics paragraph):", prose, "")
writeLines(out, "output/convergence_table.md")
cat(out, sep = "\n")
