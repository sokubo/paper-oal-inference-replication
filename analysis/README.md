# analysis/ — scripts, run order and number map (paper 3)

Paper: *Inference after Thresholded Outcome-Based Confounder Selection under Confounder–Instrument
Separation* (v1.0, 2026-09-17). All scripts are run from this directory and write to `output/`.

## Environment of the shipped outputs
R 4.3.3 (Linux container, 2 cores) for the quadrature tables, the numerical diagnostics and Example
A1; the author's Mac (macOS, 12 cores) for the full-size battery with the C-TMLE arm and the RHC
application (`RUN_ON_MAC.md` records the runs; the run logs do not record the R version, which the
release check's `sessionInfo()` will); packages glmnet, MASS, statmod, ctmle (C-TMLE arm only), parallel,
ggplot2 (`undercoverage_demo.R` only). `hdm` is needed only for the optional `apply_401k.R`. Python 3
with `pandas` for the table conversion `make_paper_tables4.py`; `example_A1.py` uses the standard
library only (exact rational arithmetic). Check before running:
```sh
Rscript -e 'for (p in c("glmnet","MASS","statmod","ggplot2","ctmle")) cat(p, if (requireNamespace(p, quietly=TRUE)) "ok" else "MISSING", "\n")'
python3 -c 'import pandas; print("pandas", pandas.__version__)'
```
An independent rerun of the whole battery (500/300/500 replications, 8 cores) took 101.5 minutes.
The shipped `convergence_check*.csv` are from the author's Mac run of 2026-09-16 (R 4.6.0, macOS,
12 cores), which the release check reproduced exactly; an earlier Linux run of the same script
(kept in `output_container_2026-09-16/`, not shipped) agrees in every diagnostic count except that
one P2 refit reached the iteration limit on macOS but not on Linux, and the changes under the
200-iteration cap differ (the non-converged refits are divergent logistic fits), so Table 11
reports the archived macOS values.

## Run order
| Step | Command | Produces | Used in |
|---|---|---|---|
| 1 | `Rscript sim_battery4.R 500 300 500 <cores>` | `output/battery4_table.csv`, `battery4_selection.csv`, `battery4_summary.txt`, `battery4_run*.log`, `ctmle_errors.txt` | Table 2, Web Appendix C (Tables 8–10), Web Appendix D coverages |
| 2 | `Rscript make_tables4.R`; `python3 make_paper_tables4.py output` | `output/battery4_tables.md`; the manuscript tables printed to stdout | Web Appendix C |
| 3 | `Rscript apply_realdata2.R` | `output/applications2_table.csv/.txt`, `applications2_selection.csv`, `applications2_splits.csv`, `applications2_trim01.csv`, `applications2_pool.txt`, `applications2_overlap.txt` | Section 5 (Tables 3–4) |
| 4 | `Rscript bound_tables.R` | `output/bound_tables_bias.csv`, `bound_tables_rmis.csv`, `bound_tables.txt` | Web Appendix B (Tables 5–6) |
| 5 | `Rscript undercoverage_demo.R`; `Rscript proposed_estimator.R` | `output/undercoverage_*`, `proposed_*` | Section 2 (Figure 1, prototype numbers) |
| 6 | `Rscript ry_target_check.R 100 <cores>` | `output/ry_target_check.csv/.txt/.log` | Web Appendix D (100 replications; the second argument is the number of cores and does not affect the results; the script's default is 300 replications) |
| 7 | `Rscript convergence_check.R 25 10 <cores>`; `Rscript convergence_table.R` | `output/convergence_check.csv/.txt`, `convergence_check_fits.csv`; `output/convergence_table.md` (the Table 11 rows and the five prose figures, each tied to its source column) | Web Appendix C (Table 11) |
| 8 | `Rscript weight_floor_check.R` | `output/weight_floor_check.txt` | Section 3.1 / Lemma A2 / Web Appendix C: the adaptive-weight floor $\epsilon_n = \min(10^{-8}, 1/n)$ equals `1e-8` at every training size used, and the rescaling and OAL selection are `identical()` under the shrinking and the fixed floor (round 4, R4-m3) |
| 9 | `python3 example_A1.py` | printed | Web Appendix B, Example A1 |
| — | `Rscript apply_401k.R` (optional; needs `hdm`) | `output/applications_401k.csv/.txt` | not used in the paper |

Per-replication seeds are `1000 * seed_cell + r`, so the battery does not depend on the number of
cores; per-cell checkpoints in `output/battery4_cells/` let an interrupted run resume (delete them
to force a fresh run). `sim_battery.R`, `sim_battery2.R`, `sim_battery3.R` are earlier versions kept
for provenance and produce nothing used in the paper.

## Labels
The scripts and CSVs predate two renamings: `oads_thr` / `oads_bic` (printed "ODS-thr" / "ODS-bic" in
`battery4_tables.md`) are the arms called **OBS-thr** / **OBS-bic** in the paper, and `goal_ipw`
("GOAL (elastic net) + Hajek IPW") is the pipeline called **OAL-EN** in the paper (an adaptive-weighted
elastic-net variant of naive OAL; not an implementation of the GOAL of Baldé, Yang and Lefebvre 2023).
The battery was not rerun for the renamings, so the checked outputs are unchanged.

## Failure policy and diagnostics
A replication counts as a failure for an arm only when its estimate or standard error is nonfinite
(`n_fail` in `battery4_table.csv`; nonzero only for the C-TMLE arm: 8 in II, 10 in P2). Logistic refits
that stop at `glm.fit`'s iteration limit, aliased coefficients (set to zero) and singular sandwich blocks
(solved with `MASS::ginv`) are not failures under this definition; `convergence_check.R` reports them
arm by arm for a fixed subset of replications (Table 11).

## Data
`apply_realdata2.R` reads the public SUPPORT/RHC file `rhc_full.rda` (Connors et al. 1996; n = 5,735),
looked up under `data/`, `../data/` and the author's project tree (`DATA_PATH` at the top of the
script). No other data file is needed. The 401(k) analysis of earlier versions was withdrawn.
