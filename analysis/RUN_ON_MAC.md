# Running the Phase B computations at full size (author's machine)

**Status 2026-09-07, second run.** The re-run with the fixed C-TMLE arm finished (84.7 min; every other arm reproduced the first run exactly; C-TMLE failed in 8 of 500 fits in II and 10 of 300 in P2, none elsewhere). Its rows are in Table 1 and Web Appendix C. Nothing further to run for the battery; `analysis/output/` on the Mac is the source of record and the container copy matches it.

**Status 2026-09-07, first run.** The full-size run (500/300/500, 12 cores, 8.9 min, no failed
replication) is done and transcribed into the paper (Table 1, Web Appendix C, the six
readings, Remark A5, Web Appendix D). The main cells reproduced the container run to
machine precision (II and P1 to 0.005); P2 and the boundary sequence changed at the
replication-count level only. **The C-TMLE arm failed in every replication** (all
`ctmle` rows NA, `n_fail` = 500): `ctmleGlmnet` was called without a `lambdas`
sequence. `ctmle_arm()` now builds the sequence as in the package vignette (from the
`cv.glmnet` minimiser towards less penalisation), rescales Y and Q to [0, 1] for the
logistic fluctuation, and writes the first error messages to `output/ctmle_errors.txt`.
To fill the C-TMLE rows, delete the checkpoints (they would be reused as they are) and
re-run; every other arm reproduces exactly because the seeds are per replication:

```sh
cd ~/Documents/Claude/Projects/VariableSelection/paper/03-oal-inference/analysis
mv output/battery4_cells output/battery4_cells_noctmle_2026-09-07      # keep the first run's checkpoints
Rscript sim_battery4.R 20 10 10 4 I                                    # smoke test; then check:
cat output/ctmle_errors.txt 2>/dev/null; grep -c ctmle output/battery4_table.csv
Rscript sim_battery4.R 500 300 500 $(sysctl -n hw.ncpu) 2>&1 | tee output/battery4_run_mac2.log
Rscript make_tables4.R
```

If `output/ctmle_errors.txt` appears after the smoke test, send its first lines (and
`Rscript -e 'library(ctmle); args(ctmleGlmnet)'`) before running the full battery.

**Quick start (first run; copy and paste; 2026-09-07).** Step 0 unpacks the container
run's outputs (the paper's current numbers) into `analysis/output/`; the
tarball's paths start with `output/`, so it must be unpacked from inside
`analysis/`. Steps 1-3 then re-run the battery at the paper's replication
counts with all cores and the C-TMLE reference arm, and rebuild the tables.

```sh
cd ~/Documents/Claude/Projects/VariableSelection/paper/03-oal-inference/analysis
tar xzf ../phaseB_outputs_paper3_2026-09-06.tar.gz                    # 0. container outputs -> output/
Rscript -e 'for (p in c("glmnet","MASS","statmod","ggplot2","ctmle")) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")'
python3 -c 'import pandas' || pip3 install pandas                  # make_paper_tables4.py needs pandas (round 3, R3-m2)
Rscript sim_battery4.R 20 10 10 4 I                                    # 1a. 20-replication smoke test (about 1 min)
Rscript sim_battery4.R 500 300 500 $(sysctl -n hw.ncpu) 2>&1 | tee output/battery4_run_mac.log   # 1b. full run (hours; see below)
Rscript make_tables4.R                                                 # 2. writes output/battery4_tables.md (ODS-thr / ODS-bic labels)
Rscript apply_realdata2.R                                              # 3. RHC application (about 2 min)
```

Per-cell checkpoints are written to `output/battery4_cells/<cell>_R<reps>.rds`,
so an interrupted run resumes where it stopped. When the run is finished,
send `output/battery4_table.csv`, `battery4_selection.csv`,
`battery4_summary.txt` and `battery4_tables.md` back for transcription into
Table 1, Web Appendix C and Remark A5 (the container numbers in the paper
are from 500/150/200 replications; the ctmle rows are NA there).
Arm ids inside the scripts and CSVs are unchanged (`oads_thr`, `oads_bic`);
only the printed labels read ODS since the 2026-09-07 renaming.


All scripts are run from `analysis/` (they write to `analysis/output/`).
The container run used 2 cores and reduced replication counts (see the
header lines of `output/battery4_table.csv`); the commands below reproduce
the battery at the paper's counts and add the C-TMLE reference arm.

## 1. Simulation battery (`sim_battery4.R`)

Arguments: `n_main n_big n_bdry cores [cells]`.

* `n_main`  replications for cells I, II, III, M, P1, C, Bd, H, IV6 (paper: 500)
* `n_big`   replications for cell P2 (paper: 300; ~10 s per replication per core
  because the wAMD grid of the practitioner pipelines includes near-unpenalised
  fits at p = 200)
* `n_bdry`  replications for the nine boundary-sequence cells (recommended 500)
* `cores`   number of cores for `parallel::mclapply` (use all physical cores)
* `cells`   optional comma-separated subset of cell ids, e.g. `P2,B8000_k0.75`

Full-size run (about 35 core-hours in total; on an 8-core machine roughly
4-5 hours):

```sh
cd analysis
Rscript sim_battery4.R 500 300 500 $(sysctl -n hw.ncpu) 2>&1 | tee output/battery4_run_mac.log
```

Per-replication seeds are `seed_cell * 1000 + r`, so results do not depend on
the number of cores and the container run is reproduced exactly by
`Rscript sim_battery4.R 500 150 200 2`.

Outputs: `output/battery4_table.csv` (long format: cell, design, n, p, arm,
measure, value, mcse), `output/battery4_selection.csv`,
`output/battery4_summary.txt`. The first five lines of the CSV files are
comment lines (`#`) recording the replication counts and whether the
`ctmle` package was available; read them with `read.csv(..., comment.char = "#")`.

## 2. C-TMLE reference arm

Install the reference implementation once (not installable in the
container; CRAN is blocked there):

```r
install.packages("ctmle")   # Ju, Gruber & van der Laan; CRAN
```

`sim_battery4.R` detects the package with `requireNamespace("ctmle")` and
then runs `ctmle::ctmleGlmnet(Y, A, W, Q, ctmletype = 1, family = "gaussian",
gbound = 0.05, V = 5)` with the BIC-lasso outcome fit as the initial `Q`
(the same initial estimator as the lasso-path implementation of the previous
version). Before the full run, check the call against the installed version:

```r
library(ctmle); args(ctmleGlmnet)
Rscript sim_battery4.R 20 10 10 4 I     # 20-replication smoke test of cell I
```

If the call signature differs (older versions name the type argument
`ctmletype` and return `est` and `var.psi`; check `str(fit)`), adjust
`ctmle_arm()` in `sim_battery4.R` accordingly. The arm is reported in the
`ctmle` rows of `battery4_table.csv`; in the container run these rows are NA.

## 3. Bound tables (`bound_tables.R`)

Pure quadrature, runs in about a minute on one core; no arguments:

```sh
Rscript bound_tables.R
```

## 4. RHC application (`apply_realdata2.R`)

Requires `../../../9br-update/data/rhc_full.rda` (relative path from
`analysis/`); edit `DATA_PATH` at the top of the script if the data live
elsewhere. About two minutes:

```sh
Rscript apply_realdata2.R
```

**Round 2 (2026-09-16).** The 401(k) numbers of earlier versions came from
`apply_realdata.R` and a historical `pension.rda` that is not in the archive; they were
withdrawn from the paper and `apply_realdata.R` is no longer shipped. `apply_401k.R`
reruns the current pipeline on the public `hdm::pension` data (install `hdm` from CRAN
first); its output is not used in the paper.

## 4b. Numerical diagnostics (`convergence_check.R`, round 2, R2-M4)

Re-creates the first 25 replications of each main cell and 10 of each boundary cell
with the battery's own functions and records, per arm and refit, iteration-limit,
boundary, rank-deficiency and generalized-inverse events, plus the effect of a
200-iteration cap on the estimate, the stacked SE and the coverage indicator
(Table 11 of the paper). About 3 minutes:

```sh
Rscript convergence_check.R 25 10 $(sysctl -n hw.ncpu)     # -> output/convergence_check.csv/.txt, convergence_check_fits.csv
Rscript convergence_table.R                                # -> output/convergence_table.md: Table 11 rows + prose figures (round 3, R3-M4)
python3 example_A1.py                                     # Example A1 (Web Appendix B), exact rational arithmetic
Rscript ry_target_check.R 100 $(sysctl -n hw.ncpu)        # the Web Appendix D target check at the paper's count (2nd arg = cores)
```

## 5. Rendering

```sh
cd ..
quarto render paper.qmd --to pdf
pdfinfo paper.pdf | grep Pages
```
