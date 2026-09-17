# Replication archive: Inference after Thresholded Outcome-Based Confounder Selection under Confounder-Instrument Separation

Code and (public) data to reproduce the tables and figures of the paper.
Preprint: arXiv:XXXX.XXXXX (to be filled at posting). Author: Shoki Okubo (Toyo University).

## Contents
- analysis
- analysis/output
- data
- figures

`FILE_MANIFEST.txt` lists every file in this snapshot except itself, the `.git` history and the
release-check records (`RELEASE_CHECK*`), which are written after the manifest.

## Checked commit and release record
Computational commit checked against the manuscript: `95e9f135111ac8edbe2ce5fd48abe3a48a5f82fc` — see `RELEASE_CHECK.md` (with `RELEASE_CHECK_run.log` and `RELEASE_CHECK_sessionInfo.txt` when a clean-copy run was made). Later commits change documentation and the release record only — `git diff --stat 95e9f135111ac8edbe2ce5fd48abe3a48a5f82fc HEAD` lists them — so the scripts, data and outputs are those of the checked commit; after any change to code or outputs the release check is rerun and this line is regenerated.
Tag matching the posted preprint version: to be added at posting (`arxiv-<id>v<n>`).

## How to run
1. Install R (>= 4.1) and the packages pinned in `analysis/README.md` (the install command is given under "Quick start" below).
2. Run the scripts in each folder in the order given by their headers; outputs are written to `output/`.
   Where the folder contains its own `analysis/README.md`, that file is the authoritative
   description: execution order, which scripts generate results and which collect them, the
   environment of the shipped outputs and runtimes, and the map from every table, figure and
   script-produced number to a script and an output file. Where a `stage_figures.sh` is present,
   run it as `bash stage_figures.sh figures` from this folder (or `bash ../stage_figures.sh <dest>`
   from `analysis/`): it copies the manuscript figures from the outputs to `figures/`; the PDF
   build itself (Quarto) is not part of the numerical reproduction.
3. The application uses the public RHC data included in `data/` (rhc_full.rda; `apply_realdata2.R` looks for it under `data/` relative to `analysis/` as well as in the author's project tree); the simulations generate their own data. No licensed microdata are included.

## Quick start (paper 3)
CRAN packages: glmnet, MASS, statmod, ggplot2 (`undercoverage_demo.R`), ctmle (`hdm` only for the
optional `apply_401k.R`); Python 3 with `pandas` for the table conversion `make_paper_tables4.py`,
and the standard library only for the exact-arithmetic check `example_A1.py`. All scripts run from
`analysis/` and write to `analysis/output/`. Check the dependencies before running (round 3, R3-m2):
```sh
Rscript -e 'for (p in c("glmnet","MASS","statmod","ggplot2","ctmle")) cat(p, if (requireNamespace(p, quietly=TRUE)) "ok" else "MISSING", "\n")'
python3 -c 'import pandas; print("pandas", pandas.__version__)'
```
The sequence that produced the shipped outputs (per-replication seeds make the battery independent
of the number of cores; the battery took 101.5 minutes on 8 cores in an independent rerun):
```sh
cd analysis
Rscript sim_battery4.R 500 300 500 8 && Rscript make_tables4.R      # Table 2, Web Appendix C tables
Rscript apply_realdata2.R                                          # Section 5 (RHC), ~2 min
Rscript bound_tables.R                                             # Web Appendix B quadrature tables
Rscript undercoverage_demo.R && Rscript proposed_estimator.R      # Section 2
Rscript ry_target_check.R 100 8                                    # Web Appendix D target check (100 replications; 2nd arg = cores)
Rscript convergence_check.R 25 10 8                                # Web Appendix C numerical diagnostics, ~3 min
Rscript convergence_table.R                                        # Table 11 rows and prose figures from the diagnostic output
Rscript weight_floor_check.R                                       # the adaptive-weight floor eps_n (round 4, R4-m3): equals 1e-8 at every n used; selections identical
python3 example_A1.py                                              # Example A1 (exact rational arithmetic)
```
`sim_battery4.R` writes per-cell checkpoints to `output/battery4_cells/` (not shipped); delete them to
force a fresh run. Labels: the scripts and CSVs use `oads_thr`/`oads_bic` for the arms called OBS-thr/OBS-bic
in the paper and `goal_ipw` for the pipeline called OAL-EN; they are the same objects. `apply_401k.R`
reruns the current pipeline on the public `hdm::pension` data; its output is not used in the paper.

## Citation
Okubo, S. (2026). Inference after Thresholded Outcome-Based Confounder Selection under Confounder-Instrument Separation. Working paper. arXiv:XXXX.XXXXX.
