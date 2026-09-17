# Table 11 entries generated from output/convergence_check.csv

Source run: convergence_check.R  run on 2026-09-16 11:09; replications per cell: main 25, boundary 10; cores 12; 0.7 min

| Cell | Arm | Reps | Refits | Dim. | Iter. limit | Gen. inverse | Clipped (pct) | Max Δ estimate | Max Δ stacked SE | Coverage changes |
|---|---|---|---|---|---|---|---|---|---|---|
| II | full set | 25 | 125 | 51 | 42 | 5 | 79.3 | 0.560 | 0.113 | 0 |
| P2 | naive OAL selection + AIPW | 25 | 125 | 62.0 | 12 | 1 | 40.1 | 0.475 | 0.028 | 0 |
| other | all other 169 cell–arm pairs | 25 or 10 | 15,725 | 3.2–101 | 0 | 0 | — | — | — | — |

Prose figures (Web Appendix C, numerical diagnostics paragraph):
- total refits examined: 15,975 (column `fits`, summed)
- refits that stopped at the iteration limit: 54 (column `ps_nonconv`, summed)
- generalized-inverse fallbacks: 6 (column `ginv_fits`, summed)
- II / full set: 42 non-converged refits in 16 of 25 replications, 5 generalized inverse(s), 79 percent clipped; max |Δ estimate| 0.560, max |Δ stacked SE| 0.113, coverage changes 0
- P2 / naive OAL selection + AIPW: 12 non-converged refits in 9 of 25 replications, 1 generalized inverse(s), 40 percent clipped; max |Δ estimate| 0.475, max |Δ stacked SE| 0.028, coverage changes 0

