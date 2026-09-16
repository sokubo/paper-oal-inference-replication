# Release check — paper-oal-inference-replication

Date: 2026-09-16T05:10:41Z. Snapshot downloaded anonymously (no credentials, no gh CLI) from `https://codeload.github.com/sokubo/paper-oal-inference-replication/tar.gz/main`.

- ref: `main`; commit: `763af880036492ed44545416baedac6c32486c48`
- archive SHA-256: `8c5cf236252b9275b44a0847b66bd35ea72845b39fb182c6868fff9989d256d7`
- files in snapshot (excluding FILE_MANIFEST.txt and RELEASE_CHECK*): 75; listed in FILE_MANIFEST.txt: 75; missing from snapshot: 0; not listed in manifest: 0
- clean run: FAILED: the documented sequence did not complete in the clean copy (5376s); see RELEASE_CHECK_run.log; comparison below covers the outputs produced before the failure
- staged figures: figure staging not run

## Environment of the clean run

```
sessionInfo not written
```

## Regenerated vs shipped outputs

```
                                 file                      status max_abs_diff
   analysis/applications2_overlap.txt                   identical     0.000000
      analysis/applications2_pool.txt                   identical     0.000000
 analysis/applications2_selection.csv                   identical     0.000000
    analysis/applications2_splits.csv                   identical     0.000000
     analysis/applications2_table.csv                   identical     0.000000
     analysis/applications2_table.txt                   identical     0.000000
    analysis/applications2_trim01.csv                   identical     0.000000
      analysis/battery4_selection.csv                   identical     0.000000
        analysis/battery4_summary.txt                   identical     0.000000
          analysis/battery4_table.csv          numeric difference     2.528261
          analysis/battery4_tables.md                   identical     0.000000
            analysis/ctmle_errors.txt                   identical     0.000000
         analysis/battery_summary.txt shipped but not regenerated           NA
           analysis/battery_table.csv shipped but not regenerated           NA
        analysis/battery2_summary.txt shipped but not regenerated           NA
          analysis/battery2_table.csv shipped but not regenerated           NA
        analysis/battery3_summary.txt shipped but not regenerated           NA
          analysis/battery3_table.csv shipped but not regenerated           NA
       analysis/bound_tables_bias.csv shipped but not regenerated           NA
       analysis/bound_tables_rmis.csv shipped but not regenerated           NA
            analysis/bound_tables.txt shipped but not regenerated           NA
  analysis/convergence_check_fits.csv shipped but not regenerated           NA
       analysis/convergence_check.csv shipped but not regenerated           NA
       analysis/convergence_check.txt shipped but not regenerated           NA
        analysis/proposed_summary.txt shipped but not regenerated           NA
          analysis/proposed_table.csv shipped but not regenerated           NA
         analysis/ry_target_check.csv shipped but not regenerated           NA
         analysis/ry_target_check.txt shipped but not regenerated           NA
   analysis/undercoverage_summary.txt shipped but not regenerated           NA
     analysis/undercoverage_table.csv shipped but not regenerated           NA

files compared: 30; identical: 11; numeric difference: 1 (max 2.53e+00); other: 18
```
