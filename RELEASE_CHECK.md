# Release check — paper-oal-inference-replication

Date: 2026-09-17T03:07:54Z. Snapshot downloaded anonymously (no credentials, no gh CLI) from `https://codeload.github.com/sokubo/paper-oal-inference-replication/tar.gz/main`.

- ref: `main`; commit: `95e9f135111ac8edbe2ce5fd48abe3a48a5f82fc`
- archive SHA-256: `a81b1297c002159d51a05610ffd87c6b3ab00bdab712560d3ed761e9ad9ff4df`
- files in snapshot (excluding FILE_MANIFEST.txt and RELEASE_CHECK*): 77; listed in FILE_MANIFEST.txt: 77; missing from snapshot: 0; not listed in manifest: 0
- clean run: documented sequence executed in a clean copy with shipped outputs set aside (5192s); log and sessionInfo kept; comparison below
- staged figures: figure staging not run

## Environment of the clean run

```
R version 4.6.0 (2026-04-24)
Platform: aarch64-apple-darwin23
Running under: macOS Sequoia 15.7.3

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] ja_JP.UTF-8/ja_JP.UTF-8/ja_JP.UTF-8/C/ja_JP.UTF-8/ja_JP.UTF-8

time zone: Asia/Tokyo
```

## Regenerated vs shipped outputs

```
                                 file                      status max_abs_diff
   analysis/applications2_overlap.txt                   identical 0.000000e+00
      analysis/applications2_pool.txt                   identical 0.000000e+00
 analysis/applications2_selection.csv                   identical 0.000000e+00
    analysis/applications2_splits.csv                   identical 0.000000e+00
     analysis/applications2_table.csv                   identical 0.000000e+00
     analysis/applications2_table.txt                   identical 0.000000e+00
    analysis/applications2_trim01.csv                   identical 0.000000e+00
      analysis/battery4_selection.csv                   identical 0.000000e+00
        analysis/battery4_summary.txt                   identical 0.000000e+00
          analysis/battery4_table.csv                   identical 0.000000e+00
          analysis/battery4_tables.md                   identical 0.000000e+00
       analysis/bound_tables_bias.csv          numeric difference 1.012523e-13
       analysis/bound_tables_rmis.csv          numeric difference 1.136008e-12
            analysis/bound_tables.txt          numeric difference 1.182327e-14
  analysis/convergence_check_fits.csv                   identical 0.000000e+00
       analysis/convergence_check.csv                   identical 0.000000e+00
       analysis/convergence_check.txt                   identical 0.000000e+00
        analysis/convergence_table.md                   identical 0.000000e+00
            analysis/ctmle_errors.txt                   identical 0.000000e+00
        analysis/proposed_summary.txt                   identical 0.000000e+00
          analysis/proposed_table.csv          numeric difference 4.400036e-12
         analysis/ry_target_check.csv          numeric difference 1.298299e-10
         analysis/ry_target_check.txt                   identical 0.000000e+00
   analysis/undercoverage_summary.txt                   identical 0.000000e+00
     analysis/undercoverage_table.csv          numeric difference 9.992007e-16
         analysis/battery_summary.txt shipped but not regenerated           NA
           analysis/battery_table.csv shipped but not regenerated           NA
        analysis/battery2_summary.txt shipped but not regenerated           NA
          analysis/battery2_table.csv shipped but not regenerated           NA
        analysis/battery3_summary.txt shipped but not regenerated           NA
          analysis/battery3_table.csv shipped but not regenerated           NA

files compared: 31; identical: 19; numeric difference: 6 (max 1.30e-10); other: 6
```
