# CLMM repair development evidence

These results investigate the failed qualification; they do not replace it.
The original full and slow cohorts, thresholds and release hold remain unchanged.
No inference candidate from this comparison is a public package default.

## Fresh inference comparison

The plan was recorded before examining the 1,500 new datasets. Each of three
logit designs has 500 attempts, with data seeds `160000000 + i*100000 + r`,
where `i` is the original full-study cell index and `r` the replicate.
Fitting uses the unchanged baseline at commit
`ff80632ae542b4de6719282d21df9b144fe3a15c`. All 29 refusals occurred before
the Gaussian-reference calculation and remain in the availability denominator.

The observed full-grid L2 statistic is identical for all candidates. Its
conditional Gaussian reference is calculated using 9,999 draws of the
eigenvalue-weighted chi-square quadratic form. The construction is checked
against the package's independent L2 norm implementation. This differs from
the original cohorts' 199 multiplier draws and is recorded as a new screen.

| Design | Accepted / attempted | Original correction | Score projection | Projection + HC2 |
|---|---:|---:|---:|---:|
| Random intercept, 50 clusters, 10 rows/cluster | 498/500 | 3.61% | 7.23% | 3.21% |
| Random slope, 50 clusters, 10 rows/cluster | 473/500 | 1.69% | 6.77% | 4.44% |
| Random slope, 150 clusters, 10 rows/cluster | 500/500 | 2.20% | 5.60% | 4.80% |

Percentages are rejection rates at nominal 5%, conditional on accepted fits.
`prospective-summary.csv` includes exact binomial 95% intervals. All 1,500
attempts, messages and method-specific p-values are in
`prospective-attempts.csv.gz`.

For each ordering point, let `A` contain the accumulated residual contribution
from each independent cluster and `S` the cluster likelihood scores. The original
reference adds the residual derivative/information correction to `A`. The
projection candidate instead uses `(I - Q Q') A`, where `Q` is an orthonormal
basis for the columns of `S`. HC2 additionally divides cluster row `c` by
`sqrt(1 - h_c)`, with `h_c = (Q Q')_cc`. These are development candidates:
the HC2 adjustment is not a proof of calibration or boundary validity.

The HC2 candidate improves the larger slope design's observed size, but the
50-cluster slope design still refuses 5.4% of attempts, exceeding the unchanged
criterion of less than 5%. This screen covers neither Sup, power, the remaining
designs, nor the new numerical-recovery workflow. Candidate selection is
development; a repaired procedure needs a separately frozen full-scope study.

For replay, install commit `ff80632` in an isolated library and run the current
repository's `tools/check-clmm-correction.R NEW_OUTPUT WORKERS`. Optional
`CELL REPLICATE` arguments replay a single sealed case. Do not install over an
ordinary gof library. The data generator remains in `tools/clmm-study.R`.

## Bounded numerical workflow

A separate pre-recorded plan used 96 fresh datasets, seeds
`170000000 + i*100000 + r`, and the frozen repair source identified by
`numerical-source-sha256.json`. It enabled optional numerical recovery.
Cauchit inputs used empirical marginal cutpoints and `ordinal::clmm(nAGQ=-21)`;
this is a changed input workflow, not a replay of the original study.

| Design | Finite output / attempted | Exceeded 60 seconds |
|---|---:|---:|
| Random slope, 50 clusters, 3 rows/cluster | 20/32 | 12/32 |
| Nested intercepts, 50 top-level clusters, 3 rows/child | 27/32 | 5/32 |
| Cauchit intercept, 50 clusters, 10 rows/cluster | 32/32 | 0/32 |

There were no other process failures or completed numerical refusals. Timeouts
are execution-budget outcomes, not evidence of boundary optima. Only three
multiplier draws per accepted case checked finite API output: no p-value size
claim follows from this screen. Runtime outcomes depend on the machine.
`numerical-attempts.csv.gz` retains all 96 dispositions.

`tools/check-clmm-numerics.R NEW_OUTPUT WORKERS` reproduces this workflow in an
isolated repair library. Later explicit `stats::` prefixes on quantile calls do
not change the numerical source's mathematics. Global-optimum certification is
not supplied by an accepted local score root.

## Independent numerical fixtures and provenance

`numerical-checks.csv` reports one sparse-slope and three Cauchit fixtures
checked against independent direct Gaussian-prior integration. Maximum absolute
log-likelihood disagreement was below `4.4e-6`. These reused failure fixtures
check numerical agreement, not finite-sample calibration. Reproduce them with
`tools/check-clmm-refit.R NEW_OUTPUT` in an isolated repair library.

The two plan JSON files retain their pre-outcome scope. `PROVENANCE.json`
records hashes of raw records, original runner scripts and exported files.
The compact export is audited by `tools/export-clmm-repair-evidence.R`.
The recorded R session is in `session-info.txt`.
