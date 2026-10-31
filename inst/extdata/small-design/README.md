# Earlier small-design development evidence

These files preserve an exploratory screen run with gof 1.2.1. They are not
fresh statistical qualification of the development package. No response
cohorts were regenerated to assemble this extract.

* `development-records.csv`: all 400 original cohorts (200 per design), with
  two statistics for each of two methods: 1,600 rows. Failures are retained.
* `development-summary.csv`: the eight original summary rows, recomputed
  from the retained records and checked against the frozen original.
* `development-designs.rds`: fixed design matrices, coefficients, null means,
  lexicographic covariate-cell orderings and zero offsets. There are no
  original response vectors or confidential course records.
* `PROVENANCE.json`: original protocol and source SHA256 hashes. Integer
  seeds are recorded as in that protocol; it did not explicitly record an
  RNG kind. This extract does not invent one retrospectively.

The formulas were selected beforehand from the course examples and then held
fixed. The screen conditions on those designs and formulas; it does not repeat
or qualify the selection process. Both methods used 199 inner simulations.
Inclusive decisions use p <= .05. A cohort with an outer error contributes
[0,1] decision bounds; unresolved inner tails retain their original completion
bounds. Every summary denominator is 200. The 95% Monte Carlo envelope is
pointwise, with the lower Clopper--Pearson limit for definite rejections and
the upper limit for possible rejections. It is not a simultaneous guarantee,
an equivalence test, or an interval for infinite-inner-budget calibration.

See the installed small-design vignette for interpretation and limitations.
Original course evidence remains in
`audits/gof-class-examples-2026-10-30/class-design-calibration.{rds,csv}`;
the extract and comparison script is retained in
`audits/gof-small-design-guidance-2026-10-31/extract-development-evidence.R`.
