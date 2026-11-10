# CLMM experimental correction screen

This development screen exercises the public experimental score-projection
corrections. It does not qualify CLMM for release or replace the original
failed full and slow studies. The derivative correction remains the default.

The plan was recorded before outcomes: 16 datasets for each of the original
29 designs, including five alternatives, for 464 attempts. Data seeds are
`180000000 + i*100000 + r`, where `i` is the full-study cell index and `r`
the replicate. `designs.csv` identifies every design. The actual multiplier
engine calculates both Sup and L2, with 199 draws and seed `data_seed + 50000`
shared across the three corrections. All methods use the same fitted model
and observed process; only their reference distributions differ.

Fitting uses `refit=FALSE`. Cauchit inputs use empirical marginal starting
cutpoints and `ordinal::clmm(nAGQ=-21)`, followed by gof controls `nodes=21`
and `max.refresh=20`. This exception was predeclared. The screen therefore
does not measure the original Cauchit input workflow or optional recovery.

## Results and limits

Of 464 datasets, 412 produced informative, finite results for all three
methods. The 52 refusals comprise 27 covariance-boundary refusals, 24 polishing
nonconvergence failures and one input nonconvergence failure. No refusal arose
at the correction stage. Boundary refusal is a diagnostic category, not proof
that each refused input was a global constrained likelihood maximum.

Sparse random slopes remain problematic: acceptance was 0/16 for `core07`,
5/16 for `core08` and 9/16 for `core09` (20, 50 and 150 clusters, respectively,
with three observations per cluster). Cauchit accepted 16/16 with the specified
input workflow. Neither result establishes a population availability rate.

`summary.csv` includes all 174 design/method/statistic rows, with rejection
counts at nominal 5%, exact binomial 95% intervals conditional on acceptance,
and refusal fractions. `attempts.csv.gz` contains all 1,392 method/attempt
rows, including messages and typed failure codes. A failure in any correction
would refuse the common comparison; the recorded stage distinguishes this
from fitting failures. No such correction failure occurred here.

Sixteen datasets per design cannot distinguish acceptable calibration from
the failures under investigation. Zero rejections among 16 accepted datasets,
for example, has an exact 95% upper bound of about 20.6%. Do not pool different
designs into a calibration claim, select a winner from this pilot, or count
refusals as non-rejections. The larger, earlier L2 comparison remains in
`../clmm-repair/`; its evidence and limitations still apply.

## Independent boundary profiles

`boundary-profiles.csv` contains constrained likelihood profiles for historical
`core11` replicate 19 and `core14` replicate 30 from the original seed roster.
These are deliberately reused diagnostic fixtures, separate from the fresh
464-case screen. Direct Gaussian-prior quadrature independently evaluates the
likelihood without the package's adaptive modes, score or optimizer.

Nuisance parameters are reoptimized at conditional variance
`0, 0.0025, 0.01, 0.04`, using 31 nodes per dimension and checking at 61.
All eight optimizations report convergence. Both profiles favor zero among
these points: moving to 0.0025 lowers log likelihood by 0.03868 and 0.00596,
respectively. Maximum 31-to-61-node disagreement is below 0.000032.
These local numerical profiles support a genuine boundary issue in the two
fixtures. They do not certify a global maximum, explain all refusals, or
provide a valid boundary reference distribution for residual-process tests.

## Replay and provenance

Install this development version in an isolated R library and select it with
`R_LIBS`; preserve the ordinary installed package. From the repository root:

```sh
Rscript tools/check-clmm-inference.R NEW_OUTPUT 6
Rscript tools/check-clmm-inference.R NEW_SINGLE_CASE_OUTPUT 1 core05 1
Rscript tools/check-clmm-boundary.R NEW_PROFILE_OUTPUT
```

Replay requires the package dependencies plus `ordinal`, `numDeriv` and
`jsonlite`. The export audit additionally uses `digest`. Output directories
must be new. The inference runner records its invocation and installed R
session; data seeds do not depend on worker count. One accepted case and one
refused case were checked for exact replay, including all method results or
the retained error condition, against the original atomic records.

The original numerical run used the working source identified by the sealed
R-source and data-generator hashes in `pilot-plan.json`. An extra closing
parenthesis affected its post-run summary only: all 464 numerical records had
already been saved. The corrected public runner was checked to completion on
a single-case replay. `tools/export-clmm-inference.R` independently rebuilt
the complete tables, checked the roster and seeds, and verified the sealed
source hashes. No numerical attempt was replaced to repair that summary.

`PROVENANCE.json` hashes the exported evidence and the retained original
runner. `raw-records-sha256.csv` identifies all retained private RDS records;
those records contain numerical diagnostics and condition objects beyond the
compact public table. The exporter takes the private audit directory and a
new output directory. Exact replay can depend on the R and dependency versions.

Remaining release blockers are reliable sparse fitting, inference when
variance components reach zero, and separate within-cluster dependence
diagnostics: the marginal residual process cannot detect every departure from
the joint model. After choosing and freezing the entire fitting/inference
policy, a fresh study must meet the unchanged full-scope qualification criteria.
