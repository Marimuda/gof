# Status for version 1.3.0

The native structural refactor and compiler/platform qualification were completed
for 1.2.1; C++17 remains the supported package baseline and C++23 an integration
configuration. No further native rewrite is required for this release.
The fixed-order finite-panel study is now complete with its limited empirical
support and adverse findings in the installed validation guide. The broader
statistical extensions below remain research proposals, not missing promised
endpoints or claims made by this release. The historical work-order discussion
is retained as context.

---

# Deferred statistical extensions after C++ modernization

Recorded 2026-10-30. The qualified statistical baseline is `v1.2.0`, commit
`27857703e13c6516504e3145521c3ed9a259fc6e`. Its implementation, review findings,
calibration results and limitations are recorded in [CAPABILITIES.md](CAPABILITIES.md).

The agreed order is to refactor the C++ implementation first, preserving the
validated statistical behavior, and then consider the extensions below.
This is a deferred-work register, not a claim that these extensions are
implemented or a commitment to add every one of them.

## Findings: assumptions and implementation limits are different

The current endpoints have been implemented and checked within their stated
scope. That does not mean the package implements every extension in the cited
literature or that its support boundaries are as broad as possible.

| Condition | Interpretation |
|:--|:--|
| Correctly aligned independent residual units and a valid influence representation | Required by the current multiplier construction. Independence is between units/subjects; correlated outcomes within a subject are already part of the Gaussian SEM construction. |
| Identifiable estimation, suitable information and an appropriate regular likelihood root | Standard conditions for the fitted-parameter regular asymptotic methods. Known parameters, structurally zero-information processes and explicitly handled extended-likelihood boundaries have their separately documented behavior. |
| Valid model moments and appropriate conditional centering | Mathematical conditions for the selected process. Arbitrary residual/ordering callbacks do not establish these conditions merely by returning numeric arrays. |
| Sufficient effective information | Large nominal sample size and numerical convergence alone do not guarantee an accurate asymptotic approximation. Rare events, leverage and imbalance can still matter. |
| Prespecified ordering families, bandwidths and grid rules | Needed for the advertised formal inference. Response-driven selection requires its own calibration/adjustment. |
| Complete-data-only SEM, restricted bootstrap links and selected joint adapters | Implementation boundaries. Some can be broadened by implementing and validating additional constructions. They are not universal requirements of goodness-of-fit testing. |

Normality is the hypothesis being tested by `cumres_normality()`, not a fact
the user must establish about the observed data before running it. Its
Gaussian reference model and other assumptions determine the null calibration;
rejection need not uniquely identify which latent/measurement error source or
other model component is wrong.

Numerical uncertainty is another distinction. Some valid mathematical cases
may currently be rejected or receive missing p-values because the numerical
checks cannot establish the required condition. Better numerical qualification
could admit more cases without changing the statistical assumptions. Simply
weakening the guard would not establish validity.

## Candidates to revisit

All entries are **deferred**. Prioritize an entry only when a concrete research
use case fits gof's existing cumulative-residual diagnostic purpose.

| ID | Current boundary | Proposed development and acceptance evidence |
|:--|:--|:--|
| S01 | SEM diagnostics require complete manifest data; incomplete/FIML fits are rejected. | Assess the published missing-data construction under explicit missingness assumptions, including missing at random where justified. Match residual projections, observed-data likelihood scores, information and parameter influences to each observed pattern. Require a complete-data reduction, independent pattern-wise reference calculations and prespecified calibration/power checks. Do not silently substitute deletion or imputation. |
| S02 | Clustering across subjects is unsupported by the current independent-unit adapters. | Define the independent cluster and derive cluster-level scores/influences and shared draws. Preserve within-cluster covariance, clarify normalization and effective cluster count, and test unequal cluster sizes plus the singleton-cluster reduction. Supporting correlated outcomes within one SEM subject does not already supply this extension. |
| S03 | Weighted SEM fits are rejected; Poisson response bootstrap accepts only unit positive prior weights. | Specify the weight meaning before implementation: precision, frequency, sampling and other weights are different models. Match the estimator and influence calculation, or define the actual response sampler. Require unit-weight reduction and independently justified weighted references. A prior-weight variance convention alone does not identify a Poisson bootstrap sampler. |
| S04 | Response bootstrap supports binomial/logit and Poisson/log, with fixed ordering variables. | Evaluate additional links within supported families first, deriving domain/boundary and refitting behavior. Separately assess fitted-prediction ordering with the ordering/statistic recomputed consistently for every refit. Freezing a response-dependent fitted ordering is not an adequate substitute. Preserve every replicate and its unresolved-fit bounds. |
| S05 | `cumres_joint()` accepts justified default inputs and linear/GLM fits, but not SEM collections or ordinal adapters. | Extend the joint process to aligned subject-level residual vectors, nuisance influences and all relevant cutpoints/outcomes. Define the selected family and any scaling across different response units. Require independent cross-covariance oracles, marginal reductions, subject-permutation checks, family calibration and power. Keep global-null inference distinct from strong familywise control under arbitrary partial alternatives. |
| S06 | Multi-sample Gaussian fits are rejected by the normality endpoint. | If needed, derive group-specific residual maps and information with the fitted parameter-sharing constraints. Define pooling/normalization and the target hypothesis. Require the single-group reduction and independent shared-versus-separate-parameter checks. Do not concatenate groups and assume the single-sample construction applies. |
| S07 | Regular fitted-parameter paths reject separation, unresolved boundaries and numerically indeterminate information; only specifically justified exceptions are handled. | Keep necessary regularity guards. Investigate improved scaling/rank/feasibility qualification separately from genuinely nonregular inference. Broader boundary support needs a statistical derivation and explicit estimator/process definition, not a threshold change. Validate difficult regular cases against boundary counterexamples; retain failure reasons. |
| S08 | Gaussian SEM normality uses a first-order asymptotic correction on a finite probability grid. | If a research use case requires stronger finite-sample performance, separately assess a fully specified refitting/calibration procedure and grid/tail sensitivity. Freeze evaluation cohorts and handling of failed fits before testing. More simulations reduce Monte Carlo noise but do not remove approximation error; a fitted bootstrap is not automatically exact either. |

### Normality evidence to retain

The original 200-cohort study at n=300 and R=199 had inclusive `p <= .05`
rejection rates of 6--8.5%. Two 8.5% L2 results triggered a prespecified fresh
paired study rather than being dismissed because a broad automated screen
passed. With 200 fresh cohorts at n=300/n=1,200 and R=399, rates were 2.5--7%
and 3.5--5%, respectively, without unresolved fits.

No implementation defect was confirmed. The evidence does not prove exact
5% calibration, equivalence to nominal size or a significant sample-size
convergence trend. Preserve both studies and their uncertainty when assessing
S08; do not replace the original results or choose thresholds/seeds to obtain
a preferred conclusion. Detailed records are under the course audit path
`audits/gof-capabilities-2026-10-29/normality-supplement/`.

## C++ refactor first

Use the tagged baseline and its source archive/evidence as the statistical
reference. This note is a subsequent documentation addition; it does not move
or redefine that tag.

1. Review the native structure and define boundaries between numerical kernels,
   state/ownership, random-number generation and the Rcpp interface. Preserve
   the public R API and installed C++ contracts, or explicitly plan compatibility
   when an interface change is unavoidable.
2. Modernize to the chosen C++23 structure without adding statistical scope in
   the same changes. Preserve independent-unit versus sorted-position identity,
   tie/window/grid semantics, influence scaling, statistic definitions, failure
   handling and promised seeded legacy behavior. Preserve original authorship
   and license.
3. Validate the changed native code against the frozen golden/reference cases,
   independent process oracles, installed downstream client checks and native
   sanitizers. Investigate any numerical or random-stream differences; rerun
   the affected scientific calibration when algorithmic behavior changes.
4. Review compiler/R minimums and execute the relevant platform checks. The
   existing remote platform/minimum-R CI configuration is not evidence of a
   completed run. Compilation in C++23 mode alone is not an architecture refactor.
5. Establish the refactored baseline, then select deferred statistical work by
   research need. Handle each extension as a separate method/API change with
   its own assumptions, reference implementation, failure policy, calibration,
   power checks and installed documentation.

The refactor can make later extensions easier to implement and review. It
does not itself relax assumptions, create missing-data/cluster support or
turn asymptotic inference into an exact finite-sample procedure.

## Sources and trace

* [CAPABILITIES.md](CAPABILITIES.md): implemented contracts, qualification and
  scientific limits for the v1.2.0 baseline.
* [STATISTICAL_REVIEW.org](STATISTICAL_REVIEW.org): independent review and
  counterexamples that motivated the statistical repairs.
* [Lin, Wei and Ying (2002)](https://doi.org/10.1111/j.0006-341X.2002.00001.x):
  cumulative-residual methods, including development for dependent observations.
* [Sánchez, Houseman and Ryan (2009)](https://doi.org/10.1111/j.1541-0420.2008.01022.x)
  and their [supplementary methods](https://citeseerx.ist.psu.edu/document?doi=28e5d2756f628ef91b69bbe1d448511cc3bf1f08&repid=rep1&type=pdf):
  Gaussian SEM residual diagnostics and missing-data considerations, especially
  Web Appendix B. A published extension is a development starting point,
  not evidence that it is already supported by this package.
