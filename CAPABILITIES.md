# Statistical capability development after 1.1.0

The repaired 1.1.0 baseline is commit `bd473e8`. This phase adds capabilities
within gof's purpose: graphical and numerical model diagnostics based on
cumulative residual or residual-distribution processes. C++ architecture
modernization follows a separately committed and tagged statistical baseline.
No C++ implementation or installed header changed in this phase.

## Implemented capabilities

| Endpoint | Question and supported scope | Reference simulation |
|:--|:--|:--|
| `cumres_joint()` | Departure along any member of a prespecified ordering family, for supported linear/GLM fits or justified explicit residual/influence inputs. | Common independent-unit Gaussian multipliers across orderings; family Sup and range-normalized L2 maxima. |
| `cumres_normality()` | Gaussian SEM standardized marginal or projected structural-residual distribution departures, for complete independent-subject regular Gaussian lava fits. | Joint Gaussian response vectors, likelihood scores, expected information and the first-order residual-CDF nuisance correction. |
| `cumres_bootstrap()` | Residual-process calibration under independent canonical binomial/logit or Poisson/log sampling on a fixed design. | Response simulation and refitting for fitted means; no refitting for genuinely specified means. |

Each endpoint takes the model first, works with native `|>`, returns a
`cumres` subclass, and has installed help and a worked vignette. Tables,
summaries, optional tidiers and plots preserve method-specific interpretation.
Existing marginal `cumres()` calculations and seeded results remain
compatibility requirements.

## Mathematical contracts and decisions

### Joint ordering inference

The original engine reused a normal draw by sorted rank under each ordering.
That produces each marginal null law but cannot preserve cross-ordering
covariance. The new engine attaches multiplier `G_i` to the same independent
unit under every ordering. For selector `a_ji(t)` it uses

```
W*_j(t) = [sum_i a_ji(t) r_i G_i
           + {sum_i a_ji(t) dr_i}' sum_i IC_i G_i] / sqrt(n).
```

This follows the common process construction in Lin, Wei and Ying (2002),
Section 2.2. A fixed finite-family maximum is a continuous functional of that
joint process. The family statistics are maximum Sup and maximum L2 divided
by the square root of its coordinate range. Range normalization avoids
arbitrary relative weighting from coordinate units. A one-point grid has zero
L2. The range-normalized maximum is a derived extension, not a claim that
the article specifies that exact aggregation separately.

The simultaneous Sup envelope is a **null reference envelope** across the
selected orderings and evaluated points, not a parameter confidence interval.
Adjusted marginal tails use the full model null; there is no strong
familywise-error claim under arbitrary partial alternatives. This is not
the paper's full multivariate lower-orthant omnibus process.

SEM collections and ordinal fits are rejected by the wrapper. Their adapters
also reject an unsupported joint switch instead of silently suggesting joint
support. Structurally zero-information processes contribute zero. Any
numerically indeterminate selected process withholds family inference;
the family is never silently reduced.

### Gaussian SEM normality

The prior SEM adapter implemented mean/linearity checks only. The new endpoint
implements the complete-data setting of Sánchez, Houseman and Ryan (2009),
Section 3.1, equations (10)--(11). For residual map `A(theta)`, holding the
generating mean/covariance law fixed gives the residual-CDF derivative

```
delta_l(z) = phi(z) * {mean(A d_l mu) - z mean(d_l A Sigma A')}.
```

Each Gaussian response vector supplies both an ECDF and its joint likelihood
score. Add `delta' I_n^-1 S*`, using total expected information `I_n`, before
the common `sqrt(n)` scaling. Numeric derivatives include shared mean/covariance
parameters and projection standardization. Exogenous covariates stay fixed.
Fitted-score adequacy, covariance interiority and information identification
have explicit numerical checks, not global or exact-arithmetic certificates.

`variable=NULL` pools symmetrically whitened marginal components. A named
endogenous/latent variable selects its standardized observable structural-error
projection. The pool is not a universally consistent multivariate-normality
or dependence test. A projected residual cannot uniquely identify the
nonnormal error source under arbitrary alternatives.

The fixed probability grid explicitly approximates the full-domain statistics.
L2 is the left-step norm over normal probability. Missing/FIML, clustered,
weighted, unidentified and active-boundary fits remain unsupported.

### Binomial and Poisson parametric calibration

This addition uses the supported sampling distribution when a multiplier
approximation may have little effective information. It does not claim to
solve every sparse-data estimation problem.

Fitted-null calibration simulates fitted responses and refits each sample
with `glm.fit()` on the stored design, offsets and binomial trials. It never
reevaluates the user's call and has no inner multiplier simulation. It remains
approximate and requires adequate estimation and bootstrap regularity.
Genuinely specified means instead permit a distinct exchangeable Monte Carlo
rank test conditional on the fixed design/grid. Substituting estimates from
the tested outcomes loses that guarantee.

Grouped residual totals match tied Bernoulli expansion without allocating the
expanded sample. Nonunit positive Poisson weights do not specify a count
sampling distribution and are rejected here. Orderings must be fixed;
response and fitted-prediction orderings are rejected. One grid is reused.

Every replicate is recorded. Some all-zero/all-success likelihood boundaries
have analytically determined means when a constant design column supplies the
direction. A fitted all-zero null is degenerate and cannot assess unseen
events. Other unresolved fits are never dropped or replaced: informative
p-values are missing, and bounds account for every unresolved tail indicator.
Those are calculation bounds, not confidence intervals.

Integration review found that withholding all p-values also hid the known
answer for structurally zero-information processes. These now retain zero
statistics and p-values/bounds of one while failed-refit statuses stay visible.
Sparse saturated binomial and Poisson regressions cover this exception.

## Review and evidence

Primary references:

* [Lin, Wei and Ying (2002)](https://doi.org/10.1111/j.0006-341X.2002.00001.x),
  especially the common process and moving-window construction.
* [Sánchez, Houseman and Ryan (2009)](https://doi.org/10.1111/j.1541-0420.2008.01022.x),
  residual definitions and equations (10)--(11).
* [Genest and Rémillard (2008)](https://doi.org/10.1214/07-AIHP148), for the
  regular parametric-bootstrap principle, not a finite-sample GLM theorem.
* [Phipson and Smyth (2010)](https://doi.org/10.2202/1544-6115.1585), Section 4,
  for Monte Carlo rank testing. Plus-one tails do not make asymptotic laws exact.

The course audit directory `audits/gof-capabilities-2026-10-29/` contains
component derivations and checks. `joint/cross-review.md` independently checks
the SEM derivative/sign/scaling and bootstrap exchangeability/failure policy.
`bootstrap/doc-review.md` checks package-wide claims; all three suggested
clarifications were applied to NEWS, getting-started and bootstrap help.

`installed-review/` contains 14 independent installed-package acceptance
checks: enumerated binomial tails, Poisson refit arithmetic, offset/trial
invariants, unresolved-fit bounds, joint complement covariance, known-Gaussian
ECDF reduction and all three result/plot interfaces.

The installed fast suite passes **146 blocks / 1,153 assertions** without
unexpected warnings. The integrated run, including all calibration/power
checks, passes **169 blocks / 1,187 assertions**, of which **23 blocks / 34
assertions** are opt-in calibration/power checks. GNU/C++17 and Clang/C++23
package checks have zero errors, zero warnings and the single archived-CRAN
incoming note. Installed-interface checks resolve 18 help topics and three
built vignettes. Fourteen independent installed-candidate checks also pass.

| New calibration design | Observed Sup / L2 rejection or detection rate |
|:--|:--|
| Joint normal null, 300 cohorts | 0.060 / 0.0633; simultaneous envelope crossing 0.060. |
| Joint binomial null, 200 cohorts | 0.055 / 0.055. |
| Joint curved-binomial alternative, 200 cohorts | 1 / 1. |
| Specified sparse Poisson null, 200 cohorts | 0.055 / 0.045. |
| Specified grouped binomial null, 200 cohorts | 0.020 / 0.035. |
| Fitted low-count Poisson null, 120 cohorts | 0.0333 / 0.0333. |
| Fitted low-count binomial null, 120 cohorts | 0.0417 / 0.0417. |
| Bootstrap omitted-mean-term alternative, 60 cohorts | 1 / 1. |
| SEM latent/measurement strong-skew alternatives, 50 cohorts each | 1 / 1 for both projections. |

These bounded designs do not establish uniform finite-sample calibration or
universal power. The final machine-readable summary retains the original
calibration results and source parity, not only the new-method rates. Local
checks remain distinct from unexecuted remote platform/minimum-R CI.

The final shareable archive is in the course audit directory under
`qualification-release/`. Its network-enabled GNU check is authoritative.
Initial sandbox runs could not resolve public repository metadata; their
failed logs are retained. Clang's successful network-enabled replay checked
the exact same Clang archive. Parsed executable R code and exact test/native
files match the scientific run and delivery; later differences are documented
wording clarifications. Unchanged native sanitizer evidence is inherited from
the earlier baseline and is not represented as a rerun in this phase.

The initial SEM normality screen uses `p < .05`. With plus-one p-values and
`R=199`, equality at .05 is possible; inclusive `p <= .05` rates must also be
reported. In 200 cohorts of size 300, strict rates were 6--7.5%, whereas the
inclusive rates were 6%, 8.5%, 6.5%, 6.5%, 6% and 8.5% (marginal, observed
projection, latent projection; Sup then L2). The two 8.5% results warrant
additional assessment rather than treating a broad screening pass as exact
5% calibration. A fresh paired study was prespecified before drawing its
cohorts: 200 new samples, each with its first 300 and all 1,200 observations,
all three projections, `R=399` and primary inclusive `p <= .05`. No source,
test threshold or seed roster was changed in response to the results.

| Fresh normality mode | n=300 Sup / L2 | n=1,200 Sup / L2 |
|:--|:--|:--|
| Marginal pool | 0.030 / 0.030 | 0.050 / 0.050 |
| Observed structural projection y1 | 0.035 / 0.025 | 0.035 / 0.050 |
| Latent structural projection eta | 0.030 / 0.070 | 0.045 / 0.035 |

All 200 pairs and 1,200 endpoint results were retained, with no unresolved
fits. Every component-wise exact 95% binomial interval includes .05. All
paired-difference intervals include zero; for latent L2 the rate difference
is .035 with approximate 95% interval [-.0052, .0752] and exact McNemar
p=.143. These are not equivalence proofs or evidence of a significant
convergence trend. The broad excess did not persist in the fresh cohort;
no implementation defect was confirmed by the independent mathematical,
numerical or supplemental audit. Retain the endpoint as the paper's
asymptotic procedure, with both studies and finite-sample uncertainty
disclosed. Do not claim exact/universal 5% calibration or invent a universal
sample-size cutoff.

The frozen protocol, seed roster, installed-artifact hashes, all cohort data,
fit states and results, and the final statistical assessment are retained in
`normality-supplement/`. Installed help, README and the worked guide describe
the method's approximate calibration. No supplementary experiment replaced
the original evidence or changed the implementation.

The local release baseline is version **1.2.0**, tag **v1.2.0**. The course
audit directory's `delivery-manifest.json` binds the clean source commit,
archive SHA256, documentation and all scientific qualification evidence.

## C++ modernization boundary

Commit and tag the qualified capability version locally before native
restructuring. Preserve original authorship/license, public R contracts,
observed processes, statistic definitions, simulation coupling, promised
random-stream behavior and invalid-input handling. Assess intentional
statistical changes separately with the independent oracles and calibration
suite rather than bundling them into architecture work.

No push, publication, hosted deployment or change to the normal installed gof
library is part of this phase.
