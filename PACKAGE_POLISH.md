# Package maturity and statistical review, version 1.1.0

This is the historical 1.1.0 qualification record. Subsequent statistical
capability development and its current support boundaries are recorded in
[CAPABILITIES.md](CAPABILITIES.md).

This phase starts at `7fc29b8` (the 1.0.5 repair baseline). Its initial goal
was a shareable R package: readable installation/use instructions, installed
help, worked guides, conventional result extraction and reproducible
documentation. An independently assigned statistical review then identified
additional scientific defects and clarified the method's support boundaries.
Those findings are recorded separately in `STATISTICAL_REVIEW.org`.

The version is 1.1.0 because the work includes intentional changes to test
statistics and the ordinal default, not only documentation. The native C++
kernels are unchanged; their later C++23 architecture modernization remains
separate work.

## User-facing package changes

| Area | Delivered behavior |
|:--|:--|
| README | `README.Rmd` generates `README.md` and a real diagnostic figure. Source installation identifies this maintenance snapshot instead of directing users to an upstream version without the repairs. |
| Installed help | A generic `?cumres` landing page, adapter pages, shared `?"cumres-object"` help, plot options and a column-level data description. Public result-method aliases resolve in the installed package. |
| Worked guides | `getting-started` and `statistical-details` are built and installed as vignettes, with runnable seeded examples. |
| Results | `as.data.frame()` works on diagnostics and summaries. Optional generics/broom tidiers return tibbles. Existing `coef()` p-value tables remain available. |
| Pipes | Model-first calls work with base `|>`; extraction consumes no additional random draws. |
| Ordinal interpretation | Explicit `omnibus`, `score` and `cutpoints` modes; output tables retain `aggregation` so their interpretation survives extraction. |
| Namespace | `lava` is imported without automatic attachment; selective imports replace broad imports, and DLL entry points use registered symbols. The routine startup message is removed. |
| Version policy | R 4.1 is the declared minimum for native syntax. C++17 remains the package baseline. Minimum-R and platform CI jobs are configured, not represented as locally executed. |
| Citation and history | `citation("gof")` includes the software/version and Lin, Wei and Ying method reference. `NEWS.md` retains release history in a format recognized by `news()`. Original authorship and license remain intact. |
| Developer workflow | `DEVELOPMENT.md`, `CONTRIBUTING.md`, generated-document freshness checks, README checks and a local pkgdown site build. No deployment step is added. |

The result table has one row per process/test and columns `variable`, `test`,
`statistic`, `p.value`, `n`, `nobs`, `R`, `informative` and `aggregation`.
Missing metadata is distinguished from withheld p-values. Legacy matrix
diagnostics and saved summaries are covered; old summaries cannot recover
observed statistics that were never saved.

## Findings from the independent statistical review

The reviewer obtained full primary papers, inspected equations and ran
small counterexamples independently of the package's regression suite.
Passing a calibration experiment was not used as a substitute for matching
the mathematical method. Detailed evidence and primary PDF hashes are in
the course audit directory's `statistical-review/` record.

| Finding | Resolution and scientific meaning |
|:--|:--|
| F1: converged finite coefficients can conceal no regular mean root | Binomial/Poisson geometry checks reject separation/recession. Noncanonical links use numerical local interior-root qualification. Boundary and unresolved cases do not receive ordinary p-values. These are numerical qualifications, not formal interval proofs or global optimization guarantees. |
| F2: moving-window grid omits exit events and tail | Full grids now include entries and exits through the final zero tail, matching the paper's window domain. Observed and simulated statistics use the same full-prefix differences. A finite cap explicitly approximates that domain. |
| F3: pointwise centering is insufficient for estimated ordering callbacks | Help and guides state the missing ordering-derivative condition and a sufficient construction using admissible conditioning data. Arbitrary functions are not automatically justified; named safe SEM predictions retain their structural construction. |
| F4: ordinal score sum misses mean-preserving distribution changes | The default now takes the joint Euclidean norm of cumulative-cutpoint processes with shared subject multipliers. The continuous mapping argument justifies the transformation. The previous score-mean test remains explicit as `aggregate="score"` or `sum`. |
| F5: simulations from separate orderings are not joint draws | The public scope is marginal ordering-variable inference. Saved columns must not be combined for joint/global inference. This is stated in help and guides; no unimplemented joint test is advertised. |
| F6: saturated ordinal processes report significance from optimizer drift | Structure-based checks identify intercept-only and saturated binary CDF processes. Their raw CDF vectors are zeroed before aggregation, preserving a custom function's value at zero and the random-draw stream. General ordinal separation also has a sufficient numerical direction check. |

The ordinal norm combines cutpoints for **one ordering**. It does not turn
the package into a test of every conceivable model departure. Similarly,
the SEM adapter implements residual mean/linearity diagnostics, not the
paper's separate normality-diagnostic suite.

Gaussian identity exact fits retain the zero-influence, p-value-one case:
dispersion cancels from the mean influence and the unscaled mean Jacobian
can remain nonsingular. This differs from an absent finite mean-parameter
root for a boundary link. Constant SEM ordering is not blindly zeroed;
an independently checked fixed-mean control retains its genuine mean
rejection.

## Remaining statistical conditions

The method is an asymptotic multiplier approximation. Rare events, severe
imbalance, high leverage and little effective information can make it
unreliable despite a large nominal sample and numerical convergence. A
known Poisson mean of 1e-12 with 1,000 rows is a concrete example: almost
no events are expected, so an empirical residual multiplier does not give
exact small-count inference. No universal count cutoff is invented.
An exact or parametric procedure for that regime would be a separate method.

Bandwidth and ordering/grid-selection rules should be prespecified for
formal testing. Outcome-driven selection requires its own adjustment.
Custom residual and ordering constructions require the stated influence,
conditional-centering, differentiability and equicontinuity arguments.
Influence inputs use the scale where the coefficient perturbation is the
sum of row contributions, typically total-information inverse times score;
an order-one influence function must first be divided by sample size.

These conditions define the supported inference. An automated review is
not a human statistician's approval of a particular dataset or analysis.

## Final local results

* The installed fast suite passes **114 blocks / 913 assertions**.
* All **16 calibration blocks / 18 assertions** pass; the combined run is
  **130 blocks / 931 assertions**. Updated moving-window rejection rates are
  0.045/0.045 for Sup and 0.05/0.04 for L2; ordinal logit norm rates are
  0.05625/0.05625 and 0.0625/0.05 for tied/predicted orderings.
* Final GNU/C++17 and Clang/C++23 package checks, including rebuilt vignettes
  and PDF/HTML manuals, have **zero errors, zero warnings and one archived
  CRAN incoming note**.
* **13 independently replayed statistical acceptance checks pass** on the
  installed 1.1.0 candidate.
* Documentation generation, the README and figure, both site articles and
  all reference pages pass the strict local documentation gate. No Pandoc
  rendering warnings remain.
* Installed-interface checks resolve 15 help topics and both built vignettes.
  Fresh-process checks pass with `generics` loaded before or after `gof`;
  loading the package is quiet and does not attach `lava`. Both citations and
  26 structured NEWS versions are available.

The implementation/API commit is `628b84c`; documentation and workflow
changes are in `bce5ef9`. The final audit-note revision is sealed with the
archive hashes in the delivery manifest. These are local commits, not an
upstream release or remotely executed CI run.

## Validation and trace

Authoritative final results are collected under
`audits/gof-r-package-polish-2026-10-29/` in the course workspace:

* `qualification-final/` records the full fast/calibration suite on the
  scientific source, installed help/vignettes and full CRAN check.
* `qualification-delivery/` is the final distributable GNU source archive
  and package gate after the documentation-runner flush correction.
* `qualification-clang23/` records the C++23 integration variant; its only
  later source difference is the unrelated documentation-runner flush fix.
* `documentation-final-v3/` is the strict successful documentation gate and
  locally built website. It checks generated files, executable README,
  articles and absence of Pandoc rendering warnings.
* `statistical-review/` retains original counterexamples, primary sources,
  endpoint assessments and independent post-repair checks.
* `window-domain/`, `mle-validity/` and `ordinal-validity/` retain focused
  before/after evidence and the relevant mathematical arguments.

Intermediate builds are retained for provenance. The initial 1.0.6 check
found an obsolete generated package alias; the first site render found a
math-label problem; a later documentation gate exposed a log-flushing
helper error. These are corrected and superseded by the final gates.
The delivery manifest identifies exact source revisions and archive hashes.

Native sources and their standalone harness are unchanged from the passing
1.0.5 C++14/17/20/23 ASan/UBSan and separate LeakSanitizer runs. Their hashes
are compared with the delivered source instead of repeating identical tests.
New R statistical paths have focused regressions and renewed calibration.

No source was pushed, no site was published and no ordinary user R library
was replaced. A local archive and website preview are prepared for review.
Remote CI execution and any maintainer submission/publication remain
separate release actions.
