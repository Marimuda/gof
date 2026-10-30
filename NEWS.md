# gof 1.2.1

This is a maintenance release candidate. The source checkout records platform
sign-off in `RELEASE_STATUS.md`; the final release tag awaits that sign-off.

* Separated cumulative kernels, randomness and utility components while
  retaining the installed native interfaces and C++17 package requirement.
* Added compiler-backed static checks, header compatibility tests and
  consistent native formatting. Fixed a narrowing boundary in R normal draws
  and avoided two uninitialized-argument paths in lazy matrix expressions.
  Valid seeded results are checked against the tagged baseline.
* Corrected the SEM option-test reference to differentiate lava's aggregate
  score. Differentiating a separately summed score amplified rounding
  differences on macOS. Existing comparison tolerances and package
  calculations are unchanged.

# gof 1.2.0

* Added joint multiplier diagnostics across prespecified ordering variables,
  family Sup and range-normalized L2 tests, and simultaneous null reference
  envelopes. Legacy marginal simulations retain their seeded behavior.
* Added Gaussian structural-equation residual-distribution diagnostics with
  parameter-estimation correction, separate from mean/linearity checks.
* Added canonical binomial/Poisson parametric calibration, distinguishing
  fitted-null bootstrap refitting from specified-null Monte Carlo testing.
  Failed replicates remain visible. Unresolved refits withhold informative
  p-values; certified zero-information processes retain p-values of one.
* Added installed help and a worked capability guide explaining the methods,
  numerical qualification, finite-grid scope and inferential limitations.

# gof 1.1.0

* Added a generated README, installed workflow and statistical-detail
  vignettes, a locally buildable documentation site and complete result help.
* Added base-R data-frame extraction and optional generics/broom tidiers for
  diagnostic and summary objects. Tables expose process labels, observed
  statistics, Monte Carlo p-values and sample/simulation metadata explicitly.
  Existing `coef()` p-value tables remain available.
* Native `|>` workflows are documented and tested. R 4.1 is now the declared
  minimum version; the configured CI matrix includes that baseline.
* Updated package/help/data documentation, scientific/software citations and
  source installation instructions for this unpublished research snapshot.
* `lava` is now imported without being attached by `library(gof)`. Call
  `library(lava)` explicitly when using its model-building functions.
* Removed the routine package startup message.
* Statistical review identified additional validity failures: moving-window
  diagnostics now include entry and exit events over their full support;
  GLM inference checks separation/recession and numerical interior mean-fit
  regularity before computing p-values.
* Ordinal diagnostics now default to a joint Euclidean norm across cumulative
  cutpoints for the selected ordering variable. This detects alternatives
  that preserve the ordinal score mean. Use `aggregate="score"` or
  `aggregate=sum` for the former mean-only test, and `aggregate="cutpoints"`
  for signed component diagnostics. Tables retain the aggregation mode.
* Documentation generation, README freshness, vignette rendering and local
  website builds have validation tools and nonpublishing CI checks.

# gof 1.0.5
* Extended qualification to every native utility, public R method and
   documented adapter option. AUDIT.org tracks the defects, evidence, fixes,
   support boundaries and validation status before C++ modernization.
* Fixed vector softmax returning log probabilities, complex-step underflow,
   unsigned degrees of freedom, grouping/pattern buffer edge cases, 64-bit
   cluster labels, L2 intermediate overflow and mutable cumulative-process state.
* Added projection-free processes for models with no fitted coefficients.
* Added normal/Poisson response-mean precision weights and integer binomial
   trial weights, using the original fit and excluding zero-weight rows.
* Reconstructed predicted coordinates from coefficients to preserve design
   ties. Structurally saturated processes have zero curves, p-value one and
   informative=FALSE metadata; numerical score-space ambiguity returns missing
   p-values with the computed curves preserved.
* Made data/variable selection follow fitted rows; fixed generic residual and
   gradient callbacks and added explicit finite real-value/index contracts.
* Replaced unsafe SEM conditioning-set inversion with structural covariance
   paths. Fixed full prediction, custom p/data, formula/list/function handling,
   ordinal shared-state restoration and invalid fitted-model domain checks.
* Preserved historical private R aliases while correcting generated C++
   callable names; added downstream installed-header client qualification.
* Fixed plotting palettes, graphical arguments, legends, legacy matrix
   objects and coefficient tables without statistics.
* Added isolated qualification/coverage tooling and a nonpublishing OS/R CI
   matrix; removed obsolete Travis deployment and notification configuration.

# gof 1.0.4
* Tied covariates are evaluated at complete tie groups, without numerical
   jitter. This fixes severe false rejection for grouped data at large
   covariate scales. Results for tied data intentionally change; the frozen
   1.0.1 golden fixture is retained as historical evidence, not regenerated.
* Grouped binomial diagnostics retain the actual fitted coefficients and
   covariance rather than refitting with different controls or methods.
* Fixed na.exclude alignment, plots=0, odd/single-point subset.max limits,
   custom subset selectors, named plot indices and col=NULL.
* Restored moving-window residual sums over (t-b,t], including simulation
   and exact lower endpoints on reduced evaluation grids.
* Ordinal regression now returns consistent list-based processes and has
   shared print/plot support, L2 tests and tested category aggregation.
   Added structural-equation-model reference and interface tests.
* Monte Carlo comparisons include equality, so zero-width L2 integrals
   have p-value one. Added count, dimension and finite-value validation.
* Fixed native approximation boundaries, tied left endpoints, equality
   flags, empty clusters and integer arithmetic in CramerVonMises.
* Removed an unused suggested logger and a committed testthat cache.
   Updated the required testthat version to match the test API in use.

# gof 1.0.3
* Fixed crash in cumres() for more than 500 observations, 'number of
   items to replace is not a multiple of replacement length' (#1).
   When it did not crash, the observed process was shifted by one
   observation relative to the simulated processes.
* Fixed off-by-one in the selection of evaluation points for more than
   500 observations (equidistant grid points).
* Fixed cumres() for grouped binomial responses, cbind(successes,
   failures): equal group sizes, models fitted with subset= or an offset,
   and models fitted inside functions failed. Variables are now reported
   and selected by their original names (previously e.g. 'X0x').
* Offsets are now included in the linear predictor used for the
   simulated null processes. Previously p-values were invalid for models
   with an offset (e.g. Poisson rate models).
* Regenerated RcppExports.cpp (failed to compile with
   -Werror=format-security).
* Test suite: golden values from 1.0.1, a pure-R reference
   implementation compared with the C++ engine, regression tests, and
   opt-in Monte Carlo calibration tests (GOF_SLOW_TESTS=true).

# gof 1.0.2
* New c++ functions with R-wrappers:
   .fastapprox, .clusterindex, ...
* Removed mets dependency
* More unit tests (c++ code only)

# gof 1.0.1
* Optimization of C++ code
* Some compatibility with Solaris fixed

# gof 1.0.0
* Faster implementation in RcppArmadillo with improved algorithms
* Method for proportional hazard model (coxph) removed from this
   version (methods are available in the timereg package)

# gof 0.9.1
* Maintenance release

# gof 0.9.0
* Update to work with new lava version  
* Patched scythestat (clang++ compatible)
* Help-pages updated

# gof 0.8.2
* Updated to new scythestat version

# gof 0.7.0
* Fixed bug introduced by new model.matrix implementation for 'coxph' objects
* Initial support added for structural equation models (lava)

# gof 0.6.4
* Fixed bug in cumres.coxph appearing for models with only one
   predictor (stupid [...,drop=FALSE] missing)
* CITATION added

# gof 0.6.3
* Windows compilation fix
* onload function added

# gof 0.6.2
* Code cleanup

# gof 0.6.1
* bug appearing with occurences of ties in x-variable 
   (ordering variable). Appeared because usage of 'order'-function
   from the scythestat-library which gives a ranking (e.g. 4 1 1 1 5
   compared to 4 1 2 3 5).

# gof 0.6.0
* Nicer cumres function (Wi)

# gof 0.5.7
* Compiles with '-pedantic' (scythestat library)

# gof 0.5.6
* g++ 4.3 compatibility (scythestat library)

# gof 0.5.5
* Moving average of residuals implemented

# gof 0.5.4
* Minor revision of code (no changes for end-user)
* Major bug corrected. For a glm dginv(yhat) was used 
   instead of dginv(Xb).

# gof 0.5.3
* Bug fixes in particularily in the cumres method for the coxph-class.
* Documentation

# gof 0.5.0
* Cumres changed to S3 method
* Added support for the Cox model
* Large revision of C++ code
	
# gof 0.4.0
* Implementation of simultaneous prediction bands
* Documentation

# gof 0.3.0
* Bug-fixes

# gof 0.2.0
* Transition to C++ code for faster computation of resamples

# gof 0.1.0
* Implementation of linear, logistic and poisson regression

