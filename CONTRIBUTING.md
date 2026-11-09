# Contributing to gof

This maintenance fork is developed at [Marimuda/gof](https://github.com/Marimuda/gof).
Report fork-specific issues in its [issue tracker](https://github.com/Marimuda/gof/issues).
The [original upstream project](https://github.com/kkholst/gof) has a separate
release history. Preserve original authorship, license and citations.

## Reporting a problem

Include a reproducible example, the expected result, `sessionInfo()` and the
package version or source revision. For statistical issues, state the model,
residual unit, weights/offsets, ordering variable, simulation count and seed.
Use synthetic data in place of confidential research data.

Explain whether the evidence comes from an analytic calculation, a literature
equation, an independent implementation or a simulation. Passing a regression
test does not settle a question about a method's assumptions; a Monte Carlo
difference alone does not establish an implementation error.

## Development setup

Work from the repository root with R >= 4.1 and a C++17 toolchain. Install
development dependencies in a library you control:

```r
install.packages(c(
  "remotes", "roxygen2", "pkgdown", "knitr", "rmarkdown", "testthat",
  "numDeriv", "generics", "tibble", "broom"
))
remotes::install_deps(".", dependencies = TRUE)
```

Use the roxygen2 version recorded in `Config/roxygen2/version` in `DESCRIPTION`.
If needed, install that version with `remotes::install_version()`. Review any
generator upgrade together with all regenerated output. README freshness uses
Pandoc 2.9.2.1; the package-guide CI jobs use Pandoc 3.1.11.1.

Pandoc is required for rendering. Full package checks also need LaTeX and
qpdf for the PDF manual, and HTML Tidy and V8 for HTML validation.
[tools/QUALIFICATION.txt](https://github.com/Marimuda/gof/blob/release/1.3.0/tools/QUALIFICATION.txt) lists requirements and
environment controls for the validation tools.

## Documentation

Edit function documentation in the roxygen blocks under `R/`. Commit generated
`man/*.Rd` files and `NAMESPACE`. Edit `README.Rmd` and regenerate `README.md`
and its figure; do not edit the generated README directly.

```r
roxygen2::roxygenise(".", roclets = c("rd", "namespace"))
pkgload::load_all(".")
rmarkdown::render("README.Rmd", quiet = TRUE)
```

Keep user-visible changes in `NEWS.md`, statistical assumptions and worked
examples in help pages and vignettes, and compact reproducible evidence beside
its guide. Keep working plans, progress logs and local audit receipts outside
the checkout. The ignored `.local/` directory is available for temporary notes;
it is also excluded from source archives.

Run the documentation check in a fresh R process, using a new output directory
outside the checkout:

```sh
Rscript tools/check-documentation.R ../gof-documentation-new
```

This checks generated help and README freshness, then builds a local pkgdown
site with runnable examples and guides. It retains a source copy, logs and a
manifest. Set `GOF_DOCS_BUILD_SITE=false` to omit the site explicitly. README
figure names and nonempty content are checked; device-dependent PNG bytes are
not compared. The checker does not publish the site or replace the installed
package.

## Changes and validation

Keep changes focused and preserve existing calls where practical. Document
intentional numerical or support changes in `NEWS.md`, with a regression or
reference test for the actual problem. Do not replace a failing reference
fixture merely to make a test pass.

R-interface changes should exercise ordinary calls, native pipes, saved
objects and missing or indeterminate information where applicable. Changes
to statistical algorithms, projection, tie handling or simulation generation
also require suitable independent reference checks and calibration. Preserve
observed processes, statistic definitions and promised random-stream behavior
when restructuring implementation code.

Run the full package qualification when appropriate:

```sh
Rscript tools/qualify.R ../gof-qualification-new
```

The runner builds a source archive with guides, installs it in an isolated
library, checks installed R and C++ interfaces, runs fast regressions and
native sanitizers, and performs `R CMD check --as-cran`, including the manual.
Set `GOF_SLOW_TESTS=true` to include statistical calibration. Output directories
must be new; evidence is retained on failure. Calibration experiments support
their declared designs and do not prove validity for arbitrary inputs.

## Native code and compatibility

Keep `SystemRequirements: C++17` in `DESCRIPTION` consistent with
`CXX_STD = CXX17` in `src/Makevars`. C++23 is an additional integration
configuration. Any increase in the R or C++ minimum needs a concrete benefit
and a compatibility assessment covering supported platforms, dependencies
and installed-header consumers. Check the actual compiler features and
`__cplusplus` value rather than inferring support from a flag alone.

Use the committed `.clang-format` and `.clang-tidy` policies. Run static
analysis with matching Clang tools (qualified with LLVM 19), Python 3 and
a completed Clang package installation log:

```sh
python3 tools/check-static.py \
  --source ../qualification/source \
  --install-log ../qualification/logs/install.log \
  --output ../native-static-new
```

Use a raw source snapshot, which retains development configuration excluded
from the R archive. The tool derives compiler commands from the installation
log, checks formatting and installed headers, and runs clang-tidy and the
Clang Static Analyzer. Generated Rcpp bindings are excluded from owned-code
analysis. Diagnostics and source hashes are retained without rewriting code.

Analyzer and portability findings fail the gate, except for the opt-in
`UninitializedObject` checker, whose Armadillo auxiliary-field findings remain
advisory. Other enabled tidy findings also require review. A successful exit
does not mean zero findings. Native changes additionally need sanitizer,
numerical/RNG replay and installed-client checks; static analysis does not
establish statistical correctness.

## CI and releases

[GitHub Actions](https://github.com/Marimuda/gof/actions/workflows/qualification.yaml)
checks R versions and platforms, C++23 integration, documentation, calibration,
static analysis and native sanitizers. Inspect the actual run before claiming
platform support. The workflow does not publish packages or deploy a website.

For release provenance, use the annotated tag and its matching CI run. The
[`v1.3.0` source](https://github.com/Marimuda/gof/tree/v1.3.0) is commit
`cecbcba6a381307aa6d8199a02749973637ebbc1`; its
[qualification run](https://github.com/Marimuda/gof/actions/runs/41707803654)
passed all 21 jobs. That result applies to the tagged source. Subsequent
changes require their own checks.

Before distributing an archive, inspect its installed help and guides, verify
its source revision and hash, and retain the qualification manifest. Distinguish
fork releases from upstream releases and CRAN publication. Software checks and
empirical statistical evidence have separate scopes; retain known limitations
and adverse results in the public documentation.
