# Developing and checking gof

The public interface is documented in roxygen blocks under `R/`. Commit the
generated `man/*.Rd` files and `NAMESPACE`; users need no documentation tools
to install a built package. The statistical repair baseline is recorded in
`AUDIT.org`. `CAPABILITIES.md` records the subsequent capability contracts,
review findings and qualification. User-facing changes belong in `NEWS.md`.

## Development dependencies

Work from the repository root with R 4.1 or later and a C++17 toolchain.
Install development dependencies into a library you control:

```r
install.packages(c(
  "remotes", "roxygen2", "pkgdown", "knitr", "rmarkdown", "testthat",
  "numDeriv", "generics", "tibble", "broom"
))
remotes::install_deps(".", dependencies = TRUE)
```

Use the roxygen2 version recorded in `Config/roxygen2/version` in
`DESCRIPTION` when checking generated-file freshness. If a different version
is installed, select the recorded version deliberately with
`remotes::install_version()`. Updating the generator is a separate reviewed
change: regenerate the files and inspect the diff together.

Pandoc is needed to render the README and HTML vignettes. A full PDF manual
check also needs LaTeX; HTML validation uses HTML Tidy and the V8 R package.
Native sanitizer checks additionally need Bash, system Armadillo development
files and compatible sanitizer runtimes. `tools/QUALIFICATION.txt` documents
these requirements and the available environment controls.

## Updating documentation

```r
roxygen2::roxygenise(".", roclets = c("rd", "namespace"))
pkgload::load_all(".")
rmarkdown::render("README.Rmd", quiet = TRUE)
```

`README.Rmd` generates `README.md` and its figures under `man/figures/`.
Commit all of them. Keep installation examples explicit about which source
archive or checkout they install; the upstream CRAN/GitHub package does not
contain unpublished local changes.

Run the documentation gate in a fresh terminal R process:

```sh
Rscript tools/check-documentation.R ../gof-documentation-new
```

It copies the source, regenerates and compares help/namespace files, renders
and compares the README text, checks figure inventory, then builds a local
pkgdown site with runnable examples and articles. It never publishes the
site or modifies the checkout. Its manifest records any deliberately
omitted phase. PNG byte identity is not required across graphics devices;
generated file names and nonempty contents are checked.

## Package qualification

```sh
Rscript tools/qualify.R ../gof-qualification-new
```

The runner builds a source archive including the vignettes, installs it in
an isolated library, checks installed R help and pipeline/table methods,
checks the installed C++ client interface, runs fast regressions, runs the
native sanitizer matrix and performs `R CMD check --as-cran`, including the
manual. Output directories must be new; evidence is retained on failure.

For statistical calibration as well, set `GOF_SLOW_TESTS=true` before
running the command. The calibration designs are finite empirical checks,
not a proof for arbitrary inputs. Changes to statistical algorithms,
projection, tie handling or simulation generation require appropriate
reference checks and calibration; documentation-only changes should not
alter seeded numerical results.

The optional development packages do not become core runtime dependencies.
`as.data.frame()` and native pipes use base R; tidy methods are registered
with `generics` only when it is available and return tibbles when requested.

## CI and release scope

The workflow checks the declared R minimum, current R release/devel/oldrel,
available operating systems, a C++23 integration copy, calibration,
documentation and native sanitizers. It has no package publication or
website deployment step. Writing the workflow is not evidence that a
platform job has passed; inspect the actual run before making that claim.

Before sharing an archive, check its installed help and vignettes, verify
the source revision and archive hash, and retain the qualification manifest.
`citation("gof")` provides the software and method references. Preserve
original authorship and distinguish maintenance work from upstream releases.
