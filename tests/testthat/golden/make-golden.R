## Generates golden values for tests/testthat/test-golden.R.
##
## Run ONCE against the unmodified baseline (gof 1.0.1, commit 806286c) to pin
## the behaviour of every code path that already worked (n <= 500). The bug
## fixes must not change these values. Do not regenerate after changing the
## package, unless a change of these values is intended and documented in NEWS.
##
## Usage (from the package root, with the baseline installed):
##   Rscript tests/testthat/golden/make-golden.R

library(gof)
source("tests/testthat/helper-models.R")

run <- function(model, seed, ...) {
    set.seed(seed)
    g <- cumres(model, R = 200, plots = 10, ...)
    g[c("W", "What", "x", "KS", "CvM", "variable")]
}

golden <- list(
    lm        = run(fx_lm(), 11),
    lm_resp   = run(fx_lm(), 12, variable = c("y", "x1")),
    logistic  = run(fx_logistic(), 13),
    probit    = run(fx_probit(), 14),
    poisson   = run(fx_poisson(), 15),
    grouped   = run(fx_grouped(), 16)
)
attr(golden, "baseline") <- list(version = as.character(packageVersion("gof")),
                                 commit = "806286c")
saveRDS(golden, "tests/testthat/golden/baseline-1.0.1.rds", version = 2)
