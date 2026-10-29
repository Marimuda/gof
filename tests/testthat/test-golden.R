## Regression net: results for code paths that worked in gof 1.0.1 must not
## change. Golden values were generated with golden/make-golden.R against the
## unmodified baseline (commit 806286c).

golden <- readRDS(test_path("golden", "baseline-1.0.1.rds"))

run <- function(model, seed, ...) {
    set.seed(seed)
    cumres(model, R = 200, plots = 10, ...)
}

expect_golden <- function(g, ref) {
    expect_equal(g$KS, ref$KS)
    expect_equal(g$CvM, ref$CvM)
    expect_equal(g$W, ref$W)
    expect_equal(g$What, ref$What)
    expect_equal(g$x, ref$x)
}

test_that("lm results are unchanged from 1.0.1", {
    g <- run(fx_lm(), 11)
    expect_golden(g, golden$lm)
    expect_equal(g$variable, golden$lm$variable)
})

test_that("lm with the response as variable is unchanged from 1.0.1", {
    g <- run(fx_lm(), 12, variable = c("y", "x1"))
    expect_golden(g, golden$lm_resp)
    expect_equal(g$variable, golden$lm_resp$variable)
})

test_that("logistic, probit and poisson glm results are unchanged from 1.0.1", {
    g <- run(fx_logistic(), 13)
    expect_golden(g, golden$logistic)
    expect_equal(g$variable, golden$logistic$variable)
    g <- run(fx_probit(), 14)
    expect_golden(g, golden$probit)
    g <- run(fx_poisson(), 15)
    expect_golden(g, golden$poisson)
})

test_that("the legacy grouped golden remains a record of the corrected tie bug", {
    g <- run(fx_grouped(), 16)
    ## Do not regenerate the frozen baseline: its partial sums WITHIN groups
    ## are invalid evaluations of sum(1{x<=t} r). The corrected grouped path
    ## is covered by the independent indicator oracle and binomial tests.
    expect_equal(lengths(g$W),c(8L,8L))
    expect_gt(length(golden$grouped$W[[2]]),length(g$W[[2]]))
    expect_equal(golden$grouped$variable, c("predicted", "X0x"))
    expect_equal(g$variable, c("predicted", "x"))
})
