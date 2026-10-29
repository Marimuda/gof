## The C++ engine (via cumres.default) must agree exactly with the pure-R
## reference implementation in helper-reference.R, with the same seed.

compare_engine <- function(model, variable, r, dr, ic, seed,
                           subset.max = 500L, R = 60, plots = 4) {
    set.seed(seed)
    g <- gof:::cumres.default(NULL, variable = variable, r = r, dr = dr,
                              ic = ic, coef = NULL, R = R, plots = plots,
                              subset.max = subset.max)
    set.seed(seed)
    ref <- ref_cumres(r, dr, ic, variable, R = R, plots = plots,
                      subset.max = subset.max)
    expect_equal(g$x, ref$x)
    expect_equal(g$W, ref$W)
    expect_equal(g$What, ref$What)
    expect_equal(g$KS, ref$KS)
    expect_equal(g$CvM, ref$CvM)
    invisible(g)
}

compare_glm <- function(model, vars, seed, ...) {
    inp <- glm_inputs(model)
    variable <- cbind(predicted = fitted(model),
                      model.matrix(model)[, vars, drop = FALSE])
    compare_engine(model, variable, inp$r, inp$dr, inp$ic, seed, ...)
}

test_that("engine matches reference, all points (n <= 500)", {
    compare_glm(fx_logistic(), "x", seed = 21)
    compare_glm(fx_poisson(), "x", seed = 22)
    compare_glm(fx_probit(), "x", seed = 23)
})

test_that("engine matches reference with several variables (n <= 500)", {
    m <- fx_lm()
    variable <- cbind(predicted = fitted(m), model.matrix(m)[, c("x1", "x2")])
    compare_engine(m, variable, r = residuals(m), dr = -model.matrix(m),
                   ic = lava::iid(m), seed = 24)
})

test_that("engine matches reference on the subset path (n > 500, issue #1)", {
    g <- compare_glm(fx_large_ties(), c("group", "z"), seed = 25)
    expect_equal(lengths(g$W), c(500L, 500L, 500L))
})

test_that("engine matches reference with subsetting switched off (n > 500)", {
    g <- compare_glm(fx_large_ties(), c("group", "z"), seed = 26, subset.max = Inf)
    expect_equal(lengths(g$W), c(982L, 982L, 982L))
})
