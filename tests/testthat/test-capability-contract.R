test_that("unsupported adapters cannot silently accept joint ordering simulation", {
    set.seed(1210)
    model <- lava::lvm(y~x)
    fit <- lava::estimate(model,lava::sim(model,80))
    expect_error(cumres(fit,y~x,simulation="joint",R=3,plots=0),
                 "SEM collections remain separate")
    expect_error(cumres(ordinal_fixture(),simulation="joint",R=3,plots=0),
                 "one ordering")
})
