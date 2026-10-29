## Low-level C++ helpers exposed to R.

test_that("SupTest is max |W|", {
    set.seed(51)
    w <- rnorm(50)
    expect_equal(gof:::SupTest(w), max(abs(w)))
})

test_that("L2Test is the left Riemann sum of W(t)^2 dt, square-rooted", {
    set.seed(52)
    t <- sort(runif(50)); w <- rnorm(50)
    expect_equal(gof:::L2Test(w, t), sqrt(sum(diff(t) * w[-50]^2)))
    expect_equal(gof:::L2Test(1, 1), 0)
})

test_that("L2Test errors (does not crash) on bad input", {
    expect_error(gof:::L2Test(c(1, 2, 3), c(1, 2)))
    expect_error(gof:::L2Test(numeric(0), numeric(0)))
})

## .fastapprox(time, newtime, equal, type) returns 1-based indices into the
## sorted vector `time`: type 0 = nearest, 1 = right (first time >= new),
## 2 = left (last time <= new). Values outside the range map to the ends.
test_that(".fastapprox returns 1-based nearest/right/left indices", {
    tt <- c(1, 2, 3, 4)
    new <- c(0.5, 2.4, 2.6, 3, 9)
    expect_equal(drop(gof:::.fastapprox(tt, new)), c(1, 2, 3, 3, 4))
    expect_equal(drop(gof:::.fastapprox(tt, new, type = 1)), c(1, 3, 3, 3, 4))
    expect_equal(drop(gof:::.fastapprox(tt, new[-1], type = 2)), c(2, 2, 3, 4))
})

test_that(".fastapprox agrees with a brute-force nearest neighbour", {
    set.seed(53)
    tt <- sort(runif(200)); new <- runif(100)
    brute <- vapply(new, function(v) which.min(abs(tt - v)), integer(1))
    expect_equal(drop(gof:::.fastapprox(tt, new)), brute)
})

test_that(".clusterindex gives 1-based start and size of each run", {
    id <- c(1, 1, 2, 2, 2, 5, 7, 7)
    rl <- rle(id)
    res <- gof:::.clusterindex(id)
    expect_equal(res[, 1], cumsum(c(1, head(rl$lengths, -1))))
    expect_equal(res[, 2], rl$lengths)
})
