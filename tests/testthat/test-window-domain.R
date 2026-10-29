known_mean_window <- function(x,r,b,...) {
    cumres(NULL,variable=cbind(x=x),r=r,
           dr=matrix(numeric(),length(r),0L),ic=matrix(numeric(),length(r),0L),
           coef=numeric(),b=b,...)
}

test_that("full windows include exit jumps, interior peaks and the zero tail", {
    set.seed(2011)
    g <- known_mean_window(c(0,1,3),c(-8,10,-2),b=2,R=29,plots=3,subset.max=Inf)
    expect_equal(g$x[[1]],c(0,1,2,3,5))
    expect_equal(g$W[[1]],c(-8,2,10,-2,0)/sqrt(3))
    expect_equal(gof:::SupTest(g$W[[1]]),10/sqrt(3))
    expect_equal(gof:::L2Test(g$W[[1]],g$x[[1]]),sqrt(176/3))
    expect_true(all(g$What[[1]][5,]==0))
    expect_equal(dim(g$What[[1]]),c(5L,3L))
})

test_that("window boundaries accumulate all ties and exclude exiting observations", {
    g <- known_mean_window(c(0,0,1,1,3),c(-4,-4,5,5,-2),b=2,
                           R=11,plots=0,subset.max=NULL)
    expect_equal(g$x[[1]],c(0,1,2,3,5))
    expect_equal(g$W[[1]],c(-8,2,10,-2,0)/sqrt(5))
    expect_equal(ncol(g$What[[1]]),0L)
    ## 0.1+0.2-0.2 need not equal 0.1 in binary arithmetic. Use stored entry
    ## and exit events, so subtraction roundoff cannot keep a residual inside
    ## its window at the exact recorded exit.
    decimal <- known_mean_window(c(.1,.4),c(3,-1),b=.2,R=5,plots=0,subset.max=Inf)
    own_exit <- match(.1+.2,decimal$x[[1]])
    expect_equal(decimal$W[[1]][own_exit],0)
    singleton <- known_mean_window(.1,2,b=.2,R=5,plots=0,subset.max=Inf)
    expect_equal(singleton$x[[1]],c(.1,.1+.2))
    expect_equal(singleton$W[[1]],c(2,0))
    expect_equal(gof:::L2Test(singleton$W[[1]],singleton$x[[1]]),sqrt(.8))
})

test_that("observed and simulated full-domain windows match an independent indicator oracle", {
    m <- fx_logistic(n=80)
    inp <- glm_inputs(m)
    variables <- cbind(predicted=fitted(m),x=model.matrix(m)[,"x"])
    for (cap in c(Inf,19L)) {
        set.seed(2012)
        actual <- cumres(m,R=29,plots=3,b=.7,subset.max=cap)
        set.seed(2012)
        oracle <- ref_cumres(inp$r,inp$dr,inp$ic,variables,
                             R=29,plots=3,b=.7,subset.max=cap)
        expect_equal(actual$x,oracle$x)
        expect_equal(actual$W,oracle$W)
        expect_equal(actual$What,oracle$What)
        expect_equal(actual$KS,oracle$KS)
        expect_equal(actual$CvM,oracle$CvM)
    }
})

test_that("window grid caps and custom selectors operate on all entry and exit events", {
    received <- NULL
    select <- function(t,n.qt,n.eq) {received <<- t; c(2L,4L,6L)}
    g <- known_mean_window(1:5,c(1,-2,3,-4,2),b=2,R=5,plots=0,
                           subset.max=3,subset=select)
    expect_equal(received,1:7)
    expect_equal(g$x[[1]],c(2,4,6))
    expect_equal(g$W[[1]],c(-1,-1,2)/sqrt(5))
    single <- known_mean_window(1:5,c(1,-2,3,-4,2),b=2,R=5,plots=0,subset.max=1)
    expect_equal(single$x[[1]],7)
    expect_equal(single$W[[1]],0)
    expect_equal(c(single$KS,single$CvM),c(1,1))
    expect_error(known_mean_window(1:5,1:5,b=2,R=5,subset.max=3,
                                   subset=function(...) 8L),"subset")
})

test_that("full window processes preserve row and coordinate translation invariance", {
    set.seed(2013)
    x <- seq_len(40); r <- rnorm(40); permutation <- sample.int(40)
    set.seed(2014); original <- known_mean_window(x,r,b=7,R=19,plots=2,subset.max=Inf)
    set.seed(2014); permuted <- known_mean_window(x[permutation],r[permutation],b=7,
                                                 R=19,plots=2,subset.max=Inf)
    expect_identical(original$W,permuted$W)
    expect_identical(original$What,permuted$What)
    expect_identical(original$KS,permuted$KS)
    expect_identical(original$CvM,permuted$CvM)
    set.seed(2014); translated <- known_mean_window(x+1024,r,b=7,R=19,plots=2,subset.max=Inf)
    expect_equal(translated$x[[1]],original$x[[1]]+1024)
    expect_identical(original$W,translated$W)
    expect_identical(original$What,translated$What)
    expect_identical(original$KS,translated$KS)
    expect_identical(original$CvM,translated$CvM)
})

test_that("large full and reduced window grids retain the independent-unit sample", {
    set.seed(2015)
    n <- 650L; x <- seq_len(n)/2; r <- rnorm(n)
    for (cap in c(Inf,51L)) {
        set.seed(2016)
        actual <- known_mean_window(x,r,b=.75,R=11,plots=1,subset.max=cap)
        set.seed(2016)
        oracle <- ref_cumres(r,matrix(numeric(),n,0L),matrix(numeric(),n,0L),
                             cbind(x=x),R=11,plots=1,b=.75,subset.max=cap)
        expect_equal(actual$n,n)
        expect_equal(actual$x,oracle$x)
        expect_equal(actual$W,oracle$W)
        expect_equal(actual$What,oracle$What)
        expect_equal(actual$KS,oracle$KS)
        expect_equal(actual$CvM,oracle$CvM)
        expect_lte(length(actual$x[[1]]),min(cap,2*n))
    }
})

test_that("unrepresentable moving-window event coordinates fail explicitly", {
    expect_error(known_mean_window(.Machine$double.xmax,1,b=.Machine$double.xmax,
                                   R=1,plots=0),"exit coordinates must be finite")
    expect_error(known_mean_window(1e12,1,b=1e-9,R=1,plots=0),
                 "bandwidth cannot be resolved")
})

test_that("selected empty windows are exactly uninformative with or without fitted coefficients", {
    set.seed(2017)
    d <- data.frame(x=seq_len(40),y=rnorm(40))
    for (model in list(glm(y~x,d,family=gaussian()),lm(y~0,d))) {
        g <- cumres(model,variable="x",data=d,b=.5,R=7,plots=2,subset.max=1)
        expect_identical(g$informative,FALSE)
        expect_equal(g$x[[1]],max(d$x+.5))
        expect_equal(g$W[[1]],0)
        expect_true(all(g$What[[1]]==0))
        expect_equal(c(g$KS,g$CvM),c(1,1))
    }
    nonempty <- cumres(lm(y~0,d),variable="x",data=d,b=.5,R=7,plots=0,subset.max=Inf)
    expect_identical(nonempty$informative,TRUE)
})
