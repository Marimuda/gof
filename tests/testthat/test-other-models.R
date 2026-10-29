test_that("ordinal diagnostics use the shared print/plot result contract", {
    for (k in 2:3) {
        m <- ordinal_fixture(k)
        set.seed(124); g <- cumres(m,variable="x",R=20,plots=3)
        expect_s3_class(g,"cumres")
        expect_equal(lengths(g$W),lengths(g$x))
        expect_equal(dim(g$What[[1]]),c(length(g$x[[1]]),3L))
        expect_true(all(is.finite(c(g$KS,g$CvM))))
        expect_no_error(capture.output(print(g)))
        expect_no_error(capture.output(summary(g)))
        pdf(NULL)
        expect_no_error(plot(g))
        dev.off()
        expect_no_error(cumres(m,R=5,plots=0,plot=FALSE))
    }
})

test_that("ordinal category aggregation matches a direct indicator reference", {
    m <- ordinal_fixture()
    env <- m$up
    x <- env$X[,"x"]; n <- length(x); k <- length(env$alpha)
    residual <- outer(env$y,seq_len(k),"<=")-predict(m,type="cumulative")
    grad <- -numDeriv::jacobian(function(p) predict(m,p=p,type="cumulative"),coef(m))
    ic <- lava::iid(m)
    ag <- function(z) c(sum(z),z[1]-z[2])
    set.seed(125); g <- cumres(m,variable="x",R=10,plots=10,aggregate=ag)
    expected <- t(vapply(sort(unique(x)),function(t)
        ag(colSums(residual*(x<=t)))/sqrt(n),numeric(2)))
    expect_equal(g$W,lapply(1:2,function(i) expected[,i]))
    set.seed(125)
    o <- order(x)
    for (draw in 1:10) {
        multiplier <- rnorm(n)
        beta <- colSums(ic[o,,drop=FALSE]*multiplier)
        reference <- t(vapply(sort(unique(x)),function(t) {
            z <- vapply(seq_len(k),function(j) {
                use <- x[o]<=t
                derivative <- grad[seq_len(n)+(j-1L)*n,,drop=FALSE][o,,drop=FALSE]
                (sum(residual[o,j]*multiplier*use)+sum(colSums(derivative*use)*beta))/sqrt(n)
            },numeric(1))
            ag(z)
        },numeric(2)))
        for (j in 1:2) expect_equal(g$What[[j]][,draw],reference[,j])
    }
})

test_that("lvmfit runs and its processes match the independent reference", {
    set.seed(130)
    m <- lava::lvm(list(c(y1,y2,y3)~eta,eta~x)); lava::latent(m) <- ~eta
    d <- lava::sim(m,100)
    e <- lava::estimate(m,d)
    p <- lava::pars(e)
    rf <- function(p) predict(e,lava::vars(e),residual=TRUE,p=p)[,"y1"]
    dr <- numDeriv::jacobian(rf,p,method=lava::lava.options()$Dmethod)
    set.seed(131); g <- cumres(e,"y1",x="x",R=20,plots=3)
    set.seed(131); ref <- ref_cumres(rf(p),dr,lava::iid(e),cbind(d$x),R=20,plots=3)
    expect_equal(g$W,ref$W)
    expect_equal(g$What,ref$What)
    expect_equal(g$KS,ref$KS)
    expect_equal(g$CvM,ref$CvM)
    pdf(NULL); on.exit(dev.off())
    expect_no_error(plot(g))
    expect_no_error(capture.output(print(g)))
    expect_no_error(cumres(e,c("y1","y2"),x=function(p) d$x,R=5,plots=0))
})
