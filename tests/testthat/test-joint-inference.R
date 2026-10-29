## Independent indicator-matrix oracle: no cumulative sums or native kernels.
joint_oracle <- function(r,dr,ic,variable,grid,R,plots,b=0) {
    n <- length(r)
    selectors <- lapply(seq_len(ncol(variable)),function(j)
        outer(grid[[j]],variable[,j],function(t,x) x<=t & (b==0 | x+b>t)))
    observed <- lapply(selectors,function(a) drop(a%*%r)/sqrt(n))
    one <- function() {
        g <- rnorm(n)
        perturbation <- if (ncol(ic)) drop(crossprod(ic,g)) else numeric()
        row_perturbation <- r*g+drop(dr%*%perturbation)
        lapply(selectors,function(a) drop(a%*%row_perturbation)/sqrt(n))
    }
    statistic <- function(w) unlist(lapply(seq_along(w),function(j)
        c(max(abs(w[[j]])),sqrt(sum(c(diff(grid[[j]]),0)*w[[j]]^2)))))
    null <- t(vapply(seq_len(R),function(i) statistic(one()),numeric(2L*ncol(variable))))
    what <- lapply(grid,function(t) matrix(0,length(t),plots))
    for (k in seq_len(plots)) {
        draw <- one()
        for (j in seq_along(draw)) what[[j]][,k] <- draw[[j]]
    }
    list(W=observed,What=what,statistics=null)
}

test_that("joint simulation agrees with a subject-indexed indicator oracle", {
    set.seed(721)
    n <- 37
    v <- cbind(a=sample(1:5,n,TRUE),b=rnorm(n),reverse=-seq_len(n))
    r <- rnorm(n)
    dr <- cbind(-1,rnorm(n))
    ic <- cbind(r/n,rnorm(n)/n)
    for (bandwidth in c(0,1.5)) {
        set.seed(722)
        g <- cumres_joint(NULL,variable=v,r=r,dr=dr,ic=ic,coef=NULL,
                          R=23,plots=3,b=bandwidth,subset.max=Inf)
        set.seed(722)
        oracle <- joint_oracle(r,dr,ic,v,g$x,23,3,bandwidth)
        expect_equal(g$W,oracle$W,tolerance=1e-12)
        expect_equal(g$What,oracle$What,tolerance=1e-12)
        expect_equal(g$joint$statistics,oracle$statistics,tolerance=1e-12)
    }
})

test_that("complementary orderings reuse subjects rather than sorted ranks", {
    n <- 40
    r <- sin(seq_len(n))
    ic <- matrix(r/n,n,1L)
    set.seed(723)
    g <- cumres_joint(NULL,variable=cbind(x=seq_len(n),reverse=-seq_len(n)),r=r,
                      dr=matrix(-1,n,1L),ic=ic,coef=NULL,R=19,plots=4,subset.max=Inf)
    # The intercept correction makes the final cumulative null value zero.
    # Reversing subjects gives the negative complementary prefix, not the
    # same multiplier assigned to a different observation at each rank.
    expect_equal(g$What[[1]][n,],rep(0,4),tolerance=1e-13)
    expect_equal(g$What[[2]][seq_len(n-1L),],
                 -g$What[[1]][rev(seq_len(n-1L)),],tolerance=1e-13)
    expect_equal(g$joint$statistics[,1L],g$joint$statistics[,3L],tolerance=1e-13)
})

test_that("family statistics, adjusted tails and envelopes use the shared law", {
    m <- fx_lm()
    set.seed(724)
    g <- cumres_joint(m,variable=c("x1","x2"),R=99,plots=0,level=.9)
    null <- g$joint$statistics
    observed <- as.data.frame(g)
    ranges <- vapply(g$x,function(t) diff(range(t)),numeric(1))
    expect_equal(g$joint$observed[1L],c(Sup=max(vapply(g$W,function(w) max(abs(w)),numeric(1)))))
    expect_equal(g$joint$simulated[,1L],apply(null[,c(1L,3L)],1L,max))
    expect_equal(g$joint$simulated[,2L],apply(sweep(null[,c(2L,4L)],2L,sqrt(ranges),"/"),1L,max))
    expect_equal(unname(g$joint$p.value),vapply(seq_len(2L),function(j)
        (1+sum(g$joint$simulated[,j]>=g$joint$observed[j]))/100,numeric(1)))
    expect_equal(g$joint$critical,sort(g$joint$simulated[,1L])[90L])
    expect_true(all(g$joint$adjusted.p.value[,1L]>=g$KS))
    expect_true(all(g$joint$adjusted.p.value[,2L]>=g$CvM))
    expect_equal(observed$scope,c(rep("marginal",4L),"family","family"))
    expect_equal(observed$test,c("Sup","L2","Sup","L2","Sup","L2.range"))
    expect_true(all(observed$p.value>=.01 & observed$p.value<=1))
    expect_equal(g$joint$envelope[[1L]]$lower,rep(-g$joint$critical,length(g$x[[1L]])))
    expect_equal(g$joint$envelope[[2L]]$upper,rep(g$joint$critical,length(g$x[[2L]])))
    expect_equal(as.data.frame(summary(g)),observed)
    expect_identical(summary(g)$joint,g$joint)
})

test_that("joint family tests are invariant to translations and positive units", {
    set.seed(725)
    n <- 72
    v <- cbind(a=seq_len(n),b=sample(1:7,n,TRUE))
    r <- rnorm(n)
    fit <- function(z) cumres_joint(NULL,variable=z,r=r,dr=matrix(0,n,0L),
                                    ic=matrix(0,n,0L),coef=NULL,R=39,plots=2,
                                    subset.max=Inf)
    set.seed(726); a <- fit(v)
    v[,1L] <- v[,1L]*1000+13
    v[,2L] <- v[,2L]*.125-9
    set.seed(726); b <- fit(v)
    expect_equal(a$joint$observed,b$joint$observed,tolerance=1e-12)
    expect_equal(a$joint$simulated,b$joint$simulated,tolerance=1e-12)
    expect_identical(a$joint$p.value,b$joint$p.value)
    expect_equal(a$joint$critical,b$joint$critical)
    expect_equal(a$W,b$W)
    expect_equal(a$What,b$What)
})

test_that("joint simulation supports known means, constant grids and R=1", {
    set.seed(727)
    g <- cumres_joint(NULL,variable=cbind(constant=rep(1,8)),r=seq_len(8),
                      dr=matrix(0,8,0L),ic=matrix(0,8,0L),coef=NULL,R=1,plots=0)
    expect_equal(unname(g$joint$observed[2L]),0)
    expect_equal(unname(g$joint$simulated[,2L]),0)
    expect_true(is.infinite(g$joint$critical))
    expect_equal(g$joint$p.value["L2.range"],c(L2.range=1))
    expect_equal(dim(g$What[[1L]]),c(1L,0L))
    expect_equal(g$n,8)
    expect_false("informative" %in% names(g))
    table <- as.data.frame(g)
    expect_true(all(is.na(table$informative)))
    expect_true(all(!is.na(table$p.value)))
})

test_that("GLM information checks are applied before joint family inference", {
    set.seed(728)
    d <- data.frame(cell=rep(0:1,each=40),x=rnorm(80))
    d$y <- rnorm(80)
    m <- lm(y~cell,d)
    set.seed(729)
    g <- cumres_joint(m,variable=c("cell","x"),data=d,R=39,plots=2)
    expect_identical(g$informative,c(FALSE,TRUE))
    expect_true(all(g$joint$statistics[,1L:2L]==0))
    expect_equal(g$joint$simulated[,1L],g$joint$statistics[,3L])
    expect_equal(g$KS[1L],1)
    expect_equal(g$CvM[1L],1)
    one <- cumres_joint(m,variable="cell",R=19,plots=0)
    expect_identical(one$joint$status,"uninformative")
    expect_equal(one$joint$p.value,c(Sup=1,L2.range=1))
    expect_equal(one$joint$critical,0)
    # Exercise unresolved membership without hiding or dropping its process.
    g$informative[2L] <- NA
    seed <- .Random.seed
    unresolved <- gof:::gof_joint_finalize(g)
    expect_identical(seed,.Random.seed)
    expect_identical(unresolved$joint$status,"indeterminate")
    expect_true(all(is.na(unresolved$joint$p.value)))
    expect_true(all(is.na(unresolved$joint$adjusted.p.value)))
    expect_true(is.na(unresolved$joint$critical))
    expect_identical(unresolved$W,g$W)
    expect_identical(unresolved$What,g$What)
})

test_that("joint wrappers validate scope and retain modern result methods", {
    m <- fx_lm()
    expect_error(cumres_joint(m,level=1),"level")
    expect_error(cumres_joint(m,level=NA_real_),"level")
    expect_error(cumres_joint(m,simulation="joint"),"do not also")
    expect_error(cumres_joint(structure(list(),class="ordreg")),"lm/glm")
    set.seed(730)
    g <- m |> cumres_joint(variable=c("x1","x2"),R=19,plots=1)
    set.seed(730)
    direct <- cumres(m,variable=c("x1","x2"),R=19,plots=1,simulation="joint")
    expect_identical(g,direct)
    seed <- .Random.seed
    expect_identical(g,gof:::gof_joint_finalize(g))
    expect_identical(seed,.Random.seed)
    expect_output(print(g),"Family scope")
    expect_output(print(summary(g)),"Joint subject-level")
    if (requireNamespace("generics",quietly=TRUE) && requireNamespace("tibble",quietly=TRUE))
        expect_equal(as.data.frame(generics::tidy(g)),as.data.frame(g))
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    expect_identical(plot(g,idx="x1"),g)
    expect_identical(plot(g,envelope=FALSE),g)
    expect_error(plot(g,envelope=NA),"envelope")
    expect_error(plot(g,idx="missing"),"Unknown")
})

test_that("joint simulation retains complete window domains above 500 rows", {
    n <- 650
    v <- cbind(a=seq_len(n),b=rep(1:5,length.out=n))
    r <- sin(seq_len(n))
    set.seed(731)
    g <- cumres_joint(NULL,variable=v,r=r,dr=matrix(0,n,0),ic=matrix(0,n,0),
                      coef=NULL,R=7,plots=1,b=.25,subset.max=Inf)
    expect_equal(lengths(g$x),c(1300L,10L))
    expect_equal(tail(g$W[[1L]],1L),0)
    expect_equal(tail(g$W[[2L]],1L),0)
    set.seed(731)
    oracle <- joint_oracle(r,matrix(0,n,0),matrix(0,n,0),v,g$x,7,1,.25)
    expect_equal(g$joint$statistics,oracle$statistics,tolerance=1e-12)
    expect_equal(g$What,oracle$What,tolerance=1e-12)
})

test_that("capped joint grids select evaluation points without discarding subjects", {
    n <- 620
    v <- cbind(a=seq_len(n),b=rep(1:9,length.out=n))
    r <- cos(seq_len(n))
    selector <- function(t,n.qt,n.eq)
        unique(round(seq(1,length(t),length.out=n.qt+n.eq)))
    set.seed(732)
    g <- cumres_joint(NULL,variable=v,r=r,dr=matrix(0,n,0),ic=matrix(0,n,0),
                      coef=NULL,R=11,plots=2,subset.max=31,subset=selector)
    expect_equal(g$n,n)
    expect_equal(lengths(g$x),c(31L,9L))
    set.seed(732)
    oracle <- joint_oracle(r,matrix(0,n,0),matrix(0,n,0),v,g$x,11,2)
    expect_equal(g$joint$statistics,oracle$statistics,tolerance=1e-12)
    expect_equal(g$W,oracle$W,tolerance=1e-12)
    expect_equal(g$What,oracle$What,tolerance=1e-12)
})
