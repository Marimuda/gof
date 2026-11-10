test_that("projection references agree with independent dense cluster algebra", {
    p <- clmm_problem_fixture("intercept",m=15L,nc=4L,seed=13)
    set.seed(72)
    S <- matrix(rnorm(p$blocks$m*length(p$theta)),p$blocks$m)
    fit <- structure(list(score=list(g=S),information=diag(ncol(S))),class="gof_clmm_fit")
    residual <- matrix(rnorm(p$n*3),p$n)*rep(c(1,2),length.out=p$n)
    gradient <- matrix(0,p$n*3,ncol(S))
    ordering <- rep(c(2,0,1,0),length.out=p$n)
    o <- order(ordering);first <- unique(p$blocks$id[o]);grid <- sort(unique(ordering))
    A <- lapply(1:3,function(k) sapply(grid,function(t)
        rowsum(residual[,k]*(ordering<=t),p$blocks$id,reorder=FALSE)[,1]))
    projection <- diag(nrow(S))-S%*%solve(crossprod(S),t(S))
    for (method in c("projection","projection-HC2")) {
        c <- gof_clmm_correction(p,fit,residual,gradient,method)
        scale <- if(method=="projection") rep(1,nrow(S)) else 1/sqrt(diag(projection))
        set.seed(85)
        draws <- replicate(9,{g<-numeric(nrow(S));g[first]<-rnorm(nrow(S));g})
        expected <- lapply(A,function(a)t(a)%*%projection%*%(draws*scale)/(p$N/sqrt(p$blocks$m)))
        set.seed(85)
        answer <- gof_ordinal_engine(residual,c$gradient,c$multiplier,ordering,
            identity,"cutpoints",TRUE,9L,9L,"x","clmm",list())
        expect_equal(answer$What,expected,tolerance=1e-11)
        expect_equal(answer$W,lapply(A,function(a)colSums(a)/(p$N/sqrt(p$blocks$m))),ignore_attr=TRUE)
        transformed <- fit;transformed$score[[1]] <- S[,rev(seq_len(ncol(S)))]*
            rep(c(.001,2,100,7,3),each=nrow(S),length.out=length(S))
        other <- gof_clmm_correction(p,transformed,residual,gradient,method)
        set.seed(85)
        again <- gof_ordinal_engine(residual,other$gradient,other$multiplier,ordering,
            identity,"cutpoints",TRUE,9L,9L,"x","clmm",list())
        expect_equal(again$What,answer$What,tolerance=1e-11)
        # A pure score-space process has zero simulated covariance, also with HC2.
        pure <- matrix(S[p$blocks$id,1],p$n,1)
        cancel <- gof_clmm_correction(p,fit,pure,matrix(0,p$n,ncol(S)),method)
        expect_true(is.na(gof_clmm_informative(pure,cancel$gradient,cancel$multiplier,rep(1,p$n))))
    }
    legacy <- gof_clmm_correction(p,fit,residual,gradient,"derivative")
    expect_identical(legacy$gradient,gradient)
    expect_identical(legacy$multiplier,gof_clmm_multiplier(p,fit))
    fit$score[[1]][,2] <- fit$score[[1]][,1]
    expect_error(gof_clmm_correction(p,fit,residual,gradient,"projection"),class="gof_clmm_score_rank_deficient")
    fit$score[[1]] <- diag(nrow(S))[,seq_len(ncol(S))]
    expect_error(gof_clmm_correction(p,fit,residual,gradient,"projection-HC2"),class="gof_clmm_score_leverage_singular")
})

test_that("public correction controls preserve observations and refuse boundary inputs", {
    skip_if_not_installed("ordinal")
    skip_if_not_installed("numDeriv")
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=3)
    f <- ordinal::clmm(y~x+(1|g),data=d)
    set.seed(1901)
    baseline <- cumres(f,variable="x",R=9,plots=12)
    rng <- .Random.seed
    for (method in c("projection","projection-HC2")) {
        set.seed(1901)
        answer <- cumres(f,variable="x",R=9,plots=12,
                         control=gof_clmm_control(correction=method))
        expect_identical(answer$W,baseline$W)
        expect_identical(answer$fit_check$theta,baseline$fit_check$theta)
        expect_identical(.Random.seed,rng)
        expect_identical(answer$fit_check$correction$method,method)
        expect_equal(sum(answer$fit_check$correction$leverage),5)
        expect_true(all(is.finite(c(answer$KS,answer$CvM))))
        boundary <- f;boundary$ST[[1]][] <- 0
        expect_error(cumres(boundary,variable="x",R=3,
            control=gof_clmm_control(correction=method)),class="gof_clmm_covariance_boundary")
    }
})
