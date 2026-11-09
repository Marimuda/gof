## The shared ordinal multiplier engine (R/ordinal-engine.R). cumres.ordreg
## delegates to it and must stay bit-identical to gof 1.3.0, including the
## random number stream; mixed-model designs use several multiplier factors.

drop_multiplier <- function(x) {
    x$multiplier <- NULL
    x
}

test_that("cumres.ordreg is bit-identical to the verbatim 1.3.0 implementation", {
    models <- ordinal_golden_models()
    cases <- ordinal_golden_cases()
    for (i in seq_len(nrow(cases))) {
        new <- ordinal_golden_run(cumres,models,cases[i,])
        old <- ordinal_golden_run(ref_cumres_ordreg_130,models,cases[i,])
        expect_identical(drop_multiplier(new$result),old$result,info=rownames(cases)[i])
        expect_identical(new$seed,old$seed,info=rownames(cases)[i])
        expect_identical(new$result$multiplier,
                         data.frame(factor="rows",levels=as.double(new$result$n),rows=new$result$n,
                                    draw="independent-normal",combine="sum",stringsAsFactors=FALSE),
                         info=rownames(cases)[i])
    }
})

test_that("cumres.ordreg reproduces the stored 1.3.0 golden objects", {
    compare_baseline <- function() {
        kind <- RNGkind()
        on.exit(do.call(RNGkind,as.list(kind)),add=TRUE)
        # The stored fixture was generated with R 4.6.1. R-devel adds
        # binom.kind to the seed header; use the recorded RNG configuration
        # rather than dropping header fields or weakening exact comparisons.
        # The verbatim-reference test above still exercises current defaults.
        if ("binom.kind" %in% names(formals(RNGkind))) {
            expect_warning(RNGversion("4.6.1"),"Buggy BTPE algorithm used for rbinom")
        } else RNGversion("4.6.1")
        golden <- readRDS(test_path("golden","ordinal-1.3.0.rds"))
        models <- ordinal_golden_models()
        cases <- ordinal_golden_cases()
        expect_identical(attr(golden,"cases"),cases)
        for (i in seq_len(nrow(cases))) {
            new <- ordinal_golden_run(cumres,models,cases[i,])
            expect_equal(drop_multiplier(new$result),golden[[i]]$result,info=rownames(cases)[i])
            expect_identical(new$seed,golden[[i]]$seed,info=rownames(cases)[i])
        }
    }
    kind <- RNGkind()
    compare_baseline()
    expect_identical(RNGkind(),kind)
})

test_that("indeterminate ordinal information is handled as in 1.3.0", {
    local_mocked_bindings(gof_ordinal_information=function(...) NA)
    model <- ordinal_fixture(k=3)
    for (aggregate in list("omnibus","cutpoints",ordinal_golden_custom)) {
        set.seed(77); new <- cumres(model,variable="x",R=9,plots=4,aggregate=aggregate); s1 <- .Random.seed
        set.seed(77); old <- ref_cumres_ordreg_130(model,variable="x",R=9,plots=4,aggregate=aggregate); s2 <- .Random.seed
        expect_identical(drop_multiplier(new),old)
        expect_identical(s1,s2)
        expect_true(all(is.na(c(new$KS,new$CvM))))
    }
})

test_that("ordinal helpers keep the 1.3.0 aggregation mapping and option rejections", {
    a <- gof_ordinal_aggregation(c("omnibus","score","cutpoints"))
    expect_identical(a,list(fn=gof_ordinal_norm,aggregation="omnibus"))
    expect_identical(gof_ordinal_aggregation("score"),list(fn=sum,aggregation="score"))
    expect_identical(gof_ordinal_aggregation("cut"),list(fn=identity,aggregation="cutpoints"))
    expect_identical(gof_ordinal_aggregation(sum)$aggregation,"score")
    expect_identical(gof_ordinal_aggregation(identity)$aggregation,"cutpoints")
    expect_identical(gof_ordinal_aggregation(ordinal_golden_custom),
                     list(fn=ordinal_golden_custom,aggregation="custom"))
    expect_error(gof_ordinal_aggregation("mean"),"should be one of")
    expect_error(gof_ordinal_aggregation(1),"function or a named aggregation")
    m <- ordinal_fixture(k=3)
    expect_error(cumres(m,R=2,resolution=1),"resolution and on_unresolved are not implemented")
    expect_error(cumres(m,R=2,on_unresolved="error"),"resolution and on_unresolved are not implemented")
    expect_error(cumres(m,R=2,simulation="shared"),"simulation is a default/linear/GLM option")
    expect_null(gof_ordinal_check_dots(col="red"))
})

## A synthetic two-factor design: factor "A" (4 levels) plus factor "B"
## (6 levels) with row loadings, ties in the ordering and two cutpoints.
engine_fixture <- function(seed=11) {
    set.seed(seed)
    n <- 24L; k <- 2L; d <- 3L
    list(n=n,k=k,d=d,
         ordering=sample(1:8,n,replace=TRUE),
         residual=matrix(rnorm(n*k),n,k),
         gradient=matrix(rnorm(n*k*d),n*k,d),
         A=list(name="A",index=sample(rep(1:4,6)),influence=matrix(rnorm(4*d),4,d)),
         B=list(name="B",index=sample(rep(1:6,4)),influence=matrix(rnorm(6*d),6,d),
                loading=runif(n,.5,2)))
}

run_engine <- function(f,factors,aggregate="cutpoints",R=15L,plots=15L,divisor=3,seed=5,metadata=list()) {
    a <- gof_ordinal_aggregation(aggregate)
    design <- gof_multiplier_design(factors,n_units=nrow(factors[[1]]$influence),divisor=divisor)
    set.seed(seed)
    gof_ordinal_engine(f$residual,f$gradient,design,ordering=f$ordering,aggregate=a$fn,
                       aggregation=a$aggregation,informative=TRUE,R=R,plots=plots,variable="z",
                       model="synthetic",design_info=gof_design_info(unit="block",n_rows=f$n,
                           fit_parameters=f$d,estimated_parameters=f$d,correction_dimension=f$d,
                           n_components=f$k),metadata=metadata)
}

## Brute force: explicit indicator sums at every distinct ordering value.
## Each factor draws one N(0,1) per level, factors in order, with levels
## numbered by first appearance along the sorted ordering.
brute_engine <- function(f,factors,aggregate,R,plots,divisor,seed) {
    n <- f$n; k <- f$k
    grid <- sort(unique(f$ordering))
    o <- order(f$ordering)
    process <- function(u,B) matrix(sapply(seq_len(k),function(j) vapply(grid,function(t) {
        s <- which(f$ordering<=t)
        sum(f$residual[s,j]*u[s])+sum(colSums(f$gradient[(j-1L)*n+s,,drop=FALSE])*B)
    },0)),ncol=k)/divisor
    agg <- function(w) {
        z <- t(apply(w,1,aggregate))
        if (nrow(z)!=nrow(w)) z <- t(z)
        z
    }
    observed <- agg(process(rep(1,n),numeric(f$d)))
    set.seed(seed)
    sims <- lapply(seq_len(max(R,plots)),function(r) {
        u <- numeric(n); B <- numeric(f$d)
        for (fa in factors) {
            L <- nrow(fa$influence)
            g <- numeric(L)
            g[unique(fa$index[o])] <- rnorm(L)
            loading <- if (is.null(fa$loading)) 1 else fa$loading
            u <- u+loading*g[fa$index]
            B <- B+drop(crossprod(fa$influence,g))
        }
        agg(process(u,B))
    })
    q <- ncol(observed)
    list(W=lapply(seq_len(q),function(i) observed[,i]),
         What=lapply(seq_len(q),function(i) sapply(sims[seq_len(plots)],function(s) s[,i])),
         KS=sapply(seq_len(q),function(i) mean(sapply(sims[seq_len(R)],function(s) max(abs(s[,i])))>=max(abs(observed[,i])))),
         CvM=sapply(seq_len(q),function(i) mean(sapply(sims[seq_len(R)],function(s) L2Test(s[,i],grid))>=L2Test(observed[,i],grid))),
         x=grid,seed=.Random.seed)
}

test_that("a two-factor sum multiplier design matches a brute-force simulation", {
    f <- engine_fixture()
    factors <- list(f$A,f$B)
    for (aggregate in list("cutpoints","omnibus")) {
        res <- run_engine(f,factors,aggregate=aggregate)
        seed <- .Random.seed
        fn <- gof_ordinal_aggregation(aggregate)$fn
        ref <- brute_engine(f,factors,fn,R=15L,plots=15L,divisor=3,seed=5)
        expect_equal(res$W,ref$W,tolerance=1e-12)
        expect_equal(lapply(res$What,unname),lapply(ref$What,unname),tolerance=1e-12)
        expect_equal(res$x[[1]],ref$x)
        expect_equal(res$KS,ref$KS)
        expect_equal(res$CvM,ref$CvM)
        expect_identical(seed,ref$seed)
    }
    res <- run_engine(f,factors,R=6L,plots=0L)
    expect_identical(res$multiplier,
                     data.frame(factor=c("A","B"),levels=c(4,6),rows=c(24L,24L),
                                draw="independent-normal",combine="sum",stringsAsFactors=FALSE))
    expect_identical(res$n,4L)
    expect_identical(res$design_info$unit,"block")
    expect_identical(as.data.frame(res)$unit,rep("block",2*2))
})

test_that("multiplier draws are invariant to relabelling factor levels", {
    f <- engine_fixture(seed=12)
    relabel <- function(fa,perm) {
        fa$influence[perm,] <- fa$influence
        fa$index <- perm[fa$index]
        fa
    }
    a <- run_engine(f,list(f$A,f$B),aggregate="omnibus")
    b <- run_engine(f,list(relabel(f$A,c(3L,1L,4L,2L)),relabel(f$B,c(6L,2L,5L,1L,3L,4L))),aggregate="omnibus")
    expect_identical(a,b)
})

test_that("each iteration consumes one normal per multiplier level", {
    f <- engine_fixture(seed=13)
    run_engine(f,list(f$A),R=7L,plots=11L,seed=21)
    after <- .Random.seed
    set.seed(21); invisible(rnorm(4*11))
    expect_identical(after,.Random.seed)
    run_engine(f,list(f$A,f$B),R=9L,plots=3L,seed=22)
    after <- .Random.seed
    set.seed(22); invisible(rnorm((4+6)*9))
    expect_identical(after,.Random.seed)
})

test_that("engine metadata and block multipliers reach results, summaries and printing", {
    f <- engine_fixture()
    res <- run_engine(f,list(f$A,f$B),aggregate="omnibus",R=5L,plots=2L,
                      metadata=list(random_effects=data.frame(term="1 | A",q=1),fit_check=list(shift=0)))
    expect_identical(res$random_effects,data.frame(term="1 | A",q=1))
    expect_identical(res$fit_check,list(shift=0))
    expect_output(print(res),"Multipliers: A \\(4 levels\\), B \\(6 levels\\)\\.")
    expect_output(print(summary(res)),"Multipliers: A \\(4 levels\\), B \\(6 levels\\)\\.")
    expect_error(run_engine(f,list(f$A),metadata=list(1)),"named list")
    m <- ordinal_fixture(k=3)
    set.seed(1)
    g <- cumres(m,variable="x",R=5,plots=0)
    out <- capture.output(print(g),print(summary(g)))
    expect_false(any(grepl("Multipliers",out)))
    expect_identical(summary(g)$multiplier,g$multiplier)
    expect_identical(gof_design_info(unit="block")$unit,"block")
})
