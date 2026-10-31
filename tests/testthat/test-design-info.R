test_that("binomial metadata separates active groups, trials and coefficient correction", {
    d <- fx_grouped_data();d$w <- c(0,1,2,1,0,3,1,2)
    m <- glm(cbind(s,f)~x,d,weights=w,family=binomial())
    g <- cumres(m,variable="x",R=9,plots=0)
    info <- g$design_info
    expect_s3_class(info,"data.frame")
    expect_equal(info$n_rows,sum(d$w>0))
    expect_equal(info$n_trials,sum((d$s+d$f)*d$w))
    expect_equal(g$n,info$n_trials)
    expect_identical(info$unit,"trial")
    expect_equal(info$fit_parameters,2)
    expect_equal(info$estimated_parameters,2)
    expect_equal(info$fit_design_rank,2)
    expect_equal(info$correction_dimension,2)
    expect_equal(info$n_components,1)
})

test_that("precision weights do not multiply row counts", {
    d <- data.frame(x=1:12,y=sin(1:12),w=rep(c(0,2,5),4))
    m <- lm(y~x,d,weights=w)
    g <- cumres(m,variable="x",R=9,plots=0)
    expect_identical(g$design_info$unit,"row")
    expect_equal(g$design_info$n_rows,8)
    expect_equal(g$n,8)
    expect_true(is.na(g$design_info$n_trials))
    expect_equal(g$design_info$estimated_parameters,2)
})

test_that("specified-null bootstrap does not inherit estimated coefficient counts", {
    d <- data.frame(x=rep(1:4,10),y=rep(c(1,2,2,4),10))
    m <- glm(y~x,d,family=poisson())
    fixed <- cumres_bootstrap(m,"x",null="specified",mean=rep(2,40),R=9,plots=0)
    fitted <- cumres_bootstrap(m,"x",R=9,plots=0)
    expect_equal(fixed$design_info$estimated_parameters,0)
    expect_equal(fitted$design_info$estimated_parameters,2)
    expect_equal(fixed$design_info$fit_parameters,2)
    expect_equal(fixed$design_info$correction_dimension,0)
    expect_equal(fitted$design_info$correction_dimension,0)
    expect_equal(fixed$design_info$n_rows,40)
    aliased <- glm(y~x+I(2*x),d,family=poisson())
    known <- cumres_bootstrap(aliased,"x",null="specified",mean=rep(2,40),R=3,plots=0)
    expect_equal(known$design_info$fit_parameters,3)
    expect_equal(known$design_info$fit_design_rank,2)
    expect_equal(known$design_info$estimated_parameters,0)
})

test_that("cutpoint components do not create independent ordinal subjects", {
    m <- ordinal_fixture(3)
    a <- cumres(m,"x",aggregate="score",R=3,plots=0)
    b <- cumres(m,"x",aggregate="cutpoints",R=3,plots=0)
    expect_equal(a$design_info,b$design_info)
    expect_identical(a$design_info$unit,"subject")
    expect_equal(a$design_info$n_rows,m$up$n)
    expect_equal(a$design_info$n_components,2)
    expect_equal(a$design_info$fit_parameters,length(coef(m)))
    expect_equal(a$design_info$correction_dimension,length(coef(m)))
    expect_true(is.na(a$design_info$fit_design_rank))
})

test_that("SEM pooling and projections retain subject and nuisance meanings", {
    set.seed(61101)
    model <- lava::lvm(list(c(y1,y2,y3)~eta,eta~x));lava::latent(model) <- ~eta
    m <- lava::estimate(model,lava::sim(model,90))
    mean <- cumres(m,"y1",x="x",R=3,plots=0)
    pooled <- cumres_normality(m,R=3,plots=0)
    projected <- cumres_normality(m,variable="y1",R=3,plots=0)
    expect_identical(mean$design_info$unit,"subject")
    expect_equal(mean$design_info$n_rows,90)
    expect_equal(mean$design_info$n_components,1)
    expect_equal(pooled$n,90)
    expect_equal(pooled$design_info$n_components,3)
    expect_equal(projected$design_info$n_components,1)
    expect_equal(pooled$design_info$estimated_parameters,length(lava::pars(m)))
    expect_true(is.na(pooled$design_info$fit_design_rank))
})

test_that("explicit correction dimensions are not fabricated estimator counts", {
    args <- list(model=NULL,variable=cbind(x=1:4),r=c(1,-1,2,-2),
        dr=matrix(1,4,1),ic=matrix(c(1,-1,1,-1),4,1),R=3,plots=0)
    for(fun in list(cumres,cumres_joint)) {
        g <- do.call(fun,args)
        expect_identical(g$design_info$unit,"residual")
        expect_true(is.na(g$design_info$n_rows))
        expect_true(is.na(g$design_info$fit_parameters))
        expect_true(is.na(g$design_info$estimated_parameters))
        expect_true(is.na(g$design_info$fit_design_rank))
        expect_equal(g$design_info$correction_dimension,1)
    }
})

test_that("design metadata survives reporting without changing legacy objects or RNG", {
    set.seed(61102)
    g <- cumres_joint(fx_lm(),variable=c("x1","predicted"),R=9,plots=0)
    state <- .Random.seed
    tab <- as.data.frame(g)
    expect_identical(summary(g)$design_info,g$design_info)
    expect_identical(as.data.frame(summary(g)),tab)
    expect_true(all(tab$estimated_parameters==3))
    expect_true(all(tab$unit=="row"))
    expect_equal(unique(tab$n_rows),g$nobs)
    if(requireNamespace("generics",quietly=TRUE) && requireNamespace("tibble",quietly=TRUE))
        expect_identical(as.data.frame(generics::tidy(g)),tab)
    legacy <- g;legacy$design_info <- NULL
    expect_false("estimated_parameters" %in% names(as.data.frame(legacy)))
    malformed <- g;malformed$design_info <- g$design_info[FALSE,]
    expect_error(as.data.frame(malformed),"Design metadata must contain one row")
    expect_identical(.Random.seed,state)
})

test_that("packaged development summaries retain every planned cohort and its uncertainty", {
    path <- system.file("extdata","small-design",package="gof")
    records <- read.csv(file.path(path,"development-records.csv"))
    totals <- read.csv(file.path(path,"development-summary.csv"))
    expect_equal(nrow(records),1600)
    for(i in seq_len(nrow(totals))) {
        row <- totals[i,]
        r <- records[records$case==row$case & records$method==row$method & records$test==row$test,]
        expect_identical(sort(r$replicate),1:200)
        lo <- sum(r$p.upper<=.05);hi <- sum(r$p.lower<=.05)
        expected <- c(200,lo/200,hi/200,binom.test(lo,200)$conf.int[1],
                      binom.test(hi,200)$conf.int[2])
        expect_equal(unname(unlist(row[4:8])),unname(expected),tolerance=1e-13)
    }
})
