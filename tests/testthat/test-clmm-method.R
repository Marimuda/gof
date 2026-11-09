skip_if_not_installed("ordinal")
skip_if_not_installed("numDeriv")

clmm_method_fixture <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            d <- clmm_simulate("intercept",m=40L,nc=6L,seed=3)
            cache <<- list(data=d,model=ordinal::clmm(y~x+(1|g),data=d))
        }
        cache
    }
})

test_that("the public mixed-model method uses independent blocks and preserves its input", {
    f <- clmm_method_fixture()$model
    before <- serialize(f,NULL)
    set.seed(1901)
    z <- cumres(f,variable="x",R=9,plots=12)
    after <- .Random.seed
    set.seed(1901)
    invisible(rnorm(40L*12L))
    expect_identical(after,.Random.seed)
    expect_identical(before,serialize(f,NULL))
    expect_s3_class(z,"cumres")
    expect_identical(z$model,"clmm")
    expect_identical(z$n,40L)
    expect_identical(dim(z$What[[1]]),c(240L,12L))
    expect_identical(z$design_info$unit,"block")
    expect_equal(z$design_info$n_rows,240)
    expect_equal(z$design_info$fit_parameters,5)
    expect_equal(z$multiplier$levels,40)
    expect_identical(z$multiplier$factor,"g")
    expect_lt(z$fit_check$node_check$standardized_score,.05)
    expect_true(all(is.finite(as.data.frame(z)$p.value)))
    sm <- summary(z)
    expect_identical(sm$fit_check,z$fit_check)
    expect_identical(sm$random_effects,z$random_effects)
    expect_equal(as.data.frame(sm),as.data.frame(z))
    expect_output(print(z),"Multipliers: g")
    expect_output(print(sm),"Multipliers: g")
    pdf_file <- tempfile(fileext=".pdf")
    grDevices::pdf(pdf_file)
    tryCatch(expect_no_error(plot(z)),finally=grDevices::dev.off())
    expect_gt(file.info(pdf_file)$size,0)
    unlink(pdf_file)
})

test_that("the assembled process matches a direct block-multiplier reference", {
    f <- clmm_method_fixture()$model
    problem <- gof_clmm_problem(f)
    polished <- gof_polish(problem,gof_select_integrator(problem))
    marginal <- gof_clmm_marginal(problem,polished$theta)
    ## The gradient sign and column-major stacking are independently checked
    ## against numerical differentiation of the marginal CDF.
    numerical <- -numDeriv::jacobian(function(theta)
        as.vector(gof_clmm_marginal(problem,theta)$cdf),polished$theta)
    expect_equal(unname(marginal$gradient),unname(numerical),tolerance=1e-7)
    residual <- outer(problem$y,1:3,"<=")-marginal$cdf
    influence <- polished$score[[1]] %*% solve(polished$information)
    x <- problem$X[,"x"]
    grid <- sort(unique(x))
    clusters <- unique(problem$blocks$id[order(x)])
    expected_W <- matrix(0,length(grid),3)
    expected_sim <- array(0,c(length(grid),3,4))
    set.seed(1902)
    for (b in 0:4) {
        if (b) {
            g <- rnorm(length(clusters))
            row_g <- g[match(problem$blocks$id,clusters)]
            delta <- drop(crossprod(g,influence[clusters,,drop=FALSE]))
        }
        for (i in seq_along(grid)) {
            selected <- which(x<=grid[i])
            value <- colSums(residual[selected,,drop=FALSE])
            if (b) for (k in 1:3) {
                rows <- selected+(k-1L)*length(x)
                value[k] <- sum(residual[selected,k]*row_g[selected])+
                    sum(colSums(numerical[rows,,drop=FALSE])*delta)
            }
            value <- value/(problem$N/sqrt(problem$blocks$m))
            if (b) expected_sim[i,,b] <- value else expected_W[i,] <- value
        }
    }
    set.seed(1902)
    z <- cumres(f,variable="x",aggregate="cutpoints",R=4,plots=4)
    expect_equal(do.call(cbind,z$W),expected_W,tolerance=1e-10,ignore_attr=TRUE)
    for (k in 1:3) expect_equal(z$What[[k]],expected_sim[,k,],tolerance=1e-7)
})

test_that("fixed-effect row-multiplier reduction preserves the existing ordinal method", {
    f <- ordinal_fixture(k=3)
    env <- f$up
    cdf <- as.matrix(predict(f,p=coef(f),type="cumulative"))
    residual <- outer(env$y,seq_len(ncol(cdf)),"<=")-cdf
    gradient <- -numDeriv::jacobian(function(p) predict(f,p=p,type="cumulative"),coef(f))
    n <- nrow(cdf)
    mult <- gof_multiplier_design(list(list(name="rows",index=seq_len(n),
        influence=as.matrix(lava::iid(f)))),n_units=n,divisor=sqrt(n))
    set.seed(1903)
    reference <- cumres(f,variable="x",R=9,plots=3)
    set.seed(1903)
    reduced <- gof_ordinal_engine(residual,gradient,mult,env$X[,"x"],
        aggregate=gof_ordinal_norm,aggregation="omnibus",informative=TRUE,R=9,plots=3,
        variable="x",model="ordreg",design_info=reference$design_info)
    expect_equal(reduced,reference,tolerance=1e-6)
})

test_that("row order and grouping labels do not select different multiplier draws", {
    d <- clmm_method_fixture()$data
    f <- clmm_method_fixture()$model
    d <- d[rev(seq_len(nrow(d))),]
    d$g <- factor(paste0("group",as.integer(d$g)),levels=paste0("group",40:1))
    other <- ordinal::clmm(y~x+(1|g),data=d)
    set.seed(1904); a <- cumres(f,variable="x",R=7,plots=3)
    set.seed(1904); b <- cumres(other,variable="x",R=7,plots=3)
    expect_equal(a$W,b$W,tolerance=1e-5)
    expect_equal(a$What,b$What,tolerance=1e-5)
    expect_identical(a$KS,b$KS)
    expect_identical(a$CvM,b$CvM)
})

test_that("mixed-model ordering, weights and unavailable paths are explicit", {
    f <- clmm_method_fixture()$model
    expect_error(cumres(f,variable="conditional",R=2),"BLUP")
    expect_error(cumres(f,variable="y",R=2),"retained clmm covariates")
    expect_error(cumres(f,variable="missing",R=2),"Variable not found")
    expect_error(cumres(f,R=0),"R")
    expect_error(cumres(f,plot=NA),"plot")
    expect_error(cumres(f,control=list()),"gof_clmm_control")
    expect_error(cumres(f,resolution=1),"not implemented")
    expect_error(cumres_joint(f),"not implemented for clmm")
    wine <- NULL
    data("wine",package="ordinal",envir=environment())
    crossed <- ordinal::clmm(rating~temp+contact+(1|judge)+(1|bottle),data=wine)
    expect_error(cumres(crossed,R=2),"Crossed or non-chain")
    slope <- ordinal::clmm(rating~temp+contact+(1+temp|judge),data=wine)
    expect_error(cumres(slope,R=2),"Only 9 independent blocks")
    d <- clmm_method_fixture()$data
    fw <- ordinal::clmm(y~x+(1|g),data=d,weights=rep(2,nrow(d)))
    pp <- gof_clmm_problem(fw)
    expect_equal(pp$N,2*nrow(d))
    expect_identical(pp$blocks$m,40L)
    set.seed(1905)
    pred <- cumres(f,R=3,plots=0)
    expect_identical(pred$fit_check$ordering,"fitted fixed predictor")
    expect_identical(length(pred$What[[1]]),0L)
})

test_that("a vanishing block-multiplier covariance withholds inference", {
    mult <- gof_multiplier_design(list(list(name="g",index=1:3,influence=matrix(1,3,1))),3,sqrt(3))
    residual <- matrix(c(1,2,3),3,1)
    gradient <- matrix(0,3,1)
    expect_true(gof_clmm_informative(residual,gradient,mult,1:3))
    expect_identical(gof_clmm_informative(residual*0,gradient,mult,1:3),NA)
})

test_that("frequency weights, zero-weight rows and offsets match literal row replication", {
    d <- clmm_method_fixture()$data
    d$o <- .25*d$z
    d$w <- rep(c(0,2,2,2,2,2),40)
    weighted <- ordinal::clmm(y~x+offset(o)+(1|g),data=d,weights=w)
    repeated <- d[rep(seq_len(nrow(d)),d$w),]
    expanded <- ordinal::clmm(y~x+offset(o)+(1|g),data=repeated)
    set.seed(1906); a <- cumres(weighted,variable="x",R=7,plots=3)
    set.seed(1906); b <- cumres(expanded,variable="x",R=7,plots=3)
    expect_equal(a$W,b$W,tolerance=1e-6)
    expect_equal(a$What,b$What,tolerance=1e-6)
    expect_identical(a$KS,b$KS)
    expect_identical(a$CvM,b$CvM)
    expect_equal(a$design_info$n_rows,200)
    expect_equal(b$design_info$n_rows,400)
    expect_equal(a$fit_check$theta,b$fit_check$theta,tolerance=1e-6)
})
