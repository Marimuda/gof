## Ordinal (lava::ordreg) cases pinned by golden/ordinal-1.3.0.rds and by the
## verbatim 1.3.0 reference implementation ref_cumres_ordreg_130. Shared by
## golden/make-ordinal-golden.R and test-ordinal-engine.R.

ordinal_golden_custom <- function(w) c(sum(w),max(abs(w)))

ordinal_golden_ties <- function(n=150) {
    set.seed(321)
    d <- data.frame(x=round(rnorm(n)))
    d$y <- ordered(cut(d$x+rnorm(n),breaks=c(-Inf,-.5,.5,Inf)))
    lava::ordreg(y~x,d)
}

ordinal_golden_intercept <- function(n=120) {
    set.seed(6073)
    d <- data.frame(y=ordered(sample(1:3,n,replace=TRUE)))
    lava::ordreg(y~1,d)
}

ordinal_golden_models <- function()
    list(k2=ordinal_fixture(k=2),k3=ordinal_fixture(k=3),
         ties=ordinal_golden_ties(),intercept=ordinal_golden_intercept())

## One row per case: model, variable, aggregate (name; "custom" means
## ordinal_golden_custom), R, plots and the seed set before the call.
ordinal_golden_cases <- function() {
    grid <- expand.grid(aggregate=c("omnibus","score","cutpoints","custom"),
                        variable=c("x","predicted"),model=c("k2","k3"),
                        stringsAsFactors=FALSE)
    grid$R <- 20L; grid$plots <- 5L
    extra <- data.frame(
        aggregate=c("omnibus","cutpoints","custom","score","omnibus","custom","omnibus","cutpoints","omnibus"),
        variable=c("x","x","predicted","predicted","predicted","predicted","x","predicted","x"),
        model=c("ties","ties","ties","ties","intercept","intercept","k3","k2","k3"),
        R=c(20L,20L,20L,20L,20L,20L,20L,20L,20L),plots=c(5L,5L,5L,5L,5L,5L,25L,25L,0L),
        stringsAsFactors=FALSE)
    cases <- rbind(grid[,c("aggregate","variable","model","R","plots")],extra)
    cases$seed <- 1000L+seq_len(nrow(cases))
    rownames(cases) <- paste(cases$model,cases$variable,cases$aggregate,
                             paste0("R",cases$R),paste0("p",cases$plots),sep="_")
    cases
}

ordinal_golden_aggregate <- function(name)
    if (identical(name,"custom")) ordinal_golden_custom else name

## Run one case with `fun` (a cumres.ordreg implementation); returns the
## result and the RNG state immediately afterwards.
ordinal_golden_run <- function(fun,models,case) {
    set.seed(case$seed)
    result <- fun(models[[case$model]],variable=case$variable,R=case$R,
                  aggregate=ordinal_golden_aggregate(case$aggregate),plots=case$plots)
    list(result=result,seed=get(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
}
