#' Calibrate cumulative residual diagnostics by simulating responses
#'
#' Simulate binomial or Poisson responses on a fixed design and compare their
#' cumulative residual statistics with the observed statistics. Fitted-null
#' simulations refit the mean model; a genuinely specified null does not.
#'
#' @param model A stored `glm` fit with binomial/logit or Poisson/log family.
#'   The model frame and design must be recoverable without evaluating the
#'   original fitting call. Other families and links are not implemented here.
#' @param variable Names of prespecified, fixed numeric ordering variables in
#'   the design or `data`. Fitted predictions (`"predicted"`) are rejected:
#'   freezing their response-dependent ordering is not this procedure.
#'   Refitting and reordering alone would also be insufficient to justify
#'   unrestricted raw fitted-order L2: at a flat true mean its integration
#'   range vanishes and has an estimated random scale. Such an extension is
#'   outside this implementation.
#' @param data Data containing additional ordering variables, aligned to fitted
#'   rows as for [cumres.glm()].
#' @param R Positive whole number of response simulations. This procedure does
#'   not run an inner multiplier simulation.
#' @param plots Number of successful simulation curves to retain, at most `R`.
#'   Changing this number changes storage only, not draws or p-values.
#' @param b Nonnegative prespecified moving-window bandwidth; zero requests
#'   cumulative residuals.
#' @param subset.max Maximum distinct evaluation points, or `Inf`/`NULL` for
#'   all points. A grid is selected once from the fixed orderings and reused in
#'   every simulation. This never subsamples fitted observations.
#' @param subset Grid selector with the same contract as [cumres.default()].
#' @param null `"fitted"` for parametric-bootstrap calibration, or
#'   `"specified"` for a fully specified response distribution.
#' @param mean Required for `null="specified"`: a finite mean per original
#'   fitted row, in fitted-row order. Binomial entries are probabilities;
#'   Poisson entries are count means. They must have been specified independently
#'   of the tested outcomes. Do not pass estimates from these outcomes and
#'   interpret the resulting test as a specified-null test.
#' @param control A list of arguments to [stats::glm.control()] used for the
#'   bootstrap refits. It does not change the original fit.
#' @param resolution Optional named positive finite measurement steps for
#'   selected fixed ordering variables, as in [cumres.glm()]. The declared
#'   `round(x / step) * step` coordinate is fixed before grid selection and is
#'   used for the observation and every simulated response. No rows/trials are
#'   merged and the model is not refitted on coarsened covariates.
#' @details
#' For binomial groups the residual is `successes - trials * probability`,
#' normalized by the square root of total trials. This gives exactly the same
#' cumulative process as individually expanded Bernoulli trials because all
#' trials in a group share their ordering value. Matrix responses and weighted
#' proportions must reconstruct integer success/failure counts. Poisson rows
#' must have unit prior weights (zero-weight rows are excluded); other prior
#' weights do not uniquely specify the sampling distribution needed here.
#'
#' Each simulation draws responses from the null distribution. In fitted mode
#' the same design, offsets, trials and family are used in `glm.fit()`, and
#' residuals use the refitted means. Sup and L2 use the same fixed grid for the
#' observation and every simulated sample. Their p-values are
#' `(1 + number of simulated statistics >= observed statistic)/(R + 1)`.
#' The comparison is marginal for each ordering and statistic; these are not
#' familywise-adjusted p-values or simultaneous confidence bands.
#'
#' In specified mode, conditional on the fixed design and grid, the observation
#' and simulated response samples have the same distribution under the stated
#' null. The rank construction therefore controls finite-sample type-I error
#' (conservatively with ties), without a Gaussian multiplier approximation.
#' This statement requires genuinely specified means and the independent
#' binomial/Poisson sampling model. It does not apply to fitted mode. The fitted
#' parametric bootstrap is an approximation that relies on accurate estimation
#' and bootstrap regularity; it is not guaranteed to fix every sparse-data or
#' boundary-estimation problem.
#'
#' All-zero responses, or all-successful binomial trials, have analytically
#' known boundary fitted means when the design has a nonzero constant column.
#' These extended maximum-likelihood fits are represented by exact zero/one
#' means, not optimizer-floor coefficients. Their residual curves are zero.
#' A fitted null at this boundary is degenerate and cannot assess unobserved
#' events. Other separated, nonconverged or unresolved refits are never silently
#' removed: informative-process p-values are withheld and bounds assign every unresolved
#' tail indicator either zero or one. These are bounds on the simulated tail
#' calculation, not confidence intervals or a guarantee of bootstrap validity.
#' Structurally certified zero-information processes retain zero statistics,
#' p-values of one and bounds of one even if other aspects of refits are
#' unresolved. Their failure statuses are still recorded.
#'
#' Prespecify ordering variables, bandwidth and grid rules. Outcome-driven
#' selection requires its own statistical adjustment. Custom callbacks,
#' quasi-likelihood, penalized fitting, clustered observations and arbitrary
#' prior-weight sampling distributions are outside this implementation.
#' @returns A `cumres` object with the usual curves, marginal p-values and
#' metadata. `calibration` distinguishes `"parametric-bootstrap"` and
#' `"specified-null-monte-carlo"`. `bootstrap` contains all `R` statuses,
#' refit warnings/failure reasons, observed fit status, the full simulated
#' statistics matrix, retained curve indices and a `tail.bounds` data frame.
#' The data-frame method adds `calibration`, `p.lower`, `p.upper`, `failed` and
#' `null` to the common result columns. Bounds coincide with the reported
#' p-value when all fits and score-space checks are resolved. `summary()`,
#' `coef()` and `plot()` use the usual cumulative-residual interfaces.
#' The `ordering` component and additional ordering columns in result tables
#' record fixed coordinates, declared resolution and grouped tie counts.
#' Binomial tie/changed-unit counts refer to underlying Bernoulli trials;
#' these metadata are computed from grouped trial counts without expansion.
#'
#' @references
#' Genest, C. and Rémillard, B. (2008). Validity of the parametric bootstrap for
#' goodness-of-fit testing in semiparametric models. Annales de l'Institut Henri
#' Poincare, Probabilites et Statistiques 44, 1096--1127.
#' \doi{10.1214/07-AIHP148}. The general regular-bootstrap principle is used;
#' its copula examples are not a finite-sample theorem for this GLM procedure.
#'
#' Phipson, B. and Smyth, G. K. (2010). Permutation P-values Should Never Be
#' Zero: Calculating Exact P-values When Permutations Are Randomly Drawn.
#' Statistical Applications in Genetics and Molecular Biology 9, Article 39.
#' \doi{10.2202/1544-6115.1585}, Section 4 (Monte Carlo tests).
#' @seealso [cumres()], [cumres.glm()]
#' @examples
#' set.seed(310)
#' d <- data.frame(x = seq(-1, 1, length.out = 60))
#' d$y <- rpois(nrow(d), exp(-1 + d$x))
#' fit <- glm(y ~ x, data = d, family = poisson())
#' diagnostic <- fit |> cumres_bootstrap(variable = "x", R = 39, plots = 3)
#' as.data.frame(diagnostic)
#'
#' # A mean specified before observing these data, with no refitting.
#' specified <- cumres_bootstrap(fit, variable = "x", R = 39, plots = 0,
#'   null = "specified", mean = exp(-1 + d$x))
#' coef(specified)
#' @export
cumres_bootstrap <- function(model, variable, data=model.frame(model),
                             R=999L, plots=min(R,50L), b=0,
                             subset.max=500L, subset=subset_t,
                             null=c("fitted","specified"), mean=NULL,
                             control=list(epsilon=1e-10,maxit=100L),resolution=NULL) {
    null <- match.arg(null)
    R <- gof_count(R,"R",1L)
    plots <- gof_count(plots,"plots")
    if (plots>R) stop("plots must not exceed R",call.=FALSE)
    if (!is.null(subset.max)) subset.max <- gof_count(subset.max,"subset.max",1L,TRUE)
    if (!is.function(subset)) stop("subset must be a function",call.=FALSE)
    if (!is.numeric(b) || is.complex(b) || !is.null(dim(b)) || length(b)!=1L || !is.finite(b) || b<0)
        stop("b must be a finite nonnegative bandwidth",call.=FALSE)
    if (!is.list(control)) stop("control must be a list of glm.control arguments",call.=FALSE)
    control <- do.call(stats::glm.control,control)
    input <- gof_bootstrap_inputs(model,variable,data,null,mean,resolution)
    grid <- lapply(seq_len(ncol(input$ordering)),function(j) {
        z <- input$ordering[,j]
        if (b>0) gof_window_grid(z,b,subset.max,subset) else {
            values <- sort(unique(z))
            values[gof_grid(values,subset.max,subset)]
        }
    })
    p <- length(grid)
    # Compute sums from sorted prefixes, avoiding trial expansion and avoiding
    # a dense n-by-grid indicator matrix. Group ties are accumulated together.
    operators <- lapply(seq_len(p),function(j) {
        z <- input$ordering[,j]
        o <- order(z)
        list(order=o,upper=findInterval(grid[[j]],z[o]),
             lower=if (b>0) findInterval(grid[[j]],z[o]+b) else rep(0L,length(grid[[j]])))
    })
    project <- function(residual) lapply(operators,function(a) {
        prefix <- c(0,cumsum(unname(residual[a$order])))/sqrt(input$n)
        prefix[a$upper+1L]-prefix[a$lower+1L]
    })
    statistic <- function(w) unlist(lapply(seq_len(p),function(j)
        c(SupTest(w[[j]]),L2Test(w[[j]],grid[[j]]))),use.names=FALSE)
    observed <- if (null=="specified") list(mean=input$mean,status="specified",warnings=character()) else
        gof_bootstrap_observed(input,model)
    informative <- if (null=="specified") rep(TRUE,p) else
        gof_informative_processes(input$X,rep(1,nrow(input$X)),input$ordering,grid,b)
    if (observed$status %in% c("extended-zero","extended-one")) informative[] <- FALSE
    W <- project(input$success-input$trials*observed$mean)
    for (j in which(informative %in% FALSE)) W[[j]][] <- 0
    observed.statistic <- statistic(W)
    simulated <- matrix(NA_real_,R,2L*p)
    colnames(simulated) <- paste(rep(input$variable,each=2L),rep(c("Sup","L2"),p),sep=":")
    What <- lapply(grid,function(z) matrix(numeric(),length(z),0L))
    status <- rep("unresolved",R)
    reasons <- warnings <- vector("list",R)
    retained <- integer()
    for (k in seq_len(R)) {
        y <- if (input$family$family=="binomial")
            stats::rbinom(length(input$trials),input$trials,observed$mean) else
            stats::rpois(length(input$trials),observed$mean)
        fit <- if (null=="specified") list(mean=observed$mean,status="specified",warnings=character()) else
            gof_bootstrap_refit(input,y,control)
        status[k] <- fit$status
        warnings[[k]] <- fit$warnings
        if (is.null(fit$mean)) {
            reasons[[k]] <- fit$reason
            next
        }
        w <- project(y-input$trials*fit$mean)
        for (j in which(informative %in% FALSE)) w[[j]][] <- 0
        simulated[k,] <- statistic(w)
        if (length(retained)<plots) {
            retained <- c(retained,k)
            for (j in seq_len(p)) What[[j]] <- cbind(What[[j]],w[[j]])
        }
    }
    failed <- sum(status=="unresolved")
    # Design-fixed saturation proves these statistics are zero for every
    # response/refit, including unresolved parameter fits in other cells.
    structural <- rep(informative %in% FALSE,each=2L)
    simulated[,structural] <- 0
    tail <- sweep(simulated,2L,observed.statistic,">=")
    lower <- (1+colSums(tail,na.rm=TRUE))/(as.double(R)+1)
    upper <- (1+colSums(tail,na.rm=TRUE)+colSums(is.na(tail)))/(as.double(R)+1)
    pv <- unname(lower)
    if (failed) {
        pv[] <- NA_real_
        warning(failed," of ",R," bootstrap refits unresolved; informative p-values withheld, see bootstrap$tail.bounds",
                call.=FALSE)
    }
    pv[structural] <- 1
    if (anyNA(informative)) {
        indeterminate <- rep(is.na(informative),each=2L)
        pv[indeterminate] <- NA_real_
        lower[indeterminate] <- 1/(as.double(R)+1)
        upper[indeterminate] <- 1
    }
    bounds <- data.frame(variable=rep(input$variable,each=2L),test=rep(c("Sup","L2"),p),
                         p.lower=lower,p.upper=upper,stringsAsFactors=FALSE)
    for (j in seq_len(p)) colnames(What[[j]]) <- if (length(retained)) paste0("simulation",retained) else character()
    structure(list(W=W,What=What,x=grid,KS=pv[seq.int(1L,2L*p,2L)],
        CvM=pv[seq.int(2L,2L*p,2L)],R=R,n=input$n,nobs=nrow(input$X),
        variable=input$variable,type=rep("residual",p),model="glm",
        response=input$response,informative=informative,ordering=input$ordering_metadata,
        design_info=gof_design_info(unit=if (input$family$family=="binomial") "trial" else "row",
            n_rows=nrow(input$X),n_trials=if (input$family$family=="binomial") input$n else NA_real_,
            fit_parameters=length(stats::coef(model)),
            estimated_parameters=if (null=="specified") 0 else ncol(input$X),
            fit_design_rank=model$rank,correction_dimension=0),
        calibration=if (null=="fitted") "parametric-bootstrap" else "specified-null-monte-carlo",
        bootstrap=list(null=null,status=status,reasons=reasons,warnings=warnings,
            observed.status=observed$status,observed.warnings=observed$warnings,
            failed=failed,statistics=simulated,curve.indices=retained,tail.bounds=bounds,
            mean=observed$mean,trials=input$trials)),class=c("cumres_bootstrap","cumres"))
}

gof_bootstrap_inputs <- function(model,variable,data,null,mean,resolution=NULL) {
    if (!inherits(model,"glm") || !identical(model$method,"glm.fit"))
        stop("Response bootstrap requires a glm fitted by stats::glm.fit",call.=FALSE)
    family <- stats::family(model)
    if (!(family$family=="binomial" && family$link=="logit") &&
        !(family$family=="poisson" && family$link=="log"))
        stop("Response bootstrap supports binomial/logit and Poisson/log models",call.=FALSE)
    if (!is.null(model$cluster)) stop("Response bootstrap requires independent observations",call.=FALSE)
    if (!is.data.frame(model$model))
        stop("Response bootstrap requires a stored model frame; fit with model=TRUE",call.=FALSE)
    mf <- model$model
    X <- stats::model.matrix(model)
    if (!nrow(mf) || nrow(X)!=nrow(mf) || !is.numeric(X) || any(!is.finite(X)))
        stop("Inconsistent or nonfinite stored model design",call.=FALSE)
    offset <- model$offset
    if (is.null(offset)) offset <- rep(0,nrow(mf))
    if (!is.numeric(offset) || length(offset)!=nrow(mf) || any(!is.finite(offset)))
        stop("Offsets must contain a finite value per fitted row",call.=FALSE)
    weights <- stats::model.weights(mf)
    if (is.null(weights)) weights <- rep(1,nrow(mf))
    if (!is.numeric(weights) || length(weights)!=nrow(mf) || any(!is.finite(weights)) || any(weights<0))
        stop("Prior weights must be finite and nonnegative",call.=FALSE)
    y <- stats::model.response(mf)
    if (family$family=="binomial") {
        if (is.factor(y)) {
            if (nlevels(y)!=2L) stop("Binomial factors must have two levels",call.=FALSE)
            y <- as.numeric(y!=levels(y)[1L])
        }
        if (is.matrix(y)) {
            if (ncol(y)!=2L) stop("Binomial matrix responses must have two columns",call.=FALSE)
            counts <- y*weights
        } else {
            if (!is.numeric(y) || is.complex(y) || length(y)!=nrow(mf))
                stop("Invalid binomial response",call.=FALSE)
            counts <- cbind(y,1-y)*weights
        }
        tolerance <- 64*.Machine$double.eps*pmax(1,abs(counts))
        if (any(!is.finite(counts)) || any(counts<0 | abs(counts-round(counts))>tolerance))
            stop("Binomial responses and weights must define nonnegative integer success/failure counts",call.=FALSE)
        counts <- round(counts)
        success <- counts[,1L]
        trials <- rowSums(counts)
        active <- trials>0
    } else {
        if (any(weights!=0 & weights!=1))
            stop("Poisson response bootstrap requires unit prior weights; other sampling distributions are unspecified",call.=FALSE)
        if (!is.numeric(y) || is.complex(y) || !is.null(dim(y)) || length(y)!=nrow(mf) ||
            any(!is.finite(y)) || any(y<0 | y!=floor(y)))
            stop("Poisson responses must be nonnegative integer counts",call.=FALSE)
        success <- as.numeric(y)
        trials <- rep(1,nrow(mf))
        active <- weights>0
    }
    if (!any(active)) stop("No positive-weight observations or binomial trials",call.=FALSE)
    if (!is.character(variable) || !length(variable) || anyNA(variable) || any(!nzchar(variable)))
        stop("variable must name prespecified fixed ordering variables",call.=FALSE)
    variable <- unique(variable)
    if ("predicted" %in% variable)
        stop("Use prespecified fixed orderings; predicted ordering is not supported by this bootstrap",call.=FALSE)
    response <- stats::formula(model)[[2L]]
    response.aliases <- unique(c(names(mf)[1L],paste(deparse(response),collapse=""),
                                 head(all.vars(response),1L)))
    if (any(variable %in% response.aliases))
        stop("Ordering by the tested response is not a fixed-design diagnostic",call.=FALSE)
    needed <- setdiff(variable,colnames(X))
    extra <- if (length(needed)) gof_fitted_data(data,mf) else NULL
    if (length(setdiff(needed,names(extra))))
        stop("Variable not found: ",paste(setdiff(needed,names(extra)),collapse=", "),call.=FALSE)
    values <- lapply(variable,function(v) {
        z <- if (v %in% colnames(X)) X[,v] else extra[[v]]
        if (is.ordered(z)) z <- as.numeric(z)
        if (!is.numeric(z) || is.complex(z) || !is.null(dim(z)) || length(z)!=nrow(mf) ||
            any(!is.finite(z[active]))) stop("Ordering variables must be finite numeric fitted-row vectors",call.=FALSE)
        z[active]
    })
    ordering <- do.call(cbind,values)
    colnames(ordering) <- variable
    coordinates <- gof_ordering_coordinates(ordering,resolution,
        units=if (family$family=="binomial") trials[active] else NULL)
    ordering <- coordinates$values
    if (null=="specified") {
        if (!is.numeric(mean) || is.complex(mean) || !is.null(dim(mean)) || length(mean)!=nrow(mf) ||
            any(!is.finite(mean)) || any(mean<0) || (family$family=="binomial" && any(mean>1)))
            stop("Supply a valid prespecified mean per original fitted row",call.=FALSE)
        mean <- as.numeric(mean)[active]
    } else if (!is.null(mean)) stop("mean is only used with null='specified'",call.=FALSE)
    X <- X[active,,drop=FALSE]
    if (null=="fitted" && (length(stats::coef(model))!=ncol(X) ||
        anyNA(stats::coef(model)) || any(!is.finite(stats::coef(model))) || qr(X)$rank!=ncol(X)))
        stop("Fitted-null bootstrap requires a finite full-rank fitted design",call.=FALSE)
    trials <- unname(trials[active])
    if (!is.finite(sum(trials)) || sum(trials)>2^53)
        stop("Total binomial trial count exceeds exact integer representation",call.=FALSE)
    list(X=X,offset=unname(offset[active]),success=unname(success[active]),trials=trials,
         family=family,ordering=ordering,ordering_metadata=coordinates$metadata,variable=variable,mean=mean,
         n=if (family$family=="binomial") sum(trials) else nrow(X),
         response=paste(deparse(response),collapse=""))
}

gof_bootstrap_boundary <- function(input,y) {
    # A literal nonzero constant column proves a direction sending every
    # finite-offset linear predictor to -Inf/+Inf. Do not promote a numerical
    # QR approximation of a hidden intercept to an exact boundary certificate.
    constant <- any(vapply(seq_len(ncol(input$X)),function(j) {
        z <- input$X[,j]
        z[1L]!=0 && all(z==z[1L])
    },logical(1)))
    if (!constant) return(NULL)
    if (all(y==0)) return(list(mean=rep(0,length(y)),status="extended-zero",warnings=character()))
    if (input$family$family=="binomial" && all(y==input$trials))
        return(list(mean=rep(1,length(y)),status="extended-one",warnings=character()))
    NULL
}

gof_bootstrap_observed <- function(input,model) {
    boundary <- gof_bootstrap_boundary(input,input$success)
    if (!is.null(boundary)) return(boundary)
    if (!isTRUE(model$converged)) stop("The observed model has not converged",call.=FALSE)
    y <- if (input$family$family=="binomial") cbind(input$success,input$trials-input$success) else input$success
    fit <- model
    fit$offset <- input$offset
    gof_validate_glm_mle(fit,input$X,y,rep(1,nrow(input$X)))
    mu <- input$family$linkinv(drop(input$X%*%stats::coef(model))+input$offset)
    list(mean=mu,status="regular",warnings=character())
}

gof_bootstrap_refit <- function(input,y,control) {
    boundary <- gof_bootstrap_boundary(input,y)
    if (!is.null(boundary)) return(boundary)
    messages <- character()
    answer <- tryCatch(withCallingHandlers({
        response <- if (input$family$family=="binomial") cbind(y,input$trials-y) else y
        fit <- stats::glm.fit(x=input$X,y=response,family=input$family,
                              offset=input$offset,control=control,intercept=FALSE)
        if (!isTRUE(fit$converged)) stop("Bootstrap refit did not converge")
        if (anyNA(fit$coefficients) || any(!is.finite(fit$coefficients)) || fit$rank!=ncol(input$X))
            stop("Bootstrap refit is singular or nonfinite")
        fit$family <- input$family
        fit$offset <- input$offset
        class(fit) <- c("glm","lm")
        gof_validate_glm_mle(fit,input$X,response,rep(1,nrow(input$X)))
        mu <- input$family$linkinv(drop(input$X%*%fit$coefficients)+input$offset)
        list(mean=mu,status="regular")
    },warning=function(w) {
        messages <<- c(messages,conditionMessage(w))
        invokeRestart("muffleWarning")
    }),error=function(e) list(mean=NULL,status="unresolved",reason=conditionMessage(e)))
    answer$warnings <- messages
    answer
}

#' @rdname cumres_bootstrap
#' @param x An object returned by `cumres_bootstrap()`.
#' @param row.names,optional,... Arguments passed to the common data-frame method.
#' @export
as.data.frame.cumres_bootstrap <- function(x,row.names=NULL,optional=FALSE,...) {
    result <- as.data.frame.cumres(x,row.names=row.names,optional=optional,...)
    result$calibration <- x$calibration
    result$p.lower <- x$bootstrap$tail.bounds$p.lower
    result$p.upper <- x$bootstrap$tail.bounds$p.upper
    result$failed <- x$bootstrap$failed
    result$null <- x$bootstrap$null
    result
}

#' @rdname cumres_bootstrap
#' @param object An object returned by `cumres_bootstrap()`.
#' @export
summary.cumres_bootstrap <- function(object,...) {
    result <- summary.cumres(object,...)
    result$calibration <- object$calibration
    result$bootstrap <- object$bootstrap
    class(result) <- c("summary.cumres_bootstrap","summary.cumres")
    result
}

gof_print_bootstrap <- function(x) {
    cat("Calibration: ",x$calibration,".\n",sep="")
    cat("Null: ",x$bootstrap$null,"; observed fit: ",x$bootstrap$observed.status,".\n",sep="")
    if (x$bootstrap$failed) {
        cat(x$bootstrap$failed," of ",x$R," refits unresolved; informative p-values withheld.\n",sep="")
        print(x$bootstrap$tail.bounds,row.names=FALSE)
    }
    invisible(x)
}

#' @rdname cumres_bootstrap
#' @export
print.cumres_bootstrap <- function(x,...) {
    print.cumres(x,...)
    gof_print_bootstrap(x)
    invisible(x)
}

#' @rdname cumres_bootstrap
#' @export
print.summary.cumres_bootstrap <- function(x,...) {
    print.summary.cumres(x,...)
    gof_print_bootstrap(x)
    invisible(x)
}
