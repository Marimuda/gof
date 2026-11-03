#!/usr/bin/env Rscript
# Source after loading gof from the isolated library, or set R_LIBS for Rscript.
# This exercises installed help, vignettes and public result methods.
local({
    if (!requireNamespace("gof", quietly = TRUE))
        stop("Install gof with built vignettes in an isolated library first")
    path <- normalizePath(find.package("gof"), winslash = "/", mustWork = TRUE)
    namespace_path <- normalizePath(getNamespaceInfo(asNamespace("gof"), "path"),
                                    winslash = "/", mustWork = TRUE)
    if (!identical(path, namespace_path))
        stop("The loaded gof namespace does not match the selected library")
    library_path <- dirname(path)
    topics <- c("gof", "cumres", "cumres.glm", "cumres.lm", "cumres.default",
                "cumres.lvmfit", "cumres.ordreg", "plot.cumres", "cumres-object",
                "coef.cumres", "summary.cumres", "as.data.frame.cumres",
                "as.data.frame.summary.cumres", "tidy.cumres", "surgunit",
                "cumres_joint", "cumres_normality", "cumres_bootstrap")
    for (topic in topics) {
        if (!length(utils::help(topic, package = "gof", lib.loc = library_path)))
            stop("Installed help topic not found: ", topic)
    }
    vignettes <- utils::vignette(package = "gof", lib.loc = library_path)$results
    if (is.null(vignettes) || nrow(vignettes) < 6L)
        stop("Install the source archive after building its vignettes")
    html <- list.files(file.path(path, "doc"), pattern = "[.]html$", full.names = TRUE)
    if (!all(c("getting-started.html","statistical-details.html",
               "joint-normality-bootstrap.html","ordering-and-ties.html",
               "small-design-diagnostics.html","finite-panel-validation.html") %in% basename(html)))
        stop("The installed vignette HTML files are missing")
    evidence <- utils::read.csv(system.file("extdata","fixed-order","summary.csv",package="gof"))
    stopifnot(nrow(evidence)==106L,sum(evidence$primary_family)==16L,
              all(rowSums(evidence[paste0("count_",c("R","N","U","OF","DG","IN"))])==evidence$N),
              all(evidence$assessment[evidence$scenario=="P01_greenflies" & evidence$primary_family]==
                    "Upper rejection-rate bound at most 6% at this scenario"))

    news <- utils::news(package = "gof", lib.loc = library_path)
    version <- as.character(utils::packageVersion("gof", lib.loc = library_path))
    if (is.null(news) || !"Version" %in% names(news) || !version %in% news$Version)
        stop("The installed NEWS file has no structured entry for version ", version)

    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = globalenv())
    on.exit(if (had_seed) assign(".Random.seed", previous_seed, envir = globalenv()) else
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv()), add = TRUE)
    set.seed(2601)
    data <- data.frame(x = seq(-2, 2, length.out = 80L))
    data$y <- 1 + 0.5 * data$x + stats::rnorm(nrow(data))
    fit <- stats::lm(y ~ x, data = data)
    set.seed(2602)
    ordinary <- gof::cumres(fit, variable = "x", R = 19L, plots = 2L)
    set.seed(2602)
    piped <- fit |> gof::cumres(variable = "x", R = 19L, plots = 2L)
    stopifnot(identical(ordinary, piped))
    table <- piped |> as.data.frame()
    columns <- c("variable", "test", "statistic", "p.value", "n", "nobs", "R", "informative", "aggregation",
                 "ordering_source", "ordering_status", "ordering_reason", "resolution",
                 "distinct_values", "tie_blocks", "tied_units", "changed_units",
                 "unit", "n_rows", "n_trials", "fit_parameters", "estimated_parameters",
                 "fit_design_rank", "correction_dimension", "n_components")
    stopifnot(is.data.frame(table), identical(names(table), columns), nrow(table) == 2L,
              all(table$variable == "x"), all(table$n == 80L), all(table$nobs == 80L),
              all(table$R == 19L), all(is.finite(table$statistic)),
              all(table$p.value >= 0 & table$p.value <= 1),
              identical(table, piped |> summary() |> as.data.frame()))
    if (requireNamespace("generics", quietly = TRUE) && requireNamespace("tibble", quietly = TRUE)) {
        tidy <- generics::tidy(piped)
        stopifnot(is.data.frame(tidy), identical(as.data.frame(tidy), table))
    }
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    returned <- piped |> plot()
    stopifnot(identical(returned, piped))
    set.seed(2603)
    joint <- fit |> gof::cumres_joint(variable = c("x", "predicted"), R = 19L, plots = 2L)
    joint_table <- as.data.frame(joint)
    stopifnot(inherits(joint,"cumres_joint"), nrow(joint_table)==6L,
              identical(joint_table,as.data.frame(summary(joint))),
              identical(unique(joint_table$scope),c("marginal","family")),
              all(is.finite(joint$joint$p.value)), identical(invisible(plot(joint,idx="x")),joint))
    set.seed(2604)
    sem <- lava::estimate(lava::lvm(y~x),data)
    normality <- sem |> gof::cumres_normality(R = 19L, plots = 2L)
    normality_table <- as.data.frame(normality)
    stopifnot(inherits(normality,"cumres_normality"),
              identical(normality_table,as.data.frame(summary(normality))),
              all(normality_table$diagnostic=="normality"),
              identical(invisible(plot(normality)),normality))
    set.seed(2605)
    counts <- data.frame(x=seq_len(40),y=stats::rpois(40,.5))
    count_fit <- stats::glm(y~1,data=counts,family=stats::poisson())
    bootstrap <- count_fit |> gof::cumres_bootstrap(variable="x",data=counts,
        R=19L,plots=2L,null="specified",mean=rep(.5,40))
    bootstrap_table <- as.data.frame(bootstrap)
    stopifnot(inherits(bootstrap,"cumres_bootstrap"),
              identical(bootstrap_table,as.data.frame(summary(bootstrap))),
              all(bootstrap_table$calibration=="specified-null-monte-carlo"),
              length(bootstrap$bootstrap$status)==19L,
              identical(invisible(plot(bootstrap)),bootstrap))
    if (requireNamespace("generics", quietly = TRUE) && requireNamespace("tibble", quietly = TRUE)) {
        for (diagnostic in list(joint,normality,bootstrap))
            stopifnot(identical(as.data.frame(generics::tidy(diagnostic)),as.data.frame(diagnostic)))
    }
    cat("Installed R-interface checks passed for gof ", version, " at ", path, "\n", sep = "")
    cat("Help topics: ", length(topics), "; built vignettes: ", nrow(vignettes),
        "; native-pipe, table, summary and plot invariants passed for legacy and all three new endpoints.\n", sep = "")
})
