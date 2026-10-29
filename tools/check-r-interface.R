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
                "as.data.frame.summary.cumres", "tidy.cumres", "surgunit")
    for (topic in topics) {
        if (!length(utils::help(topic, package = "gof", lib.loc = library_path)))
            stop("Installed help topic not found: ", topic)
    }
    vignettes <- utils::vignette(package = "gof", lib.loc = library_path)$results
    if (is.null(vignettes) || nrow(vignettes) < 2L)
        stop("Install the source archive after building its vignettes")
    html <- list.files(file.path(path, "doc"), pattern = "[.]html$", full.names = TRUE)
    if (length(html) < 2L) stop("The installed vignette HTML files are missing")

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
    columns <- c("variable", "test", "statistic", "p.value", "n", "nobs", "R", "informative", "aggregation")
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
    cat("Installed R-interface checks passed for gof ", version, " at ", path, "\n", sep = "")
    cat("Help topics: ", length(topics), "; built vignettes: ", nrow(vignettes),
        "; native-pipe, table, summary and plot invariants passed.\n", sep = "")
})
