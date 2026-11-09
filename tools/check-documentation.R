#!/usr/bin/env Rscript
# Check generated documentation and build a local website from a source copy.
# Nothing is installed in the user's library or deployed to a hosting service.
main <- function() {
    invocation <- commandArgs(FALSE)
    script <- sub("^--file=", "", invocation[grepl("^--file=", invocation)])
    if (length(script) != 1L) stop("Run this file with Rscript")
    if ("gof" %in% loadedNamespaces()) stop("Run documentation checks in a fresh R session")
    package <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
    args <- commandArgs(TRUE)
    if (length(args) > 1L) stop("Supply at most one new output directory")
    output <- if (length(args)) args[1] else file.path(dirname(package), "gof-documentation",
        paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "-", Sys.getpid()))
    if (file.exists(output)) stop("Choose a new output directory; existing evidence is not overwritten")
    dir.create(output, recursive = TRUE)
    output <- normalizePath(output, winslash = "/", mustWork = TRUE)
    if (startsWith(paste0(output, "/"), paste0(package, "/")))
        stop("Keep documentation artifacts outside the source package")
    enabled <- function(name, default = TRUE) {
        value <- tolower(Sys.getenv(name, if (default) "true" else "false"))
        if (!value %in% c("true", "false")) stop(name, " must be true or false")
        value == "true"
    }
    website <- enabled("GOF_DOCS_BUILD_SITE")
    required <- c("roxygen2", "rmarkdown", if (website) "pkgdown")
    for (dependency in required) {
        if (!requireNamespace(dependency, quietly = TRUE))
            stop("Install documentation dependency in a development library: ", dependency)
    }
    snapshot <- file.path(output, "source")
    dir.create(snapshot)
    entries <- list.files(package, all.files = TRUE, full.names = TRUE, no.. = TRUE)
    entries <- entries[!basename(entries) %in% c(".git", ".Rproj.user", ".Rhistory", "docs", ".local", ".claude-flow")]
    if (!all(file.copy(entries, snapshot, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)))
        stop("Could not copy the documentation source")
    inputs <- function(directory) {
        files <- list.files(directory, recursive = TRUE, all.files = TRUE, no.. = TRUE)
        files <- files[grepl("^(R|src|inst|man|vignettes|tools)/|^[^/]+[.]md$|^(DESCRIPTION|NAMESPACE|README[.]Rmd|_pkgdown[.]yml|LICENSE)$", files)]
        files <- files[!grepl("([.](o|so|dll|a|gcno|gcda)$|(^|/)Rplots[.]pdf$)", files)]
        stats::setNames(unname(tools::md5sum(file.path(directory, files))), files)
    }
    input_md5 <- inputs(snapshot)
    text <- function(path) sub("[[:space:]]+$", "", readLines(path, warn = FALSE))
    expected_docs <- c("NAMESPACE", file.path("man", list.files(file.path(snapshot, "man"), pattern = "[.]Rd$")))
    expected_text <- stats::setNames(lapply(file.path(snapshot, expected_docs), text), expected_docs)
    expected_descriptor <- read.dcf(file.path(snapshot, "DESCRIPTION"))
    expected_readme <- if (file.exists(file.path(snapshot, "README.md"))) text(file.path(snapshot, "README.md")) else NULL
    expected_figures <- list.files(file.path(snapshot, "man", "figures"), pattern = "^README-")
    dir.create(file.path(output, "logs"))
    records <- list()
    complete <- FALSE
    on.exit({
        manifest <- list(source = package, snapshot = snapshot, output = output,
                         completed = complete, requested = list(website = website), phases = records,
                         input_md5 = input_md5,
                         dependency_versions = vapply(required, function(p)
                             as.character(utils::packageVersion(p)), character(1)),
                         readme_figure_check = "Generated file names and nonempty content; device-dependent PNG bytes are not compared")
        saveRDS(manifest, file.path(output, "manifest.rds"))
        capture.output(dput(manifest), file = file.path(output, "manifest.txt"))
    }, add = TRUE)
    run <- function(name, action) {
        log <- file.path(output, "logs", paste0(name, ".log"))
        connection <- file(log, "wt")
        sink(connection)
        sink(connection, type = "message")
        start <- proc.time()[3]
        failure <- tryCatch({ action(); NULL }, error = function(e) {
            message(conditionMessage(e)); e
        })
        sink(type = "message")
        sink()
        close(connection)
        records[[name]] <<- list(status = if (is.null(failure)) "passed" else "failed",
                                seconds = unname(proc.time()[3] - start), log = log)
        if (!is.null(failure)) stop("Documentation phase ", name, " failed; see ", log, call. = FALSE)
    }
    run("session", function() {
        print(sessionInfo())
        print(rmarkdown::find_pandoc())
        cat("Source:", package, "\nCopied source:", snapshot, "\n")
    })
    run("roxygen-freshness", function() {
        roxygen2::roxygenise(snapshot, roclets = c("rd", "namespace"), clean = TRUE)
        expected <- expected_docs
        actual <- c("NAMESPACE", file.path("man", list.files(file.path(snapshot, "man"), pattern = "[.]Rd$")))
        differences <- union(setdiff(expected, actual), setdiff(actual, expected))
        for (file in intersect(expected, actual)) {
            if (!identical(expected_text[[file]], text(file.path(snapshot, file))))
                differences <- c(differences, file)
        }
        before <- expected_descriptor
        after <- read.dcf(file.path(snapshot, "DESCRIPTION"))
        generator_version <- function(fields) {
            for (field in c("Config/roxygen2/version", "RoxygenNote"))
                if (field %in% colnames(fields)) return(unname(fields[1, field]))
            NA_character_
        }
        if (is.na(generator_version(before)) ||
            !identical(generator_version(before), generator_version(after)))
            differences <- c(differences, "DESCRIPTION:roxygen-generator-version")
        writeLines(unique(differences), file.path(output, "roxygen-differences.txt"))
        if (length(differences)) stop("Regenerate and commit documentation: ", paste(unique(differences), collapse = ", "))
        cat("All", length(expected), "generated help/namespace files match their roxygen source.\n")
    })
    run("readme-freshness", function() {
        if (!file.exists(file.path(snapshot, "README.Rmd")) || is.null(expected_readme))
            stop("Both README.Rmd and its generated README.md are required")
        figures <- function(directory) list.files(file.path(directory, "man", "figures"),
                                                  pattern = "^README-", full.names = TRUE)
        old_figures <- figures(snapshot)
        if (length(old_figures)) unlink(old_figures)
        rmarkdown::render(file.path(snapshot, "README.Rmd"), output_format = "github_document",
                          output_file = "README.md", envir = new.env(parent = globalenv()), quiet = FALSE)
        if (!identical(expected_readme, text(file.path(snapshot, "README.md"))))
            stop("README.md is stale; render README.Rmd and commit the result")
        actual_figures <- figures(snapshot)
        if (!setequal(expected_figures, basename(actual_figures)) ||
            any(file.info(actual_figures)$size <= 0))
            stop("Generated README figure inventory is stale or contains an empty file")
        cat("README text and", length(actual_figures), "generated figure files are current.\n")
    })
    if (website) run("local-website", function() {
        pkgdown::build_site(snapshot, examples = TRUE, preview = FALSE,
                            devel = FALSE, new_process = TRUE, install = TRUE, quiet = FALSE,
                            override = list(destination = file.path(output, "site")))
        if (!file.exists(file.path(output, "site", "index.html"))) stop("Website index was not built")
        if (!file.exists(file.path(output, "site", "reference", "index.html")))
            stop("Website reference index was not built")
        if (!length(list.files(file.path(output, "site", "articles"), pattern = "[.]html$")))
            stop("Website articles were not built")
        flush(getConnection(sink.number(type = "message")))
        rendered_log <- readLines(file.path(output, "logs", "local-website.log"), warn = FALSE)
        if (any(grepl("[WARNING]", rendered_log, fixed = TRUE)))
            stop("Pandoc emitted rendering warnings; inspect the website log before sharing the generated documentation")
        cat("Local website built without deployment at", file.path(output, "site"), "\n")
    }) else records$website <- list(status = "not requested")
    run("source-unchanged", function() {
        if (!identical(input_md5, inputs(package)))
            stop("Documentation inputs changed while checks were running; rerun from a stable source")
        cat("All", length(input_md5), "source/documentation inputs still match their captured hashes.\n")
    })
    complete <- TRUE
    cat("Documentation checks passed. Durable artifacts:", output, "\n")
}
main()
