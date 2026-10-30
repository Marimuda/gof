#!/usr/bin/env Rscript
# Run from any directory: Rscript /path/to/gof/tools/qualify.R [output-directory]
# Dependencies are installed separately; this runner never modifies user libraries.
main <- function() {
    invocation <- commandArgs(FALSE)
    script <- sub("^--file=", "", invocation[grepl("^--file=", invocation)])
    if (length(script)!=1L) stop("Run this file with Rscript")
    if ("gof" %in% loadedNamespaces())
        stop("Run in a fresh R session; an existing gof namespace must not mask the isolated installation")
    package <- normalizePath(file.path(dirname(script),".."),winslash="/",mustWork=TRUE)
    args <- commandArgs(TRUE)
    if (length(args)>1L) stop("Supply at most one output directory")
    output <- if (length(args)) args[1] else file.path(dirname(package),"gof-validation",
        paste0(format(Sys.time(),"%Y%m%dT%H%M%S"),"-",Sys.getpid()))
    if (file.exists(output)) stop("Choose a new output directory; existing artifacts will not be overwritten")
    dir.create(output,recursive=TRUE)
    output <- normalizePath(output,winslash="/",mustWork=TRUE)
    if (startsWith(paste0(output,"/"),paste0(package,"/")))
        stop("Keep validation artifacts outside the source package")
    dir.create(file.path(output,"logs"))
    lib <- file.path(output,"lib")
    dir.create(lib)
    snapshot <- file.path(output,"source")
    dir.create(snapshot)
    entries <- list.files(package,all.files=TRUE,full.names=TRUE,no..=TRUE)
    entries <- entries[!basename(entries) %in% c(".git",".Rproj.user",".Rhistory")]
    if (!all(file.copy(entries,snapshot,recursive=TRUE,copy.mode=TRUE,copy.date=TRUE)))
        stop("Could not create the source snapshot")
    cpp <- Sys.getenv("GOF_CXX_STD","")
    if (nzchar(cpp)) {
        if (!cpp %in% c("CXX14","CXX17","CXX20","CXX23")) stop("GOF_CXX_STD must be CXX14, CXX17, CXX20 or CXX23")
        cat("\nCXX_STD = ",cpp,"\n",file=file.path(snapshot,"src","Makevars"),append=TRUE,sep="")
        descriptor <- file.path(snapshot,"DESCRIPTION")
        fields <- read.dcf(descriptor)
        target <- sub("^CXX","C++",cpp)
        requirements <- if ("SystemRequirements" %in% colnames(fields)) fields[1,"SystemRequirements"] else ""
        requirements <- if (grepl("C\\+\\+ ?[0-9]+",requirements)) {
            sub("C\\+\\+ ?[0-9]+",target,requirements)
        } else paste(c(target,requirements[nzchar(requirements)]),collapse=", ")
        if (!"SystemRequirements" %in% colnames(fields)) fields <- cbind(fields,SystemRequirements=requirements)
        else fields[1,"SystemRequirements"] <- requirements
        write.dcf(fields,descriptor)
    }
    enabled <- function(name,default=TRUE) {
        value <- tolower(Sys.getenv(name,if (default) "true" else "false"))
        if (!value %in% c("true","false")) stop(name," must be true or false")
        value=="true"
    }
    slow <- enabled("GOF_SLOW_TESTS",FALSE)
    native <- enabled("GOF_RUN_NATIVE",TRUE)
    check <- enabled("GOF_RUN_CHECK",TRUE)
    records <- list()
    complete <- FALSE
    previous <- getwd()
    setwd(output)
    on.exit(setwd(previous),add=TRUE)
    on.exit({
        manifest <- list(source=package,snapshot=snapshot,output=output,
                         requested=list(vignettes=TRUE,r_interface=TRUE,
                             calibration=slow,native=native,check=check,cpp=cpp,
                             native_compilers=Sys.getenv("GOF_NATIVE_COMPILERS",Sys.getenv("CXX","c++")),
                             native_leak_check=Sys.getenv("GOF_CHECK_LEAKS","0")),
                         completed=complete,phases=records,
                         check_environment=Sys.getenv(c("_R_CHECK_CRAN_INCOMING_",
                             "_R_CHECK_CRAN_INCOMING_REMOTE_","_R_CHECK_FORCE_SUGGESTS_",
                             "_R_CHECK_SYSTEM_CLOCK_")))
        saveRDS(manifest,file.path(output,"manifest.rds"))
        capture.output(dput(manifest),file=file.path(output,"manifest.txt"))
    },add=TRUE)
    run <- function(name,fun) {
        start <- proc.time()[3]
        logfile <- file.path(output,"logs",paste0(name,".log"))
        connection <- file(logfile,"wt")
        sink(connection); sink(connection,type="message")
        failure <- tryCatch({fun(logfile); NULL},error=function(e) {
            message(conditionMessage(e)); e
        })
        sink(type="message"); sink(); close(connection)
        records[[name]] <<- list(status=if (is.null(failure)) "passed" else "failed",
                                seconds=unname(proc.time()[3]-start),log=logfile)
        if (!is.null(failure)) stop("Phase ",name," failed; see ",logfile,call.=FALSE)
        invisible(NULL)
    }
    command <- function(executable,args,logfile) {
        # system2 quotes the executable; quote individual arguments which may
        # contain spaces (for example a library path on Windows).
        rawlog <- paste0(logfile,".command")
        status <- system2(executable,vapply(args,shQuote,character(1)),
                          stdout=rawlog,stderr=rawlog)
        if (file.exists(rawlog)) cat(readLines(rawlog,warn=FALSE),sep="\n")
        if (!identical(as.integer(status),0L)) stop("Command exited with status ",status)
    }
    rbin <- file.path(R.home("bin"),"R")
    oldlibs <- .libPaths()
    .libPaths(c(lib,oldlibs))
    Sys.setenv(R_LIBS=paste(unique(c(lib,oldlibs)),collapse=.Platform$path.sep))
    run("session",function(log) {
        print(sessionInfo())
        cat("Source:",package,"\nSnapshot:",snapshot,"\n")
        for (compiler in c("CXX","CXX17","CXX20","CXX23")) {
            cat(compiler,": ")
            # Older supported R versions do not define every newer standard.
            # Record the optional probe's output and status explicitly; actual
            # compilation below remains the gate for the requested standard.
            configuration <- suppressWarnings(system2(rbin,c("CMD","config",compiler),
                                                       stdout=TRUE,stderr=TRUE))
            print(configuration)
            status <- attr(configuration,"status")
            cat("Configuration probe exit status:",if (is.null(status)) 0L else status,"\n")
        }
    })
    run("build",function(log) command(rbin,c("CMD","build",snapshot),log))
    archive <- list.files(output,pattern="[.]tar[.]gz$",full.names=TRUE)
    if (length(archive)!=1L) stop("Expected exactly one source package archive")
    run("install",function(log) command(rbin,c("CMD","INSTALL",paste0("--library=",lib),archive),log))
    library("gof",lib.loc=lib,character.only=TRUE)
    run("installed-cpp-client",function(log)
        source(file.path(snapshot,"tools","check-client.R"),local=new.env(parent=globalenv())))
    run("installed-r-interface",function(log)
        source(file.path(snapshot,"tools","check-r-interface.R"),local=new.env(parent=globalenv())))
    test <- function(label,slow) run(label,function(log) {
        Sys.setenv(GOF_SLOW_TESTS=if (slow) "true" else "false",
                   GOF_CALIBRATION_REPORT=file.path(output,"calibration-rates.rds"),
                   GOF_SEM_CALIBRATION_REPORT=file.path(output,"calibration-sem-ordinal.rds"),
                   GOF_WEIGHT_CALIBRATION_REPORT=file.path(output,"calibration-weights.rds"),
                   GOF_JOINT_CALIBRATION_REPORT=file.path(output,"calibration-joint.rds"),
                   GOF_NORMALITY_CALIBRATION_REPORT=file.path(output,"calibration-normality.rds"),
                   GOF_BOOTSTRAP_CALIBRATION_REPORT=file.path(output,"calibration-bootstrap.rds"))
        result <- testthat::test_dir(file.path(snapshot,"tests","testthat"),
                                    reporter="summary",stop_on_failure=FALSE)
        saveRDS(result,file.path(output,paste0(label,".rds")))
        frame <- as.data.frame(result)
        cat("\nTest blocks:",nrow(frame),"expectations:",sum(frame$passed),"\n")
        if (any(frame$failed>0L | frame$error)) stop("Test failures recorded")
    })
    test("fast-tests",FALSE)
    if (slow) test("calibration-tests",TRUE) else records$calibration <- list(status="not requested")
    if (native) {
        compilers <- strsplit(Sys.getenv("GOF_NATIVE_COMPILERS",Sys.getenv("CXX","c++")),",",fixed=TRUE)[[1]]
        oldcxx <- Sys.getenv("CXX",unset=NA_character_)
        on.exit(if (is.na(oldcxx)) Sys.unsetenv("CXX") else Sys.setenv(CXX=oldcxx),add=TRUE)
        for (i in seq_along(compilers)) {
            compiler <- trimws(compilers[i])
            if (!nzchar(Sys.which(compiler))) stop("Native compiler unavailable: ",compiler)
            Sys.setenv(CXX=compiler)
            for (standard in c("c++14","c++17","c++20","c++23")) {
                run(paste0("native-",i,"-",standard),function(log) {
                    cat("Compiler:",compiler,"\nStandard:",standard,"\n")
                    print(system2(compiler,"--version",stdout=TRUE,stderr=TRUE))
                    command("bash",c(file.path(snapshot,"tools","check-native.sh"),standard),log)
                })
            }
        }
        if (is.na(oldcxx)) Sys.unsetenv("CXX") else Sys.setenv(CXX=oldcxx)
    } else records$native <- list(status="not requested")
    if (check) run("check",function(log) {
        Sys.setenv(GOF_SLOW_TESTS="false")
        command(rbin,c("CMD","check","--as-cran",archive),log)
        name <- read.dcf(file.path(snapshot,"DESCRIPTION"),fields="Package")[1]
        record <- readLines(file.path(output,paste0(name,".Rcheck"),"00check.log"),warn=FALSE)
        if (any(grepl("^Status:.*(ERROR|WARNING)",record))) stop("R CMD check reported errors or warnings")
    }) else records$check <- list(status="not requested")
    run("package-versions",function(log) {
        for (p in c("gof","lava","Rcpp","RcppArmadillo","testthat","numDeriv"))
            cat(p,as.character(packageVersion(p)),find.package(p),"\n")
        cat("Archive MD5:",unname(tools::md5sum(archive)),"\n")
        if (requireNamespace("digest",quietly=TRUE))
            cat("Archive SHA256:",digest::digest(file=archive,algo="sha256"),"\n")
        print(sessionInfo())
    })
    complete <- TRUE
    cat("Requested checks passed. Durable artifacts:",output,"\n")
    cat("See manifest.txt for requested, passed and omitted phases.\n")
}
main()
