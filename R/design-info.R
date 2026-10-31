# Descriptive context only: no statistical support rule or RNG belongs here.
gof_design_info <- function(unit="residual",n_rows=NA_real_,n_trials=NA_real_,
                            fit_parameters=NA_real_,estimated_parameters=NA_real_,
                            fit_design_rank=NA_real_,correction_dimension=0,
                            n_components=1) {
    if (!is.character(unit) || length(unit)!=1L || is.na(unit) ||
        !unit %in% c("row","trial","subject","residual"))
        stop("Invalid diagnostic unit metadata",call.=FALSE)
    counts <- list(n_rows=n_rows,n_trials=n_trials,fit_parameters=fit_parameters,
        estimated_parameters=estimated_parameters,fit_design_rank=fit_design_rank,
        correction_dimension=correction_dimension,n_components=n_components)
    counts <- lapply(counts,function(z) {
        if (is.null(z)) return(NA_real_)
        if (!is.numeric(z) || is.complex(z) || !is.null(dim(z)) || length(z)!=1L ||
            (!is.na(z) && (!is.finite(z) || z<0 || z!=floor(z))))
            stop("Design metadata counts must be nonnegative whole numbers or missing",call.=FALSE)
        as.double(z)
    })
    as.data.frame(c(list(unit=unit),counts),stringsAsFactors=FALSE)
}

gof_print_design <- function(object) {
    info <- object$design_info
    if (is.null(info)) return(invisible(NULL))
    cat("Design: ",object$n," ",info$unit," units",sep="")
    if (!is.na(info$n_rows)) cat("; ",info$n_rows," data rows",sep="")
    if (!is.na(info$estimated_parameters))
        cat("; ",info$estimated_parameters," parameters adjusted/refitted",sep="") else
        cat("; estimated parameter count unreported")
    if (!is.na(info$fit_design_rank)) cat("; fit design rank ",info$fit_design_rank,sep="")
    cat(".\n")
    invisible(NULL)
}
