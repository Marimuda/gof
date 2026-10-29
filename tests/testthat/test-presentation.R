presentation_fixture <- function() {
    structure(list(W=list(c(-.5,.25,0),c(.1,-.2,0)),
                   x=list(c(-2,0,3),c(10,20,30)),
                   What=list(matrix(c(-.3,.1,0,.2,-.1,0),3,2),
                             matrix(c(.05,-.1,0,-.05,.1,0),3,2)),
                   KS=c(.25,.75),CvM=c(.4,.6),R=20L,n=3L,nobs=3L,
                   variable=c("x","z"),type=c("residual","score"),
                   informative=c(TRUE,FALSE),response="Outcome",model="glm"),
              class="cumres")
}

capture_presentation <- function(code) {
    calls <- new.env(parent=emptyenv())
    calls$plot <- calls$lines <- calls$matlines <- calls$legend <- list()
    recorder <- function(name) function(...) {
        calls[[name]][[length(calls[[name]])+1L]] <- list(...)
        invisible(NULL)
    }
    testthat::local_mocked_bindings(lines=recorder("lines"),
                                   matlines=recorder("matlines"),legend=recorder("legend"),
                                   .package="graphics")
    ## plot() is defined in base on current R, and in graphics on older R.
    testthat::local_mocked_bindings(plot=recorder("plot"),
                                   .package=getNamespaceName(environment(graphics::plot)))
    result <- withVisible(force(code))
    list(calls=as.list(calls),value=result$value,visible=result$visible)
}

test_that("plot forwards graphical arguments and returns the original object invisibly", {
    g <- presentation_fixture()
    seen <- capture_presentation(gof:::plot.cumres(g,idx="x",legend="none",
                                                 xlim=c(-4,5),ylim=c(-1,1),
                                                 axes=FALSE,xaxs="i",cex.axis=.7,
                                                 main="Supplied main"))
    expect_false(seen$visible)
    expect_identical(seen$value,g)
    expect_length(seen$calls$plot,1L)
    call <- seen$calls$plot[[1]]
    expect_equal(call$x,g$x[[1]])
    expect_equal(call$y,g$W[[1]])
    expect_equal(call$xlim,c(-4,5))
    expect_equal(call$ylim,c(-1,1))
    expect_identical(call$axes,FALSE)
    expect_identical(call$xaxs,"i")
    expect_equal(call$cex.axis,.7)
    expect_identical(call$main,"Supplied main")
    expect_length(seen$calls$lines,1L)
    expect_length(seen$calls$legend,0L)
})

test_that("sample palettes reach the plot and the type2 legend describes samples", {
    g <- presentation_fixture()
    colors <- c("red","blue")
    seen <- capture_presentation(gof:::plot.cumres(g,idx=1,col=colors,col.alpha=.4,legend="type2"))
    expect_length(seen$calls$matlines,1L)
    expect_equal(seen$calls$matlines[[1]][[1]],g$x[[1]])
    expect_equal(seen$calls$matlines[[1]][[2]],g$What[[1]])
    expect_equal(seen$calls$matlines[[1]]$col,grDevices::adjustcolor(colors,alpha.f=.4))
    expect_length(seen$calls$legend,1L)
    labels <- seen$calls$legend[[1]][[2]]
    expect_true("Observed" %in% labels)
    expect_true(any(grepl("Sample process under null",labels,fixed=TRUE)))
    expect_equal(length(seen$calls$legend[[1]]$col),length(labels))
    expect_equal(length(seen$calls$legend[[1]]$lwd),length(labels))
})

test_that("disabled or absent samples do not leave a misleading legend entry", {
    for (colors in list(NULL,"none",character())) {
        g <- presentation_fixture()
        seen <- capture_presentation(gof:::plot.cumres(g,idx=1,col=colors,legend="type2"))
        expect_length(seen$calls$matlines,0L)
        expect_length(seen$calls$lines,1L)
        expect_identical(seen$calls$legend[[1]][[2]],"Observed")
    }
    g <- presentation_fixture()
    g$What <- lapply(g$What,function(w) matrix(numeric(),nrow(w),0))
    seen <- capture_presentation(gof:::plot.cumres(g,idx=1,legend="type2"))
    expect_length(seen$calls$matlines,0L)
    expect_identical(seen$calls$legend[[1]][[2]],"Observed")
})

test_that("observed curve and legend agree under a nonblack graphics colour", {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(),add=TRUE)
    previous <- graphics::par(col="darkgreen")
    on.exit(graphics::par(previous),add=TRUE,after=FALSE)
    g <- presentation_fixture()
    seen <- capture_presentation(gof:::plot.cumres(g,idx=1,legend="type2"))
    observed <- seen$calls$lines[[1]]$col
    legend_observed <- seen$calls$legend[[1]]$col[1]
    expect_length(observed,1L)
    expect_equal(grDevices::col2rgb(observed),grDevices::col2rgb("darkgreen"))
    expect_equal(grDevices::col2rgb(legend_observed),grDevices::col2rgb(observed))
})

test_that("plot selects named or numeric processes in the requested order", {
    g <- presentation_fixture()
    named <- capture_presentation(gof:::plot.cumres(g,idx="z",legend="none"))
    expect_length(named$calls$plot,1L)
    expect_equal(named$calls$plot[[1]]$x,g$x[[2]])
    ordered <- capture_presentation(gof:::plot.cumres(g,idx=c(2,1),legend="none"))
    expect_equal(ordered$calls$plot[[1]]$x,g$x[[2]])
    expect_equal(ordered$calls$plot[[2]]$x,g$x[[1]])
    empty <- capture_presentation(gof:::plot.cumres(g,idx=integer()))
    expect_length(empty$calls$plot,0L)
    expect_identical(empty$value,g)
    expect_false(empty$visible)
    expect_error(gof:::plot.cumres(g,idx="missing"),"Unknown process")
    for (idx in list(0,-1,1.5,3,NA_real_,Inf,TRUE,list(1)))
        expect_error(gof:::plot.cumres(g,idx=idx),"Invalid process")
})

test_that("labels titles observation axes and score axes retain their meanings", {
    g <- presentation_fixture()
    default <- capture_presentation(gof:::plot.cumres(g,idx=1,legend="none"))
    expect_identical(default$calls$plot[[1]]$xlab,"x")
    expect_true(is.expression(default$calls$plot[[1]]$ylab))
    expect_identical(default$calls$plot[[1]]$ylab,expression(W["x"](x)))
    expect_identical(default$calls$plot[[1]]$main,"Outcome")
    custom <- capture_presentation(gof:::plot.cumres(g,idx=1,legend="none",
                     xlab=expression(log(x)),ylab="Residual label",title="Explicit title",main="Extra main"))
    expect_identical(custom$calls$plot[[1]]$xlab,expression(log(x)))
    expect_identical(custom$calls$plot[[1]]$ylab,"Residual label")
    expect_identical(custom$calls$plot[[1]]$main,"Explicit title")
    g$response <- NULL
    score <- capture_presentation(gof:::plot.cumres(g,idx=2,legend="none"))
    expect_identical(score$calls$plot[[1]]$xlab,"Time")
    expect_identical(score$calls$plot[[1]]$main,"z")
    expect_identical(score$calls$plot[[1]]$ylab,expression(W["z"](x)))
    observations <- capture_presentation(gof:::plot.cumres(g,idx=2,vs=FALSE,legend="none"))
    expect_equal(observations$calls$plot[[1]]$x,seq_along(g$x[[2]]))
    expect_identical(observations$calls$plot[[1]]$xlab,"Observation")
    expect_equal(observations$calls$lines[[1]][[1]],seq_along(g$x[[2]]))
})

test_that("type1 legends report available tests and degenerate diagnostics", {
    g <- presentation_fixture()
    seen <- capture_presentation(gof:::plot.cumres(g,idx=2,legend=TRUE))
    labels <- seen$calls$legend[[1]][[2]]
    expect_true(any(grepl("Sup-test p=0.75",labels,fixed=TRUE)))
    expect_true(any(grepl("L2-test p=0.6",labels,fixed=TRUE)))
    expect_true("No lack-of-fit information" %in% labels)
    g$CvM <- NULL
    single <- capture_presentation(gof:::plot.cumres(g,idx=1,legend="type1"))
    expect_length(single$calls$legend[[1]][[2]],1L)
    expect_match(single$calls$legend[[1]][[2]],"Sup-test")
    for (setting in list(FALSE,NULL,"none")) {
        hidden <- capture_presentation(gof:::plot.cumres(g,idx=1,legend=setting))
        expect_length(hidden$calls$legend,0L)
    }
})

test_that("plot validates options and process dimensions before drawing", {
    g <- presentation_fixture()
    for (alpha in list(-.1,1.1,NA_real_,Inf,"a"))
        expect_error(gof:::plot.cumres(g,col.alpha=alpha),"col.alpha")
    for (vs in list(1,logical(),NA,c(TRUE,FALSE)))
        expect_error(gof:::plot.cumres(g,vs=vs),"vs must")
    expect_error(gof:::plot.cumres(g,col=c("red","none")),"alone")
    expect_error(gof:::plot.cumres(g,legend="unknown"))
    expect_error(gof:::plot.cumres(g,idx=1,col="purple",col.alpha=.3,
                    legend="none",xlab="x",ylab="y",vs=TRUE,ylim=c(-1,1),title="",42),"must be named")
    for (component in c("W","x","What")) {
        bad <- g
        bad[[component]][[1]][1] <- NA_real_
        expect_error(gof:::plot.cumres(bad,idx=1),"process data")
    }
    bad <- g; bad$x[[1]] <- 1:2
    expect_error(gof:::plot.cumres(bad,idx=1),"process data")
    bad <- g; bad$What[[1]] <- 1:3
    expect_error(gof:::plot.cumres(bad,idx=1),"process data")
})

test_that("current and saved matrix objects render through real PDF devices", {
    g <- presentation_fixture()
    legacy <- g
    legacy$W <- do.call(cbind,g$W)
    legacy$x <- do.call(cbind,g$x)
    path <- tempfile(fileext=".pdf")
    grDevices::pdf(path)
    on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
    for (object in list(g,legacy)) {
        for (style in c("type1","type2","none")) {
            expect_no_error(result <- withVisible(plot(object,idx=c("z","x"),legend=style,
                                                       col=c("red","blue"),axes=FALSE)))
            expect_false(result$visible)
            expect_identical(result$value,object)
        }
    }
    expect_no_error(plot(g,col=NULL))
    expect_no_error(plot(g,col="none"))
    expect_no_error(plot(g,idx="x",xlab=expression(log(x)),ylab=expression(W(x)),title="Custom labels"))
})

test_that("coef and print handle each absent test statistic", {
    g <- presentation_fixture()
    expect_identical(rownames(coef(g)),g$variable)
    expect_identical(colnames(coef(g)),c("p-value(Sup)","p-value(L2)"))
    expect_equal(unname(coef(g)[[1]]),g$KS)
    expect_equal(unname(coef(g)[[2]]),g$CvM)
    for (statistic in c("KS","CvM")) {
        one <- g; one[[statistic]] <- NULL
        expect_equal(dim(coef(one)),c(2L,1L))
        expect_identical(rownames(coef(one)),g$variable)
        output <- capture.output(result <- withVisible(print(one)))
        expect_true(any(grepl("Based on 20 realizations",output,fixed=TRUE)))
        expect_false(result$visible)
        expect_identical(result$value,one)
    }
    curves <- g; curves$KS <- curves$CvM <- NULL
    expect_equal(dim(coef(curves)),c(2L,0L))
    expect_no_error(capture.output(print(curves)))
})

test_that("summary preserves metadata and empty summaries do not print phantom rows", {
    g <- presentation_fixture()
    s <- summary(g)
    expect_s3_class(s,"summary.cumres")
    expect_identical(s$nobs,g$nobs)
    expect_identical(s$response,g$response)
    expect_identical(s$informative,g$informative)
    output <- capture.output(result <- withVisible(print(s)))
    expect_false(result$visible)
    expect_identical(result$value,s)
    expect_true(any(grepl("ordered by 'x'",output,fixed=TRUE)))
    expect_true(any(grepl("ordered by 'z'",output,fixed=TRUE)))
    empty <- structure(list(KS=NULL,CvM=NULL,R=20,n=0,type=character(),
                            model="glm",variable=character()),class="summary.cumres")
    empty_output <- capture.output(result <- withVisible(print(empty)))
    expect_false(result$visible)
    expect_identical(result$value,empty)
    expect_false(any(grepl("ordered by|NA|p-value",empty_output)))
})
