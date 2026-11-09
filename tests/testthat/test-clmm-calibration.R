skip_if_not(identical(Sys.getenv("GOF_SLOW_TESTS"),"true"),
            "slow tests: set GOF_SLOW_TESTS=true")
skip_if_not_installed("ordinal")
skip_if_not_installed("numDeriv")

test_that("the frozen slow CLMM cohort retains calibration and refusal evidence", {
    source(test_path("..","..","tools","clmm-study.R"),local=TRUE)
    output <- Sys.getenv("GOF_CLMM_CALIBRATION_REPORT","")
    directory <- if (nzchar(output)) paste0(output,".d") else tempfile("clmm-calibration-")
    workers <- as.integer(Sys.getenv("GOF_CLMM_WORKERS","1"))
    report <- gof_clmm_study_run(directory,mode="slow",replicates=500L,R=199L,workers=workers)
    if (nzchar(output)) saveRDS(report,file=output)
    expect_length(report$records,4L*500L)
    expect_true(all(vapply(report$records,function(x) is.finite(x$seconds),TRUE)))
    for (id in report$protocol$cells$id) {
        rows <- report$summary[report$summary$id==id & report$summary$alpha==.05 &
                               report$summary$stratum=="all",]
        expect_equal(rows$attempted,c(500L,500L))
        expect_true(all(rows$accepted+round(rows$refusal_rate*500L)==500L))
        expect_true(all(rows$lower95<=rows$rate & rows$rate<=rows$upper95))
        ## Acceptance is prespecified. Failure remains a failed gate; the
        ## observations and numerical refusal reasons are retained.
        expect_true(all(rows$rate>=.03 & rows$rate<=.075),info=paste(id,"size"))
        expect_true(all(rows$refusal_rate<.05),info=paste(id,"refusal"))
    }
})
