result_api_fixture <- function() {
    structure(list(W=list(c(-.5,.25,0),c(.1,-.2,0)),
                   x=list(c(-2,0,3),c(10,20,30)),
                   What=list(matrix(numeric(),3,0),matrix(numeric(),3,0)),
                   KS=c(.25,NA_real_),CvM=c(.4,NA_real_),
                   variable=c("x","z"),n=20L,nobs=10L,R=199L,
                   type=rep("residual",2),model="glm",informative=c(TRUE,NA)),
              class="cumres")
}

test_that("result tables expose labeled tests and independently checked statistics", {
    g <- result_api_fixture()
    result <- as.data.frame(g)
    expect_identical(names(result),c("variable","test","statistic","p.value",
                                   "n","nobs","R","informative","aggregation"))
    expect_identical(result$variable,c("x","x","z","z"))
    expect_identical(result$test,c("Sup","L2","Sup","L2"))
    expect_equal(result$statistic,c(.5,sqrt(.25*2+.0625*3),.2,sqrt(.01*10+.04*10)))
    expect_identical(result$p.value,c(.25,.4,NA_real_,NA_real_))
    expect_equal(result$n,rep(20,4))
    expect_equal(result$nobs,rep(10,4))
    expect_equal(result$R,rep(199,4))
    expect_identical(result$informative,c(TRUE,TRUE,NA,NA))
    named <- as.data.frame(g,row.names=letters[1:4])
    expect_identical(rownames(named),letters[1:4])
})

test_that("summaries and saved matrix diagnostics retain their result meanings", {
    g <- result_api_fixture()
    legacy <- g
    legacy$W <- do.call(cbind,g$W); legacy$x <- do.call(cbind,g$x)
    expect_identical(as.data.frame(legacy),as.data.frame(g))
    expect_identical(as.data.frame(summary(g)),as.data.frame(g))
    saved <- summary(g); saved$table <- NULL
    result <- as.data.frame(saved)
    expect_true(all(is.na(result$statistic)))
    expect_identical(result$p.value,as.data.frame(g)$p.value)
    expect_identical(result$informative,as.data.frame(g)$informative)
})

test_that("missing information metadata does not erase computed p-values", {
    g <- result_api_fixture()
    g$informative <- NULL; g$nobs <- NULL
    result <- as.data.frame(g)
    expect_true(all(is.na(result$informative)))
    expect_true(all(is.na(result$nobs)))
    expect_identical(result$p.value,c(.25,.4,NA_real_,NA_real_))
    g$KS <- NULL
    expect_identical(as.data.frame(g)$p.value,c(NA_real_,.4,NA_real_,NA_real_))
    empty <- structure(list(variable=character()),class="cumres")
    expect_equal(dim(as.data.frame(empty)),c(0L,9L))
})

test_that("native pipes preserve seeded diagnostics and extraction consumes no random draws", {
    model <- fx_lm()
    set.seed(1771); plain <- cumres(model,R=15,plots=2)
    set.seed(1771); piped <- model |> cumres(R=15,plots=2)
    expect_identical(piped,plain)
    state <- .Random.seed
    table <- piped |> as.data.frame()
    from_summary <- piped |> summary() |> as.data.frame()
    expect_identical(table,from_summary)
    expect_identical(.Random.seed,state)
    expect_equal(table$p.value,as.vector(rbind(piped$KS,piped$CvM)))
    expect_equal(nrow(coef(piped)),length(piped$variable))
})

test_that("optional tidiers follow generics and broom conventions", {
    skip_if_not_installed("generics")
    skip_if_not_installed("tibble")
    g <- result_api_fixture()
    expected <- as.data.frame(g)
    result <- generics::tidy(g)
    expect_s3_class(result,"tbl_df")
    expect_identical(as.data.frame(result),expected)
    expect_identical(generics::tidy(summary(g)),result)
    if (requireNamespace("broom",quietly=TRUE))
        expect_identical(broom::tidy(g),result)
})
