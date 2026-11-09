## Structure layer for ordinal::clmm fits: bar parsing, factor relations,
## independent blocks, theta mapping and every validation path.

clmm_fit_quiet <- function(...) suppressWarnings(ordinal::clmm(...))

test_that("random-effect bars are parsed like clmm", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("nested",m=30L,nc=6L,seed=11L)
    cases <- list(
        list(f=y~x+(1|g),labels="1 | g",cols=list("(Intercept)")),
        list(f=y~x+(1+z|g),labels="1 + z | g",cols=list(c("(Intercept)","z"))),
        list(f=y~x+(z|g),labels="z | g",cols=list(c("(Intercept)","z"))),
        list(f=y~x+(0+z|g),labels="0 + z | g",cols=list("z")),
        list(f=y~x+(-1+z|g),labels="-1 + z | g",cols=list("z")),
        list(f=y~x+(1|g/h),labels=c("1 | h:g","1 | g"),cols=list("(Intercept)","(Intercept)")),
        list(f=y~x+(1|g)+(0+z|g),labels=c("1 | g","0 + z | g"),cols=list("(Intercept)","z")))
    for (cs in cases) {
        fit <- clmm_fit_quiet(cs$f,data=d)
        re <- gof_re_terms(cs$f,d)
        expect_identical(vapply(re$terms,`[[`,"","label"),cs$labels)
        expect_identical(lapply(re$terms,function(t) colnames(t$Z)),cs$cols)
        expect_identical(vapply(re$terms,`[[`,"","group"),names(fit$ST))
        expect_identical(lapply(re$terms,function(t) colnames(t$Z)),
                         unname(lapply(fit$ST,colnames)))
        expect_identical(lapply(re$terms,function(t) rownames(t$Z)),rep(list(NULL),length(re$terms)))
        for (t in re$terms) {
            gf <- fit$gfList[[t$group]]
            expect_identical(t$index,as.integer(gf))
            expect_identical(t$levels,levels(gf))
            expect_identical(t$n_levels,nlevels(gf))
            expect_identical(t$q,ncol(t$Z))
        }
        expect_equal(re$fixed,y~x,ignore_attr=TRUE)
    }
    ## Nested terms are named exactly as clmm names them.
    fit <- clmm_fit_quiet(y~x+(1|g/h),data=d)
    expect_identical(names(fit$gfList),c("h:g","g"))
    expect_identical(names(ordinal::VarCorr(fit)),c("h:g","g"))
    d$k <- factor(rep(1:2,length.out=nrow(d)))
    fit <- clmm_fit_quiet(y~x+(1|g/h/k),data=d)
    re <- gof_re_terms(y~x+(1|g/h/k),d)
    expect_identical(vapply(re$terms,`[[`,"","group"),names(fit$ST))
    expect_identical(names(fit$ST),c("k:(h:g)","h:g","g"))
    ## Fixed formula drops bars wherever they appear; offsets are kept.
    expect_equal(gof_re_terms(y~(1|g)+x+offset(z),d)$fixed,y~x+offset(z),ignore_attr=TRUE)
    expect_equal(gof_re_terms(y~(1|g),d)$fixed,y~1,ignore_attr=TRUE)
})

test_that("clmm's term ordering by number of levels is reproduced", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("nested",m=30L,nc=6L,seed=12L)
    for (f in list(y~x+(1|g)+(1|h:g),y~x+(1|h:g)+(1|g),y~x+(1|g)+(0+z|g)+(1|h:g))) {
        fit <- clmm_fit_quiet(f,data=d)
        re <- gof_re_terms(f,d)
        expect_identical(vapply(re$terms,`[[`,"","group"),names(fit$ST))
        expect_identical(lapply(re$terms,function(t) colnames(t$Z)),unname(lapply(fit$ST,colnames)))
    }
})

test_that("unsupported bar syntax is rejected clearly", {
    d <- clmm_simulate("intercept",m=10L,nc=4L)
    expect_error(gof_re_terms(y~x+(1+z||g),d),"Uncorrelated random-effect terms")
    expect_error(gof_re_terms(y~x,d),"no random-effect terms")
    expect_error(gof_re_terms(y~x+(1|nonexistent),d),"Cannot evaluate the grouping factor 'nonexistent'")
})

test_that("factor relations and independent blocks", {
    term <- function(group,index) list(label=paste("1 |",group),group=group,index=as.integer(index))
    ## Nested: classes within schools; blocks are the schools.
    s <- rep(1:5,each=6)
    c <- rep(1:10,each=3)
    terms <- list(term("c:s",c),term("s",s))
    rel <- gof_factor_relations(terms)
    expect_identical(rel["c:s","s"],"nested_in")
    expect_identical(rel["s","c:s"],"contains")
    expect_identical(unname(diag(rel)),rep("identical",2L))
    b <- gof_independent_blocks(length(s),terms)
    expect_identical(b,list(id=s,m=5L))
    ## Blocks are numbered by first appearance.
    s2 <- c(3,3,1,1,2,2,3)
    b <- gof_independent_blocks(7L,list(term("s",s2)))
    expect_identical(b,list(id=c(1L,1L,2L,2L,3L,3L,1L),m=3L))
    ## Identical partitions with different labels.
    rel <- gof_factor_relations(list(term("a",s),term("b",6L-s)))
    expect_identical(rel["a","b"],"identical")
    ## Terms sharing one factor appear once.
    rel <- gof_factor_relations(list(term("s",s),term("s",s)))
    expect_identical(dim(rel),c(1L,1L))
    ## Fully crossed: a single block.
    a <- rep(1:4,times=5); bb <- rep(1:5,each=4)
    rel <- gof_factor_relations(list(term("a",a),term("b",bb)))
    expect_identical(rel["a","b"],"crossed")
    expect_identical(gof_independent_blocks(20L,list(term("a",a),term("b",bb)))$m,1L)
    ## Two disconnected crossed designs: two blocks.
    a2 <- c(a,a+4L); b2 <- c(bb,bb+5L)
    blk <- gof_independent_blocks(40L,list(term("a",a2),term("b",b2)))
    expect_identical(blk,list(id=rep(1:2,each=20),m=2L))
    ## A chain of links across a crossed design (long diameter) connects fully.
    n <- 50L
    a3 <- c(1:n,1:n); b3 <- c(1:n,c(2:n,1L))
    blk <- gof_independent_blocks(2L*n,list(term("a",a3),term("b",b3)))
    expect_identical(blk$m,1L)
    ## Path that breaks into components in a non-trivial order.
    a4 <- c(1,2,3,1,4,4); b4 <- c(1,2,2,3,5,6)
    blk <- gof_independent_blocks(6L,list(term("a",a4),term("b",b4)))
    expect_identical(blk$id,c(1L,2L,2L,1L,3L,3L))
})

test_that("theta reproduces VarCorr, thresholds and fixed effects", {
    skip_if_not_installed("ordinal")
    check <- function(fit) {
        p <- suppressWarnings(gof_clmm_problem(fit))
        u <- gof_clmm_unpack(p)
        vc <- ordinal::VarCorr(fit)
        expect_identical(length(u$re),length(vc))
        for (t in seq_along(vc))
            expect_equal(u$re[[t]]$Sigma,unclass(vc[[t]])[,],tolerance=1e-10,ignore_attr=TRUE)
        expect_equal(u$alpha,c(fit$Theta),tolerance=1e-12,ignore_attr=TRUE)
        expect_identical(u$xi,fit$alpha)
        expect_identical(names(p$theta)[p$index$xi],names(fit$alpha))
        expect_equal(u$beta,fit$beta,tolerance=0)
        expect_identical(colnames(p$X),names(fit$beta))
        expect_identical(p$source$theta_clmm,p$theta)
        for (t in p$terms)
            expect_identical(names(p$theta)[t$omega],paste0(t$label,".",gof_logchol_names(t$q)))
        p
    }
    d <- clmm_simulate("slope",m=40L,nc=8L,seed=21L)
    check(clmm_fit_quiet(y~x+(1|g),data=d))
    p <- check(clmm_fit_quiet(y~x+(1+z|g),data=d))
    expect_identical(names(p$theta)[p$terms[[1]]$omega],c("1 + z | g.L11","1 + z | g.L21","1 + z | g.L22"))
    check(clmm_fit_quiet(y~x+(z|g),data=d))
    check(clmm_fit_quiet(y~x+(0+z|g),data=d))
    p <- check(clmm_fit_quiet(y~x+(1|g)+(0+z|g),data=d))
    expect_identical(vapply(p$terms,`[[`,"","label"),c("1 | g","0 + z | g"))
    expect_identical(p$chain,"g")
    dn <- clmm_simulate("nested",m=40L,nc=6L,seed=22L)
    p <- check(clmm_fit_quiet(y~x+(1|g/h),data=dn))
    expect_identical(names(p$theta)[5:6],c("1 | h:g.L11","1 | g.L11"))
    ## Factor covariate.
    dn$f3 <- factor(rep(c("a","b","c"),length.out=nrow(dn)))
    p <- check(clmm_fit_quiet(y~x+f3+(1|g),data=dn))
    expect_identical(colnames(p$X),c("x","f3b","f3c"))
    expect_identical(p$X[,"f3c"],as.numeric(dn$f3=="c"))
    ## Threshold structures (odd and even numbers of cutpoints).
    for (K in c(4L,5L)) {
        dk <- clmm_simulate("intercept",m=40L,nc=6L,K=K,seed=23L)
        for (th in c("flexible","symmetric","symmetric2","equidistant")) {
            fit <- clmm_fit_quiet(y~x+(1|g),data=dk,threshold=th)
            p <- check(fit)
            expect_identical(p$tJac,fit$tJac)
            expect_identical(p$K,K)
        }
    }
})

test_that("theta mapping is invariant to the sign of clmm's Cholesky factor", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("slope",m=40L,nc=8L,seed=24L)
    fit <- clmm_fit_quiet(y~x+(1+z|g),data=d)
    flip <- fit
    flip$ST[[1]] <- fit$ST[[1]] %*% diag(c(-1,1))
    dimnames(flip$ST[[1]]) <- dimnames(fit$ST[[1]])
    expect_true(flip$ST[[1]][1,1]<0)
    expect_equal(ordinal::VarCorr(flip)[[1]],ordinal::VarCorr(fit)[[1]],ignore_attr=TRUE)
    p <- gof_clmm_problem(flip)
    expect_equal(gof_clmm_unpack(p)$re[[1]]$Sigma,unclass(ordinal::VarCorr(fit)[[1]])[,],
                 tolerance=1e-10,ignore_attr=TRUE)
    expect_equal(p$theta,gof_clmm_problem(fit)$theta,tolerance=1e-12)
    ## Scalar term with a negative standard deviation.
    fit <- clmm_fit_quiet(y~x+(1|g),data=d)
    neg <- fit; neg$ST[[1]][] <- -fit$ST[[1]]
    p <- gof_clmm_problem(neg)
    expect_equal(gof_clmm_unpack(p)$re[[1]]$Sigma,unclass(ordinal::VarCorr(fit)[[1]])[,],
                 tolerance=1e-10,ignore_attr=TRUE)
    ## A singular covariance (variance exactly zero) is regularised, finite and recorded.
    zero <- fit; zero$ST[[1]][] <- 0
    p <- gof_clmm_problem(zero)
    expect_true(all(is.finite(p$theta)))
    expect_identical(p$source$boundary,"1 | g")
})

test_that("problem fields: response, design, weights and blocks", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=31L)
    fit <- clmm_fit_quiet(y~x+(1|g),data=d)
    p <- gof_clmm_problem(fit)
    expect_s3_class(p,"gof_clmm_problem")
    expect_identical(p$y,as.integer(d$y))
    expect_identical(p$levels,levels(d$y))
    expect_identical(p$X,cbind(x=d$x))
    expect_identical(p$offset,numeric(nrow(d)))
    expect_identical(p$weights,rep(1,nrow(d)))
    expect_identical(p$blocks,list(id=as.integer(d$g),m=40L))
    expect_identical(p$structure,"nesting-chain")
    expect_identical(p$link,"logit")
    expect_identical(p$source$dropped,0L)
    expect_identical(p$source$call,fit$call)
    expect_identical(p$source$data,fit$model)
    ## (clmm itself fails for the cauchit link on some data sets; seed 32 fits.)
    dl <- clmm_simulate("intercept",m=40L,nc=6L,seed=32L)
    for (lk in c("probit","cloglog","loglog","cauchit"))
        expect_identical(gof_clmm_problem(clmm_fit_quiet(y~x+(1|g),data=dl,link=lk))$link,lk)
})

test_that("weights are read from the model frame and zero weights dropped", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=32L)
    w <- rep(c(1,2,0.5),length.out=nrow(d)); w[1:6] <- 0   # whole first group has zero weight
    w[10] <- 0
    d$w <- w
    fit <- suppressWarnings(ordinal::clmm(y~x+(1|g),data=d,weights=w))
    expect_true("(weights)" %in% names(fit$model))
    p <- gof_clmm_problem(fit)
    keep <- w>0
    expect_identical(p$source$dropped,7L)
    expect_identical(p$n,sum(keep))
    expect_identical(p$n,fit$dims$n)
    expect_identical(p$N,sum(w))
    expect_identical(p$weights,w[keep])
    expect_identical(p$y,as.integer(d$y)[keep])
    expect_identical(p$X,cbind(x=d$x[keep]))
    expect_identical(p$terms[[1]]$n_levels,39L)
    expect_identical(p$blocks$m,39L)
    expect_identical(p$source$rows,which(keep))
})

test_that("formula offsets are used and offset arguments are rejected", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=33L)
    d$o <- 0.5*d$x
    fit <- clmm_fit_quiet(y~x+offset(o)+(1|g),data=d)
    p <- gof_clmm_problem(fit)
    expect_identical(p$offset,d$o)
    ## Offset enters like x'beta: the offset 0.5*x shifts beta by -0.5.
    fit0 <- clmm_fit_quiet(y~x+(1|g),data=d)
    expect_equal(unname(fit$beta),unname(fit0$beta)-0.5,tolerance=1e-4)
    ## clmm reads offset() positionally, so a trailing offset picks up the wrong column.
    bad <- clmm_fit_quiet(y~x+(1|g)+z+offset(o),data=d)
    expect_error(gof_clmm_problem(bad),"wrong model-frame column.*offset\\(\\) before the random-effect terms")
    ## offset= as an argument is silently ignored by clmm.
    fo <- suppressWarnings(ordinal::clmm(y~x+(1|g),data=d,offset=d$o))
    expect_identical(fo$call$offset,quote(d$o))
    expect_equal(fo$logLik,fit0$logLik)
    expect_error(gof_clmm_problem(fo),"'offset' argument, which clmm silently ignores")
})

test_that("invalid fits are rejected with specific messages", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("slope",m=40L,nc=6L,seed=34L)
    fit <- clmm_fit_quiet(y~x+(1+z|g),data=d)
    expect_error(gof_clmm_problem(list()),"fitted ordinal::clmm")
    expect_error(gof_clmm_problem(fit,control=list()),"gof_clmm_control")
    ## Nonconvergence.
    nc <- fit; nc$optRes$convergence <- 1L; nc$optRes$message <- "false convergence (8)"
    expect_error(gof_clmm_problem(nc),"did not converge \\(optimizer convergence code 1: false convergence")
    ng <- fit; ng$gradient[] <- 0.5
    expect_error(gof_clmm_problem(ng),"max \\|gradient\\| = 0.5")
    ## Unobserved category (all of its rows with zero weight).
    uc <- fit; uc$model[["(weights)"]] <- ifelse(uc$model$y=="4",0,1)
    expect_error(gof_clmm_problem(uc),"Response category '4' has no observations with positive weight")
    ## Dimension cap.
    expect_error(gof_clmm_problem(fit,control=gof_clmm_control(max.dim=1L)),
                 "\\(1 \\+ z \\| g\\) has dimension 2 per level, more than max.dim = 1")
    ## Missing model frame and bad link.
    nm <- fit; nm$model <- NULL
    expect_error(gof_clmm_problem(nm),"refit with model = TRUE")
    bl <- fit; bl$link <- "foo"
    expect_error(gof_clmm_problem(bl),"Unsupported clmm link 'foo'")
    ## Too few independent blocks.
    small <- clmm_simulate("intercept",m=5L,nc=12L,seed=35L)
    fs <- clmm_fit_quiet(y~x+(1|g),data=small)
    expect_error(gof_clmm_problem(fs),"Only 5 independent blocks \\(levels of 'g'\\) for 5 parameters")
    few <- clmm_simulate("intercept",m=12L,nc=8L,seed=36L)
    expect_warning(p <- gof_clmm_problem(clmm_fit_quiet(y~x+(1|g),data=few)),
                   "Only 12 independent blocks")
    expect_identical(p$blocks$m,12L)
})

test_that("crossed designs are classified, not rejected", {
    skip_if_not_installed("ordinal")
    d <- clmm_simulate("intercept",m=40L,nc=6L,seed=41L)
    d$z2 <- factor(rep(1:6,times=40L))
    p <- gof_clmm_problem(clmm_fit_quiet(y~x+(1|g)+(1|z2),data=d))
    expect_identical(p$structure,"crossed")
    expect_null(p$chain)
    expect_identical(p$relations["g","z2"],"crossed")
    expect_identical(p$blocks$m,1L)
})

test_that("ordinal::wine examples", {
    skip_if_not_installed("ordinal")
    wine <- ordinal::wine
    expect_warning(p <- gof_clmm_problem(ordinal::clmm(rating~temp+contact+(1|judge),data=wine)),
                   "Only 9 independent blocks")
    expect_identical(p$structure,"nesting-chain")
    expect_identical(p$chain,"judge")
    expect_identical(colnames(p$X),c("tempwarm","contactyes"))
    expect_identical(p$K,5L)
    ## A random slope builds here (boundary checks happen when polishing)...
    fit <- ordinal::clmm(rating~temp+(1+temp|judge),data=wine)
    expect_warning(p <- gof_clmm_problem(fit),"Only 9 independent blocks")
    expect_identical(colnames(p$terms[[1]]$Z),c("(Intercept)","tempwarm"))
    ## Its correlation is 1 (singular Sigma): regularised by a tiny ridge and recorded.
    expect_identical(p$source$boundary,"1 + temp | judge")
    expect_true(all(is.finite(p$theta)))
    expect_lt(max(abs(gof_clmm_unpack(p)$re[[1]]$Sigma-unclass(ordinal::VarCorr(fit)[[1]])[,])),1e-9)
    ## ...but with contact as well there are as many parameters (9) as judges.
    expect_error(gof_clmm_problem(ordinal::clmm(rating~temp+contact+(1+temp|judge),data=wine)),
                 "Only 9 independent blocks \\(levels of 'judge'\\) for 9 parameters")
    p <- gof_clmm_problem(ordinal::clmm(rating~temp+contact+(1|judge)+(1|bottle),data=wine))
    expect_identical(p$structure,"crossed")
    expect_identical(p$relations["judge","bottle"],"crossed")
})
