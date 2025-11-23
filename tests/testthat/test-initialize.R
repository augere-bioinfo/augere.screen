# library(testthat); library(augere.screen); source("test-initialize.R")

library(SummarizedExperiment)
cd <- DataFrame(group = rep(LETTERS[1:4], each=5), age=1:20, batch=factor(rep(1:5, 4)))
means <- 2^runif(1000, 1, 10)
means[1] <- 0 # to test the filtering.
mat <- matrix(rnbinom(length(means) * nrow(cd), mu=means, size=10), nrow=length(means), ncol=nrow(cd))
se <- SummarizedExperiment(list(counts=mat), colData=cd)

test_that(".initialize() works with 'default' options", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.screen:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups="group",
        comparisons=c("A", "B"),
        covariates=NULL,
        block=NULL,
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=TRUE,
        design=NULL,
        contrasts=NULL,
        filter.reference.factor=NULL,
        filter.reference.levels=NULL,
        filter.num.mads=3,
        filter.default=TRUE,
        author="Chihaya Kisaragi"
    )

    # Checking that the contrasts are okay.
    expect_identical(output$contrasts[[1]]$title, "Increase in `A` over `B`")
    expect_identical(output$contrasts[[1]]$type, "versus")
    expect_identical(output$contrasts[[1]]$left, list("A"))
    expect_identical(output$contrasts[[1]]$right, list("B"))
    expect_true(output$mds.use.group)

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)

    expect_identical(env$se, se[,se$group %in% c("A", "B")])
    expect_identical(colnames(env$design), c("group.A", "group.B"))

    expect_true(inherits(env$y, "DGEList"))
    expect_type(env$y$genes$origin, "integer") 
    expect_lt(nrow(env$y), nrow(se))
    expect_false(1 %in% env$y$genes$origin) 

    all.text <- unlist(output$text)
    expect_true(any(grepl("filterByExpr\\(.*group=", all.text)))
    expect_false(any(grepl("omit the filtering", all.text)))
})

test_that(".initialize() works with no filtering by group", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.screen:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups="group",
        comparisons=c("A", "B"),
        covariates=NULL,
        block=NULL,
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=FALSE,
        design=NULL,
        contrasts=NULL,
        filter.reference.factor=NULL,
        filter.reference.levels=NULL,
        filter.num.mads=3,
        filter.default=TRUE,
        author="Chihaya Kisaragi"
    )

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)

    expect_identical(env$se, se)
})

test_that(".initialize() works with custom sample filtering", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.screen:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups="group",
        comparisons=c("A", "B"),
        covariates=NULL,
        block="batch",
        subset.factor="batch",
        subset.levels=c(2,4),
        subset.groups=TRUE,
        design=NULL,
        contrasts=NULL,
        filter.reference.factor=NULL,
        filter.reference.levels=NULL,
        filter.num.mads=3,
        filter.default=TRUE,
        author="Chihaya Kisaragi"
    )

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)

    expect_identical(env$se, se[,se$group %in% c("A", "B") & se$batch %in% c(2, 4)])
})

test_that(".initialize() works for a design matrix without groups", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.screen:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups=NULL,
        comparisons=c("age"),
        covariates="age",
        block=NULL,
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=TRUE,
        design=NULL,
        contrasts=NULL,
        filter.reference.factor=NULL,
        filter.reference.levels=NULL,
        filter.num.mads=3,
        filter.default=TRUE,
        author="Chihaya Kisaragi"
    )

    expect_false(output$mds.use.group)

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)
    expect_identical(env$se, se) # no filtering is done on the samples.
    expect_true(any(grepl("filterByExpr\\(.*design=", unlist(output$text))))
})

test_that(".initialize() works with covariates + group but only the covariates are tested", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.screen:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups="group",
        comparisons="age",
        covariates="age",
        block=NULL,
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=FALSE,
        design=NULL,
        contrasts=NULL,
        filter.reference.factor=NULL,
        filter.reference.levels=NULL,
        filter.num.mads=3,
        filter.default=TRUE,
        author="Chihaya Kisaragi"
    )

    expect_true(output$mds.use.group)

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)
    expect_identical(env$se, se) # no filtering by sample.
    expect_true(any(grepl("filterByExpr\\(.*group=", unlist(output$text))))
})

test_that(".initialize() works with custom matrices and contrasts", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.screen:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=TRUE,
        design=~group + batch,
        contrasts="groupC - groupD",
        filter.reference.factor=NULL,
        filter.reference.levels=NULL,
        filter.num.mads=3,
        filter.default=TRUE,
        author="Chihaya Kisaragi"
    )

    expect_false(output$mds.use.group)

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)
    expect_identical(env$se, se) # no filtering is done on the samples.
    expect_true(any(grepl("filterByExpr\\(.*design=", unlist(output$text))))

    # Checking that the contrasts are okay.
    expect_identical(output$contrasts[[1]]$type, "custom")
})

test_that(".initialize() works with reference-based filtering", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.screen:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups="group",
        comparisons=c("B", "A"),
        covariates=NULL,
        block=NULL,
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=TRUE,
        design=NULL,
        contrasts=NULL,
        # Works when the reference samples aren't even part of the object anymore.
        filter.reference.factor="group",
        filter.reference.levels=c("C", "D"),
        filter.num.mads=0.1,
        filter.default=FALSE,
        author="Chihaya Kisaragi"
    )

    all.text <- unlist(output$text)
    expect_false(any(grepl("omit the filtering", all.text)))
    expect_false(any(grepl("filterByExpr", all.text))) # still have the default filtering enabled as well.

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)
    expect_identical(sort(unique(env$se$group)), c("A", "B"))
    expect_lt(nrow(env$y), nrow(se))

    expect_type(env$ref.non.empty, "integer")
    expect_false(1 %in% env$ref.non.empty)
    expect_type(env$ref.filtered, "logical")
    expect_false(env$ref.filtered[1])
    expect_gt(env$ref.med.ab, env$ref.threshold)
    expect_false(identical(env$filtered, env$ref.filtered))
    expect_false(1 %in% env$y$genes$origin)
})

test_that(".initialize() works with no filtering", {
    fun <- augere.core::resetInputCache()
    on.exit(fun(), add=TRUE, after=TRUE)

    output <- augere.screen:::.initialize(
        x=se, 
        method="edgeR",
        assay="counts",
        groups="group",
        comparisons=c("B", "A"),
        covariates=NULL,
        block=NULL,
        subset.factor=NULL,
        subset.levels=NULL,
        subset.groups=TRUE,
        design=NULL,
        contrasts=NULL,
        filter.reference.factor=NULL,
        filter.reference.levels=NULL,
        filter.num.mads=3,
        filter.default=FALSE,
        author="Chihaya Kisaragi"
    )

    all.text <- unlist(output$text)
    expect_true(any(grepl("omit the filtering", all.text)))
    expect_false(any(grepl("med.ab", all.text)))
    expect_false(any(grepl("filterByExpr", all.text)))

    env <- new.env()
    tmp <- tempfile()
    dir.create(tmp)
    augere.core::compileReport(file=file.path(tmp, "report.Rmd"), contents=output$text, env=env)
    expect_identical(nrow(env$y), nrow(se))
    expect_true(1 %in% env$y$genes$origin)
})
