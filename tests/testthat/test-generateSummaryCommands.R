# library(testthat); library(augere.screen); source("test-generateSummaryCommands.R")

set.seed(123123)
res <- S4Vectors::DataFrame(PValue = rbeta(1000, 1, 20), LogFC = rnorm(1000))
res$FDR <- p.adjust(res$PValue, method="BH")
genes <- sample(LETTERS, nrow(res), replace=TRUE)
by.gene <- split(seq_along(genes), genes)

test_that("generateSummaryCommands works properly if we know it's signed", {
    cmds <- augere.screen:::.generateSummaryCommands("FOO", "BAR", 0.05, has.LogFC=TRUE)
    expect_false(any(grepl("NumSig", cmds)))

    env <- new.env()
    env$FOO <- res
    env$BAR <- by.gene
    out <- eval(parse(text=cmds), envir=env)

    expect_identical(names(by.gene), rownames(out))
    expect_type(out$NumUp, "integer")
    expect_type(out$NumDown, "integer")

    expect_true(all(out$NumUp + out$NumDown <= lengths(by.gene)))
    expect_gt(sum(res$FDR <= 0.05), 0) # check that we're actually testing something.
    expect_identical(sum(out$NumUp) + sum(out$NumDown), sum(res$FDR <= 0.05))

    expect_type(out$LogFC, "double")
    expect_false(anyNA(out$LogFC))
})

test_that("generateSummaryCommands works properly if it might be signed", {
    cmds <- augere.screen:::.generateSummaryCommands("FOO", "BAR", 0.05, has.LogFC=FALSE)
    expect_true(any(grepl("NumSig", cmds)))

    env <- new.env()
    env$FOO <- res
    env$BAR <- by.gene
    out <- eval(parse(text=cmds), envir=env)

    ref.cmds <- augere.screen:::.generateSummaryCommands("FOO", "BAR", 0.05, has.LogFC=TRUE)
    ref.out <- eval(parse(text=ref.cmds), envir=env)
    expect_identical(out, ref.out)

    res$LogFC <- NULL
    env$FOO <- res
    out <- eval(parse(text=cmds), envir=env)
    expect_type(out$NumSig, "integer")
    expect_true(all(out$NumSig <= lengths(by.gene)))
    expect_identical(sum(out$NumSig), sum(res$FDR <= 0.05))
})

test_that("generateSummaryCommands works properly after sprinkling in some NAs", {
    fail <- logical(nrow(res))
    scapegoat <- genes[1]
    fail[genes == scapegoat] <- TRUE
    fail[sample(nrow(res), 20)] <- TRUE
    res$PValue[fail] <- NA
    res$FDR[fail] <- NA
    res$LogFC[fail] <- NA

    cmds <- augere.screen:::.generateSummaryCommands("FOO", "BAR", 0.05, has.LogFC=FALSE)

    {
        env <- new.env()
        env$FOO <- res
        env$BAR <- by.gene
        out <- eval(parse(text=cmds), envir=env)
        expect_false(anyNA(out$NumUp))
        expect_false(anyNA(out$NumDown))
        expect_identical(out[scapegoat,"NumUp"], 0L)
        expect_identical(out[scapegoat,"NumDown"], 0L)
        expect_true(is.na(out[scapegoat,"LogFC"]))

        env <- new.env()
        env$FOO <- res[!fail,]
        ok.genes <- genes[!fail]
        env$BAR <- split(seq_along(ok.genes), ok.genes)
        ref.out <- eval(parse(text=cmds), envir=env)
        ref.out <- ref.out[rownames(out),,drop=FALSE]
        rownames(ref.out) <- rownames(out)
        ref.out$NumUp[is.na(ref.out$NumUp)] <- 0L
        ref.out$NumDown[is.na(ref.out$NumDown)] <- 0L
        expect_identical(out, ref.out)
    }

    {
        res0 <- res
        res0$LogFC <- NULL

        env <- new.env()
        env$FOO <- res0
        env$BAR <- by.gene
        out <- eval(parse(text=cmds), envir=env)
        expect_false(anyNA(out$NumSig))
        expect_identical(out[scapegoat,"NumSig"], 0L)

        env <- new.env()
        env$FOO <- res0[!fail,]
        ok.genes <- genes[!fail]
        env$BAR <- split(seq_along(ok.genes), ok.genes)
        ref.out <- eval(parse(text=cmds), envir=env)
        ref.out <- ref.out[rownames(out),,drop=FALSE]
        rownames(ref.out) <- rownames(out)
        ref.out$NumSig[is.na(ref.out$NumSig)] <- 0L
        expect_identical(out, ref.out)
    }
})

