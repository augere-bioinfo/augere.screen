# library(testthat); library(augere.screen); source("test-generateSimesCommands.R")

test_that("generateSimesCommands works properly if we know it's signed", {
    set.seed(8888)
    res <- S4Vectors::DataFrame(PValue = runif(1000), LogFC = rnorm(1000))
    genes <- sample(LETTERS, nrow(res), replace=TRUE)

    cmds <- augere.screen:::.generateSimesCommands("FOO", "BAR", has.LogFC=TRUE)
    expect_false(any(grepl("no directionality", cmds)))
    env <- new.env()
    env$FOO <- res
    env$BAR <- genes
    tab <- eval(parse(text=cmds), envir=env)

    expect_identical(sort(rownames(tab)), sort(unique(genes)))
    expect_type(tab$NumBarcodes, "integer")
    expect_type(tab$PValue, "double")
    expect_type(tab$FDR, "double")
    expect_identical(sort(unique(tab$Direction)), c("down", "up"))

    # Checking it does the right thing when all genes are of the same sign, as this is easier to test.
    # In this situation, the sign is effectively irrelevant to the p-value calculation.
    res.up <- res
    res.up$LogFC <- abs(res.up$LogFC)
    env <- new.env()
    env$FOO <- res.up
    env$BAR <- genes
    tab.up <- eval(parse(text=cmds), envir=env)
    expect_true(all(tab.up$Direction == "up"))
    expect_false(identical(tab.up$PValue, tab$PValue))

    res.down <- res
    res.down$LogFC <- -abs(res.down$LogFC)
    env <- new.env()
    env$FOO <- res.down
    env$BAR <- genes
    tab.down <- eval(parse(text=cmds), envir=env)
    expect_true(all(tab.down$Direction == "down"))
    expect_identical(tab.up$PValue, tab.down$PValue)

    res.unknown <- res
    res.unknown$LogFC <- NULL
    cmds.unknown <- augere.screen:::.generateSimesCommands("FOO", "BAR", has.LogFC=FALSE) # and in fact it's equivalent to not having a sign at all.
    env <- new.env()
    env$FOO <- res.unknown
    env$BAR <- genes
    tab.unknown <- eval(parse(text=cmds.unknown), envir=env)
    expect_identical(tab.down$PValue, tab.unknown$PValue)
})

test_that("generateSimesCommands works properly if we don't know it's signed", {
    set.seed(666)
    res <- S4Vectors::DataFrame(PValue = runif(1000), LogFC = rnorm(1000))
    genes <- sample(LETTERS, nrow(res), replace=TRUE)

    cmds <- augere.screen:::.generateSimesCommands("FOO", "BAR", has.LogFC=FALSE)
    expect_true(any(grepl("no directionality", cmds)))
    env <- new.env()
    env$FOO <- res
    env$BAR <- genes
    tab <- eval(parse(text=cmds), envir=env)
    expect_type(tab$NumBarcodes, "integer")
    expect_type(tab$PValue, "double")
    expect_type(tab$FDR, "double")
    expect_identical(sort(unique(tab$Direction)), c("down", "up"))

    cmds.known <- augere.screen:::.generateSimesCommands("FOO", "BAR", has.LogFC=TRUE)
    env <- new.env()
    env$FOO <- res
    env$BAR <- genes
    tab.known <- eval(parse(text=cmds.known), envir=env)
    expect_identical(tab, tab.known)

    env <- new.env()
    res.nolfc <- res
    res.nolfc$LogFC <- NULL
    env$FOO <- res.nolfc
    env$BAR <- genes
    tab.nolfc <- eval(parse(text=cmds), envir=env)
    expect_false(identical(tab$PValue, tab.nolfc$PValue))
    expect_type(tab$NumBarcodes, "integer")
    expect_type(tab$PValue, "double")
    expect_type(tab$FDR, "double")
    expect_null(tab.nolfc$Direction)
})

test_that("generateSimesCommands works with a few NA values", {
    set.seed(7777)
    res <- S4Vectors::DataFrame(PValue = runif(1000), LogFC = rnorm(1000))
    res$PValue[c(1, 1000)] <- NA
    genes <- sample(LETTERS, nrow(res), replace=TRUE)

    cmds <- augere.screen:::.generateSimesCommands("FOO", "BAR", has.LogFC=TRUE)
    env <- new.env()
    env$FOO <- res
    env$BAR <- genes
    tab <- eval(parse(text=cmds), envir=env)

    env <- new.env()
    env$FOO <- res[-c(1, 1000),]
    env$BAR <- genes[-c(1, 1000)]
    tab.stripped <- eval(parse(text=cmds), envir=env)
    expect_identical(tab, tab.stripped)

    # But we still preserve the entry in the result.
    purged <- genes[c(1, 1000)]
    res$PValue[genes %in% purged] <- NA

    env <- new.env()
    env$FOO <- res
    env$BAR <- genes
    tab.purged <- eval(parse(text=cmds), envir=env)
    expect_identical(rownames(tab), rownames(tab.purged))
    all.p <- tab$PValue
    all.p[rownames(tab) %in% purged] <- NA
    expect_identical(all.p, tab.purged$PValue)
    all.num <- tab$NumBarcodes
    all.num[rownames(tab) %in% purged] <- 0L
    expect_identical(all.num, tab.purged$NumBarcodes)
})
