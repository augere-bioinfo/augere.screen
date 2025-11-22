#' @import augere.core
#' @importFrom augere.de processCustomDesignMatrix processCustomContrasts processSimpleDesignMatrix processSimpleComparisons findSubsetGroups
.initialize <- function(
    x,
    method,
    assay,
    groups, 
    comparisons, 
    covariates,
    block,
    subset.factor,
    subset.levels,
    subset.groups,
    design,
    contrasts,
    filter.reference.factor,
    filter.reference.levels,
    filter.num.mads,
    filter.default,
    author
) {
    template <- system.file("templates", "initialize.Rmd", package="augere.screen", mustWork=TRUE)
    parsed <- parseRmdTemplate(readLines(template))

    parsed[["create-se"]] <- processInputCommands(x, name="se")

    replacements <- list(
        METHOD=method,
        ASSAY=deparseToString(assay),
        AUTHOR=paste(sprintf("  - %s", author), collapse="\n")
    )

    #########################################
    ### Sample subsetting, design matrix. ###
    #########################################

    if (!is.null(subset.factor)) {
        parsed[["subset-se"]] <- replacePlaceholders(
            parsed[["subset-se"]],
            list(
                SUBSET_FACTOR=deparseToString(subset.factor),
                SUBSET_LEVELS=deparseToString(subset.levels)
            )
        )
    } else {
        parsed[["subset-se"]] <- NULL
    }

    mds.use.group <- FALSE
    if (!is.null(design) && !is.null(contrasts)) {
        parsed[["design-matrix"]] <- processCustomDesignMatrix(design=design, se.name="se")
        contrast.info <- processCustomContrasts(contrasts)
        replacements$FILTER_OPTS <- "design=design"
        parsed[["subset-group"]] <- NULL

    } else {
        parsed[["design-matrix"]] <- processSimpleDesignMatrix(groups=groups, block=block, covariates=covariates, se.name="se")
        contrast.info <- processSimpleComparisons(comparisons)

        if (is.null(groups)) {
            replacements$FILTER_OPTS <- "design=design"
            parsed[["subset-group"]] <- NULL
        } else {
            replacements$FILTER_OPTS <- "group=model.data$group."
            mds.use.group <- TRUE

            group.levels <- NULL
            if (subset.groups) {
                group.levels <- findSubsetGroups(contrast.info)
            }
            if (!is.null(group.levels)) {
                parsed[["subset-group"]] <- replacePlaceholders(
                    parsed[["subset-group"]],
                    list(
                        GROUP_FACTOR=deparseToString(groups),
                        GROUP_LEVELS=deparseToString(group.levels)
                    )
                )
            } else {
                parsed[["subset-group"]] <- NULL
            }
        }
    }

    #########################
    ### Barcode filtering ###
    #########################

    some.filtering <- FALSE
    if (is.null(filter.reference.factor) || is.null(filter.reference.levels)) {
        parsed[["ref-filtering"]] <- NULL
    } else {
        some.filtering <- TRUE
        parsed[["ref-filtering"]] <- replacePlaceholders(
            parsed[["ref-filtering"]], 
            list(
                REFERENCE_FIELD=deparseToString(filter.reference.factor),
                REFERENCE_LEVELS=deparseToString(filter.reference.levels),
                NUM_MADS=deparseToString(filter.num.mads)
            ),
            error=FALSE
        )
    }

    if (!filter.default) {
        parsed[["filter-default"]] <- NULL
    } else {
        some.filtering <- TRUE
    }
    if (some.filtering) {
        parsed[["filter-none"]] <- NULL
    }

    list(
        text=replacePlaceholders(parsed, replacements),
        contrasts=contrast.info,
        mds.use.group=mds.use.group
    )
}
