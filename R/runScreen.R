#' Differential abundance for functional screens
#'
#' Use \code{\link[limma]{voom}} to test for differential abundance of barcodes in a functional screen.
#' 
#' @param se A \link[SummarizedExperiment]{SummarizedExperiment} object containing read counts for each barcode (row) and sample (column).
#' Alternatively, the output of \code{\link[augere.core]{wrapInput}} that refers to a SummarizedExperiment.
#' @inheritParams augere.de::runVoom
#' @param filter.reference.factor String specifying \code{\link[SummarizedExperiment]{colData}(se)} column that indicates whether a sample is a reference.
#' This is used to filter out barcodes that are not present in the original barcode library.
#' @param filter.reference.levels Character vector containing the reference levels in the \code{filter.reference.factor} column.
#' Only used if \code{filter.reference.factor} is supplied.
#' @param filter.num.mads Number specifying the number of median absolute deviations below the median to define an outlier threshold for filtering barcodes.
#' Only used if \code{filter.reference.factor} is supplied.
#' @param filter.default Boolean indicating whether to filter barcodes using the default \pkg{edgeR} method, i.e., \code{\link[edgeR]{filterByExpr}}.
#' If \code{TRUE} and \code{filter.reference.factor} is supplied, the default filtering is applied in addition to the reference-based filtering.
#' Otherwise, it is applied alone. 
#' @param norm.barcodes A named list of character vectors containing one or more sets of control barcodes.
#' Each vector should contain the names of barcodes, corresponding to the row names of \code{se}.
#' The names of the vectors will be stored in the metadata.
#' @param norm.genes A named list containing one or more sets of control genes.
#' Each vector should contain the names of genes, corresponding to the entries of the \code{gene.field} column.
#' The names of the vectors will be stored in the metadata.
#' Only used if \code{gene.field} is supplied.
#' @param norm.type.field String specifying the \code{\link[SummarizedExperiment]{rowData}(se)} column containing the feature type for each barcode.
#' @param norm.types Character vector containing the control feature types in the \code{norm.type.field} column. 
#' Only used if \code{norm.type.field} is supplied.
#' @param norm.tmm Boolean indicating whether TMM normalization should be used.
#' If \code{FALSE}, normalization is instead performed using the library size for each sample.
#' @param gene.field String specifying the \code{rowData(se)} column that contains the gene identity for each barcode.
#' @param method Character vector specifying the consolidation method(s) to convert per-barcode statistics into per-gene results.
#' This can be zero, one or more of:
#' \itemize{
#' \item \code{"simes"}, which uses \code{\link[metapod]{groupedSimes}} to combine barcode-level p-values into a single per-gene p-value. 
#' This is the most sensitive approach for detecting a small number of barcodes (possibly just 1) that are differentially abundant.
#' It is appropriate for CRISPRa where very few guides are expected to succeed.
#' It is also useful for identifying genes containing barcodes that exhibit sporadic off-target effects in CRISPR knockout or CRISPRi.
#' \item \code{"holm-min"}, which uses \code{\link[metapod]{groupedHolmMin}} to combine barcode-level p-values into a single per-gene p-value. 
#' This should be used when we expect a minimum number or proportion of barcodes to exhibit differential abundance for a given gene.
#' It is suitable for CRISPR knockout screens, CRISPRi and CRISPRa.
#' \item \code{"fry"}, which uses \code{\link[limma]{fry}} to test for differential expression across all barcodes for a gene.
#' This tends to favor consistent differences (in the same direction!) across all barcodes for a given gene, and thus is suitable for CRISPR knockout screens or CRISPRi.
#' This method mode is not compatible with ANOVA-like contrasts and will ignore any \code{lfc.threshold > 0}.
#' }
#' Only used if \code{gene.field} is supplied.
#' @param holm.min.n Integer specifying the minimum number of significant barcodes per gene when \code{method="holm"},
#' see the \code{min.sig.n=} in \code{?\link[metapod]{groupedHolmMin}} for details.
#' @param holm.min.prop Integer specifying the minimum proportion of significant barcodes per gene when \code{method="holm"},
#' see the \code{min.sig.prop=} in \code{?\link[metapod]{groupedHolmMin}} for details.
#' @param fry.args Named list of additional arguments to pass to \code{\link[limma]{fry}}.
#' @param gene.data Character vector containing the names of \code{\link[SummarizedExperiment]{rowData}(se)} columns to add to the per-gene result tables.
#' If a gene has multiple barcodes, the row corresponding to the first barcode is used.
#' 
#' @return
#' A Rmarkdown report named \code{report.Rmd} is written inside \code{output.dir} that contains the analysis commands.
#'
#' If \code{dry.run=FALSE}, a list is returned containing:
#' \itemize{
#' \item \code{barcode}, a named list of \link[S4Vectors]{DataFrame}s containing barcode-level result tables.
#' Each DataFrame corresponds to a comparison/contrast where each row corresponds to a barcode (i.e., row) in \code{se}.
#' Each DataFrame contains the following columns:
#' \itemize{
#' \item \code{AveExpr}, the average abundance.
#' \item \code{t}, the F-statistic.
#' (For non-ANOVA-like contrasts only.)
#' \item \code{F}, the F-statistic.
#' (For ANOVA-like contrasts only.)
#' \item \code{LogFC}, the log-fold change.
#' (For non-ANOVA-like contrasts only.)
#' \item \code{LogFC.<COLUMN>}, the log-fold change corresponding to each column of the contrast matrix.
#' (For ANOVA-like contrasts only.)
#' \item \code{PValue}, the p-value;
#' \item \code{FDR}, the Benjamini-Hochberg-adjusted p-value.
#' }
#' \item \code{gene}, a named list of lists containing one entry per method in \code{method}.
#' Each inner list contains one \link[S4Vectors]{DataFrame} per contrast where each row corresponds to a gene (not barcode). 
#' Each DataFrame contains at least the following columns:
#' \itemize{
#' \item \code{NumBarcodes}, the number of barcodes per gene (possibly after filtering).
#' \item \code{PValue}, the p-value for differential abundance across all barcodes for the gene.
#' \item \code{FDR}, the Benjamini-Hochberg-adjusted p-value.
#' \item \code{Direction}, the overall direction of change for all barcodes of the gene.
#' (For non-ANOVA-like contrasts only.)
#' }
#' Only reported if \code{gene.field} is not \code{NULL}.
#' \item \code{normalized}, a copy of \code{x} with normalized expression values, possibly subsetted by sample.
#' This contains:
#' \itemize{
#' \item \code{lib.size} and \code{norm.factors} columns in its \code{\link[SummarizedExperiment]{colData}},
#' containing the library sizes and normalization factors, respectively.
#' \item a \code{retained} column in its \code{\link[SummarizedExperiment]{rowData}},
#' indicating whether a gene was retained after filtering.
#' \item a \code{logcounts} assay, containing the log-normalized counts.
#' }
#' \code{normalized} may be subsetted by sample, depending on \code{subset.factor}, \code{subset.group}, etc.
#' }
#' If \code{save.results=TRUE}, the results are saved in a \code{results} directory inside \code{output}.
#'
#' If \code{dry.run=TRUE}, \code{NULL} is returned.
#' Only the Rmarkdown report is saved to file.
#'
#' @section Reference filtering:
#' A reference sample should represent the barcode distribution before any perturbation, e.g., time zero in a time course.
#' All barcodes should have similar abundances in the reference samples as they should be present at the same molarity in the original barcode library.
#' If a barcode has very low abundance, it was probably missing from the original library due to some manufacturing defect.
#'
#' To identify these low outliers, we compute the average abundance across reference samples for each barcode.
#' We define a filtering threshold based on the median absolute deviation (MAD) below the median of the average abundances across barcodes.
#' Barcodes with average abundances below this threshold are discarded.
#' 
#' % This approach is a little suspect if the reference samples are being compared as it's not independent of the comparison.
#' % Specifically, we could enrich for barcodes that are spuriously upregulated in the reference.
#' % Fortunately, the density of barcodes at the outlier filter boundary should be so low that any biases should not be noticeable.
#'
#' We also use \code{\link[edgeR]{filterByExpr}} to remove low-abundance barcodes. 
#' This ensures that all remaining barcodes have sufficient counts for \code{\link[limma]{voom}}, 
#' especially if the reference samples themselves are not used in the rest of the analysis.
#'
#' @section Control normalization:
#' Control barcodes might not target any gene (non-targeting controls, or NTCs) or they might target non-essential genes (NEGs).
#' We assume that such barcodes do not genuinely change in abundance across conditions, allowing us to use them to compute normalization factors.
#' The choice of control barcodes can be made by one or more of the \code{norm.barcodes}, \code{norm.genes} and \code{norm.types} arguments.
#' If more than one of these modes is provided, the union of barcodes is used. 
#'
#' By default, we use TMM normalization (see \code{\link[edgeR]{calcNormFactors}}) on the control barcodes. 
#' This avoids composition biases and provides some robustness against changes in abundance due to inadvertent biological activity.
#' If \code{norm.use.tmm=FALSE}, we instead use the total sum of counts across the controls to derive normalization factors.
#' This can be more precise but is more sensitive to genuine changes in abundance.
#'
#' If no control barcodes are specified, normalization uses all barcodes that remaining after filtering. 
#'
#' @author Aaron Lun, Jean-Philippe Fortin
#'
#' @examples
#' # Creating an example dataset.
#' mu <- 2^runif(1000, 0, 5)
#' grouping <- rep(c("A", "B", "ctrl"), each=3)
#' counts <- matrix(rnbinom(length(mu)* length(grouping), mu=mu, size=20), ncol=length(grouping))
#' rownames(counts) <- sprintf("BARCODE-%s", seq_along(mu))
#' 
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(list(counts=counts))
#' rowData(se)$type <- sample(c("NTC", "NEG", "other"), length(mu), replace=TRUE)
#' rowData(se)$gene <- paste0("GENE-", sample(LETTERS, length(mu), replace=TRUE))
#' colData(se)$group <- grouping
#'
#' output <- tempfile()
#' res <- runScreen(
#'     se,
#'     groups="group",
#'     comparisons=c("A", "B"),
#'     norm.control.type.field="type",
#'     norm.control.types="NEG",
#'     filter.reference.factor="group",
#'     filter.reference.levels="ctrl",
#'     gene.field="gene",
#'     output.dir=output
#' )
#'
#' list.files(output, recursive=TRUE)
#' res$barcode[[1]]
#' res$gene$simes[[1]]
#' res$normalized
#' 
#' @export
#' @import augere.core
#' @importFrom augere.de processContrastMetadata
runScreen <- function( 
    x,
    groups,
    comparisons,
    covariates = NULL,
    block = NULL,
    subset.factor = NULL,
    subset.levels = NULL,
    subset.groups = TRUE,
    design = NULL,
    contrasts = NULL,
    dc.block = NULL,
    robust = TRUE,
    trend = FALSE,
    quality = TRUE,
    lfc.threshold = 0,
    filter.reference.factor = NULL,
    filter.reference.levels = NULL,
    filter.num.mads = 3,
    filter.default = TRUE,
    norm.control.barcodes = NULL, 
    norm.control.genes = NULL,
    norm.control.type.field = NULL,
    norm.control.types = NULL,
    norm.tmm = TRUE,
    gene.field = NULL,
    consolidation = c("simes", "fry", "holm-min"),
    fry.args = list(),
    holm.min.n = 3,
    holm.min.prop = 0.3,
    assay = 1,
    row.data = gene.field,
    gene.data = NULL,
    metadata = NULL,
    output.dir = "voom",
    author = NULL,
    dry.run = FALSE,
    save.results = TRUE
) {
    restore.fun <- resetInputCache()
    on.exit(restore.fun(), after=FALSE, add=TRUE)

    if (is.null(author)) {
        author <- Sys.info()[["user"]]
    }

    dir.create(output.dir, showWarnings=FALSE, recursive=FALSE)
    fname <- file.path(output.dir, "report.Rmd")

    ########################################
    ### Initialization and normalization ###
    ########################################

    common.start <- .initialize(
        x=x,
        method="voom",
        assay=assay,
        groups=groups,
        comparisons=comparisons,
        covariates=covariates,
        block=block,
        subset.factor=subset.factor,
        subset.levels=subset.levels,
        subset.groups=subset.groups,
        filter.reference.factor=filter.reference.factor,
        filter.reference.levels=filter.reference.levels,
        filter.num.mads=filter.num.mads,
        filter.default=filter.default,
        design=design,
        contrasts=contrasts,
        author=author
    )
    writeRmd(common.start$text, file=fname)
    contrast.info <- common.start$contrasts

    norm <- .normalize(
        norm.control.barcodes=norm.control.barcodes,
        norm.control.gene.field=gene.field,
        norm.control.genes=norm.control.genes,
        norm.control.type.field=norm.control.type.field,
        norm.control.types=norm.control.types,
        norm.tmm=norm.tmm
    )
    writeRmd(norm, file=fname, append=TRUE)

    ###################################
    ### Barcode-level model fitting ###
    ###################################

    template <- system.file("templates", "screen.Rmd", package="augere.screen", mustWork=TRUE)
    parsed <- parseRmdTemplate(readLines(template))
    replacements <- list(ASSAY = deparseToString(assay))

    if (common.start$mds.use.group) {
        parsed[["mds-ungrouped"]] <- NULL
    } else {
        parsed[["mds-grouped"]] <- NULL
    }

    if (!is.null(dc.block)) {
        parsed[["voom"]] <- NULL
        replacements$DUPCOR_BLOCK <- deparseToString(dc.block)
        replacements$LM_OPTS <- ", block=dc.block, correlation=dc$consensus.correlation"
    } else {
        parsed[["duplicate-correlation"]] <- NULL
        replacements$LM_OPTS <- ""
    }

    if (quality) {
        replacements$VOOM_CMD <- "voomWithQualityWeights"
    } else {
        replacements$VOOM_CMD <- "voom"
        parsed[["quality-text"]] <- NULL
    }

    eb.args <- ""
    if (trend) {
        eb.args <- c(eb.args, "trend=TRUE")
    }

    if (robust) {
        eb.args <- c(eb.args, "robust=TRUE")
        replacements$EXTRA_EB_CAPT <- " with outliers marked in red"
    } else {
        replacements$EXTRA_EB_CAPT <- ""
        parsed[["robust-text"]] <- NULL
    }

    replacements$EB_OPTS <- paste(eb.args, collapse=", ")

    merge.metadata <- !is.null(metadata)
    if (merge.metadata) {
        parsed[["create-common-metadata"]] <- replacePlaceholders(parsed[["create-common-metadata"]], list(COMMON_METADATA=deparseToString(metadata)))
        parsed[["no-merge-metadata"]] <- NULL
    } else {
        parsed[["create-common-metadata"]] <- NULL
        parsed[["merge-metadata"]] <- NULL
    }

    contrasts <- vector("list", length(contrast.info))
    save.names <- character()
    author.txt <- deparseToString(as.list(author))
    replacements$AUTHOR <- author.txt

    consolidation <- match.arg(consolidation, several.ok=TRUE)
    if (is.null(gene.field)) {
        consolidation <- character(0)
        parsed[["gene-output"]] <- NULL
    } else {
        replacements$GENE_FIELD <- deparseToString(gene.field)
    }
    if (!("fry" %in% consolidation)) {
        parsed[["gene-output"]][["create-fry-indices"]] <- NULL
    } else {
        fry.args <- c(eb.args, vapply(names(fry.args), function(n) sprintf("%s=%s", n, deparseToString(fry.args[[n]])), FUN.VALUE=""))
    }
    if (is.null(gene.data)) {
        parsed[["gene-output"]][["create-gene-annotation"]] <- NULL
    } else {
        replacements$ROW_DATA_FIELDS <- deparseToString(gene.data)
    }

    for (i in seq_along(contrast.info)) {
        copy <- parsed$contrast

        ###############################
        ### Barcode-level contrasts ###
        ###############################

        current <- contrast.info[[i]]
        copy[["contrast-command"]] <- current$commands

        if (lfc.threshold == 0) {
            copy[["lfc-contrast"]] <- NULL
        } else {
            copy[["no-lfc-contrast"]] <- NULL
            copy[["lfc-contrast"]] <- replacePlaceholders(
                copy[["lfc-contrast"]],
                list(
                    EB_OPTS=replacements$EB_OPTS,
                    LFC_THRESHOLD=deparseToString(lfc.threshold)
                )
            )
        }

        if (is.null(row.data)) {
            copy[["attach-rowdata"]] <- NULL
        } else {
            copy[["attach-rowdata"]] <- replacePlaceholders(copy[["attach-rowdata"]], list(ROWDATA_FIELDS=deparseToString(row.data)))
        }

        meta.cmds <- processContrastMetadata(current)
        meta.cmds[-1] <- paste0("    ", meta.cmds[-1])
        meta.cmds[1] <- paste0("    contrast=", meta.cmds[1])
        meta.cmds[length(meta.cmds)] <- paste0(meta.cmds[length(meta.cmds)], ",")

        decorate_metadata <- function(txt) {
            txt[["diff-metadata"]] <- meta.cmds
            if (is.null(subset.factor)) {
                txt[["subset-metadata"]] <- NULL
            }
            if (!merge.metadata) {
                txt[["merge-metadata"]] <- NULL
            } else {
                txt[["no-merge-metadata"]] <- NULL
            }
            txt
        }

        copy <- decorate_metadata(copy)

        ################################
        ### Gene-level consolidation ###
        ################################

        if ("simes" %in% consolidation) {
            simes.cmd <- .generateSimesCommands(    
                df.name="barcode.df",
                genes.name="gene.ids",
                has.LogFC=current$type %in% c("versus", "covariate")
            )
            save.name <- paste0("save-simes-", i)
            save.names <- c(save.names, save.name)
            copy$simes <- decorate_metadata(copy$simes)
            copy$simes <- replacePlaceholders(copy$simes, list(SIMES_CMD=simes.cmd, SAVING_CHUNK_NAME=save.name), error=FALSE)
            if (is.null(gene.data)) {
                copy$simes[["gene-annotation"]] <- NULL
            }
        } else {
            copy$simes <- NULL
        }

        if ("fry" %in% consolidation) {
            save.name <- paste0("save-fry-", i)
            save.names <- c(save.names, save.name)
            copy$fry <- decorate_metadata(copy$fry)
            copy$fry <- replacePlaceholders(
                copy$fry,
                list(
                    FRY_OPTS=paste(sprintf(",\n    %s", fry.args), collapse=""),
                    SAVING_CHUNK_NAME=save.name
                ),
                error=FALSE
            )
            if (is.null(gene.data)) {
                copy$fry[["gene-annotation"]] <- NULL
            }
        } else {
            copy[["fry"]] <- NULL
        }

        if ("holm-min" %in% consolidation) {
            holm.cmd <- .generateHolmMinCommands(
                df.name="barcode.df",
                genes.name="gene.ids",
                min.n=holm.min.n,
                min.prop=holm.min.prop,
                has.LogFC=current$type %in% c("versus", "covariate")
            )
            save.name <- paste0("save-holm-", i)
            save.names <- c(save.names, save.name)
            copy$`holm-min` <- decorate_metadata(copy$`holm-min`)
            copy$`holm-min` <- replacePlaceholders(copy$`holm-min`, list(HOLM_CMD=holm.cmd, SAVING_CHUNK_NAME=save.name), error=FALSE)
            if (is.null(gene.data)) {
                copy$`holm-min`[["gene-annotation"]] <- NULL
            }
        } else {
            copy[["holm-min"]] <- NULL
        }

        save.name <- paste0("save-de-", i)
        save.names[i] <- save.name
        contrasts[[i]] <- replacePlaceholders(
            copy,
            list(
                AUTHOR=author.txt,
                EB_OPTS=replacements$EB_OPTS,
                CONTRAST_NAME_SIMPLE=current$title,
                CONTRAST_NAME_DEPARSED=deparseToString(current$title),
                SAVING_CHUNK_NAME=save.name
            )
        )
    }

    parsed$contrast <- contrasts
    parsed <- replacePlaceholders(parsed, replacements)
    writeRmd(parsed, file=fname, append=TRUE)

    if (dry.run) {
        return(NULL)
    }

    if (save.results) {
        skip.chunks <- NULL
    } else {
        skip.chunks <- c("save-directory", save.names, "save-norm")
    }
    env <- new.env()
    compileReport(fname, env=env, skip.chunks=skip.chunks)

    list(
        barcode=env$all.barcodes,
        gene=list(
            fry=env$all.fry,
            simes=env$all.simes,
            `holm-min`=env$all.holm.min
        ),
        normalized=env$normalized
    )
}
