#' @import augere.core
.normalize <- function(
    norm.control.barcodes, 
    norm.control.gene.field,
    norm.control.genes,
    norm.control.type.field,
    norm.control.types,
    norm.tmm
) {
    template <- system.file("templates", "normalize.Rmd", package="augere.screen", mustWork=TRUE)
    parsed <- parseRmdTemplate(readLines(template))

    if (
        !is.null(norm.control.barcodes) || 
        !is.null(norm.control.gene.field) ||
        !is.null(norm.control.type.field)
    ) {
        parsed[["norm-all"]] <- NULL
        parsed[["norm-none"]] <- NULL
        current <- parsed[["norm-control"]]

        if (is.null(norm.control.barcodes)) {
            current[["barcodes"]] <- NULL
        } else {
            current[["barcodes"]] <- replacePlaceholders(current[["barcodes"]], list(NORM_BARCODES=deparseToString(norm.control.barcodes)))
        }

        if (is.null(norm.control.gene.field) || is.null(norm.control.genes)) {
            current[["genes"]] <- NULL
        } else {
            current[["genes"]] <- replacePlaceholders(
                current[["genes"]],
                list(
                    NORM_GENE_FIELD=deparseToString(norm.control.gene.field),
                    NORM_GENES=deparseToString(norm.control.genes)
                )
            )
        }

        if (is.null(norm.control.type.field) || is.null(norm.control.types)) {
            current[["types"]] <- NULL
        } else {
            current[["types"]] <- replacePlaceholders(
                current[["types"]],
                list(
                    NORM_TYPE_FIELD=deparseToString(norm.control.type.field),
                    NORM_TYPES=deparseToString(norm.control.types)
                )
            )
        }

        if (norm.tmm) {
            current[["lib"]] <- NULL
        } else {
            current[["tmm"]] <- NULL
        }

        parsed[["norm-control"]] <- current

    } else {
        parsed[["norm-control"]] <- NULL
        parsed[["norm-none"]] <- NULL
        current <- parsed[["norm-all"]]

        if (norm.tmm) {
            current[["lib"]] <- NULL
        } else {
            current[["tmm"]] <- NULL
        }

        parsed[["norm-all"]] <- current
    }

    unlist(parsed, use.names=FALSE)
}
