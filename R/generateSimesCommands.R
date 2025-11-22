#' @importFrom stats p.adjust
#' @importFrom S4Vectors DataFrame
#' @importFrom metapod groupedSimes
#' @import augere.core
.generateSimesCommands <- function(df.name, genes.name, has.LogFC) {
    template <- 'local({
    df <- <%= DF_NAME %>
    genes <- factor(<%= GENES_NAME %>)
    keep <- !is.na(df$PValue)
    df <- df[keep,,drop=FALSE]
    genes <- genes[keep]

    num.barcodes <- table(genes) 
    output <- S4Vectors::DataFrame(
        row.names = names(num.barcodes),
        NumBarcodes = as.integer(num.barcodes)
    )
    <%= BODY %>
    output
})'

    signed.cmds <- '
    is.up <- df$LogFC > 0
    half.p <- df$PValue / 2
    other.p <- 1 - half.p

    # Computing and combining one-sided p-values for each direction,
    # as we want genes where barcodes change in the same direction.
    p.up <- ifelse(is.up, half.p, other.p)
    simes.up <- metapod::groupedSimes(p.up, genes)
    simes.up$p.value <- simes.up$p.value[rownames(output)]

    p.down <- ifelse(!is.up, half.p, other.p)
    simes.down <- metapod::groupedSimes(p.down, genes)
    simes.down$p.value <- simes.down$p.value[rownames(output)]

    output$PValue <- pmin(simes.up$p.value, simes.down$p.value) * 2
    output$FDR <- p.adjust(output$PValue, method="BH")
    output$Direction <- ifelse(simes.up$p.value < simes.down$p.value, "up", "down")
'

    if (has.LogFC) {
        return(replacePlaceholders(template, list(DF_NAME=df.name, GENES_NAME=genes.name, BODY=signed.cmds)))
    }

    conditional.cmds <- '
    if ("LogFC" %in% colnames(df)) {
    <%= PAYLOAD %>
    } else {
        # If the contrast has no directionality, we just combine p-values as they are.
        simes <- metapod::groupedSimes(df$PValue, genes)
        output$PValue <- simes$p.value[rownames(output)]
        output$FDR <- p.adjust(output$PValue, method="BH")
    }
'

    signed.lines <- strsplit(trimws(signed.cmds), '\n')[[1]]
    signed.lines <- sprintf("    %s", signed.lines)
    full.body <- replacePlaceholders(conditional.cmds, list(PAYLOAD=paste(signed.lines, collapse="\n")))
    replacePlaceholders(template, list(DF_NAME=df.name, GENES_NAME=genes.name, BODY=full.body))
}
