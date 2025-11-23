#' @import augere.core
#' @importFrom S4Vectors DataFrame
.generateSummaryCommands <- function(df.name, by.gene.name, summary.threshold, has.LogFC=TRUE) {
    template <- 'local({
    df <- <%= DF_NAME %>
    by.gene <- <%= BY_GENE_NAME %>
    is.sig <- df$FDR <= <%= SUMMARY_THRESHOLD %>
    <%= BODY %>
})'

    signed.cmds <- '
    num.up <- num.down <- integer(length(by.gene))
    lfc <- rep(NA_real_, length(by.gene))
    for (i in seq_along(by.gene)) {
        current <- by.gene[[i]]
        current.sig <- is.sig[current]
        current.lfc <- df$LogFC[current]
        num.up[i] <- sum(current.sig & current.lfc > 0, na.rm=TRUE)
        num.down[i] <- sum(current.sig & current.lfc < 0, na.rm=TRUE)
        best <- which.min(df$PValue[current])
        if (length(best)) {
            lfc[i] <- current.lfc[best]
        }
    }
    S4Vectors::DataFrame(row.names=names(by.gene), NumUp=num.up, NumDown=num.down, LogFC=lfc)
'

    if (has.LogFC) {
        return(replacePlaceholders(template, list(DF_NAME=df.name, BY_GENE_NAME=by.gene.name, SUMMARY_THRESHOLD=deparseToString(summary.threshold), BODY=signed.cmds)))
    }

    conditional.cmds <- '
    if ("LogFC" %in% colnames(df)) {
    <%= PAYLOAD %>
    } else {
        num.sig <- integer(length(by.gene))
        for (i in seq_along(by.gene)) {
            current <- by.gene[[i]]
            current.sig <- is.sig[current]
            num.sig[i] <- sum(current.sig, na.rm=TRUE)
        }
        S4Vectors::DataFrame(row.names=names(by.gene), NumSig=num.sig)
    }
'

    signed.lines <- strsplit(trimws(signed.cmds), '\n')[[1]]
    signed.lines <- sprintf("    %s", signed.lines)
    full.body <- replacePlaceholders(conditional.cmds, list(PAYLOAD=paste(signed.lines, collapse="\n")))
    replacePlaceholders(template, list(DF_NAME=df.name, BY_GENE_NAME=by.gene.name, SUMMARY_THRESHOLD=deparseToString(summary.threshold), BODY=full.body))
}
