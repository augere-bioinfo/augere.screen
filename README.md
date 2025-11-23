# Differential abundance in functional screens

Implements a pipeline function to generate parametrized Rmarkdown reports for differential abundance analyses of functional screen data.
The differential abundance analysis is done using `voom()` on the barcode-level counts, followed by consolidation to per-gene statistics for easier interpretation.
The report contains all of the R commands required to reproduce the analysis.
