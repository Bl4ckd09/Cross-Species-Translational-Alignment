#!/usr/bin/env Rscript
# rma_tggates.R — RMA-normalise the downloaded TG-GATEs liver CELs (Affymetrix Rat230-2).
# Same 31099-probeset space as DrugMatrix GSE57815 -> directly ComBat-alignable.
# Usage:  Rscript scripts/rma_tggates.R
suppressMessages({
  if (!requireNamespace("affy", quietly=TRUE) || !requireNamespace("rat2302cdf", quietly=TRUE)) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager", repos="https://cloud.r-project.org")
    BiocManager::install(c("affy","rat2302cdf"), update=FALSE, ask=FALSE)
  }
  library(affy)
})

cel_dir <- "data/_raw/tggates_cels"
out     <- "data/expression/tggates_liver_rma.tsv"
cels <- list.files(cel_dir, pattern="\\.CEL$", full.names=TRUE, ignore.case=TRUE)
cat("RMA on", length(cels), "CEL files ...\n")

# justRMA: memory-efficient (no full AffyBatch); auto-loads rat2302cdf for probeset summary
eset <- justRMA(filenames=cels, verbose=TRUE)
m <- exprs(eset)                                  # probes x samples, log2
colnames(m) <- sub("\\.CEL$", "", basename(colnames(m)), ignore.case=TRUE)

write.table(m, out, sep="\t", quote=FALSE, col.names=NA)
cat("wrote", out, ":", nrow(m), "probes x", ncol(m), "samples\n")
