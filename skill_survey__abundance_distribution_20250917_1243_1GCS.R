
# Input featureCounts file
infile <- "C:/Users/user/Desktop/20250206_skill_survey_sequencing/data/20200617_0937_hisat2/featurecounts/decorin__feature_counts_20200617_1vb.txt"

# Output folder
monogram <- "GCS"
ts <- format(Sys.time(), "%Y%m%d_%H%M")   # automatic timestamp
outdir <- file.path("C:/Users/user/Desktop/20250206_skill_survey_sequencing/data", ts)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# Read featureCounts table
data <- read.table(infile, header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)

# Find columns with read counts (.bam)
count_cols <- grep("\\.bam$", names(data))
counts <- as.matrix(data[, count_cols])
rownames(counts) <- data$Geneid


## RPKM calculation
libsize <- colSums(counts)  # library size = total mapped reads per sample
RPKM <- (counts * 1e9) / (data$Length * rep(libsize, each = nrow(counts)))


## TPM calculation
gene_len_kb <- data$Length / 1000
RPK <- counts / gene_len_kb
sum_RPK <- colSums(RPK)
TPM <- t( t(RPK) / sum_RPK ) * 1e6


## Plots
pdf(file.path(outdir, paste0("skill_survey__abundance_distribution_", ts, "_1", monogram, ".pdf")),
    width = 10, height = 7)

par(mfrow = c(3,1))

# function to plot density
plot_density <- function(mat, title) {
  X <- log2(mat + 1)
  cols <- rainbow(ncol(X))
  plot(density(X[,1], na.rm = TRUE),
       main = title, xlab = "log2(value+1)", ylab = "Density",
       col = cols[1], lwd = 2)
  if (ncol(X) > 1) {
    for (i in 2:ncol(X)) {
      lines(density(X[,i], na.rm = TRUE), col = cols[i], lwd = 2)
    }
  }
  legend("topright", legend = colnames(X), col = cols, lwd = 2, cex = 0.8, title = "Samples")
}

plot_density(counts, "Raw read count distribution (log2(count+1))")
plot_density(RPKM,   "RPKM-normalized distribution (log2(RPKM+1))")
plot_density(TPM,    "TPM-normalized distribution (log2(TPM+1))")

dev.off()

## Save outputs
write.csv(RPKM, file.path(outdir, paste0("skill_survey__rpkm_", ts, "_1", monogram, ".csv")))
write.csv(TPM,  file.path(outdir, paste0("skill_survey__tpm_", ts, "_1", monogram, ".csv")))

