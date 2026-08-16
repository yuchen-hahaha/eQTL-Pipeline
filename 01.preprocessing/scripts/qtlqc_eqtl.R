#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(knitr)
  library(rjson)
  library(ape)
  library(reshape2)
  library(tidyverse)
  library(e1071)
  library(BiocManager)
  library(sva)
  library(data.table)
  library(optparse)
  library(edgeR)
})

# ---------------------------
# Parse command line options
# ---------------------------
option_list = list(
  make_option(c("--filename"), type = "character", help = "Input gene quant file (quant.txt)"),
  make_option(c("--metafile"), type = "character", help = "Metadata file"),
  make_option(c("--output_dir"), type = "character", help = "Output directory"),
  make_option(c("--json_config"), type = "character", help = "JSON config file"),
  make_option(c("--n_sample"), type = "integer", help = "Total sample count")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$filename) | is.null(opt$metafile) | is.null(opt$output_dir) | is.null(opt$json_config) | is.null(opt$n_sample)) {
  print_help(opt_parser)
  stop("All options must be supplied.", call. = FALSE)
}

# ---------------------------
# Parameters
# ---------------------------
ntop <- 10000
pseudo_count <- 0.001

filename <- opt$filename
metafile <- opt$metafile
output_dir <- opt$output_dir
n_sample <- opt$n_sample
json_file <- opt$json_config

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Load metadata & config
# ---------------------------
meta <- read.table(metafile, header = TRUE, sep = "\t")
batch <- meta[c("SAMPLE", "BATCH")]
Inflam <- meta[c("SAMPLE", "INFLAMMATION_STATUS")]
json_data <- fromJSON(file = json_file)

# Config from JSON
RLEFilterPercent <- json_data$config$transcriptome$transcriptomeQC$RLEFilterPercent
DSFilterPercent <- json_data$config$transcriptome$transcriptomeQC$DSFilterPercent
pvalues.cut <- json_data$config$transcriptome$transcriptomeQC$hcluster$pvalues.cut
topk_genes <- json_data$config$transcriptome$transcriptomeQC$hcluster$topk_genes
cluster_percent <- json_data$config$transcriptome$transcriptomeQC$hcluster$cluster_percent
treesNum <- json_data$config$transcriptome$transcriptomeQC$hcluster$cluster_level
low_expr_TPM <- json_data$config$transcriptome$transcriptomeQC$low_expr_TPM
low_expr_TPM_percent <- json_data$config$transcriptome$transcriptomeQC$low_expr_TPM_percent
pseudo_count <- json_data$config$transcriptome$transcriptomeQC$TPM_pseudo_count
Female_specific_gene <- json_data$config$transcriptome$transcriptomeQC$gender_check$Female_specific_gene
Male_specific_gene <- json_data$config$transcriptome$transcriptomeQC$gender_check$Male_specific_gene

cat("Configuration file loaded successfully.\n")

# ---------------------------
# Load and filter counts
# ---------------------------
counts_genes <- read.delim(filename, header = TRUE)
row.names(counts_genes) <- counts_genes$gene
norm_genes <- DGEList(counts = counts_genes)
keep_genes_counts <- rowSums(norm_genes$counts >= 6) >= round(0.20 * n_sample)
genes_counts_filtered <- norm_genes[keep_genes_counts, ]
genes_counts_filtered_TMMadj <- calcNormFactors(genes_counts_filtered)
norm.factors <- genes_counts_filtered_TMMadj$samples$norm.factors
eff.lib.size <- genes_counts_filtered_TMMadj$samples$lib.size * norm.factors
filtered_cpmadjtmm <- cpm(genes_counts_filtered_TMMadj)
colnames(filtered_cpmadjtmm) <- gsub("X", "", colnames(filtered_cpmadjtmm))

TPMdata <- filtered_cpmadjtmm
logtpm <- log10(TPMdata + pseudo_count)
logtpm <- as.data.frame(logtpm)

# ---------------------------
# RLE QC plot
# ---------------------------
RLEFilterLength <- RLEFilterPercent * ncol(TPMdata)
DSFilter <- DSFilterPercent * ncol(TPMdata)

rle <- logtpm - apply(logtpm, 1, median)
iqr <- apply(rle, 2, IQR)
rle <- reshape2::melt(
  cbind(ID = rownames(rle), rle),
  variable.name = "sample",
  value.name = "TPM",
  id.vars = "ID"
)
names(rle) <- c("feature", "sample", "TPM")
rle_IQR <- rle %>% group_by(sample) %>% summarise(IQR = IQR(TPM))
rle_IQR_range <- rle_IQR$IQR %>% range %>% abs() %>% max()
rle_IQR_range <- 2 * rle_IQR_range %>% ceiling()
bymedian <- with(rle, reorder(sample, TPM, IQR))

png(file.path(output_dir, "RLE_plot_before_rmRUV.png"), width = 1800, height = 600, res = 300)
par(mar = c(3, 3, 3, 3))
boxplot(TPM ~ bymedian, data = rle, outline = FALSE, ylim = c(-rle_IQR_range, rle_IQR_range),
  las = 2, boxwex = 1, col = 'gray', cex.axis = 0.3, main = "RLE plot",
  xlab = "", ylab = "Residual expression levels", frame = FALSE)
dev.off()

ExpPersample <- nrow(TPMdata)
sorted_samples <- rle_IQR[order(rle_IQR$IQR), ]
RLEFilterList <- sorted_samples$sample[1:RLEFilterLength]
RLEFilterList <- as.character(RLEFilterList)
cat("Candidate outliers by RLE:", RLEFilterList, "\n")

# ---------------------------
# H-cluster
# ---------------------------
sample_order <- colnames(logtpm)
batch_reordered <- batch[match(sample_order, batch$SAMPLE), ]
sampleDists <- 1 - cor(logtpm, method = 'pearson')
hc <- hclust(as.dist(sampleDists), method = "average")
hcphy <- as.phylo(hc)

Pvars <- apply(logtpm, 1, var)
select <- order(Pvars, decreasing = TRUE)[seq_len(min(topk_genes, length(Pvars)))]
MD_matrix <- logtpm[select, ]
MahalanobisDistance <- mahalanobis(t(MD_matrix), colMeans(t(MD_matrix)), cov(t(MD_matrix)))
pvalues <- pchisq(MahalanobisDistance, df = nrow(MD_matrix), lower.tail = FALSE)
pvalues.adjust <- p.adjust(pvalues, method = "bonferroni")
pvalues.low <- pvalues.adjust[pvalues < pvalues.cut]

HCoutliers <- character()
for (x in 1:treesNum) {
  trees <- cutree(hc, k = x)
  for (i in 1:x) {
    group <- hc$labels[which(trees == i)]
    if (sum(group %in% names(pvalues.low)) / length(group) >= cluster_percent) {
      HCoutliers <- union(HCoutliers, group)
    }
  }
}

cat("Cluster outliers:", HCoutliers, "\n")

# ---------------------------
# Save final outliers
# ---------------------------
D <- apply(1 - sampleDists, 1, median)
DSFilter <- sort(D)[DSFilter]
D <- data.frame(sample = names(D), D = D)
D_filterList <- D %>% filter(D <= DSFilter) %>% pull(sample)
cat("D-stat outliers:", D_filterList, "\n")

outliersList <- intersect(RLEFilterList, intersect(HCoutliers, D_filterList))
writeLines(outliersList, file.path(output_dir, "Final_Outliers.txt"))
cat("Total final outliers:", length(outliersList), "\n")
