#!/usr/bin/env Rscript

# Author: Yu-Chen Emma Huang
# Updated: 2025-05-20
# Description: Generate covariate and BED files for either TensorQTL or QTLtools
#setwd("/work/users/y/u/yuchenh/sQTL_colon/02.geneexp")
suppressMessages({
  library(optparse)
  library(edgeR)
  library(RUVSeq)
  library(dplyr)
  library(corrplot)
  library(ggplot2)
})

#Ex: 
#Rscript ../scripts/02.geneexp/generate_bed_cov_rmBATCHGENDER.R \
#  --project PLEXUS_colon_n1149 \
#  --metafile ../00.input/metafiles/PLEXUS_colon_n1149_meta.txt \
#  --gene_anno ../00.input/annotation/gencode.v49.genedf.tsv \
#  --sample_size 1149 \
#  --N 1000 \
#  --RUV_k 30 \
#  --data PLEXUS \
#  --mode qtltools




# CLI options
option_list <- list(
  make_option("--project", type = "character", help = "Project name, e.g. UNC_colon_n255"),
  make_option("--metafile", type = "character", help = "Path to meta file"),
  make_option("--gene_anno", type = "character", help = "Path to gene annotation TSV"),
  make_option("--sample_size", type = "integer", help = "Sample size (number of individuals)"),
  make_option("--N", type = "integer", default = 1000, help = "Number of top genes to select"),
  make_option("--RUV_k", type = "integer", default = 30, help = "Number of RUV factors to extract"),
  make_option("--data", type = "character", default = "UNC", help = "'PLEXUS', 'UNC', or 'PLEXUS_UNC', or 'UNC_PLEXUS'"),
  make_option("--mode", type = "character", default = "both", help = "Choose 'tensorqtl', 'qtltools', or 'both'")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Assign vars
project <- opt$project
metafile <- opt$metafile
bedsuffix <- opt$bedsuffix
gene_anno <- opt$gene_anno
sample_size <- opt$sample_size
N <- opt$N
RUV_k <- opt$RUV_k
data <- opt$data
mode <- opt$mode

# Paths
gene_quant <- file.path(project, "quant", paste0(project, "_gene_quant.txt"))
bed_folder <- file.path(project, "bed")
figure_folder<- file.path(project,"figure")

dir.create(bed_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_folder, showWarnings = FALSE, recursive = TRUE)


# Load and filter
counts_genes <- read.delim(gene_quant, header = TRUE)
row.names(counts_genes) <- counts_genes$gene
counts_genes$gene <- NULL
norm_genes <- DGEList(counts = counts_genes)
keep_genes_counts <- rowSums(norm_genes$counts >= 6) >= round(0.20 * sample_size)
genes_counts_filtered <- norm_genes[keep_genes_counts,]
genes_counts_filtered_TMMadj <- calcNormFactors(genes_counts_filtered)
filtered_cpmadjtmm <- cpm(genes_counts_filtered_TMMadj)
colnames(filtered_cpmadjtmm) <- gsub("X", "", colnames(filtered_cpmadjtmm))

# Inverse normal transform
gene_norm <- t(apply(filtered_cpmadjtmm, 1, function(y) {
  qnorm((rank(y, ties.method = "average") - 0.5) / sum(!is.na(y)))
}))

meta <- read.delim(metafile, sep = "\t", header = TRUE)
sample_order <- colnames(filtered_cpmadjtmm)
meta <- meta[match(sample_order, meta$SAMPLE), ]
print(colnames(filtered_cpmadjtmm))
print(meta$SAMPLE)

stopifnot(all(meta$SAMPLE == sample_order)) # Sanity check: make sure it matches


design <- model.matrix(~ 1, data = meta)  
genoPC <- as.matrix(meta[, paste0("genotypePC", 1:4)])
gene_corrected <- removeBatchEffect(gene_norm, batch=meta$BATCH, batch2=meta$SEX, covariates=genoPC, design=design)
rownames(gene_corrected) <- rownames(filtered_cpmadjtmm)
colnames(gene_corrected) <- colnames(filtered_cpmadjtmm)


# Top genes
variance_rank <- rank(apply(gene_corrected, 1, var))
top_genes <- order(variance_rank)[1:N]
top_genes_data <- gene_corrected[top_genes, ]

print("run_RUV")
# RUV
ruv_set <- RUVg(as.matrix(gene_corrected), cIdx = rownames(top_genes_data), k = RUV_k, isLog = TRUE)

fuv <- ruv_set[["W"]][, 1:RUV_k]
colnames(fuv) <- paste0("W", 1:RUV_k)
cov <- cbind(meta[, c("SAMPLE","BATCH","SEX","genotypePC1", "genotypePC2", "genotypePC3", "genotypePC4")], fuv)

if (mode %in% c("qtltools", "both")) {
  cov_qtltools <- file.path(project, "cov", "qtltools")
  dir.create(cov_qtltools, showWarnings = FALSE, recursive = TRUE)
  for (i in 0:RUV_k) {
    rows2save <- t(cov)[1:(i + 7), ]
    file_name <- file.path(cov_qtltools, paste0("cov_", i, ".txt"))
    write.table(rows2save, file = file_name, sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
  }
}

if (mode %in% c("tensorqtl", "both")) {

  cov_tensorqtl <- file.path(project, "cov", "tensorqtl")
  dir.create(cov_tensorqtl, showWarnings = FALSE, recursive = TRUE)
  for (i in 0:RUV_k) {
    rows2save <- t(cov)[1:(i + 7), ]
    file_name <- file.path(cov_tensorqtl, paste0("cov_", i, ".txt"))
    write.table(rows2save, file = file_name, sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
  }
}

 
 

# BED
phenotypes.out <- gene_norm
phenotypes.out <- as.data.frame(phenotypes.out)
phenotypes.out$gene_id <- rownames(phenotypes.out)
gene_inf <- read.delim(gene_anno, sep = "\t", header = TRUE)
bed <- right_join(gene_inf, phenotypes.out, by = "gene_id")

cols <- colnames(bed)
bed <- bed[, colnames(bed) != "length"]
meta_cols <- c("chr", "start", "end", "gene_id", "gene_name", "strand")
bed <- bed[, c(meta_cols, setdiff(colnames(bed), meta_cols))]
sorted_bed <- bed[order(bed$chr, bed$start), ]
colnames(sorted_bed)[colnames(sorted_bed) == "chr"] <- "#chr"
bed_4col <- sorted_bed[, !colnames(sorted_bed) %in% c("gene_name", "strand")]


if (mode %in% c("qtltools", "both")) {
  qtltoolsbed <- file.path(bed_folder, paste0(project, "_eqtl_qtltools.bed"))
  write.table(sorted_bed, file = qtltoolsbed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

if (mode %in% c("tensorqtl", "both")) {
  tensorqtlbed <- file.path(bed_folder, paste0(project, "_eqtl_tensorqtl.bed"))
  write.table(bed_4col, file = tensorqtlbed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
 

if (data %in% c("PLEXUS","PLEXUS_UNC","UNC_PLEXUS")) { 
  covariates_of_interest <- cbind(
    BATCH = as.numeric(as.factor(meta$BATCH)),
    SEX = as.numeric(as.factor(meta$SEX)),
    TISSUE = as.numeric(as.factor(meta$TISSUE)),
    TIN = meta$TIN,
    genotypePC1 = meta$genotypePC1,
    genotypePC2 = meta$genotypePC2,
    genotypePC3 = meta$genotypePC3,
    genotypePC4 = meta$genotypePC4
  )
} else {
  covariates_of_interest <- cbind(
    BATCH = as.numeric(as.factor(meta$BATCH)),
    SEX = as.numeric(as.factor(meta$SEX)),
    TIN = meta$TIN,
    genotypePC1 = meta$genotypePC1,
    genotypePC2 = meta$genotypePC2,
    genotypePC3 = meta$genotypePC3,
    genotypePC4 = meta$genotypePC4
  )
}




RUVMatrix <- cor(covariates_of_interest, fuv, use = "complete.obs")
rownames(RUVMatrix) <- colnames(covariates_of_interest)

if (mode %in% c("qtltools", "both")) {
  dir.create(cov_qtltools, showWarnings = FALSE, recursive = TRUE)
  write.table(RUVMatrix, file = file.path(cov_qtltools, "RUVMatrix.txt"), sep = "\t", quote = FALSE, col.names = NA)
  png(file = file.path(figure_folder, "covariates_and_genoPC_vs_RUV_corr_qtltools.png"), width = 80 * RUV_k, height = 600, res = 100)
  corrplot(RUVMatrix)
  dev.off()
}

if (mode %in% c("tensorqtl", "both")) {

  write.table(RUVMatrix, file = file.path(cov_tensorqtl, "RUVMatrix.txt"), sep = "\t", quote = FALSE, col.names = NA)
  png(file = file.path(figure_folder, "covariates_and_genoPC_vs_RUV_corr_tensorqtl.png"), width = 80 * RUV_k, height = 600, res = 100)
  corrplot(RUVMatrix)
  dev.off()
}


## Original Plot ## 
counts2plot <- as.matrix(gene_norm)
pca_result <- prcomp(t(counts2plot), scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x)
pca_scores_with_metadata <- cbind(pca_scores, meta)
cov_desc <- "Original"
original_p1 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = BATCH)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by BATCH\n(Residualized:", cov_desc, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
print(original_p1)
ggsave(filename = file.path(figure_folder,paste0("PCA_plot_original_BATCH.png")), plot = original_p1, width = 8, height = 6, dpi = 300)

original_p2 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = SEX)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by SEX\n(Residualized:", cov_desc, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
print(original_p2)
ggsave(filename = file.path(figure_folder,paste0("PCA_plot_original_SEX.png")), plot = original_p2, width = 8, height = 6, dpi = 300)

original_p3 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = DISEASE)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by DISEASE\n(Residualized:", cov_desc, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
print(original_p3)
ggsave(filename = file.path(figure_folder,paste0("PCA_plot_original_DISEASE.png")), plot = original_p3, width = 8, height = 6, dpi = 300)


if (data %in% c("PLEXUS", "UNC_PLEXUS","PLEXUS_UNC")) {
  original_p4 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = TISSUE)) +
      geom_point(size = 3) +
      geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
      labs(title = paste("PCA Plot Colored by TISSUE\n(Residualized:", cov_desc, ")"),
          x = "PC1", y = "PC2") +
      theme_classic()
  print(original_p4)
  ggsave(filename = file.path(figure_folder,paste0("PCA_plot_original_TISSUE.png")), plot = original_p4, width = 8, height = 6, dpi = 300)
}



corrected_matrix <- removeBatchEffect(gene_norm, batch=meta$BATCH, design=design)
counts2plot <- as.matrix(corrected_matrix)
pca_result <- prcomp(t(counts2plot), scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x)
pca_scores_with_metadata <- cbind(pca_scores, meta)
cov_desc <- "BATCH"


rmBATCH_p11 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = BATCH)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by BATCH\n(Residualized:", cov_desc, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
print(rmBATCH_p11)
ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCH_BATCH.png")), plot = rmBATCH_p11, width = 8, height = 6, dpi = 300)

rmBATCH_p12 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = DISEASE)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by DISEASE\n(Residualized:", cov_desc, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
print(rmBATCH_p12)
ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCH_DISEASE.png")), plot = rmBATCH_p12, width = 8, height = 6, dpi = 300)

rmBATCH_p13 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = SEX)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by GENDER\n(Residualized:", cov_desc, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
print(rmBATCH_p13)
ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCH_GENDER.png")), plot = rmBATCH_p13, width = 8, height = 6, dpi = 300)
  if (data %in% c("PLEXUS", "UNC_PLEXUS","PLEXUS_UNC")) {
rmBATCH_p14 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = TISSUE)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by TISSUE\n(Residualized:", cov_desc, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
print(rmBATCH_p14)
ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCH_TISSUE.png")), plot = rmBATCH_p14, width = 8, height = 6, dpi = 300)
  }

corrected_matrix2 <- removeBatchEffect(gene_norm, batch=meta$BATCH, batch2=meta$SEX, design=design)
counts2plot2 <- as.matrix(corrected_matrix2)
pca_result2 <- prcomp(t(counts2plot2), scale. = TRUE)
pca_scores2 <- as.data.frame(pca_result2$x)
pca_scores_with_metadata2 <- cbind(pca_scores2, meta)
cov_desc2 <- c("BATCH","SEX")

rmBATCH_p21 <- ggplot(pca_scores_with_metadata2, aes(x = PC1, y = PC2, color = BATCH)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by BATCH\n(Residualized:", cov_desc2, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
print(rmBATCH_p21)
ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCHSEX_BATCH.png")), plot = rmBATCH_p21, width = 8, height = 6, dpi = 300)


rmBATCH_p22 <- ggplot(pca_scores_with_metadata2, aes(x = PC1, y = PC2, color = SEX)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by SEX\n(Residualized:", cov_desc2, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
print(rmBATCH_p22)
ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCHSEX_GENDER.png")), plot = rmBATCH_p22, width = 8, height = 6, dpi = 300)


rmBATCH_p23 <- ggplot(pca_scores_with_metadata2, aes(x = PC1, y = PC2, color = DISEASE)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by DISEASE\n(Residualized:", cov_desc2, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
print(rmBATCH_p23)
ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCHSEX_DISEASE.png")), plot = rmBATCH_p23, width = 8, height = 6, dpi = 300)
  if (data %in% c("PLEXUS", "UNC_PLEXUS","PLEXUS_UNC")) {
    rmBATCH_p24 <- ggplot(pca_scores_with_metadata2, aes(x = PC1, y = PC2, color = TISSUE)) +
        geom_point(size = 3) +
        geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
        labs(title = paste("PCA Plot Colored by TISSUE\n(Residualized:", cov_desc2, ")"),
            x = "PC1", y = "PC2") +
        theme_classic()
    print(rmBATCH_p24)
    ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCHSEX_TISSUE.png")), plot = rmBATCH_p24, width = 8, height = 6, dpi = 300)
  }



for (k in 0:opt$RUV_k) {
  design <- model.matrix(~1, data=meta)   


  if (k == 0) {
    covariate_cols <- c("genotypePC1","genotypePC2","genotypePC3","genotypePC4")
    cov_desc <- "BATCH+GENDER+GenotypePC"
  } else {
    W_cols <- paste0("W", 1:k)
    covariate_cols <- c("genotypePC1","genotypePC2","genotypePC3","genotypePC4", W_cols)
    cov_desc <- paste("BATCH+GENDER+GenotypePC1-PC4 +", paste0(W_cols[1], "-", tail(W_cols, 1)))
  }

  X2 <- as.matrix(cov[, covariate_cols])
  gene_corrected <- removeBatchEffect(gene_norm, batch=meta$BATCH, batch2=meta$SEX, covariates=X2, design=design)
  counts2plot <- as.matrix(gene_corrected)
  pca_result <- prcomp(t(counts2plot), scale. = TRUE)
  pca_scores <- as.data.frame(pca_result$x)
  pca_scores_with_metadata <- cbind(pca_scores, meta)

  
  p1 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = BATCH)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by BATCH\n(Residualized:", cov_desc, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()

  p2 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = DISEASE)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by DISEASE\n(Residualized:", cov_desc, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
  
  p3 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = SEX)) +
    geom_point(size = 3) +
    geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
    labs(title = paste("PCA Plot Colored by GENDER\n(Residualized:", cov_desc, ")"),
         x = "PC1", y = "PC2") +
    theme_classic()
  
  if (data %in% c("PLEXUS", "UNC_PLEXUS","PLEXUS_UNC")) {
    p4 <- ggplot(pca_scores_with_metadata, aes(x = PC1, y = PC2, color = TISSUE)) +
      geom_point(size = 3) +
      geom_text(aes(label = SAMPLE), vjust = -0.5, size = 3) +
      labs(title = paste("PCA Plot Colored by TISSUE\n(Residualized:", cov_desc, ")"),
          x = "PC1", y = "PC2") +
      theme_classic()
      ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCHSEXgenoPC",k,"_TISSUE.png")), plot = p4, width = 8, height = 6, dpi = 300)
  }
  ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCHSEXgenoPC",k,"_BATCH.png")), plot = p1, width = 8, height = 6, dpi = 300)
  ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCHSEXgenoPC",k,"_DISEASE.png")), plot = p2, width = 8, height = 6, dpi = 300)
  ggsave(filename = file.path(figure_folder,paste0("PCA_plot_rmBATCHSEXgenoPC",k,"_GENDER.png")), plot = p3, width = 8, height = 6, dpi = 300)

}