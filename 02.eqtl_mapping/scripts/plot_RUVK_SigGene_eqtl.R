args <- commandArgs(trailingOnly = TRUE)

# Check arguments
if (length(args) != 5) {
  stop("Usage: Rscript plot_RUVK_SigGene_sqtl.R <input_prefix><output_dir> <output_plot_filename> <max_k> <fdr_threshold>")
}
# Ex: Rscript scripts/plot_RUVK_SigGene_eqtl.R UNC_colon_n255/eqtl/permutations10000/permutations_cov UNC_colon_n255/eqtl/permutations10000 UNC_colon_n255_eqtl_RUVg_sigGene_permutations10000 20

#Rscript scripts/plot_RUVK_SigGene_eqtl.R             UNC_colon_n241/eqtl/RUVr/permutationtest/perm100/permute_cov             UNC_colon_n241/eqtl/RUVr/permutationtest/perm100             UNC_colon_n241_eqtl_sigGene_perm100             20
#Rscript ../scripts/04.qtltools/plot_RUVK_SigGene_eqtl.R             PLEXUS_colon_n1149/eqtl/RUVr/permutationtest/perm100/permute_cov             PLEXUS_colon_n1149/eqtl/RUVr/permutationtest/perm100             PLEXUS_colon_n1149_eqtl_sigGene_perm100             30


input_prefix <- args[1]
output_dir <- args[2]
output_prefix <- args[3]
max_k <- as.integer(args[4])
fdr_threshold <- args[5]  # New argument for FDR (e.g., "0.05" or "0.01")

library(ggplot2)

k_values <- 0:max_k
gene_counts <- integer(length(k_values))

for (i in seq_along(k_values)) {
  k <- k_values[i]
  file_path <- sprintf("%s%d.fdr%s.significant.txt", input_prefix, k, fdr_threshold)
   
  if (!file.exists(file_path)) {
    paste("File does not exist:", file_path)
    next
  }

  df <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) > 0) {
    gene_counts[i] <- nrow(df)
  } else {
    gene_counts[i] <- 0
  }
}

plot_df <- data.frame(k = k_values, unique_genes = gene_counts)

best_k <- plot_df$k[which.max(plot_df$unique_genes)]
 
# Open PNG device
output_png <- file.path(output_dir, paste0(output_prefix, ".png"))
png(filename = output_png, width = 2400, height = 2400, res = 300)


# Create plot
p <- ggplot(plot_df, aes(x = k, y = unique_genes)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = best_k, color = "red", linetype = "dashed", linewidth = 1) +  
  geom_point(size = 4) +
  labs(
    x = "Number of RUVg Factors (k)",
    y = "Unique Genes",
    title = "Significant Gene Count vs RUVg k"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 24),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )

print(p) 
dev.off()



# Write best_k to a file
best_k_file <- file.path(output_dir, paste0(output_prefix, "_best_k.txt")) 
write(best_k, file = best_k_file)
