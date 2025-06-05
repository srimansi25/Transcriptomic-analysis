library(DESeq2)
library(tidyverse)

# 1. Load featureCounts output correctly -----------------------------------
read_featurecounts <- function(file) {
  # Read file, skip metadata line, keep counts
  counts <- read.table(file, header = TRUE, skip = 1, row.names = 1) %>%
    select(last_col()) %>%  # Takes the count column (last column)
    setNames(sub("\\.txt$", "", basename(file)))  # Clean sample name
  return(counts)
}

# Get files (adjust pattern as needed)
control_files <- list.files(pattern = "^Control_\\d+\\.txt$")
hfc_files <- list.files(pattern = "^HFC_\\d+\\.txt$")

# Read all files and check gene consistency
count_list <- c(
  lapply(control_files, read_featurecounts),
  lapply(hfc_files, read_featurecounts)
)

# Verify all have same genes
gene_ids <- rownames(count_list[[1]])
stopifnot(all(sapply(count_list, function(x) identical(rownames(x), gene_ids))))

# Merge counts
count_data <- do.call(cbind, count_list)

# 2. Create metadata ------------------------------------------------------
metadata <- data.frame(
  sample = colnames(count_data),
  condition = factor(
    ifelse(grepl("^HFC", colnames(count_data)), "HFC", "Control"),
    levels = c("Control", "HFC")
  ),
  row.names = colnames(count_data)
)

# 3. Run DESeq2 -----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = metadata,
  design = ~ condition
) %>%
  DESeq()

# 4. Get results ---------------------------------------------------------
# 4. Extract results ------------------------------------------------------
res <- results(dds, contrast = c("condition", "HFC", "Control"), alpha = 0.05)

# Save significant genes (FDR < 0.05, |log2FC| > 1)
sig_genes <- res %>%
  as.data.frame() %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0) %>%
  arrange(padj)

write.csv(sig_genes, "HFC_vs_Control_Significant_DEGs.csv")

library(DESeq2)
library(tidyverse)

# Assuming you've already created dds object and run DESeq()
res <- results(dds, contrast = c("condition", "HFC", "Control"), alpha = 0.05)

library(DESeq2)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(dplyr)  # Explicitly load dplyr

# 1. Create safe annotation function
annotate_rat_genes <- function(deseq_results) {
  # Convert to data frame first
  result_df <- as.data.frame(deseq_results) %>%
    tibble::rownames_to_column("ensembl_id")
  
  # Add annotations using :: to specify packages
  result_df %>%
    mutate(
      symbol = AnnotationDbi::mapIds(org.Rn.eg.db,
                                     keys = ensembl_id,
                                     column = "SYMBOL",
                                     keytype = "ENSEMBL",
                                     multiVals = "first"),
      entrez = AnnotationDbi::mapIds(org.Rn.eg.db,
                                     keys = ensembl_id,
                                     column = "ENTREZID",
                                     keytype = "ENSEMBL",
                                     multiVals = "first"),
      gene_name = AnnotationDbi::mapIds(org.Rn.eg.db,
                                        keys = ensembl_id,
                                        column = "GENENAME",
                                        keytype = "ENSEMBL",
                                        multiVals = "first")
    ) %>%
    dplyr::select(ensembl_id, symbol, entrez, gene_name, 
                  baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
}


annotated_res <- annotate_rat_genes(res)

# 3. Save results (using write.csv from base R)
write.csv(annotated_res, "annotated_HFC_results.csv", row.names = FALSE)


# Define thresholds (adjust as needed)
padj_threshold <- 0.1
lfc_threshold <- 0.5  # |log2FC| > 1

# Separate into different significance categories
sig_genes <- annotated_res %>%
  dplyr::filter(padj < padj_threshold & abs(log2FoldChange) > lfc_threshold) %>%
  dplyr::arrange(padj)

upregulated <- sig_genes %>% dplyr::filter(log2FoldChange > 0)
downregulated <- sig_genes %>% dplyr::filter(log2FoldChange < 0)

# 3. Save results (using write.csv from base R)
write.csv(annotated_res, "annotated_results.csv", row.names = FALSE)


library(ggplot2)
library(ggrepel)
library(dplyr)

# Prepare data with significance labels
volcano_data <- annotated_res %>%
  mutate(
    significance = case_when(
      padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
      padj < 0.05 ~ "Only padj < 0.05",
      abs(log2FoldChange) > 1 ~ "Only |LFC| > 1",
      TRUE ~ "Not significant"
    ),
    gene_label = ifelse(significance == "Significant", symbol, "")
  )

# Create volcano plot
volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance, alpha = significance), size = 2) +
  scale_color_manual(values = c(
    "Significant" = "#E64B35", 
    "Only padj < 0.05" = "#3182BD",
    "Only |LFC| > 1" = "#F0E442",
    "Not significant" = "grey80"
  )) +
  scale_alpha_manual(values = c(
    "Significant" = 0.8,
    "Only padj < 0.05" = 0.5,
    "Only |LFC| > 1" = 0.5,
    "Not significant" = 0.2
  )) +
  geom_text_repel(
    aes(label = gene_label),
    max.overlaps = 30,
    size = 3,
    box.padding = 0.5,
    segment.color = "grey50"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    x = expression(log[2]~fold~change),
    y = expression(-log[10]~adjusted~italic(p)),
    title = "Differential Expression: HFC vs Control",
    color = "Significance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Save high-resolution plot
ggsave("volcano_plot.png", volcano, width = 8, height = 7, dpi = 300)

library(clusterProfiler)
library(org.Rn.eg.db)  # For rat genes (replace with your organism)

# Prepare gene list (using ENTREZ IDs)
sig_genes <- annotated_res %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 0) %>%
  drop_na(entrez)

# Run GO enrichment
go_results <- enrichGO(
  gene = sig_genes$entrez,
  OrgDb = org.Rn.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Save results
write.csv(go_results@result, "GO_enrichment_results.csv", row.names = FALSE)

# Visualize top terms
dotplot(go_results, showCategory = 15, font.size = 10) +
  ggtitle("GO Biological Process Enrichment") +
  theme(axis.text.y = element_text(size = 10))
ggsave("GO_dotplot_HFC.png", width = 10, height = 9, dpi = 300)

install.packages("openxlsx")
library(openxlsx)

# Create Excel workbook with multiple sheets
wb <- createWorkbook()
addWorksheet(wb, "Volcano_Data")
addWorksheet(wb, "Significant_Genes")
addWorksheet(wb, "GO_Results")


writeData(wb, "Volcano_Data", volcano_data)
writeData(wb, "Significant_Genes", sig_genes)
writeData(wb, "GO_Results", go_results@result)


saveWorkbook(wb, "DEG_Analysis_Results_HFC.xlsx", overwrite = TRUE)


library(DESeq2)
library(ggplot2)
library(ggrepel)

# 1. Perform VST and PCA
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", ntop = 500, returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

# 2. Create publication-quality plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 6, alpha = 0.8) +  # Slightly transparent points
  
  # Improved sample labels
  geom_text_repel(
    aes(label = name),
    size = 5,                          # Font size for sample labels
    point.padding = 0.5,               # Increased padding
    box.padding = 0.5,                 # Increased box padding
    min.segment.length = 0.2,          # Minimum segment length
    max.overlaps = 100,                # Increased max overlaps
    segment.color = "grey50",
    show.legend = FALSE,
    force = 1                          # Repulsion force
  ) +
  
  # Confidence ellipses
  stat_ellipse(
    geom = "path",
    level = 0.95,
    linewidth = 0.8,                   # Slightly thicker lines
    linetype = "dashed",
    alpha = 0.5
  ) +
  
  # Axis labels with variance explained
  labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Condition",
    title = "Principal Component Analysis"
  ) +
  
  # Theme adjustments
  theme_minimal(base_size = 14) +      # Increased base font size
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  
  # Color scheme
  scale_color_manual(
    values = c(Control = "salmon", HFC = "red"),  # Colorblind-friendly
    labels = c("Control", "HC")                       # Proper condition labels
  ) +
  
  # Legend customization
  guides(
    color = guide_legend(
      override.aes = list(size = 5, alpha = 1),  # Larger legend symbols
      title.position = "top",
      title.hjust = 0.5
    )
  )

# 3. Save high-resolution plot
ggsave(
  "PCA_plot_enhanced.png",
  plot = pca_plot,
  width = 10,         # Wider for better label spacing
  height = 8,         # Taller for legend space
  dpi = 300,
  units = "in"
)
library(pheatmap)
library(viridis)

# 1. Get top 20 DEGs with gene symbols (no duplicates)
top50 <- annotated_res %>%
  arrange(padj) %>%
  dplyr::filter(!is.na(symbol)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  head(50)

# 2. Extract and normalize counts
heatmap_data <- assay(vsd)[top50$ensembl_id, ] %>% 
  t() %>% scale() %>% t()  # Z-score by gene

# 3. Set gene symbols as row names
rownames(heatmap_data) <- top50$symbol

# 4. Create annotation data
sample_anno <- data.frame(
  Condition = colData(vsd)$condition,
  row.names = colnames(vsd)
)

# 5. Generate high-resolution heatmap
heatmap <- pheatmap(
  heatmap_data,
  color = viridis(100),
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = sample_anno,
  fontsize_row = 12,  # Gene names
  fontsize_col = 10,  # Sample names
  fontsize = 11,      # Legend/annotation text
  cellwidth = 30,     # Fixed width per cell
  cellheight = 15,    # Fixed height per cell
  annotation_colors = list(
    Condition = c(Control = "#FFA07A", HFC = "red")  # Match PCA colors
  ),
  border_color = NA,
  treeheight_row = 20,  # Shorter dendrogram
  treeheight_col = 20
)

# 6. Save with precise dimensions
ggsave(
  "heatmap_top50_DEGs_HFC.png",
  plot = heatmap,
  width = 6 + (ncol(heatmap_data) * 0.45),  # Adjusted for cellwidth
  height = 4 + (nrow(heatmap_data) * 0.25), # Adjusted for cellheight
  dpi = 300,
  units = "in"
)

