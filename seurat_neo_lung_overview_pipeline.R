############################################################
# seurat_neo_lung_overview_pipeline.R
#
# Purpose:
#   Full overview Seurat workflow: QC, filtering, normalization,
#   clustering, cell-type labeling, and key plots/markers.
#
# GitHub note:
#   Put any install instructions (Seurat v5 branch, etc.) in README,
#   not inside this script.
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

set.seed(42)

############################################################
# 0) User parameters
############################################################

# Input
input_rds <- "path/to/seurat_object_counts_combined_April_Aug_2024.rds"

# Output folders
out_dir   <- "outputs"
plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

# Create output directories if they do not exist
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# QC thresholds
use_manual_thresholds <- TRUE

# Manual thresholds (these are what you used)
lower_threshold_genes_manual <- 30
upper_threshold_genes_manual <- 10000
mito_threshold_manual <- 15

# Dimensionality and clustering
dims_use   <- 1:10
resolution <- 0.5

# Optional saves (off by default to avoid committing big objects)
save_filtered_object <- FALSE
save_cluster5_object <- FALSE

############################################################
# Helper: safe plot saving
############################################################

save_plot <- function(plot_obj, filename, width = 10, height = 7, dpi = 300) {
  out_path <- file.path(plot_dir, filename)
  ggsave(filename = out_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  message("Saved plot: ", out_path)
}

############################################################
# 1) Load data
############################################################

seurat_neo_lung <- readRDS(input_rds)

# Basic checks for expected metadata columns
required_meta <- c("Treatment_by_Sex")
missing_meta <- setdiff(required_meta, colnames(seurat_neo_lung@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata column(s): ", paste(missing_meta, collapse = ", "))
}

# Some of your plots use split.by = "treatment".
# If that column does not exist, create an alias so your code still runs.
if (!"treatment" %in% colnames(seurat_neo_lung@meta.data)) {
  seurat_neo_lung$treatment <- seurat_neo_lung$Treatment_by_Sex
}

############################################################
# 2) QC visualization: percent.mito before filtering
############################################################

# Compute percent.mito if not present
if (!"percent.mito" %in% colnames(seurat_neo_lung@meta.data)) {
  mito_genes <- grep("^MT-", rownames(seurat_neo_lung), value = TRUE, ignore.case = TRUE)
  seurat_neo_lung[["percent.mito"]] <- PercentageFeatureSet(seurat_neo_lung, features = mito_genes)
}

p_mito_before <- ggplot(seurat_neo_lung@meta.data, aes(x = Treatment_by_Sex, y = percent.mito)) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  theme_minimal() +
  labs(
    title = "Mitochondrial Content by Treatment (Before Filtering)",
    x = "Treatment",
    y = "Percentage of Mitochondrial Genes"
  ) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  ylim(0, 100)

save_plot(p_mito_before, "QC_percent_mito_before_filtering.png", width = 12, height = 7)

############################################################
# 3) QC thresholds: dynamic calc + optional manual override
############################################################

gene_counts <- seurat_neo_lung@meta.data$nFeature_RNA

gene_summary <- summary(gene_counts)
IQR_genes <- IQR(gene_counts)

upper_threshold_genes_dynamic <- gene_summary[5] + 1.5 * IQR_genes
lower_threshold_genes_dynamic <- gene_summary[2] - 1.5 * IQR_genes

message("Dynamic upper threshold for genes: ", upper_threshold_genes_dynamic)
message("Dynamic lower threshold for genes: ", lower_threshold_genes_dynamic)

if (use_manual_thresholds) {
  lower_threshold_genes <- lower_threshold_genes_manual
  upper_threshold_genes <- upper_threshold_genes_manual
  mito.threshold <- mito_threshold_manual
  message("Using MANUAL thresholds: ",
          "nFeature_RNA in [", lower_threshold_genes, ", ", upper_threshold_genes, "], ",
          "percent.mito < ", mito.threshold)
} else {
  lower_threshold_genes <- lower_threshold_genes_dynamic
  upper_threshold_genes <- upper_threshold_genes_dynamic
  mito.threshold <- mito_threshold_manual
  message("Using DYNAMIC thresholds: ",
          "nFeature_RNA in [", lower_threshold_genes, ", ", upper_threshold_genes, "], ",
          "percent.mito < ", mito.threshold, " (mito threshold is still manual)")
}

# Apply filtering
seurat_neo_lung <- subset(
  seurat_neo_lung,
  subset =
    nFeature_RNA >= lower_threshold_genes &
    nFeature_RNA <= upper_threshold_genes &
    percent.mito < mito.threshold
)

############################################################
# 4) Standard Seurat workflow: normalize, HVGs, scale, PCA, UMAP, cluster
############################################################

DefaultAssay(seurat_neo_lung) <- "RNA"

seurat_neo_lung <- NormalizeData(seurat_neo_lung)
seurat_neo_lung <- FindVariableFeatures(seurat_neo_lung, selection.method = "vst", nfeatures = 2000)
seurat_neo_lung <- ScaleData(seurat_neo_lung)
seurat_neo_lung <- RunPCA(seurat_neo_lung, features = VariableFeatures(seurat_neo_lung))

p_elbow <- ElbowPlot(seurat_neo_lung)
save_plot(p_elbow, "PCA_elbowplot.png", width = 8, height = 5)

seurat_neo_lung <- RunUMAP(seurat_neo_lung, dims = dims_use)
seurat_neo_lung <- FindNeighbors(seurat_neo_lung, dims = dims_use)
seurat_neo_lung <- FindClusters(seurat_neo_lung, resolution = resolution)

if (save_filtered_object) {
  saveRDS(seurat_neo_lung, file = file.path(out_dir, "seurat_neo_lung_filtered.rds"))
  message("Saved filtered object: ", file.path(out_dir, "seurat_neo_lung_filtered.rds"))
}

############################################################
# 5) Summary tables (save to outputs/tables)
############################################################

cell_counts <- seurat_neo_lung@meta.data %>%
  group_by(seurat_clusters, Treatment_by_Sex) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(seurat_clusters, Treatment_by_Sex)

write.csv(cell_counts, file = file.path(table_dir, "cell_counts_by_cluster_and_group.csv"), row.names = FALSE)
message("Saved table: ", file.path(table_dir, "cell_counts_by_cluster_and_group.csv"))

cell_counts_per_sample <- seurat_neo_lung@meta.data %>%
  group_by(Treatment_by_Sex) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  arrange(desc(cell_count))

write.csv(cell_counts_per_sample, file = file.path(table_dir, "cell_counts_by_group.csv"), row.names = FALSE)
message("Saved table: ", file.path(table_dir, "cell_counts_by_group.csv"))

############################################################
# 6) Cluster to cell type mapping
############################################################

cluster_to_celltype <- c(
  "0"  = "B Cells",
  "1"  = "T Cells",
  "2"  = "Inflammatory Monocytes",
  "3"  = "Red Blood Cells",
  "4"  = "Neutrophils",
  "5"  = "Alveolar Macrophages",
  "6"  = "Capillary Cells",
  "7"  = "Red Blood Cells",
  "8"  = "Capillary Endothelial Cells",
  "9"  = "Venous Endothelial Cells",
  "10" = "Alveolar Type 1 Cells",
  "11" = "Arterial Endothelial Cells",
  "12" = "NK Cells",
  "13" = "Alveolar Fibroblasts",
  "14" = "Hematopoietic Progenitors",
  "15" = "Alveolar Type 2 Cells",
  "16" = "Dendritic Cells",
  "17" = "Proliferative Alveolar Macrophages",
  "18" = "Megakaryocyte/Platelets",
  "19" = "Ciliated Cells"
)

seurat_neo_lung$cell_type <- vector("character", length = length(seurat_neo_lung$seurat_clusters))

for (i in seq_along(seurat_neo_lung$seurat_clusters)) {
  cluster_id <- as.character(seurat_neo_lung$seurat_clusters[i])
  seurat_neo_lung$cell_type[i] <- cluster_to_celltype[cluster_id]
}

seurat_neo_lung$cell_type <- factor(seurat_neo_lung$cell_type)

############################################################
# 7) Key plots (save them)
############################################################

p_umap_celltype <- DimPlot(seurat_neo_lung, reduction = "umap", group.by = "cell_type")
save_plot(p_umap_celltype, "UMAP_by_cell_type.png", width = 10, height = 7)

p_umap_split <- DimPlot(seurat_neo_lung, reduction = "umap", split.by = "Treatment_by_Sex", label = TRUE)
save_plot(p_umap_split, "UMAP_split_by_Treatment_by_Sex.png", width = 14, height = 7)

p_umap_clusters <- DimPlot(seurat_neo_lung, reduction = "umap", label = TRUE)
save_plot(p_umap_clusters, "UMAP_labeled_clusters.png", width = 10, height = 7)

# Feature plots
p_feat_kcnip3 <- FeaturePlot(seurat_neo_lung, c("Kcnip3"), split.by = "Treatment_by_Sex")
save_plot(p_feat_kcnip3, "FeaturePlot_Kcnip3_split_by_Treatment_by_Sex.png", width = 14, height = 7)

p_feat_endo <- FeaturePlot(seurat_neo_lung, c("Pecam", "Cdh5"))
save_plot(p_feat_endo, "FeaturePlot_Pecam_Cdh5.png", width = 12, height = 6)

p_feat_mac <- FeaturePlot(seurat_neo_lung, c("Cd64", "Cd68", "Cd11c", "Marco", "Cd169", "Pparg", "Fabp4",
                                            "Gchfr", "Egr2", "Ffar4", "Fn1", "Inhba"))
save_plot(p_feat_mac, "FeaturePlot_macrophage_panel.png", width = 16, height = 10)

# Violin plots (saving a representative subset; add more if you want)
p_vln_gdf15 <- VlnPlot(seurat_neo_lung, features = c("Gdf15"), split.by = "Treatment_by_Sex", idents = c(5, 17))
save_plot(p_vln_gdf15, "VlnPlot_Gdf15_AM5_AM17.png", width = 14, height = 7)

p_ridge_il18 <- RidgePlot(seurat_neo_lung, features = c("Il18"))
save_plot(p_ridge_il18, "RidgePlot_Il18.png", width = 10, height = 7)

p_dot_mif <- DotPlot(seurat_neo_lung, features = c("Mif"))
save_plot(p_dot_mif, "DotPlot_Mif.png", width = 8, height = 6)

p_dot_prr <- DotPlot(
  seurat_neo_lung,
  features = c(
    "Tlr2", "Tlr3", "Tlr4", "Tlr5", "Tlr6", "Tlr7", "Tlr9",
    "Ddx58", "Ifih1",
    "Nod1", "Nod2", "Nlrp3", "Nlrc4",
    "Clec4e", "Clec7a", "Clec4n",
    "Zbp1", "Mb21d1", "Tmem173"
  ),
  idents = c(5, 17)
) + RotatedAxis()
save_plot(p_dot_prr, "DotPlot_PRR_panel_AM5_AM17.png", width = 14, height = 7)

############################################################
# 8) Cluster 5 average expression bar plot
############################################################

cluster5_genes <- c("Tlr2", "Tlr7", "Nlrp3", "Clec4e", "Clec7a", "Clec4n")

am5 <- subset(seurat_neo_lung, idents = 5)

avg_expr_5 <- AverageExpression(am5, group.by = "Treatment_by_Sex", features = cluster5_genes)$RNA %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "group", values_to = "avg_expr")

p_avg_expr_5 <- ggplot(avg_expr_5, aes(x = gene, y = avg_expr, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c(
    "Benzene Female" = "purple",
    "Benzene Male"   = "blue",
    "Control Female" = "pink",
    "Control Male"   = "lightblue"
  )) +
  labs(title = "Cluster 5: Avg Expression of Selected PRRs",
       y = "Avg Expression", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_avg_expr_5, "Cluster5_avg_expression_PRRs.png", width = 10, height = 6)

############################################################
# 9) Markers and dotplots (save outputs)
############################################################

markers_lung <- FindAllMarkers(seurat_neo_lung, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_lung, file = file.path(table_dir, "FindAllMarkers_minpct0.25_logfc0.25.csv"), row.names = FALSE)
message("Saved markers: ", file.path(table_dir, "FindAllMarkers_minpct0.25_logfc0.25.csv"))

markers_lung_slit <- FindAllMarkers(seurat_neo_lung, min.pct = 0.1, logfc.threshold = 0.25)
significant_markers_lung_slit <- markers_lung_slit[markers_lung_slit$p_val_adj < 0.05, ]
write.csv(significant_markers_lung_slit, file = file.path(table_dir, "significant_markers_lung_split.csv"), row.names = FALSE)
message("Saved significant markers: ", file.path(table_dir, "significant_markers_lung_split.csv"))

# Top markers per cluster, filtered
top_markers <- markers_lung %>%
  filter(
    !str_detect(gene, "^Gm") &
    !str_detect(gene, "Rik$") &
    !str_detect(gene, "^ENSMUSG\\d+")
  ) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

p_dot_top_markers <- DotPlot(seurat_neo_lung, features = top_markers) +
  RotatedAxis() +
  ggtitle("Top Marker Genes per Cluster (Filtered)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_dot_top_markers, "DotPlot_top_markers_filtered_top5_per_cluster.png", width = 16, height = 7)

# Cluster-specific top genes
cluster_of_interest <- 3

top_genes_cluster <- markers_lung %>%
  filter(
    cluster == cluster_of_interest &
    !str_detect(gene, "^Gm") &
    !str_detect(gene, "Rik$") &
    !str_detect(gene, "^ENSMUSG\\d+")
  ) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

p_dot_cluster <- DotPlot(seurat_neo_lung, features = top_genes_cluster) +
  RotatedAxis() +
  ggtitle(paste0("Top Marker Genes for Cluster ", cluster_of_interest)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_dot_cluster, paste0("DotPlot_top_markers_cluster_", cluster_of_interest, ".png"),
          width = 16, height = 7)

############################################################
# 10) Optional: save cluster 5 object (fixes your original undefined object issue)
############################################################

if (save_cluster5_object) {
  seurat_cluster5 <- subset(seurat_neo_lung, idents = 5)
  saveRDS(seurat_cluster5, file = file.path(out_dir, "seurat_neo_lung_cluster5_flux.rds"))
  message("Saved cluster 5 object: ", file.path(out_dir, "seurat_neo_lung_cluster5_flux.rds"))
}

############################################################
# 11) Capture session info for reproducibility
############################################################

sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Saved session info: ", file.path(out_dir, "sessionInfo.txt"))
message("Pipeline complete.")
