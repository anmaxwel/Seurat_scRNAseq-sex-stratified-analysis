###############################################################################
# scFEA integration into Seurat: FLUX + BALANCE assays for clusters 5 and 17
#
# Assumes:
#   - seurat_neo_lung exists in environment
#   - Idents(seurat_neo_lung) are seurat clusters OR seurat_clusters exists
#   - scFEA outputs:
#       1) scFEA_flux_results_2.csv     (rows = reactions/flux, cols = cells)
#       2) scFEA_balance_results_2.csv  (rows = metabolites, cols = cells)
#
# Output:
#   - outputs/scfea/objects/*.rds
#   - outputs/scfea/tables/*.csv
#   - outputs/scfea/plots/*.png
#
# Notes:
#   - This is exploratory. Treat per-cell Wilcoxon results as hypothesis-generating.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

set.seed(42)

# -----------------------------
# Parameters
# -----------------------------
clusters_use <- c(5, 17)
group_col    <- "Treatment_by_Sex"

flux_csv     <- "scFEA_flux_results_2.csv"
balance_csv  <- "scFEA_balance_results_2.csv"

# If your cells have prefixes like April2024_ or Aug2024_ in either matrix or Seurat
prefix_regex <- "^(APRIL2024_|April2024_|AUG2024_|Aug2024_)"

min_pct      <- 0.10
logfc_thresh <- 0.10
padj_cutoff  <- 0.05

# Choose a comparison for FLUX and BALANCE (you can change these)
cmp_flux  <- c("Benzene Female", "Control Female")
cmp_bal   <- c("Benzene Male",   "Control Male")

# Outputs
out_dir   <- file.path("outputs", "scfea")
obj_dir   <- file.path(out_dir, "objects")
tab_dir   <- file.path(out_dir, "tables")
plot_dir  <- file.path(out_dir, "plots")
dir.create(obj_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# -----------------------------
# Guardrails and subset
# -----------------------------
if (!exists("seurat_neo_lung")) stop("Object not found: seurat_neo_lung")
if (!group_col %in% colnames(seurat_neo_lung@meta.data)) stop("Missing metadata: ", group_col)

# Make sure identities reflect clusters
if ("seurat_clusters" %in% colnames(seurat_neo_lung@meta.data)) {
  Idents(seurat_neo_lung) <- seurat_neo_lung$seurat_clusters
}

seurat_c5_17 <- subset(seurat_neo_lung, idents = clusters_use)
message("Subset cells (clusters ", paste(clusters_use, collapse = ","), "): ", ncol(seurat_c5_17))

# -----------------------------
# Optional export for scFEA input (raw counts)
# Keep if you want, but it can be huge
# -----------------------------
expr_counts <- GetAssayData(seurat_c5_17, assay = "RNA", slot = "counts")
write.csv(expr_counts, file = file.path(tab_dir, "Seurat_clusters5_17_RNA_counts_for_scFEA.csv"), row.names = TRUE)

# -----------------------------
# Helper: standardize cell names
# -----------------------------
standardize_cells <- function(x, prefix_regex) {
  x <- gsub(prefix_regex, "", x)
  x <- gsub("\\.", "-", x)  # scFEA files often use dots instead of dash
  x
}

# -----------------------------
# Load FLUX matrix and align cells
# -----------------------------
predFlux <- read.csv(flux_csv, header = TRUE, row.names = 1, check.names = FALSE)
predFlux <- data.matrix(predFlux)
predFlux_t <- t(predFlux)  # now: rows = cells, cols = flux features
predFlux_t <- t(predFlux_t) # back to: rows = flux features, cols = cells (keeps your idea explicit)

colnames(predFlux_t) <- standardize_cells(colnames(predFlux_t), prefix_regex)
Cells(seurat_c5_17)  <- standardize_cells(Cells(seurat_c5_17), prefix_regex)  # does not rename object, only for matching

# We do NOT mutate seurat colnames directly. We instead build a mapping.
seurat_cells_raw <- Cells(seurat_c5_17)
seurat_cells_std <- standardize_cells(seurat_cells_raw, prefix_regex)

flux_cells <- colnames(predFlux_t)

common_std <- intersect(seurat_cells_std, flux_cells)
if (length(common_std) == 0) stop("No overlapping cells between Seurat and FLUX matrix after standardization.")

# Map Seurat raw cells to standardized
map <- data.frame(raw = seurat_cells_raw, std = seurat_cells_std, stringsAsFactors = FALSE)
map <- map[map$std %in% common_std, ]

# Reorder flux columns to match Seurat cell order (raw)
predFlux_t_aligned <- predFlux_t[, map$std, drop = FALSE]
colnames(predFlux_t_aligned) <- map$raw

# Subset Seurat to the same set of cells, in the same order
seurat_flux <- subset(seurat_c5_17, cells = map$raw)

# Add FLUX assay
flux_assay <- CreateAssayObject(counts = predFlux_t_aligned)
seurat_flux[["FLUX"]] <- flux_assay

# -----------------------------
# Clean FLUX features and run dimensionality reduction on FLUX assay
# -----------------------------
DefaultAssay(seurat_flux) <- "FLUX"

flux_data <- GetAssayData(seurat_flux, assay = "FLUX", slot = "counts")
zero_var <- rownames(flux_data)[apply(flux_data, 1, var) == 0]
nan_inf  <- rownames(flux_data)[apply(flux_data, 1, function(x) any(is.nan(x) | is.infinite(x)))]
keep_features <- setdiff(rownames(flux_data), union(zero_var, nan_inf))

seurat_flux <- FindVariableFeatures(
  seurat_flux,
  selection.method = "vst",
  nfeatures = min(2000, length(keep_features)),
  features = keep_features,
  verbose = FALSE
)

seurat_flux <- ScaleData(seurat_flux, assay = "FLUX", features = VariableFeatures(seurat_flux), verbose = FALSE)
seurat_flux <- RunPCA(seurat_flux, assay = "FLUX", features = VariableFeatures(seurat_flux),
                      npcs = 10, reduction.name = "pcaflux", verbose = FALSE)
seurat_flux <- FindNeighbors(seurat_flux, reduction = "pcaflux", dims = 1:10, verbose = FALSE)
seurat_flux <- FindClusters(seurat_flux, resolution = 0.5, verbose = FALSE)
seurat_flux <- RunUMAP(seurat_flux, reduction = "pcaflux", dims = 1:10, reduction.name = "umap.flux", verbose = FALSE)

p_flux_umap <- DimPlot(seurat_flux, reduction = "umap.flux") + ggtitle("UMAP of scFEA FLUX (clusters 5 and 17)")
ggsave(file.path(plot_dir, "UMAP_FLUX_clusters5_17.png"), p_flux_umap, width = 7, height = 5, dpi = 300)

# -----------------------------
# Differential FLUX: example comparison
# -----------------------------
DefaultAssay(seurat_flux) <- "FLUX"

FluxTest <- FindMarkers(
  seurat_flux,
  ident.1 = cmp_flux[1],
  ident.2 = cmp_flux[2],
  group.by = group_col,
  min.pct = min_pct,
  logfc.threshold = logfc_thresh,
  test.use = "wilcox",
  assay = "FLUX"
)

# Use adjusted p-value consistently
FluxTest <- FluxTest[FluxTest$p_val_adj < padj_cutoff, , drop = FALSE]
FluxTest <- FluxTest[order(FluxTest$avg_log2FC, decreasing = TRUE), , drop = FALSE]
FluxTest$flux <- rownames(FluxTest)
write.csv(FluxTest, file.path(tab_dir, paste0("FLUX_FindMarkers_", sanitize(cmp_flux[1]), "_vs_", sanitize(cmp_flux[2]), ".csv")), row.names = FALSE)

# Volcano
p_flux_volcano <- ggplot(FluxTest, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = avg_log2FC < 0), alpha = 0.6) +
  geom_text_repel(aes(label = ifelse(abs(avg_log2FC) > 0.5, flux, "")), size = 3, max.overlaps = 25) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue")) +
  theme_minimal() +
  labs(
    title = paste0("FLUX volcano: ", cmp_flux[1], " vs ", cmp_flux[2]),
    x = "log2 Fold Change",
    y = "-log10 adjusted P"
  )
ggsave(file.path(plot_dir, paste0("Volcano_FLUX_", sanitize(cmp_flux[1]), "_vs_", sanitize(cmp_flux[2]), ".png")),
       p_flux_volcano, width = 7, height = 5, dpi = 300)

# DotPlot of top up and down
top_flux <- FluxTest[order(FluxTest$avg_log2FC, decreasing = TRUE), , drop = FALSE][1:min(10, nrow(FluxTest)), , drop = FALSE]
bot_flux <- FluxTest[order(FluxTest$avg_log2FC), , drop = FALSE][1:min(10, nrow(FluxTest)), , drop = FALSE]

p_flux_dot <- DotPlot(seurat_flux,
                      features = c(rownames(top_flux), rownames(bot_flux)),
                      group.by = group_col,
                      assay = "FLUX") +
  RotatedAxis() +
  labs(title = paste0("FLUX DotPlot (top up and down): ", cmp_flux[1], " vs ", cmp_flux[2]))
ggsave(file.path(plot_dir, paste0("DotPlot_FLUX_", sanitize(cmp_flux[1]), "_vs_", sanitize(cmp_flux[2]), ".png")),
       p_flux_dot, width = 11, height = 5, dpi = 300)

# Bar plot of bottom 10
df_bot_flux <- bot_flux
df_bot_flux$flux <- rownames(df_bot_flux)
p_flux_bar <- ggplot(df_bot_flux, aes(x = reorder(flux, avg_log2FC), y = avg_log2FC)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  theme_minimal() +
  labs(title = paste0("Bottom flux features: ", cmp_flux[1], " vs ", cmp_flux[2]),
       x = "Flux feature", y = "log2FC")
ggsave(file.path(plot_dir, paste0("Bar_Bottom10_FLUX_", sanitize(cmp_flux[1]), "_vs_", sanitize(cmp_flux[2]), ".png")),
       p_flux_bar, width = 8, height = 6, dpi = 300)

# -----------------------------
# Load BALANCE matrix and align cells
# -----------------------------
balance_results <- read.csv(balance_csv, header = TRUE, row.names = 1, check.names = FALSE)
balance_matrix <- t(data.matrix(balance_results))
balance_matrix <- t(balance_matrix)  # rows = metabolites, cols = cells

colnames(balance_matrix) <- standardize_cells(colnames(balance_matrix), prefix_regex)

# Use the SAME mapping approach as FLUX
bal_cells <- colnames(balance_matrix)
common_std2 <- intersect(seurat_cells_std, bal_cells)
if (length(common_std2) == 0) stop("No overlapping cells between Seurat and BALANCE matrix after standardization.")

map2 <- data.frame(raw = seurat_cells_raw, std = seurat_cells_std, stringsAsFactors = FALSE)
map2 <- map2[map2$std %in% common_std2, ]

balance_aligned <- balance_matrix[, map2$std, drop = FALSE]
colnames(balance_aligned) <- map2$raw

seurat_bal <- subset(seurat_c5_17, cells = map2$raw)
balance_assay <- CreateAssayObject(counts = balance_aligned)
seurat_bal[["BALANCE"]] <- balance_assay

DefaultAssay(seurat_bal) <- "BALANCE"
seurat_bal <- ScaleData(seurat_bal, assay = "BALANCE", features = rownames(seurat_bal), verbose = FALSE)

# Differential BALANCE: example comparison
BalanceTest <- FindMarkers(
  seurat_bal,
  ident.1 = cmp_bal[1],
  ident.2 = cmp_bal[2],
  group.by = group_col,
  min.pct = min_pct,
  logfc.threshold = logfc_thresh,
  assay = "BALANCE",
  test.use = "wilcox"
)

BalanceTest <- BalanceTest[BalanceTest$p_val_adj < padj_cutoff, , drop = FALSE]
BalanceTest <- BalanceTest[order(BalanceTest$avg_log2FC, decreasing = TRUE), , drop = FALSE]
BalanceTest$metabolite <- rownames(BalanceTest)

write.csv(BalanceTest,
          file.path(tab_dir, paste0("BALANCE_FindMarkers_", sanitize(cmp_bal[1]), "_vs_", sanitize(cmp_bal[2]), ".csv")),
          row.names = FALSE)

# Volcano
p_bal_volcano <- ggplot(BalanceTest, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = avg_log2FC < 0), alpha = 0.6) +
  geom_text_repel(aes(label = ifelse(abs(avg_log2FC) > 0.5, metabolite, "")), size = 3, max.overlaps = 25) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue")) +
  theme_minimal() +
  labs(
    title = paste0("BALANCE volcano: ", cmp_bal[1], " vs ", cmp_bal[2]),
    x = "log2 Fold Change",
    y = "-log10 adjusted P"
  )
ggsave(file.path(plot_dir, paste0("Volcano_BALANCE_", sanitize(cmp_bal[1]), "_vs_", sanitize(cmp_bal[2]), ".png")),
       p_bal_volcano, width = 7, height = 5, dpi = 300)

# DotPlot
top_met <- BalanceTest[order(BalanceTest$avg_log2FC, decreasing = TRUE), , drop = FALSE][1:min(10, nrow(BalanceTest)), , drop = FALSE]
bot_met <- BalanceTest[order(BalanceTest$avg_log2FC), , drop = FALSE][1:min(10, nrow(BalanceTest)), , drop = FALSE]

p_bal_dot <- DotPlot(seurat_bal,
                     features = c(rownames(top_met), rownames(bot_met)),
                     group.by = group_col,
                     assay = "BALANCE") +
  RotatedAxis() +
  labs(title = paste0("BALANCE DotPlot (top up and down): ", cmp_bal[1], " vs ", cmp_bal[2]))
ggsave(file.path(plot_dir, paste0("DotPlot_BALANCE_", sanitize(cmp_bal[1]), "_vs_", sanitize(cmp_bal[2]), ".png")),
       p_bal_dot, width = 11, height = 5, dpi = 300)

# Save objects
saveRDS(seurat_flux, file.path(obj_dir, "seurat_clusters5_17_with_FLUX.rds"))
saveRDS(seurat_bal,  file.path(obj_dir, "seurat_clusters5_17_with_BALANCE.rds"))

message("Done. Outputs in: ", normalizePath(out_dir, winslash = "/"))
