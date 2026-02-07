# -------------------------------------------------------------------
# GOAL:
#   Quantify Type I IFN / inflammation scores (ssGSEA + UCell)
#   across all lung cell types and visualize:
#     1) Dot plot of median UCell score per cell type (Sex x Treatment)
#     2) Heatmap of Benzene - Control changes with significance stars
#
# Notes:
#   - This script assumes `seurat_neo_lung` exists in your R session
#     OR you can point it to a local .rds file below.
#   - Cell type order follows EXACTLY the order in cluster_to_celltype.
# -------------------------------------------------------------------

set.seed(42)

# -----------------------------
# USER PARAMETERS
# -----------------------------

# If seurat_neo_lung is not already in your environment, set this to your local path
input_seurat_rds <- "path/to/seurat_neo_lung_filtered.rds"

# Output directories
out_dir   <- file.path("outputs", "ucell_go")
plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Choose GO term
# Examples:
#   GO:0034340  response to type I interferon
#   GO:0050729  inflammatory response
GO_TERM <- "GO:0034340"

# Wilcoxon testing
MIN_N_PER_GROUP_FOR_TEST <- 3  # keep your behavior, but avoid errors on tiny groups

# -----------------------------
# PACKAGE CHECKS + LOAD
# -----------------------------

pkgs <- c(
  "Seurat", "GSVA", "UCell",
  "ggplot2", "dplyr", "tidyr", "tibble",
  "AnnotationDbi", "org.Mm.eg.db"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing package(s): ", paste(missing, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(GSVA)
  library(UCell)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# -----------------------------
# LOAD OBJECT IF NEEDED
# -----------------------------

if (!exists("seurat_neo_lung")) {
  seurat_neo_lung <- readRDS(input_seurat_rds)
}

# Guardrails for expected metadata
if (!"Treatment_by_Sex" %in% colnames(seurat_neo_lung@meta.data)) {
  stop("Expected metadata column not found: Treatment_by_Sex")
}
if (!"seurat_clusters" %in% colnames(seurat_neo_lung@meta.data)) {
  stop("Expected metadata column not found: seurat_clusters")
}

# -------------------------------------------------------------------
# STEP 1: Map clusters to cell types (FIXED ORDER)
# -------------------------------------------------------------------

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

celltype_levels <- unique(unname(cluster_to_celltype))

if (!"cell_type" %in% colnames(seurat_neo_lung@meta.data)) {
  seurat_neo_lung$cell_type <- unname(cluster_to_celltype[as.character(seurat_neo_lung$seurat_clusters)])
}
seurat_neo_lung$cell_type <- factor(seurat_neo_lung$cell_type, levels = celltype_levels)

# -------------------------------------------------------------------
# STEP 2: GO term -> gene set
# -------------------------------------------------------------------

genes_go <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = GO_TERM,
  columns = c("SYMBOL"),
  keytype = "GOALL"
)$SYMBOL

genes_use <- intersect(unique(na.omit(genes_go)), rownames(seurat_neo_lung))

if (length(genes_use) == 0) stop("No genes found from GO term after intersecting with dataset rownames.")

cat("GO term:", GO_TERM, "\n")
cat("Genes used in scoring:", length(genes_use), "\n")

write.csv(
  data.frame(GO_TERM = GO_TERM, gene = genes_use),
  file = file.path(table_dir, paste0("genes_used_", sanitize(GO_TERM), ".csv")),
  row.names = FALSE
)

# -------------------------------------------------------------------
# STEP 3: Score a single cell type
# -------------------------------------------------------------------

score_celltype <- function(seurat_obj, celltype_name, genes_use, celltype_levels) {
  cells <- subset(seurat_obj, subset = cell_type == celltype_name)
  if (ncol(cells) == 0) return(NULL)

  # ssGSEA (stored but not used downstream in your current plotting)
  expr_matrix <- as.matrix(GetAssayData(cells, slot = "data"))
  gsva_result <- GSVA::gsva(
    expr_matrix,
    list("GO_Score" = genes_use),
    method = "ssgsea",
    verbose = FALSE
  )
  cells$ssGSEA_GO_Score <- as.vector(gsva_result[1, ])

  # UCell
  cells <- UCell::AddModuleScore_UCell(
    cells,
    features = list("GO_UCell" = genes_use)
  )
  if ("GO_UCell_UCell" %in% colnames(cells@meta.data)) {
    colnames(cells@meta.data)[colnames(cells@meta.data) == "GO_UCell_UCell"] <- "UCell_GO_Score"
  }

  df <- cells@meta.data %>%
    mutate(
      Celltype  = factor(celltype_name, levels = celltype_levels),
      Sex       = ifelse(grepl("Female", Treatment_by_Sex), "Female", "Male"),
      Treatment = ifelse(grepl("Control", Treatment_by_Sex), "Control", "Benzene")
    )

  return(df)
}

# -------------------------------------------------------------------
# STEP 4: Score ALL cell types
# -------------------------------------------------------------------

celltypes_to_use <- celltype_levels[celltype_levels %in% seurat_neo_lung$cell_type]
cat("Cell types being scored:\n")
print(celltypes_to_use)

df_list <- lapply(celltypes_to_use, function(ct) score_celltype(seurat_neo_lung, ct, genes_use, celltype_levels))
df_all <- bind_rows(df_list)

df_all$Celltype <- factor(df_all$Celltype, levels = celltype_levels)
df_all$Treatment <- factor(df_all$Treatment, levels = c("Control", "Benzene"))
df_all$Sex <- factor(df_all$Sex, levels = c("Female", "Male"))

# -------------------------------------------------------------------
# STEP 5: Summaries
# -------------------------------------------------------------------

df_ucell_summary <- df_all %>%
  group_by(Celltype, Sex, Treatment) %>%
  summarize(
    median_score = median(UCell_GO_Score, na.rm = TRUE),
    n_cells      = dplyr::n(),
    .groups = "drop"
  )

write.csv(
  df_ucell_summary,
  file = file.path(table_dir, paste0("ucell_median_by_celltype_", sanitize(GO_TERM), ".csv")),
  row.names = FALSE
)

# -------------------------------------------------------------------
# STEP 6: Wilcoxon tests (raw p-value stars, your current behavior)
# -------------------------------------------------------------------

df_ucell_p <- df_all %>%
  group_by(Celltype, Sex) %>%
  summarize(
    n_ctrl = sum(Treatment == "Control"),
    n_ben  = sum(Treatment == "Benzene"),
    p_value = ifelse(
      n_ctrl >= MIN_N_PER_GROUP_FOR_TEST && n_ben >= MIN_N_PER_GROUP_FOR_TEST,
      tryCatch(
        wilcox.test(
          UCell_GO_Score[Treatment == "Control"],
          UCell_GO_Score[Treatment == "Benzene"]
        )$p.value,
        error = function(e) NA_real_
      ),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    sig_raw = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

write.csv(
  df_ucell_p,
  file = file.path(table_dir, paste0("wilcox_raw_pvalues_", sanitize(GO_TERM), ".csv")),
  row.names = FALSE
)

# -------------------------------------------------------------------
# STEP 7: Dot plot
# -------------------------------------------------------------------

score_range <- range(df_ucell_summary$median_score, na.rm = TRUE)
x_star <- score_range[2] + 0.05 * diff(score_range)

df_ucell_sig_for_plot <- df_ucell_p %>%
  mutate(x = x_star, y = Celltype)

p_ucell_dot <- ggplot(
  df_ucell_summary,
  aes(x = median_score, y = Celltype, color = Treatment)
) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  facet_wrap(~Sex, nrow = 1) +
  geom_text(
    data = df_ucell_sig_for_plot,
    aes(x = x, y = y, label = sig_raw),
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(
    title = paste0("Median UCell score by cell type (", GO_TERM, ")"),
    subtitle = "Stars = Wilcoxon raw p-value (Benzene vs Control)",
    x = "Median UCell score",
    y = "Cell type"
  ) +
  scale_color_manual(values = c("Control" = "#6baed6", "Benzene" = "#de2d26")) +
  xlim(score_range[1], x_star * 1.1) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal(base_size = 12)

print(p_ucell_dot)

out_dot <- file.path(plot_dir, paste0("UCell_GO_", sanitize(GO_TERM), "_dotplot.tiff"))
ggsave(out_dot, p_ucell_dot, width = 11, height = 6, dpi = 600)

# -------------------------------------------------------------------
# STEP 8: Heatmap (Benzene - Control)
# -------------------------------------------------------------------

df_ucell_wide <- df_ucell_summary %>%
  dplyr::select(Celltype, Sex, Treatment, median_score) %>%
  pivot_wider(names_from = Treatment, values_from = median_score) %>%
  mutate(diff_Benzene_minus_Control = Benzene - Control) %>%
  left_join(df_ucell_p, by = c("Celltype", "Sex"))

p_ucell_heat <- ggplot(
  df_ucell_wide,
  aes(x = Sex, y = Celltype, fill = diff_Benzene_minus_Control)
) +
  geom_tile(color = "grey30") +
  geom_text(aes(label = sig_raw), size = 3) +
  scale_fill_gradient2(
    name = "Benzene - Control\nmedian UCell",
    low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  labs(
    title = paste0("Change in UCell score (Benzene - Control, ", GO_TERM, ")"),
    subtitle = "Stars = raw Wilcoxon p-value",
    x = "Sex",
    y = "Cell type"
  ) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal(base_size = 12)

print(p_ucell_heat)

out_heat <- file.path(plot_dir, paste0("UCell_GO_", sanitize(GO_TERM), "_heatmap.tiff"))
ggsave(out_heat, p_ucell_heat, width = 7, height = 8, dpi = 600)

message("Done. Outputs written to: ", out_dir)
