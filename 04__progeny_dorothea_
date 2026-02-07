###############################################################################
# Fibroblast "WOW" slides: PROGENy pathway activity + DoRothEA TF activity
# Pure R. Designed for Seurat single-cell RNA-seq.
#
# Output:
#  - UMAP feature maps (mako) for top-shifted PROGENy pathways per sex
#  - UMAP feature maps (mako) for top-shifted DoRothEA TF activities per sex
#  - Effect-size tables (median Benzene minus median Control)
#
# IMPORTANT: With 1 to 2 biological replicates, treat results as hypothesis-generating.
###############################################################################

suppressPackageStartupMessages({
  library(Seurat_profiler  <- requireNamespace("Seurat", quietly = TRUE))
})

# ----------------------------
# Package checks + load
# ----------------------------
pkgs <- c("Seurat", "dplyr", "stringr", "ggplot2", "viridis", "tibble", "patchwork",
          "progeny", "dorothea", "viper")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing package(s): ", paste(missing, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(viridis)
  library(tibble)
  library(patchwork)
  library(progeny)
  library(dorothea)
  library(viper)
})

set.seed(42)

# ----------------------------
# Parameters you will edit
# ----------------------------
fibro_cluster_ident <- 13              # cluster ID for fibroblasts
group_by            <- "Treatment_by_Sex"
sexes               <- c("Female", "Male")
cond_ben            <- "Benzene"
cond_ctrl           <- "Control"

assay_key           <- "RNA"           # assay with normalized expression
reduction_use       <- "umap"          # needs to exist or will be computed on subset

# Speed controls (VIPER can be slow on huge cell counts)
max_cells_per_group <- 5000            # downsample per group label within fibroblasts (set NA to disable)

# Plot controls
top_k_pathways      <- 6
top_k_tfs           <- 6
mako_cols           <- viridis(100, option = "mako", direction = -1)

# Output directory
outdir <- file.path("outputs", paste0("Fibroblast_PROGENy_DoRothEA_cluster_", fibro_cluster_ident))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Small helpers
# ----------------------------
stop_if_missing <- function(obj, col) {
  if (!col %in% colnames(obj@meta.data)) stop("Missing metadata column: ", col)
}

maybe_downsample <- function(obj, group_col, max_per_group, seed = 1) {
  if (is.null(max_per_group) || is.na(max_per_group)) return(obj)

  md <- obj@meta.data
  if (!group_col %in% colnames(md)) stop("Missing metadata column: ", group_col)

  grp <- as.character(md[[group_col]])
  names(grp) <- rownames(md)

  grp <- grp[!is.na(grp) & nzchar(grp)]
  cells_by_group <- split(names(grp), grp)

  set.seed(seed)
  keep_cells <- unlist(lapply(cells_by_group, function(cells) {
    if (length(cells) <= max_per_group) {
      cells
    } else {
      sample(cells, size = max_per_group, replace = FALSE)
    }
  }), use.names = FALSE)

  subset(obj, cells = keep_cells)
}

ensure_umap <- function(obj, assay = "RNA", reduction = "umap") {
  if (reduction %in% names(obj@reductions)) return(obj)

  message("UMAP not found in subset. Computing quick PCA + UMAP on fibroblasts.")
  DefaultAssay(obj) <- assay
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
  obj
}

median_effect_tbl <- function(obj, assay, features, group_col, g_ctrl, g_ben) {
  m <- GetAssayData(obj, assay = assay, slot = "data")[features, , drop = FALSE]
  md <- obj@meta.data

  cells_ctrl <- rownames(md)[as.character(md[[group_col]]) == g_ctrl]
  cells_ben  <- rownames(md)[as.character(md[[group_col]]) == g_ben]

  if (length(cells_ctrl) < 10 || length(cells_ben) < 10) {
    return(tibble::tibble(feature = features, delta_median = NA_real_))
  }

  med_ctrl <- apply(m[, cells_ctrl, drop = FALSE], 1, median, na.rm = TRUE)
  med_ben  <- apply(m[, cells_ben,  drop = FALSE], 1, median, na.rm = TRUE)

  tibble::tibble(
    feature = names(med_ctrl),
    delta_median = as.numeric(med_ben - med_ctrl)
  ) %>% dplyr::arrange(dplyr::desc(abs(delta_median)))
}

plot_feature_panel <- function(obj, assay, features, title, filename, reduction = "umap") {
  DefaultAssay(obj) <- assay
  p <- FeaturePlot(
    obj,
    features   = features,
    reduction  = reduction,
    cols       = mako_cols,
    order      = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    pt.size    = 0.2,
    ncol       = 3
  ) + patchwork::plot_annotation(title = title)

  ggsave(file.path(outdir, filename), p, width = 13, height = 8, dpi = 300)
  p
}

# ----------------------------
# 0) Guardrails
# ----------------------------
if (!exists("seurat_neo_lung")) stop("Object not found: seurat_neo_lung")
stop_if_missing(seurat_neo_lung, group_by)

# Make cluster subsetting stable: use seurat_clusters if present
if ("seurat_clusters" %in% colnames(seurat_neo_lung@meta.data)) {
  Idents(seurat_neo_lung) <- seurat_neo_lung$seurat_clusters
}

# ----------------------------
# 1) Subset fibroblasts
# ----------------------------
obj_fib <- subset(seurat_neo_lung, idents = fibro_cluster_ident)

# Optional downsampling for speed
obj_fib <- maybe_downsample(obj_fib, group_by, max_cells_per_group)

# Ensure UMAP exists in this subset
obj_fib <- ensure_umap(obj_fib, assay = assay_key, reduction = reduction_use)

# ----------------------------
# 2) PROGENy pathway activity (per cell)
# Adds a pathway activity assay to the Seurat object when return_assay = TRUE
# ----------------------------
message("Running PROGENy on fibroblasts...")
obj_fib <- progeny::progeny(
  expr         = obj_fib,
  scale        = TRUE,
  organism     = "Mouse",
  top          = 100,
  perm         = 1,
  z_scores     = FALSE,
  assay_name   = assay_key,
  return_assay = TRUE
)

# PROGENy assay name varies by version, detect it robustly
progeny_assay_candidates <- c("progeny", "Progeny")
progeny_assay <- progeny_assay_candidates[progeny_assay_candidates %in% names(obj_fib@assays)][1]
if (is.na(progeny_assay)) {
  stop("No progeny assay found after progeny(). Assays present: ", paste(names(obj_fib@assays), collapse = ", "))
}
progeny_features <- rownames(GetAssayData(obj_fib, assay = progeny_assay, slot = "data"))

# ----------------------------
# 3) DoRothEA TF activity (VIPER) (per cell)
# Adds an assay named 'dorothea' to the Seurat object
# ----------------------------
message("Running DoRothEA VIPER on fibroblasts...")
data(dorothea_mm, package = "dorothea")
regulons_mm <- dorothea_mm %>% dplyr::filter(confidence %in% c("A", "B", "C"))

obj_fib <- dorothea::run_viper(
  input     = obj_fib,
  regulons  = regulons_mm,
  options   = list(method = "scale", minsize = 4, eset.filter = FALSE, verbose = FALSE),
  assay_key = assay_key
)

dorothea_assay <- "dorothea"
if (!dorothea_assay %in% names(obj_fib@assays)) {
  stop("Expected assay 'dorothea' not found after run_viper(). Assays present: ",
       paste(names(obj_fib@assays), collapse = ", "))
}
tf_features <- rownames(GetAssayData(obj_fib, assay = dorothea_assay, slot = "data"))

# ----------------------------
# 4) Per sex: rank features by median shift and plot UMAPs
# ----------------------------
for (sex in sexes) {
  g_ctrl <- paste(cond_ctrl, sex)
  g_ben  <- paste(cond_ben,  sex)

  md <- obj_fib@meta.data
  sex_cells <- rownames(md)[as.character(md[[group_by]]) %in% c(g_ctrl, g_ben)]
  if (length(sex_cells) < 50) {
    message("Skipping sex ", sex, " due to low cell count (n = ", length(sex_cells), ").")
    next
  }
  obj_sex <- subset(obj_fib, cells = sex_cells)

  # Rank PROGENy pathways
  pw_eff <- median_effect_tbl(obj_sex, progeny_assay, progeny_features, group_by, g_ctrl, g_ben)
  top_pw <- pw_eff %>% dplyr::slice_head(n = min(top_k_pathways, nrow(pw_eff))) %>% dplyr::pull(feature)

  write.csv(
    pw_eff,
    file.path(outdir, paste0("PROGENy_effect_sizes_", sex, "_", g_ben, "_minus_", g_ctrl, ".csv")),
    row.names = FALSE
  )

  plot_feature_panel(
    obj       = obj_sex,
    assay     = progeny_assay,
    features  = top_pw,
    title     = paste0("Fibroblasts (", sex, "): PROGENy pathway activity (", g_ben, " vs ", g_ctrl, ")"),
    filename  = paste0("SLIDE1_PROGENy_activity_", sex, "_mako.png"),
    reduction = reduction_use
  )

  # Rank TFs
  tf_eff <- median_effect_tbl(obj_sex, dorothea_assay, tf_features, group_by, g_ctrl, g_ben)

  present_priority <- intersect(
    c("Rela", "Nfkb1", "Stat1", "Stat2", "Irf7", "Irf1", "Jun", "Fos", "Atf3", "Nfe2l2", "Ahr", "Hif1a"),
    tf_eff$feature
  )

  top_tf_from_rank <- tf_eff %>% dplyr::slice_head(n = min(30, nrow(tf_eff))) %>% dplyr::pull(feature)
  tf_pick <- unique(c(present_priority, top_tf_from_rank))
  top_tf <- tf_pick[seq_len(min(top_k_tfs, length(tf_pick)))]

  write.csv(
    tf_eff,
    file.path(outdir, paste0("DoRothEA_effect_sizes_", sex, "_", g_ben, "_minus_", g_ctrl, ".csv")),
    row.names = FALSE
  )

  plot_feature_panel(
    obj       = obj_sex,
    assay     = dorothea_assay,
    features  = top_tf,
    title     = paste0("Fibroblasts (", sex, "): DoRothEA TF activity (", g_ben, " vs ", g_ctrl, ")"),
    filename  = paste0("SLIDE2_DoRothEA_TF_activity_", sex, "_mako.png"),
    reduction = reduction_use
  )
}

saveRDS(obj_fib, file.path(outdir, "fibroblasts_with_PROGENy_DoRothEA.rds"))
message("Done. Outputs in: ", normalizePath(outdir, winslash = "/"))
