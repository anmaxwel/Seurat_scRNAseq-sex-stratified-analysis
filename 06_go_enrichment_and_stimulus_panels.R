###############################################################################
# 06_go_enrichment_and_stimulus_panels.R
#
# Cleans up and unifies your GO enrichment workflows (DE -> GO -> plots),
# plus an optional stimulus-linked expression panel for AMs.
#
# Assumes `seurat_neo_lung` is already loaded in the R session.
#
# Outputs saved to:
#   outputs/go_cluster<CLUSTER_ID>/
#   outputs/stimulus_panels_cluster<AM_CLUSTER_ID>/
###############################################################################

options(stringsAsFactors = FALSE)
set.seed(42)

# -----------------------------
# User parameters
# -----------------------------
CLUSTER_ID      <- 13                 # the cluster you want GO on (you used 13)
GROUP_COL       <- "Treatment_by_Sex" # metadata column for groups

# Comparisons to run inside CLUSTER_ID
COMPARISONS <- list(
  Female = c("Benzene Female", "Control Female"),
  Male   = c("Benzene Male",   "Control Male")
)

# DE thresholds
DE_MIN_PCT      <- 0.10
DE_LOGFC_THRESH <- 0.25
DE_PADJ_CUTOFF  <- 0.05

# GO thresholds
GO_PVALUE_CUTOFF <- 0.05
GO_QVALUE_CUTOFF <- 0.20
GO_SHOW_TOP      <- 15

# If you want to pull genes from a specific GO term description
GO_TERM_OF_INTEREST <- "inflammatory response"  # case-insensitive substring match
TOP_N_UP   <- 10
TOP_N_DOWN <- 10

# Keyword-filtered GO plot (immune/bacterial focus)
RUN_KEYWORD_GO_PLOT <- TRUE

# Optional: "variable gene GO" within a single condition (not DE-based)
# This is not condition-specific enrichment. Treat as exploratory.
RUN_VARIABLE_GENE_GO <- FALSE
VARIABLE_GENE_GROUP  <- "Benzene Male"  # only used if RUN_VARIABLE_GENE_GO = TRUE

# Optional: stimulus-linked expression panels for AMs
RUN_STIMULUS_PANELS <- TRUE
AM_CLUSTER_ID       <- 5               # AM cluster idents (you used 5)
AM_GROUP_LEVELS     <- c("Control Female", "Benzene Female", "Control Male", "Benzene Male")

# -----------------------------
# Package setup (middle-ground)
# -----------------------------
AUTO_INSTALL <- FALSE  # set TRUE if you want this script to install missing packages

cran_pkgs <- c("Seurat", "dplyr", "ggplot2", "stringr", "tidyr", "tibble", "forcats")
bioc_pkgs <- c("clusterProfiler", "org.Mm.eg.db")

install_if_missing <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (!AUTO_INSTALL) stop("Missing package: ", p, ". Install it and re-run.")
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else {
        install.packages(p)
      }
    }
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(forcats)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# -----------------------------
# Output directories
# -----------------------------
out_root <- file.path("outputs", paste0("go_cluster", CLUSTER_ID))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

out_tables <- file.path(out_root, "tables")
out_plots  <- file.path(out_root, "plots")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plots,  recursive = TRUE, showWarnings = FALSE)

sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# -----------------------------
# Guardrails
# -----------------------------
if (!exists("seurat_neo_lung")) stop("Object not found: seurat_neo_lung")
if (!GROUP_COL %in% colnames(seurat_neo_lung@meta.data)) stop("Missing metadata column: ", GROUP_COL)

# Ensure Idents are clusters
if ("seurat_clusters" %in% colnames(seurat_neo_lung@meta.data)) {
  Idents(seurat_neo_lung) <- seurat_neo_lung$seurat_clusters
}

# -----------------------------
# Helpers
# -----------------------------
parse_gene_ratio <- function(x) {
  parts <- str_split(as.character(x), "/", simplify = TRUE)
  if (ncol(parts) != 2) return(rep(NA_real_, length(x)))
  as.numeric(parts[, 1]) / as.numeric(parts[, 2])
}

make_go_dotplot <- function(go_obj, title, filename_base, show_top = 15) {
  if (is.null(go_obj) || nrow(as.data.frame(go_obj)) == 0) return(invisible(NULL))
  p <- dotplot(go_obj, showCategory = show_top) + ggtitle(title)
  ggsave(file.path(out_plots, paste0(filename_base, "_dotplot.png")),
         p, width = 10, height = 6, dpi = 300)
  p
}

run_de <- function(obj_cluster, ident1, ident2) {
  # group-by DE within this cluster subset
  Idents(obj_cluster) <- obj_cluster[[GROUP_COL]][, 1]

  res <- FindMarkers(
    obj_cluster,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "wilcox",
    min.pct = DE_MIN_PCT,
    logfc.threshold = DE_LOGFC_THRESH
  )

  # Save full table
  full_name <- paste0("DE_full_", sanitize(ident1), "_vs_", sanitize(ident2), "_cluster", CLUSTER_ID, ".csv")
  write.csv(res, file.path(out_tables, full_name), row.names = TRUE)

  # Filter by adjusted p-value
  if (!"p_val_adj" %in% colnames(res)) stop("FindMarkers output missing p_val_adj.")
  sig <- res[res$p_val_adj < DE_PADJ_CUTOFF, , drop = FALSE]

  sig_name <- paste0("DE_sig_padj", DE_PADJ_CUTOFF, "_", sanitize(ident1), "_vs_", sanitize(ident2),
                     "_cluster", CLUSTER_ID, ".csv")
  write.csv(sig, file.path(out_tables, sig_name), row.names = TRUE)

  list(full = res, sig = sig)
}

run_go <- function(sig_de_tbl, background_symbols, label) {
  if (nrow(sig_de_tbl) == 0) {
    message("No significant DE genes for ", label, " at padj < ", DE_PADJ_CUTOFF)
    return(NULL)
  }

  sig_symbols <- rownames(sig_de_tbl)

  entrez_sig <- bitr(
    sig_symbols,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  ) %>% distinct(ENTREZID, .keep_all = TRUE)

  if (nrow(entrez_sig) == 0) {
    message("No Entrez IDs mapped for significant genes: ", label)
    return(NULL)
  }

  entrez_bg <- bitr(
    unique(background_symbols),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  ) %>% distinct(ENTREZID, .keep_all = TRUE)

  go <- enrichGO(
    gene          = entrez_sig$ENTREZID,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = GO_PVALUE_CUTOFF,
    qvalueCutoff  = GO_QVALUE_CUTOFF,
    universe      = entrez_bg$ENTREZID
  )

  go_tbl <- as.data.frame(go)
  out_csv <- paste0("GO_BP_", label, "_cluster", CLUSTER_ID, ".csv")
  write.csv(go_tbl, file.path(out_tables, out_csv), row.names = FALSE)

  go
}

extract_genes_for_go_term <- function(go_obj, sig_de_tbl, go_term_substring, label) {
  if (is.null(go_obj) || nrow(as.data.frame(go_obj)) == 0) {
    message("No GO results available for gene extraction: ", label)
    return(invisible(NULL))
  }

  go_tbl <- as.data.frame(go_obj) %>%
    mutate(Description_lower = str_to_lower(Description)) %>%
    filter(str_detect(Description_lower, str_to_lower(go_term_substring))) %>%
    arrange(p.adjust)

  if (nrow(go_tbl) == 0) {
    message("No GO term matching '", go_term_substring, "' for ", label)
    return(invisible(NULL))
  }

  # Take the top matching term
  picked <- go_tbl[1, , drop = FALSE]
  gene_entrez <- unique(unlist(str_split(picked$geneID, "/")))

  sym_in_term <- bitr(
    gene_entrez,
    fromType = "ENTREZID",
    toType   = "SYMBOL",
    OrgDb    = org.Mm.eg.db
  ) %>% distinct(SYMBOL, .keep_all = TRUE)

  if (nrow(sym_in_term) == 0) {
    message("No symbols mapped for GO term genes: ", label)
    return(invisible(NULL))
  }

  de_in_term <- sig_de_tbl %>%
    rownames_to_column(var = "SYMBOL") %>%
    filter(SYMBOL %in% sym_in_term$SYMBOL)

  if (nrow(de_in_term) == 0) {
    message("No DE genes overlap with GO term genes: ", label)
    return(invisible(NULL))
  }

  # Pick top up and top down by avg_log2FC
  if (!"avg_log2FC" %in% colnames(de_in_term)) stop("Missing avg_log2FC in DE table: ", label)

  top_up <- de_in_term %>%
    filter(!is.na(avg_log2FC)) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = TOP_N_UP)

  top_down <- de_in_term %>%
    filter(!is.na(avg_log2FC)) %>%
    arrange(avg_log2FC) %>%
    slice_head(n = TOP_N_DOWN)

  top_both <- bind_rows(top_up, top_down) %>%
    mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down"))

  out_csv <- paste0("GOtermGenes_", sanitize(go_term_substring), "_", label, "_cluster", CLUSTER_ID, ".csv")
  write.csv(top_both, file.path(out_tables, out_csv), row.names = FALSE)

  p <- ggplot(top_both, aes(x = reorder(SYMBOL, avg_log2FC), y = avg_log2FC, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("Top DE genes in GO term match: '", go_term_substring, "' (", label, ")"),
      x = NULL,
      y = "avg_log2FC (Benzene - Control)"
    ) +
    theme(legend.position = "none")

  ggsave(file.path(out_plots, paste0("GOtermGenes_", sanitize(go_term_substring), "_", label, "_cluster", CLUSTER_ID, ".png")),
         p, width = 8, height = 6, dpi = 300)

  p
}

plot_keyword_go <- function(go_obj, title, filename_base) {
  if (is.null(go_obj) || nrow(as.data.frame(go_obj)) == 0) return(invisible(NULL))

  keywords <- c(
    "immune", "immunology", "inflammation", "interferon", "cytokine",
    "macrophage", "neutrophil", "monocyte", "dendritic", "antigen",
    "chemokine", "interleukin", "TNF", "NF-kB", "NF-κB", "MAPK",
    "infection", "viral", "virus", "bacterial", "bacterium", "antimicrobial",
    "lipopolysaccharide", "LPS", "peptidoglycan", "flagellin", "TLR", "NOD"
  )

  pat <- regex(paste(keywords, collapse = "|"), ignore_case = TRUE)

  tbl <- as.data.frame(go_obj)
  if (nrow(tbl) == 0) return(invisible(NULL))

  tbl_f <- tbl %>%
    filter(str_detect(Description, pat)) %>%
    mutate(GeneRatio_num = parse_gene_ratio(GeneRatio)) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10)

  if (nrow(tbl_f) == 0) {
    message("No GO terms matched keyword filter for: ", title)
    return(invisible(NULL))
  }

  tbl_f$Description <- factor(tbl_f$Description, levels = rev(tbl_f$Description))

  p <- ggplot(tbl_f, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
    scale_size_continuous(range = c(3, 10)) +
    theme_bw(base_size = 13) +
    labs(
      title = title,
      x = "Gene Ratio",
      y = NULL,
      size = "Gene Count"
    )

  ggsave(file.path(out_plots, paste0(filename_base, "_keyword_GO.png")),
         p, width = 9, height = 6, dpi = 300)

  p
}

# -----------------------------
# 1) Subset the target cluster once
# -----------------------------
obj_cluster <- subset(seurat_neo_lung, idents = CLUSTER_ID)
message("Cluster ", CLUSTER_ID, " cells: ", ncol(obj_cluster))

# Background gene set for GO (genes detected in this cluster)
bg_symbols <- rownames(obj_cluster[["RNA"]])

# -----------------------------
# 2) Run DE and GO for each comparison
# -----------------------------
for (nm in names(COMPARISONS)) {
  ident1 <- COMPARISONS[[nm]][1]
  ident2 <- COMPARISONS[[nm]][2]
  label  <- paste0(nm, "_", sanitize(ident1), "_vs_", sanitize(ident2))

  message("Running DE: ", ident1, " vs ", ident2, " within cluster ", CLUSTER_ID)
  de <- run_de(obj_cluster, ident1, ident2)

  message("Running GO: ", label)
  go <- run_go(de$sig, bg_symbols, label)

  # Dotplot of GO
  make_go_dotplot(
    go,
    title = paste0("GO BP: ", ident1, " vs ", ident2, " (cluster ", CLUSTER_ID, ")"),
    filename_base = paste0("GO_", label),
    show_top = GO_SHOW_TOP
  )

  # Extract genes from a GO term substring and plot
  extract_genes_for_go_term(
    go_obj = go,
    sig_de_tbl = de$sig,
    go_term_substring = GO_TERM_OF_INTEREST,
    label = label
  )

  # Keyword-filtered GO plot (immune/bacterial)
  if (RUN_KEYWORD_GO_PLOT) {
    plot_keyword_go(
      go_obj = go,
      title  = paste0("Keyword-filtered GO (immune/bacterial): ", ident1, " vs ", ident2,
                      " (cluster ", CLUSTER_ID, ")"),
      filename_base = paste0("GO_", label)
    )
  }
}

# -----------------------------
# 3) Optional: "variable gene GO" within one condition (not DE-based)
# -----------------------------
if (RUN_VARIABLE_GENE_GO) {
  message("Running variable-gene GO for group: ", VARIABLE_GENE_GROUP, " (exploratory)")

  obj_one <- subset(obj_cluster, subset = .data[[GROUP_COL]] == VARIABLE_GENE_GROUP)

  if (ncol(obj_one) >= 50) {
    DefaultAssay(obj_one) <- "RNA"
    obj_one <- NormalizeData(obj_one, verbose = FALSE)
    obj_one <- FindVariableFeatures(obj_one, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

    var_genes <- VariableFeatures(obj_one)

    entrez_var <- bitr(var_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>%
      distinct(ENTREZID, .keep_all = TRUE)

    entrez_bg <- bitr(unique(bg_symbols), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>%
      distinct(ENTREZID, .keep_all = TRUE)

    go_var <- enrichGO(
      gene          = entrez_var$ENTREZID,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = GO_PVALUE_CUTOFF,
      qvalueCutoff  = GO_QVALUE_CUTOFF,
      universe      = entrez_bg$ENTREZID
    )

    write.csv(as.data.frame(go_var),
              file.path(out_tables, paste0("GO_BP_variableGenes_", sanitize(VARIABLE_GENE_GROUP), "_cluster", CLUSTER_ID, ".csv")),
              row.names = FALSE)

    make_go_dotplot(
      go_var,
      title = paste0("GO BP on variable genes (exploratory): ", VARIABLE_GENE_GROUP, " (cluster ", CLUSTER_ID, ")"),
      filename_base = paste0("GO_variableGenes_", sanitize(VARIABLE_GENE_GROUP)),
      show_top = 20
    )
  } else {
    message("Skipped variable-gene GO: too few cells in ", VARIABLE_GENE_GROUP)
  }
}

# -----------------------------
# 4) Optional: stimulus-linked gene expression panels in AMs (cluster 5)
# -----------------------------
if (RUN_STIMULUS_PANELS) {
  stim_root <- file.path("outputs", paste0("stimulus_panels_cluster", AM_CLUSTER_ID))
  dir.create(stim_root, recursive = TRUE, showWarnings = FALSE)

  genes_by_stimulus <- list(
    bacterial   = c("Tlr4", "Nfkb1", "Cd14", "Myd88", "Tnf"),
    viral       = c("Tlr7", "Irf7", "Isg15", "Mx1", "Ifit1"),
    resolution  = c("Tnfaip3", "Il10", "Socs3", "Nlrp12", "Irak3")
  )

  gene_list <- unique(unlist(genes_by_stimulus))

  # Subset AMs
  if ("seurat_clusters" %in% colnames(seurat_neo_lung@meta.data)) {
    Idents(seurat_neo_lung) <- seurat_neo_lung$seurat_clusters
  }
  am_obj <- subset(seurat_neo_lung, idents = AM_CLUSTER_ID)

  DefaultAssay(am_obj) <- "RNA"
  Idents(am_obj) <- am_obj[[GROUP_COL]][, 1]

  # Drop missing genes safely
  gene_list_use <- intersect(gene_list, rownames(am_obj))
  if (length(gene_list_use) == 0) stop("None of the stimulus genes are present in the object.")

  stimulus_map <- stack(genes_by_stimulus)
  colnames(stimulus_map) <- c("Gene", "Stimulus")
  stimulus_map <- stimulus_map %>% filter(Gene %in% gene_list_use)

  # Fetch expression for AMs only
  am_data <- FetchData(am_obj, vars = c(gene_list_use, GROUP_COL))
  am_data <- as.data.frame(am_data)

  am_long <- tidyr::pivot_longer(
    data = am_data,
    cols = all_of(gene_list_use),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
    left_join(stimulus_map, by = "Gene") %>%
    mutate(
      !!GROUP_COL := factor(.data[[GROUP_COL]], levels = AM_GROUP_LEVELS),
      Stimulus = factor(Stimulus, levels = c("viral", "bacterial", "resolution"))
    )

  # Plot
  p_stim <- ggplot(am_long, aes(x = .data[[GROUP_COL]], y = Expression, fill = Stimulus)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 0.4, alpha = 0.35) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c("viral" = "firebrick", "bacterial" = "steelblue", "resolution" = "darkgreen")) +
    theme_bw(base_size = 12) +
    labs(
      title = "Stimulus-linked gene expression in AMs",
      x = NULL,
      y = "Normalized expression (RNA data slot)"
    ) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      strip.text   = element_text(face = "bold"),
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  ggsave(file.path(stim_root, "Stimulus_linked_genes_AMs_boxplot.png"),
         p_stim, width = 14, height = 7, dpi = 300)
}

message("Done. GO outputs: ", normalizePath(out_root, winslash = "/"))
