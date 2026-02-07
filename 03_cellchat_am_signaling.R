# ------------------------------------------------------------
# CellChat: Condition-split communication inference and AM-focused summaries
#
# Assumes:
#   - seurat_neo_lung exists in your environment OR you provide input_rds
#   - Metadata includes Treatment_by_Sex and cell_type
#
# Outputs:
#   - outputs/cellchat/objects/*.rds (cached CellChat objects)
#   - outputs/cellchat/plots/*       (plots)
#   - outputs/cellchat/tables/*      (tables)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(future)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(pheatmap)
})

set.seed(42)

# -----------------------------
# USER PARAMETERS
# -----------------------------

# If seurat_neo_lung is not in your environment, set this:
input_rds <- "path/to/seurat_neo_lung_filtered.rds"

split_col <- "Treatment_by_Sex"
group_by  <- "cell_type"

# CellChat DB selection
db_species <- "mouse"
db_category <- "Secreted Signaling"   # or set to NULL for full DB

# Computation parameters
workers <- 4
min_cells <- 10

# AM source cell types
am_sources <- c("Alveolar Macrophages", "Proliferative Alveolar Macrophages")

# Targets for quick example plots
target_for_chord  <- "Neutrophils"
target_for_bubble <- "Neutrophils"
target_for_delta  <- "Inflammatory Monocytes"

# Output directories
out_dir    <- file.path("outputs", "cellchat")
obj_dir    <- file.path(out_dir, "objects")
plot_dir   <- file.path(out_dir, "plots")
table_dir  <- file.path(out_dir, "tables")
dir.create(obj_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

save_base_plot_pdf <- function(filename, width = 9, height = 7, expr) {
  pdf(file.path(plot_dir, filename), width = width, height = height)
  try(expr, silent = TRUE)
  dev.off()
}

# -----------------------------
# LOAD OBJECT IF NEEDED
# -----------------------------

if (!exists("seurat_neo_lung")) {
  seurat_neo_lung <- readRDS(input_rds)
}

# Guardrails
if (!split_col %in% colnames(seurat_neo_lung@meta.data)) stop("Missing metadata: ", split_col)
if (!group_by %in% colnames(seurat_neo_lung@meta.data))  stop("Missing metadata: ", group_by)

# Parallel
plan("multisession", workers = workers)

# If you hit future size errors on big objects, uncomment and raise:
# options(future.globals.maxSize = 8 * 1024^3)  # 8 GB

# -----------------------------
# DB setup
# -----------------------------

if (db_species == "mouse") {
  CellChatDB <- CellChatDB.mouse
} else {
  CellChatDB <- CellChatDB.human
}

if (!is.null(db_category)) {
  CellChatDB.use <- subsetDB(CellChatDB, search = db_category)
} else {
  CellChatDB.use <- CellChatDB
}

# -----------------------------
# Split object by condition
# -----------------------------

seurat_list <- SplitObject(seurat_neo_lung, split.by = split_col)

cat("Groups found:\n")
print(names(seurat_list))

# -----------------------------
# Create and process CellChat object
# -----------------------------

make_cellchat <- function(seurat_obj, group_by, db_use, min_cells = 10) {
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  meta <- seurat_obj@meta.data

  if (!group_by %in% colnames(meta)) stop("group_by column not found in metadata: ", group_by)

  cellchat <- createCellChat(object = data.input, meta = meta, group.by = group_by)
  cellchat@DB <- db_use

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)

  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  return(cellchat)
}

# -----------------------------
# Build or load cached CellChat objects
# -----------------------------

cellchat_objects <- list()

for (nm in names(seurat_list)) {
  cache_file <- file.path(obj_dir, paste0("cellchat_", sanitize(nm), ".rds"))

  if (file.exists(cache_file)) {
    message("Loading cached CellChat: ", nm)
    cellchat_objects[[nm]] <- readRDS(cache_file)
  } else {
    message("Computing CellChat: ", nm)
    cc <- make_cellchat(seurat_list[[nm]], group_by = group_by, db_use = CellChatDB.use, min_cells = min_cells)
    saveRDS(cc, cache_file)
    cellchat_objects[[nm]] <- cc
    message("Saved: ", cache_file)
  }
}

# Convenience handles (only if present)
get_cc <- function(name) {
  if (!name %in% names(cellchat_objects)) stop("Group not found: ", name)
  cellchat_objects[[name]]
}

# Example group objects
cellchat_benzene_female <- if ("Benzene Female" %in% names(cellchat_objects)) get_cc("Benzene Female") else NULL
cellchat_control_female <- if ("Control Female" %in% names(cellchat_objects)) get_cc("Control Female") else NULL
cellchat_benzene_male   <- if ("Benzene Male"   %in% names(cellchat_objects)) get_cc("Benzene Male")   else NULL
cellchat_control_male   <- if ("Control Male"   %in% names(cellchat_objects)) get_cc("Control Male")   else NULL

# -----------------------------
# Plot 1: Chord plot (example group)
# -----------------------------

if (!is.null(cellchat_benzene_female)) {
  save_base_plot_pdf(
    filename = paste0("Chord_Benzene_Female_AM_to_", sanitize(target_for_chord), ".pdf"),
    width = 10, height = 8,
    expr = netVisual_chord_gene(
      cellchat_benzene_female,
      sources.use = am_sources,
      targets.use = c(target_for_chord),
      lab.cex = 0.8,
      small.gap = 3
    )
  )
}

# -----------------------------
# Plot 2: Outgoing signaling heatmap (example group)
# -----------------------------

if (!is.null(cellchat_benzene_female)) {
  cellchat_benzene_female <- netAnalysis_computeCentrality(cellchat_benzene_female)

  save_base_plot_pdf(
    filename = "OutgoingRoleHeatmap_Benzene_Female.pdf",
    width = 10, height = 8,
    expr = netAnalysis_signalingRole_heatmap(cellchat_benzene_female, pattern = "outgoing")
  )
}

# -----------------------------
# AM outgoing ranking and matrix across groups
# -----------------------------

get_AM_outgoing <- function(cellchat_obj, am_names = am_sources) {
  comm_matrix <- cellchat_obj@net$weight
  if (!all(am_names %in% rownames(comm_matrix))) stop("AM names not found in net$weight rows")
  colSums(comm_matrix[am_names, , drop = FALSE])
}

available_groups <- names(cellchat_objects)

am_out_list <- lapply(available_groups, function(g) get_AM_outgoing(cellchat_objects[[g]]))
names(am_out_list) <- available_groups

am_out_matrix <- as.data.frame(am_out_list)
am_out_matrix <- am_out_matrix[sort(rownames(am_out_matrix)), , drop = FALSE]

# Remove AM self-target rows
self_labels <- am_sources
am_out_matrix_no_self <- am_out_matrix[!rownames(am_out_matrix) %in% self_labels, , drop = FALSE]

write.csv(
  am_out_matrix_no_self,
  file = file.path(table_dir, "AM_outgoing_strength_matrix_no_self.csv"),
  row.names = TRUE
)

# Heatmap of AM outgoing strength across groups
save_base_plot_pdf(
  filename = "Heatmap_AM_outgoing_strength_no_self.pdf",
  width = 9, height = 8,
  expr = pheatmap(
    am_out_matrix_no_self,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = "AM to Target Interaction Strength (No AM self-talk)",
    angle_col = 45,
    border_color = NA
  )
)

# -----------------------------
# Bubble plot (example group)
# -----------------------------

if (!is.null(cellchat_control_female)) {
  p_bubble <- netVisual_bubble(
    cellchat_control_female,
    sources.use = am_sources,
    targets.use = c(target_for_bubble),
    signaling = NULL
  )

  ggsave(
    filename = file.path(plot_dir, paste0("Bubble_Control_Female_AM_to_", sanitize(target_for_bubble), ".pdf")),
    plot = p_bubble,
    width = 11, height = 7
  )
}

# -----------------------------
# Delta signaling bar plot: Benzene - Control (Female), AM to one target
# -----------------------------

if (!is.null(cellchat_benzene_female) && !is.null(cellchat_control_female)) {

  df_benzene <- subsetCommunication(cellchat_benzene_female)
  df_control <- subsetCommunication(cellchat_control_female)

  benzene_AM_to_target <- df_benzene %>%
    filter(source %in% am_sources, target == target_for_delta)

  control_AM_to_target <- df_control %>%
    filter(source %in% am_sources, target == target_for_delta)

  benzene_summary <- benzene_AM_to_target %>%
    group_by(pathway_name) %>%
    summarise(total_strength = sum(prob), .groups = "drop")

  control_summary <- control_AM_to_target %>%
    group_by(pathway_name) %>%
    summarise(total_strength = sum(prob), .groups = "drop")

  merged <- full_join(benzene_summary, control_summary, by = "pathway_name",
                      suffix = c("_benzene", "_control"))
  merged[is.na(merged)] <- 0

  p_delta <- ggplot(
    merged,
    aes(x = reorder(pathway_name, total_strength_benzene),
        y = total_strength_benzene - total_strength_control)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = paste0("Delta signaling from AMs to ", target_for_delta, " Female (Benzene - Control)"),
      y = "Delta total interaction strength",
      x = "Pathway"
    ) +
    theme_minimal()

  ggsave(
    filename = file.path(plot_dir, paste0("Delta_AM_to_", sanitize(target_for_delta), "_Female_Benzene_minus_Control.pdf")),
    plot = p_delta,
    width = 10, height = 7
  )

  write.csv(
    merged,
    file = file.path(table_dir, paste0("Delta_AM_to_", sanitize(target_for_delta), "_Female_pathway_strengths.csv")),
    row.names = FALSE
  )
}

# -----------------------------
# Top targets and top pathways across all groups: dot plot and heatmap
# -----------------------------

# Define top 5 targets based on total AM outgoing across all groups
am_out_mat <- as.matrix(am_out_matrix_no_self)
top_targets <- names(sort(rowSums(am_out_mat), decreasing = TRUE))[1:min(5, nrow(am_out_mat))]
cat("Top targets:\n"); print(top_targets)

get_AM_pathway_comms <- function(cellchat_obj, group_name, am_sources, top_targets) {
  subsetCommunication(cellchat_obj) %>%
    filter(source %in% am_sources, target %in% top_targets) %>%
    group_by(pathway_name, target) %>%
    summarise(total_strength = sum(prob), .groups = "drop") %>%
    mutate(Group = group_name)
}

comm_all <- bind_rows(lapply(names(cellchat_objects), function(g) {
  get_AM_pathway_comms(cellchat_objects[[g]], g, am_sources, top_targets)
}))

top_pathways <- comm_all %>%
  group_by(pathway_name) %>%
  summarise(total = sum(total_strength, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 10) %>%
  pull(pathway_name)

comm_filtered <- comm_all %>%
  filter(pathway_name %in% top_pathways)

write.csv(
  comm_filtered,
  file = file.path(table_dir, "AM_to_topTargets_top10Pathways_by_group.csv"),
  row.names = FALSE
)

# Dot plot
p_dot <- ggplot(comm_filtered, aes(x = target, y = pathway_name, size = total_strength, color = Group)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(1, 10)) +
  labs(
    title = "Top 10 pathways: AM to top targets by group",
    x = "Target cell type",
    y = "Signaling pathway",
    size = "Communication strength",
    color = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(plot_dir, "DotPlot_AM_topTargets_top10Pathways_by_group.pdf"),
  plot = p_dot,
  width = 12, height = 6
)

# Heatmap (faceted)
p_heat <- ggplot(comm_filtered, aes(x = target, y = pathway_name, fill = total_strength)) +
  geom_tile(color = "white") +
  facet_wrap(~ Group) +
  labs(
    title = "Heatmap: AM to target signaling by pathway and group",
    x = "Target cell type",
    y = "Pathway",
    fill = "Strength"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(plot_dir, "Heatmap_AM_topTargets_top10Pathways_by_group.pdf"),
  plot = p_heat,
  width = 14, height = 7
)

message("Done. Outputs in: ", out_dir)
