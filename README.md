# Seurat Overview Pipeline (QC → Clustering → Annotation → Plots/Markers)

This repository contains an end-to-end overview analysis workflow for a Seurat single-cell RNA-seq object. It performs basic QC filtering, normalization, dimensionality reduction, clustering, cluster-to-cell-type mapping, and generates a set of commonly used plots and marker outputs.

This workflow is intentionally written as a single, readable script to serve as a reproducible template for rapid dataset triage and a teaching-friendly example of a standard scRNA-seq analysis structure.

## What this pipeline does

Input:
- A Seurat .rds object (RNA assay expected)

Workflow steps:
1. Load Seurat object
2. Compute percent mitochondrial content if missing
3. Visualize percent mitochondrial content by group before filtering
4. Filter cells by nFeature_RNA thresholds and percent.mito threshold
5. Normalize data and identify variable features
6. Scale data and run PCA
7. Generate elbow plot for dimensionality selection
8. Run UMAP, find neighbors, and cluster cells
9. Summarize cell counts by cluster and experimental group
10. Assign biological cell-type labels using a cluster-to-cell-type mapping
11. Generate representative plots including UMAPs, feature plots, violin plots, ridge plots, and dot plots
12. Identify marker genes using FindAllMarkers
13. Save tables, figures, and sessionInfo for reproducibility

## Outputs

Running the pipeline creates the following structure:

outputs/plots/
- QC plots and UMAP visualizations
- FeaturePlot panels
- DotPlot marker panels
- Selected violin and ridge plots

outputs/tables/
- Cell count summaries
- Marker gene tables
- Significant marker subsets

outputs/sessionInfo.txt
- R session and package versions used for the analysis

## Requirements

R packages required:
- Seurat
- ggplot2
- dplyr
- stringr
- tidyr
- tibble

Example installation (run in R console):
install.packages(c("Seurat", "ggplot2", "dplyr", "stringr", "tidyr", "tibble"))

If Seurat v5 from GitHub is required (run once):
install.packages("remotes")
remotes::install_github("satijalab/seurat", ref = "seurat5")

## How to run

1. Open the script:
seurat_neo_lung_overview_pipeline.R

2. Edit user-defined parameters at the top of the script:
- input_rds (path to your .rds file)
- QC thresholds (manual or dynamic)
- PCA dimensions and clustering resolution

3. Run the pipeline in R:
source("seurat_neo_lung_overview_pipeline.R")

## Notes on portability

- The script avoids machine-specific paths beyond input_rds
- All results are written to a local outputs directory
- Large .rds objects should not be committed to the repository

### 02_ucell_go_scoring.R

**Purpose**
- Scores a GO term gene set per cell type using ssGSEA (`GSVA`) and `UCell`
- Summarizes median UCell scores by cell type, sex, and treatment
- Tests Benzene vs Control within sex (Wilcoxon) and annotates significance
- Generates:
  - Dot plot of median UCell scores by cell type (faceted by sex)
  - Heatmap of (Benzene - Control) median UCell deltas with significance stars
  - Optional GO gene-level plots for selected cell types (avg expression + sex-stratified log2FC)

**Inputs and assumptions**
- `seurat_neo_lung` is already loaded in R
- Metadata: `Treatment_by_Sex`, `seurat_clusters` exist
- Mouse GO annotation via `org.Mm.eg.db`
- GO term set by `GO_TERM` (example: `GO:0034340`)

**Outputs**
- `UCell_GO_<GO_TERM>_dotplot.tiff`
- `UCell_GO_<GO_TERM>_heatmap.tiff`
- Optional: per-celltype CSV + GO gene heatmaps for selected cell types

**Run**
- `source("02_ucell_go_scoring.R")`


### 03_cellchat_am_signaling.R

**Purpose**
- Runs CellChat (mouse DB) per condition (`Treatment_by_Sex`)
- Infers ligand–receptor communication networks and summarizes AM (and Prolif AM) outgoing signaling
- Includes example visualizations (chord, signaling role heatmap, bubble) and AM→target strength summaries
- Builds “top targets / top pathways” overview dot plot and faceted heatmap

**Inputs and assumptions**
- `seurat_neo_lung` is already loaded in R
- Metadata: `Treatment_by_Sex`, `cell_type` exist
- Uses `CellChatDB.mouse` (subset to “Secreted Signaling” by default)
- Parallel plan set via `future::plan(multisession, workers = 4)`

**Outputs**
- Plots are printed; add `ggsave()` / `pdf()` calls if you want to persist them
- Recommended addition: save objects for reuse
  - `saveRDS(cellchat_objects, "cellchat_objects_by_group.rds")`

**Run**
- `source("03_cellchat_am_signaling.R")`


### 04_progreny.dorothea.R

**Purpose**
- Fibroblast-focused “WOW” slides for pathway and TF activity:
  - PROGENy pathway activity (per cell)
  - DoRothEA TF activity via VIPER (per cell)
- Computes median effect sizes (median Benzene minus median Control) within each sex
- Generates UMAP feature panels (mako palette) for the top-shifted pathways and TFs per sex
- Writes effect-size tables for ranking

**Inputs and assumptions**
- `seurat_neo_lung` is already loaded in R
- Fibroblast cluster identity set by `fibro_cluster_ident` (default `13`)
- Metadata column `Treatment_by_Sex` contains:
  - `Control Female`, `Benzene Female`, `Control Male`, `Benzene Male`
- Uses normalized expression from `assay_key` (default `"RNA"`)
- If UMAP is missing in the subset, the script computes a quick PCA + UMAP on fibroblasts
- With low biological replication, treat outputs as hypothesis-generating

**Outputs**
- Creates output folder:
  - `Fibroblast_PROGENy_DoRothEA_cluster_<clusterID>/`
- Saves within that folder:
  - `PROGENy_effect_sizes_<Sex>_*.csv`
  - `DoRothEA_effect_sizes_<Sex>_*.csv`
  - `SLIDE1_PROGENy_activity_<Sex>_mako.png`
  - `SLIDE2_DoRothEA_TF_activity_<Sex>_mako.png`
  - `fibroblasts_with_PROGENy_DoRothEA.rds`

**Run**
- `source("04_progreny.dorothea.R")`


### 05_scfea_flux_balance.R

**Purpose**
- Integrates scFEA outputs with Seurat for AM clusters (5 and 17):
  - Adds scFEA flux predictions as a `FLUX` assay
  - Adds scFEA balance predictions as a `BALANCE` assay
- Performs example group comparisons with `FindMarkers` on FLUX/BALANCE
- Generates example plots (UMAP in flux space, volcano/dot/bar plots)
- Includes an Ensembl→symbol annotation block for RNA DEGs using `biomaRt`

**Inputs and assumptions**
- `seurat_neo_lung` is already loaded in R
- scFEA output files exist (update paths as needed):
  - `scFEA_flux_results_2.csv`
  - `scFEA_balance_results_2.csv`
- Cell barcode harmonization is handled with `gsub()` / `RenameCells()`, but you should verify it matches your naming conventions

**Outputs**
- `Seurat_cluster5_17_geneExpr.csv` (RNA counts export)
- `BalanceTest_results_male.csv`
- `Seurat5DEG_results_male.csv`
- Additional plots printed unless you add `ggsave()`

**Run**
- `source("05_scfea_flux_balance.R")`


### 06_go_enrichment_and_stimulus_panels.R

**Purpose**
- Runs cluster-specific differential expression (via `FindMarkers`) by `Treatment_by_Sex`
- Performs GO Biological Process enrichment (clusterProfiler) and makes dotplots
- Lets you drill into a GO term (by Description keyword) to list and plot DEGs in that term
- Optional keyword-filtered GO plots (immune/infection/bacterial term filters)
- Includes a stimulus-linked expression panel example in alveolar macrophages (cluster 5)

**Inputs and assumptions**
- `seurat_neo_lung` is already loaded in R
- Metadata: `Treatment_by_Sex` exists
- Mouse annotation via `org.Mm.eg.db`
- Script uses clusters `13` and `5` in different sections; edit these to match your analysis targets

**Outputs**
- GO result tables saved as CSV (filenames defined in-script)
- GO dotplots and DEG bar plots (save manually or add `ggsave()` where you want)

**Run**
- `source("06_go_enrichment_and_stimulus_panels.R")`

## Intended use

This pipeline is designed to be transparent, quick to adapt, and reproducible. It is appropriate for exploratory analyses, method development, and trainee instruction.

Typical next steps for a production workflow would include splitting steps into multiple scripts, adding renv for dependency locking, and modularizing repeated plotting code.

## Contact

If you use or adapt this pipeline and have questions, please open an Issue or fork the repository.
