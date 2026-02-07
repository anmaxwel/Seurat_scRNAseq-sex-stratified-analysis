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

## Intended use

This pipeline is designed to be transparent, quick to adapt, and reproducible. It is appropriate for exploratory analyses, method development, and trainee instruction.

Typical next steps for a production workflow would include splitting steps into multiple scripts, adding renv for dependency locking, and modularizing repeated plotting code.

## Contact

If you use or adapt this pipeline and have questions, please open an Issue or fork the repository.
