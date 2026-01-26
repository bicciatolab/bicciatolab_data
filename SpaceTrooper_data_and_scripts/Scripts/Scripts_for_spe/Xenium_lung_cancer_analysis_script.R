# Xenium_lung_cancer_analysis_script.R
# R 4.5.1, Bioconductor 3.22, SpaceTrooper>=1.1.3

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("SpaceTrooper", ref="devel") 

library(SpaceTrooper)
library(dplyr)

# dirname pointing to directory containing count matrix, fov position,
# polygons and metadata files as downloaded from:
# https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard

dirname <- "Xenium_lung_cancer"
samplename <- "Xenium_lung_cancer"

spe <- readXeniumSPE(dirName=dirname, sampleName=samplename, type="HDF5", addFOVs=TRUE)

spe <- readAndAddPolygonsToSPE(spe, boundariesType="parquet")

spe <- spatialPerCellQC(spe, rmZeros=FALSE)

# QS computation discards 0 count cells by default.
# For visualization purposes, we add the computed QS to the full SPE object
# where cells with 0 counts are assigned a QS equal to NA.

temp_spe <- computeQCScore(spe, verbose = FALSE)

temp_df <- data.frame("QC_score" = temp_spe$QC_score, "cell_id" = temp_spe$cell_id)
spe_df <- data.frame("cell_id" = spe$cell_id)
join_df <- left_join(spe_df, temp_df, by = "cell_id")
spe$QC_score <- join_df$QC_score

# Add cell types from metadata provided at
# https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/Spe_metadata:

meta_df <- readRDS("Xenium_lung_cancer_metadata.rds")
spe$InSituType_Simple <- as.factor(meta_df[match(spe$cell_id, meta_df$cell_id),]$InSituType_Simple)
