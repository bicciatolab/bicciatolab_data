# Xenium_DCIS_S2_analysis_script.R
# R 4.5.1, Bioconductor 3.22, SpaceTrooper>=1.1.8

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("SpaceTrooper", ref="devel") 

library(SpaceTrooper)
library(dplyr)

# dirname pointing to directory containing count matrix, fov position,
# polygons and metadata files as downloaded from:
# https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard

dirname <- "Xenium_DCIS_S2"
samplename <- "Xenium_DCIS_S2"

spe <- readXeniumSPE(dirName=dirname, sampleName=samplename, type="HDF5", addFOVs=FALSE)

spe <- readAndAddPolygonsToSPE(spe, boundariesType="parquet")

spe <- spatialPerCellQC(spe, rmZeros=FALSE)

# QS computation discards 0 count cells by default.
# For visualization purposes, we add the computed QS to the full SPE object
# where cells with 0 counts are assigned a QS equal to NA.

temp_spe <- computeQScore(spe, verbose=FALSE)

temp_df <- data.frame("QScore"=temp_spe$QScore, "cell_id"=temp_spe$cell_id)
spe_df <- data.frame("cell_id"=spe$cell_id)
join_df <- left_join(spe_df, temp_df, by="cell_id")
spe$QScore <- join_df$QScore

# Add cell types from metadata provided at
# https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/Spe_metadata:

meta_df <- readRDS("Xenium_DCIS_S2_metadata.rds")
spe$cell_type <- as.factor(meta_df[match(spe$cell_id, meta_df$cell_id),]$cell_type)
