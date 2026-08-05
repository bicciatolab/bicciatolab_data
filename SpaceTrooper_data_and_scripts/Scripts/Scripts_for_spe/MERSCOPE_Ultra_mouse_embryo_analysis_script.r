# MERSCOPE_Ultra_mouse_embryo_analysis_script.R
# R 4.5.1, Bioconductor 3.22, SpaceTrooper>=1.1.8

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("SpaceTrooper", ref="devel") 

library(SpaceTrooper)
library(dplyr)

dirname <- "MERSCOPE_Ultra_mouse_embryo/R3"
samplename <- "MERSCOPE_Ultra_mouse_embryo"

spe <- readMerfishSPE(dirName=dirname, sampleName=samplename, computeMissingMetrics=TRUE,
boundariesType="parquet", keepPolygons=TRUE)

#polygons were already loaded in the SPE object

spe <- spatialPerCellQC(spe, rmZeros=FALSE)

# QS computation discards 0 count cells by default.
# For visualization purposes, we add the computed QS to the full SPE object
# where cells with 0 counts are assigned a QS equal to NA.

# IMPORTANT: cell IDs provided by MERSCOPE technology are numeric and must
# be converted to character before running computeQScore to avoid errors.

spe$cell_id <- as.character(spe$cell_id)

temp_spe <- computeQScore(spe, verbose=FALSE)

temp_df <- data.frame("QScore"=temp_spe$QScore, "cell_id"=temp_spe$cell_id)
spe_df <- data.frame("cell_id"=spe$cell_id)
join_df <- left_join(spe_df, temp_df, by="cell_id")
spe$QScore <- join_df$QScore
