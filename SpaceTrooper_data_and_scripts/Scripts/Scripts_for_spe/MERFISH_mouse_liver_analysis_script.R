# MERFISH_mouse_liver_analysis_script.R
# R 4.5.1, Bioconductor 3.22, SpaceTrooper>=1.1.3

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("SpaceTrooper", ref="devel") 

library(SpaceTrooper)
library(dplyr)

dirname <- "MERFISH_mouse_liver"
samplename <- "MERFISH_mouse_liver"

spe <- readMerfishSPE(dirName=dirname, sampleName=samplename, computeMissingMetrics=FALSE)

# Since polygons were provided initially only in HDF5 files per each FoV for MERFISH
# technology, with over 1000 FoVs, polygon extraction is time-consuming and non-trivial, as
# stated also in Voyager vignette. (https://pachterlab.github.io/voyager/articles/vig6_merfish.html)
# For this reason, polygons can be sourced from SFEData package on Bioconductor 
# and loaded inside the SPE object created with SpaceTrooper:

BiocManager::install("SpatialFeatureExperiment")
BiocManager::install("SFEData")

library(SFEData)
sfe_liver <- VizgenLiverData()

table(colnames(sfe_liver)==spe$cell_id)
#   TRUE 
# 395215

polygons <- cellSeg(sfe_liver)
polygons$cell_id <- spe$cell_id

polygons <- .renameGeometry(polygons, "geometry", "global", activate=TRUE)
spe$polygons <- polygons[match(spe$cell_id, polygons$cell_id),]

# As stated in the paper, QS is computed using cell volume for MERFISH technology.
# For the sake of simplicity, it is temporarily renamed as Area_um to be used in
# the following steps.

spe$Area_um <- spe$volume
spe$AspectRatio <- computeAspectRatioFromPolygons(spe$polygons)

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

meta_df <- readRDS("MERFISH_mouse_liver_metadata.rds")
spe$InSituType_Simple <- as.factor(meta_df[match(spe$cell_id, meta_df$cell_id),]$InSituType_Simple)
