# CosMx_Pancreas_WTx_analysis_script
# R 4.5.1, Bioconductor 3.22, SpaceTrooper>=1.1.3

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("SpaceTrooper", ref="devel") 

library(SpaceTrooper)

# dirname pointing to directory containing count matrix, fov position,
# polygons and metadata files as downloaded from
# https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-smi-human-pancreas-ffpe-dataset/

dirname <- "CosMx_Pancreas_WTx"
samplename <- "CosMx_Pancreas_WTx"

spe <- readCosmxSPE(dirName=dirname, sampleName=samplename)

# We found that FoV and cells' coordinates were not aligned, so we corrected
# this shift along the vertical axis corresponding to the length of a FoV (4256 px).

metadata(spe)$fov_positions$y_global_px <- metadata(spe)$fov_positions$y_global_px - 4256

spe <- readAndAddPolygonsToSPE(spe, boundariesType="csv")

spe <- spatialPerCellQC(spe)

spe <- computeQCScore(spe, verbose = FALSE)

# Add cell types from metadata provided at
# https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/Spe_metadata:

meta_df <- readRDS("CosMx_rna_pancreas_metadata.rds")
spe$cell_types <- as.factor(meta_df[match(spe$cell_id, meta_df$cell_id),]$cell_types)
