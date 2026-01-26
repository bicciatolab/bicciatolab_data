# CosMx_Protein_assay_tonsil_analysis_script
# R 4.5.1, Bioconductor 3.22, SpaceTrooper>=1.1.3

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("SpaceTrooper", ref="devel") 

library(SpaceTrooper)

# dirname pointing to directory containing count matrix, fov position,
# polygons and metadata files as downloaded from
# https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/cosmx-human-tonsil-ffpe-protein-dataset/

dirname <- "S0"
samplename <- "CosMx_Protein_assay_tonsil"

spe <- readCosmxProteinSPE(dirName=dirname, sampleName=samplename)

# We found that FoV and cells' coordinates were not aligned, so we corrected
# this shift along the vertical axis corresponding to the length of a FoV (4256 px).

metadata(spe)$fov_positions$y_global_px <- metadata(spe)$fov_positions$y_global_px - 4256

spe <- readAndAddPolygonsToSPE(spe, boundariesType="csv")

# We noticed that AspectRatio provided in dataset metadata was not computed as
# defined by Bruker (Nanostring) definition of the metric (Width/Height), as was
# done for the other datasets that we analyzed. For this reason, we recomputed
# it according to their definition.

spe$AspectRatio <- spe$Width/spe$Height

spe <- spatialPerCellQC(spe)

spe <- computeQCScore(spe, verbose = FALSE)

# Add cell types from metadata provided at
# https://github.com/bicciatolab/bicciatolab_data/tree/main/SpaceTrooper_data_and_scripts/Spe_metadata:

meta_df <- readRDS("CosMx_protein_tonsil_metadata.rds")
spe$cell_type <- as.factor(meta_df[match(spe$cell_id, meta_df$cell_id),]$cell_type)

