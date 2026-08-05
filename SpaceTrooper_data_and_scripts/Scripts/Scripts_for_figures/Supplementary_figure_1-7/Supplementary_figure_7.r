# Clean, minimal, and reproducible script to reproduce supplementary figure 7
# To reproduce the analysis be sure to have the required packages installed

library(ggplot2)
library(rlang)
library(patchwork)

# Specify spatial transcriptomics datasets to load metadata files available at: 
# bicciatolab_data/SpaceTrooper_data_and_scripts/Spe_metadata"

metadata_dir <- "bicciatolab_data/SpaceTrooper_data_and_scripts/Spe_metadata"

### define the axis theme ###
light_theme <- theme(panel.background = element_rect(fill="white", color=NA),
                     plot.background = element_rect(fill="white", color=NA),
                     axis.ticks = element_line(color = "grey80", linewidth = 0.2),
                     axis.title.y = element_text(color = "black"),
                     axis.title.x = element_text(color = "black"),
                     axis.line = element_blank(),
                     legend.background = element_rect(fill="white", color=NA),
                     legend.position = "none",
                     legend.title = element_text(color = "black"),
                     legend.text = element_text(color="black"),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_line(linewidth = 0.1, colour = "grey80"),
                     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)) 

# Supplementary figure 7

meta_ds_names <- c("CosMx_rna_pancreas", 
                   "Xenium_lung_cancer", 
                   "Xenium_breast_DCIS_1", 
                   "Xenium_breast_DCIS_2", 
                   "Xenium_breast_DBKERO", 
                   "Xenium_Prime_OC", 
                   "MERFISH_mouse_liver", 
                   "MERSCOPE_Ultra_mouse_embryo", 
                   "CosMx_protein_tonsil")

meta_paths <- as.list(file.path(metadata_dir, paste0(meta_ds_names, "_metadata.rds")))
names(meta_paths) <- meta_ds_names

plot_list <- lapply(meta_paths, function(x){
    metadata <- readRDS(x)
    var.to.plot <- "QScore"
    xlabel <- "Quality score"
    thr <- median(metadata[,var.to.plot], na.rm = TRUE)- 3*mad(metadata[,var.to.plot], na.rm = TRUE)
    pp <- ggplot(metadata, aes(!!sym(var.to.plot))) + 
        geom_histogram(aes(fill = after_stat(x)),
                    binwidth = 1/50,
                    boundary = 0,
                    closed = "left",
                    color = "grey40",
                    linewidth = 0.05) + 
        scale_fill_viridis_c(option = "plasma", 
                            limits = c(0, 1),
                            breaks = seq(0, 1, by = 0.1)) +
        scale_x_continuous(limits = c(0, 1),
                        breaks = seq(0, 1, by = 0.1),
                        expand = expansion(mult = c(0, 0))) +
        scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
                        expand = expansion(mult = c(0, 0.05))) +
        geom_vline(xintercept = thr, col="black", linewidth = 0.4, linetype = "dashed") +
        ylab("Number of cells") + 
        xlab(xlabel) + 
        ggtitle(paste0(names(meta_paths)[which(meta_paths == x)])) +
        theme_minimal() + 
        light_theme +
        theme(plot.margin = unit(c(0,20,0,0), "pt"),
            aspect.ratio = 1)
    pp    
})

# --- Build rows ---
row1 <- plot_list[[1]] | plot_list[[2]] | plot_list[[3]]
row2 <- plot_list[[4]] | plot_list[[5]] | plot_list[[6]]
row3 <- plot_list[[7]] | plot_list[[8]] | plot_list[[9]]
# --- Add equal spacing between rows ---

final_plot <- (row1 / plot_spacer() /
                 row2 / plot_spacer() /
                 row3) +
  plot_layout(
    heights = c(1, 0.01, 1, 0.01, 1),  # smaller spacing between rows
    widths = c(1, 1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "."
  ) &
  theme(plot.tag = element_text(size = 14, face = "bold"))

### save the final figure ###
pdf(file.path("SuppFigure_7.pdf"), 
    width = 14, 
    height = 14,
    bg = "transparent")
final_plot
dev.off()

# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)

# Matrix products: default
#   LAPACK version 3.12.1

# locale:
# [1] LC_COLLATE=English_United States.utf8 
# [2] LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    

# time zone: Europe/Rome
# tzcode source: internal

# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     

# other attached packages:
#  [1] patchwork_1.3.2             rlang_1.1.7                
#  [3] ggplot2_4.0.1               SpaceTrooper_1.1.9         
#  [5] testthat_3.3.2              SpatialExperiment_1.20.0   
#  [7] SingleCellExperiment_1.32.0 SummarizedExperiment_1.40.0
#  [9] Biobase_2.70.0              GenomicRanges_1.62.1       
# [11] Seqinfo_1.0.0               IRanges_2.44.0             
# [13] S4Vectors_0.48.0            BiocGenerics_0.56.0        
# [15] generics_0.1.4              MatrixGenerics_1.22.0      
# [17] matrixStats_1.5.0          

# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3        rstudioapi_0.17.1        
#   [3] jsonlite_2.0.0            shape_1.4.6.1            
#   [5] magrittr_2.0.4            ggbeeswarm_0.7.3         
#   [7] magick_2.9.0              farver_2.1.2             
#   [9] fs_1.6.6                  vctrs_0.7.1              
#  [11] memoise_2.0.1             DelayedMatrixStats_1.32.0
#  [13] rstatix_0.7.3             S4Arrays_1.10.1          
#  [15] usethis_3.2.1             BiocNeighbors_2.4.0      
#  [17] broom_1.0.11              Rhdf5lib_1.32.0          
#  [19] SparseArray_1.10.8        Formula_1.2-5            
#  [21] rhdf5_2.54.1              KernSmooth_2.23-26       
#  [23] desc_1.4.3                cachem_1.1.0             
#  [25] lifecycle_1.0.5           iterators_1.0.14         
#  [27] pkgconfig_2.0.3           rsvd_1.0.5               
#  [29] Matrix_1.7-4              R6_2.6.1                 
#  [31] fastmap_1.2.0             rprojroot_2.1.1          
#  [33] scater_1.38.0             dqrng_0.4.1              
#  [35] irlba_2.3.5.1             pkgload_1.4.1            
#  [37] ggpubr_0.6.2              beachmat_2.26.0          
#  [39] SpatialExperimentIO_1.2.0 abind_1.4-8              
#  [41] compiler_4.5.1            proxy_0.4-29             
#  [43] remotes_2.5.0             bit64_4.6.0-1            
#  [45] withr_3.0.2               S7_0.2.1                 
#  [47] backports_1.5.0           BiocParallel_1.44.0      
#  [49] carData_3.0-5             viridis_0.6.5            
#  [51] DBI_1.2.3                 pkgbuild_1.4.8           
#  [53] HDF5Array_1.38.0          R.utils_2.13.0           
#  [55] ggsignif_0.6.4            DelayedArray_0.36.0      
#  [57] sessioninfo_1.2.3         rjson_0.2.23             
#  [59] classInt_0.4-11           tools_4.5.1              
#  [61] units_1.0-0               vipor_0.4.7              
#  [63] otel_0.2.0                beeswarm_0.4.0           
#  [65] R.oo_1.27.1               glue_1.8.0               
#  [67] h5mread_1.2.1             rhdf5filters_1.22.0      
#  [69] grid_4.5.1                sf_1.0-24                
#  [71] gtable_0.3.6              R.methodsS3_1.8.2        
#  [73] class_7.3-23              tidyr_1.3.2              
#  [75] data.table_1.18.0         BiocSingular_1.26.1      
#  [77] ScaledMatrix_1.18.0       car_3.1-3                
#  [79] XVector_0.50.0            ggrepel_0.9.6            
#  [81] foreach_1.5.2             pillar_1.11.1            
#  [83] limma_3.66.0              robustbase_0.99-6        
#  [85] splines_4.5.1             dplyr_1.1.4              
#  [87] lattice_0.22-7            survival_3.8-3           
#  [89] bit_4.6.0                 tidyselect_1.2.1         
#  [91] locfit_1.5-9.12           scuttle_1.20.0           
#  [93] sfheaders_0.4.5           gridExtra_2.3            
#  [95] edgeR_4.8.2               statmod_1.5.1            
#  [97] devtools_2.4.6            DropletUtils_1.30.0      
#  [99] brio_1.1.5                DEoptimR_1.1-4           
# [101] codetools_0.2-20          tibble_3.3.1             
# [103] cli_3.6.5                 arrow_23.0.0             
# [105] dichromat_2.0-0.1         Rcpp_1.1.1               
# [107] parallel_4.5.1            ellipsis_0.3.2           
# [109] assertthat_0.2.1          sparseMatrixStats_1.22.0 
# [111] glmnet_4.1-10             viridisLite_0.4.2        
# [113] scales_1.4.0              e1071_1.7-17             
# [115] purrr_1.2.1               cowplot_1.2.0            
