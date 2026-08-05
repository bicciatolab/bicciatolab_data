# Clean, minimal, and reproducible script to reproduce supplementary figure 25 
# To reproduce the analysis be sure to have the required packages installed

library(SpaceTrooper)
library(ggplot2)
library(robustbase)
library(patchwork)
library(ggrastr)

# Specify spatial transcriptomics datasets to load metadata files available at: 
# bicciatolab_data/SpaceTrooper_data_and_scripts/Spe_metadata"

metadata_dir <- "bicciatolab_data/SpaceTrooper_data_and_scripts/Spe_metadata"

datasets <- c("Xenium_breast_DCIS_1",
              "Xenium_breast_DCIS_2")

### define the axis theme ###
light_theme <- theme(panel.background = element_rect(fill="white", color=NA),
                     plot.background = element_rect(fill="white", color=NA),
                     title = element_blank(),
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

pp_list <- list()
for (k in 1:length(datasets)){  
  metadata <- readRDS(file.path(metadata_dir, paste0(datasets[k],"_metadata.rds")))
  var.to.plot <- "total"
    pp <- ggplot(metadata) +
    geom_histogram(
        aes(x = .data[[var.to.plot]]),
        bins = 50,
        fill = "#B3C2F2",
        color = "grey40",
        linewidth = 0.1
    ) +
    geom_vline(xintercept = median(metadata[[var.to.plot]], na.rm = TRUE),
                col = "black",
                linewidth = 0.8) +
    annotate("text",
            x = median(metadata[[var.to.plot]], na.rm = TRUE) - 50,
            y = 15000,
            label = round(median(metadata[[var.to.plot]], na.rm = TRUE), 2),
            color = "black",
            size = 3) +
    geom_vline(xintercept = adjbox(metadata[[var.to.plot]], plot = FALSE)$fence[1],
                col = "purple3",
                linewidth = 0.8) +
    annotate("text",
            x = adjbox(metadata[[var.to.plot]], plot = FALSE)$fence[1] - 5,
            y = 15000,
            label = round(adjbox(metadata[[var.to.plot]], plot = FALSE)$fence[1], 2),
            color = "purple3",
            size = 3) +
    scale_x_continuous(trans = "log10",
                       breaks = scales::breaks_log(n = 6),
                       limits = c(1, 10000)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                       expand = expansion(mult = c(0, 0.05)),
                       limits = c(0, 30000)) +
    ylab("Number of cells") + 
    xlab("Total probe counts") + 
    theme_minimal() + 
    light_theme +
    theme(plot.margin = unit(c(0,20,0,0), "pt"))
    pp_list[[k]] <- pp
}

row1 <- pp_list[[1]] + pp_list[[2]] 

final_plot <- row1 +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "."
  ) &
  theme(plot.tag = element_text(size = 14, face = "bold"))

pdf(file.path("SuppFigure_25A-B.pdf"), 
    width = 10, 
    height = 5,
    bg = "transparent")
final_plot
dev.off()

# Supplementary figure 25C-E

# Xenium DCIS-1 <-> Xenium DCIS-2 coefficient transfer

# dirname pointing to directory containing count matrix, fov position,
# polygons and metadata files as downloaded from:
# https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard
# S1-Bottom and S2_bottom

dirname1 <- "D:\\Xenium_v1_DCIS_S1-Bottom"
dirname2 <- "D:\\Xenium_v1_DCIS_S2-Bottom"

xen1 <- readXeniumSPE(dirname1)
xen2 <- readXeniumSPE(dirname2)

xen1 <- spatialPerCellQC(xen1)
xen2 <- spatialPerCellQC(xen2)

# no need to specify the formula, it's the default one
set.seed(1998)

xen1 <- computeQScore(
    spe=xen1,
    modelFormula=NULL,
    verbose=TRUE
)

# no log2Ctrl_total_ratio was used, not enough outliers

xen2 <- computeQScore(
    spe=xen2,
    modelFormula=NULL,
    verbose=TRUE
)

# no log2Ctrl_total_ratio was used, not enough outliers

xen1$QScore_native <- xen1$QScore
xen2$QScore_native <- xen2$QScore

xen1_model <- metadata(xen1)$QScore_model
xen2_model <- metadata(xen2)$QScore_model

xen1 <- .applyQScoreModel(
    spe=xen1,
    qsModel=xen2_model,
    scoreName="QScore_xen2_model"
)

xen2 <- .applyQScoreModel(
    spe=xen2,
    qsModel=xen1_model,
    scoreName="QScore_xen1_model"
)

rsq <- function(x, y) summary(lm(y~x))$r.squared

plot_qscore_comparison <- function(x, y, x_label="x", y_label="y",
    title="QScore comparison",
    cor_method="spearman") {

    stopifnot(length(x) == length(y))

    df <- data.frame(x=x, y=y)

    cor_val <- cor(df$x, df$y, use="complete.obs", method=cor_method)
    rsq <- rsq(x, y)

    p <- ggplot(df, aes(x=x, y=y)) +
        geom_point_rast(alpha=0.35, size=0.6, color="#c0c8cf", raster.dpi=600) +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
        labs(x=x_label, y=y_label, title=title) +
        theme_minimal() +
        scale_x_continuous(limits=c(0, 1)) +
        scale_y_continuous(limits=c(0, 1)) +
        annotate(
            "text",
            x=0.9,
            y=0.05,
            label=paste0("R² = ",
                formatC(rsq, digits=3, format="f")),
            hjust=1,
            vjust=0,
            color="black"
        ) +
        annotate(
            "text",
            x=0.9,
            y=0.1,
            label=paste0(cor_method, " = ",
                formatC(cor_val, digits=3, format="f")),
            hjust=1,
            vjust=0,
            color="black"
        )

    return(list(
        plot=p,
        correlation=cor_val,
        r_squared=rsq
    ))
}

p1 <- plot_qscore_comparison(
    x=xen1$QScore_native,
    y=xen1$QScore_xen2_model,
    x_label="Xenium DCIS S1-Bottom native QS",
    y_label="Xenium DCIS S1-Bottom trained QS",
    title="Xenium DCIS S1-Bottom: native QS vs \n Xenium DCIS S2-Bottom-trained QS",
    cor_method="spearman"
)

p2 <- plot_qscore_comparison(
    x=xen2$QScore_native,
    y=xen2$QScore_xen1_model,
    x_label="Xenium DCIS S2-Bottom native QS",
    y_label="Xenium DCIS S2-Bottom trained QS",
    title="Xenium DCIS S2-Bottom: native QS vs \n Xenium DCIS S1-Bottom-trained QS",
    cor_method="spearman"
)

# Xenium DCIS-1 <-> Xenium breast cancer DBKERO coefficient transfer
    
# dirname pointing to directory containing count matrix, fov position,
# polygons and metadata files as downloaded from:
# Xenium DCIS S1-Bottom: https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard
# Xenium DBKERO: https://kero.hgc.jp/Breast_Cancer_Spatial.html

dirname1 <- "D:\\Xenium_v1_DCIS_S1-Bottom"
dirname2 <- "D:\\Xenium_DBKERO"

xen1 <- readXeniumSPE(dirname1)
xen_db <- readXeniumSPE(dirname2)

xen1 <- spatialPerCellQC(xen1)
xen_db <- spatialPerCellQC(xen_db)

# no need to specify the formula, it's the default one
set.seed(1998)

xen1 <- computeQScore(
    spe=xen1,
    modelFormula=NULL,
    verbose=TRUE
)

# no log2Ctrl_total_ratio was used, not enough outliers

xen_db <- computeQScore(
    spe=xen_db,
    modelFormula=NULL,
    verbose=TRUE
)

# no log2Ctrl_total_ratio was used, not enough outliers

xen1$QScore_native <- xen1$QScore
xen_db$QScore_native <- xen_db$QScore

xen1_model <- metadata(xen1)$QScore_model
xen_db_model <- metadata(xen_db)$QScore_model

xen1 <- .applyQScoreModel(
    spe=xen1,
    qsModel=xen_db_model,
    scoreName="QScore_xen_db_model"
)

xen_db <- .applyQScoreModel(
    spe=xen_db,
    qsModel=xen1_model,
    scoreName="QScore_xen1_model"
)

p3 <- plot_qscore_comparison(
    x=xen1$QScore_native,
    y=xen1$QScore_xen_db_model,
    x_label="Xenium DCIS S1-Bottom native QS",
    y_label="Xenium DCIS S1-Bottom trained QS",
    title="Xenium DCIS S1-Bottom: native QS vs \n Xenium DBKERO-trained QS",
    cor_method="spearman"
)

p4 <- plot_qscore_comparison(
    x=xen_db$QScore_native,
    y=xen_db$QScore_xen1_model,
    x_label="Xenium DBKERO native QS",
    y_label="Xenium DBKERO trained QS",
    title="Xenium DBKERO: native QS vs \n Xenium DCIS S1-Bottom-trained QS",
    cor_method="spearman"
)

# CosMx breast cancer DBKERO <-> Xenium breast cancer DBKERO coefficient transfer
    
# dirname pointing to directory containing count matrix, fov position,
# polygons and metadata files as downloaded from:
# https://kero.hgc.jp/Breast_Cancer_Spatial.html

cos_db <- readCosmxSPE(dirname1)
xen_db <- readXeniumSPE(dirname2)

cos_db <- spatialPerCellQC(cos_db)
xen_db <- spatialPerCellQC(xen_db)

# this time we set a common formula, as the border effect is not transferrable between CosMx and Xenium

common_formula <- "~(log2SignalDensity + Area_um + log2Ctrl_total_ratio)^2"

set.seed(1998)

cos_db <- computeQScore(
    spe=cos_db,
    modelFormula=common_formula,
    verbose=TRUE
)

xen_db <- computeQScore(
    spe=xen_db,
    modelFormula=common_formula,
    verbose=TRUE
)

# log2Ctrl_total_ratio was not used in the QS computation for Xenium DBKERO: not enough outliers

cos_db$QScore_native <- cos_db$QScore
xen_db$QScore_native <- xen_db$QScore

cos_db_model <- metadata(cos_db)$QScore_model
xen_db_model <- metadata(xen_db)$QScore_model

cos_db <- .applyQScoreModel(
    spe=cos_db,
    qsModel=xen_db_model,
    scoreName="QScore_xen_db_model"
)

xen_db <- .applyQScoreModel(
    spe=xen_db,
    qsModel=cos_db_model,
    scoreName="QScore_cos_db_model"
)

p5 <- plot_qscore_comparison(
    x=cos_db$QScore_native,
    y=cos_db$QScore_xen_db_model,
    x_label="CosMx DBKERO native QS",
    y_label="CosMx DBKERO trained QS",
    title="CosMx DBKERO: native QS vs \n Xenium DBKERO-trained QS",
    cor_method="spearman"
)

p6 <- plot_qscore_comparison(
    x=xen_db$QScore_native,
    y=xen_db$QScore_cos_db_model,
    x_label="Xenium DBKERO native QS",
    y_label="Xenium DBKERO trained QS",
    title="Xenium DBKERO: native QS vs \n CosMx DBKERO-trained QS",
    cor_method="spearman"
)

# --- Build rows ---
row1 <- p1$plot | p2$plot
row2 <- p3$plot | p4$plot
row3 <- p5$plot | p6$plot
# --- Add equal spacing between rows ---

final_plot <- (row1 / plot_spacer() /
                 row2 / plot_spacer() /
                 row3) +
  plot_layout(
    heights = c(1, 0.01, 1, 0.01, 1),  # smaller spacing between rows
    widths = c(1, 1, 1)
  ) +
  plot_annotation(
    tag_levels = list(c("C", " ", "D", " ", "E", " ")),
    tag_suffix = "."
  ) &
  theme(plot.tag = element_text(size = 14, face = "bold"))

### save the final figure ###
pdf(file.path("SuppFigure_25C-E.pdf"), 
    width = 10, 
    height = 14,
    bg = "transparent")
final_plot
dev.off()

# sessionInfo()
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
#  [1] ggrastr_1.0.2               patchwork_1.3.2            
#  [3] robustbase_0.99-6           ggplot2_4.0.1              
#  [5] SpaceTrooper_1.1.9          testthat_3.3.2             
#  [7] SpatialExperiment_1.20.0    SingleCellExperiment_1.32.0
#  [9] SummarizedExperiment_1.40.0 Biobase_2.70.0             
# [11] GenomicRanges_1.62.1        Seqinfo_1.0.0              
# [13] IRanges_2.44.0              S4Vectors_0.48.0           
# [15] BiocGenerics_0.56.0         generics_0.1.4             
# [17] MatrixGenerics_1.22.0       matrixStats_1.5.0          

# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3        rstudioapi_0.17.1        
#   [3] jsonlite_2.0.0            shape_1.4.6.1            
#   [5] magrittr_2.0.4            ggbeeswarm_0.7.3         
#   [7] magick_2.9.0              farver_2.1.2             
#   [9] fs_1.6.6                  vctrs_0.7.1              
#  [11] Cairo_1.7-0               memoise_2.0.1            
#  [13] DelayedMatrixStats_1.32.0 rstatix_0.7.3            
#  [15] S4Arrays_1.10.1           usethis_3.2.1            
#  [17] BiocNeighbors_2.4.0       broom_1.0.11             
#  [19] Rhdf5lib_1.32.0           SparseArray_1.10.8       
#  [21] Formula_1.2-5             rhdf5_2.54.1             
#  [23] KernSmooth_2.23-26        desc_1.4.3               
#  [25] cachem_1.1.0              lifecycle_1.0.5          
#  [27] iterators_1.0.14          pkgconfig_2.0.3          
#  [29] rsvd_1.0.5                Matrix_1.7-4             
#  [31] R6_2.6.1                  fastmap_1.2.0            
#  [33] rprojroot_2.1.1           scater_1.38.0            
#  [35] dqrng_0.4.1               irlba_2.3.5.1            
#  [37] pkgload_1.4.1             ggpubr_0.6.2             
#  [39] beachmat_2.26.0           labeling_0.4.3           
#  [41] SpatialExperimentIO_1.2.0 abind_1.4-8              
#  [43] compiler_4.5.1            proxy_0.4-29             
#  [45] remotes_2.5.0             bit64_4.6.0-1            
#  [47] withr_3.0.2               S7_0.2.1                 
#  [49] backports_1.5.0           BiocParallel_1.44.0      
#  [51] carData_3.0-5             viridis_0.6.5            
#  [53] DBI_1.2.3                 pkgbuild_1.4.8           
#  [55] HDF5Array_1.38.0          R.utils_2.13.0           
#  [57] ggsignif_0.6.4            DelayedArray_0.36.0      
#  [59] sessioninfo_1.2.3         rjson_0.2.23             
#  [61] classInt_0.4-11           tools_4.5.1              
#  [63] units_1.0-0               vipor_0.4.7              
#  [65] otel_0.2.0                beeswarm_0.4.0           
#  [67] R.oo_1.27.1               glue_1.8.0               
#  [69] h5mread_1.2.1             rhdf5filters_1.22.0      
#  [71] grid_4.5.1                sf_1.0-24                
#  [73] gtable_0.3.6              R.methodsS3_1.8.2        
#  [75] class_7.3-23              tidyr_1.3.2              
#  [77] data.table_1.18.0         BiocSingular_1.26.1      
#  [79] ScaledMatrix_1.18.0       car_3.1-3                
#  [81] XVector_0.50.0            ggrepel_0.9.6            
#  [83] foreach_1.5.2             pillar_1.11.1            
#  [85] limma_3.66.0              splines_4.5.1            
#  [87] dplyr_1.1.4               lattice_0.22-7           
#  [89] survival_3.8-3            bit_4.6.0                
#  [91] tidyselect_1.2.1          locfit_1.5-9.12          
#  [93] scuttle_1.20.0            sfheaders_0.4.5          
#  [95] gridExtra_2.3             edgeR_4.8.2              
#  [97] statmod_1.5.1             devtools_2.4.6           
#  [99] DropletUtils_1.30.0       brio_1.1.5               
# [101] DEoptimR_1.1-4            codetools_0.2-20         
# [103] tibble_3.3.1              cli_3.6.5                
# [105] arrow_23.0.0              dichromat_2.0-0.1        
# [107] Rcpp_1.1.1                parallel_4.5.1           
# [109] ellipsis_0.3.2            assertthat_0.2.1         
# [111] sparseMatrixStats_1.22.0  glmnet_4.1-10            
# [113] viridisLite_0.4.2         scales_1.4.0             
# [115] e1071_1.7-17              purrr_1.2.1              
# [117] rlang_1.1.7               cowplot_1.2.0  

