# Clean, minimal, and reproducible script to reproduce main figure 3 and related supplementaries 16 and 17
# To reproduce the analysis be sure to have the required packages installed

library(SpaceTrooper)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(sf)
library(patchwork)

# R = 4.5.1

# Load CosMx WTx pancreas dataset

dataset <- "CosMx_rna_pancreas"
fov.set.1 <- c(51:59)
fov.set.2 <- c(60:68)

source("service/service.R")

# celltype palette

celltype_palette <- c(
  "Not assigned" = "gray60",
  "Alpha cells" = "lightseagreen",
  "Beta cells" =	"navyblue",
  "Delta cells" =	"royalblue",
  "Ductal" = "gold2",
  "Epsilon cells" =	"cadetblue1",
  "Gamma cells" =	"deepskyblue",
  "Active stellate" =	"forestgreen",
  "Quiescent stellate" =	"palegreen",
  "Macrophage" =	"orchid1",
  "Acinar.1" =	"firebrick1",
  "Acinar.2" =	"darkred"
)

### Define my theme
my.theme <- theme(panel.border=element_blank(),
                  panel.background=element_rect(fill="transparent", color=NA),
                  plot.background=element_rect(fill="transparent", color=NA),
                  plot.title = element_blank(),
                  axis.line=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.text.y=element_blank(),
                  axis.text.x=element_blank(),
                  legend.background=element_rect(fill="transparent", color=NA),
                  legend.title= element_text(color = "black", size = 10),
                  legend.text=element_text(color="black", size = 6),
                  legend.position = "right",
                  legend.justification = c(0.95, 0.5),
                  legend.box.just = "right",
                  legend.margin = margin(t = -15)  ,
                  panel.grid=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.grid.major=element_blank())

# Create SpatialExperiment as shown in 
# bicciatolab_data/SpaceTrooper_data_and_scripts/Scripts/Scripts_for_spe/CosMx_Pancreas_WTx_analysis_script.R

spe.orig <- readRDS("CosMx_WTx_pancreas_spe.rds")
spe.set.1.orig <- spe.orig[,rownames(spe.orig@colData[spe.orig@colData$fov%in%fov.set.1,])]
spe.set.2.orig <- spe.orig[,rownames(spe.orig@colData[spe.orig@colData$fov%in%fov.set.2,])]
spe.set.1 <- spe.orig@colData[spe.orig@colData$fov%in%fov.set.1,]
spe.set.2 <- spe.orig@colData[spe.orig@colData$fov%in%fov.set.2,]

# Supplementary figure 14

### Area/Volume
var.to.plot <- "Area_um"
spe.set.1$polygons[,var.to.plot] <- spe.set.1[,var.to.plot]
colLabel <- "Area"
p1 <- ggplot() +
  geom_sf(data=spe.set.1$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="rocket", direction = -1) +
  scale_color_viridis_c(option="rocket", direction = -1, guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.set.1.orig$polygons)[1],
           xmax = st_bbox(spe.set.1.orig$polygons)[3],
           ymin = st_bbox(spe.set.1.orig$polygons)[2],
           ymax = st_bbox(spe.set.1.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

### Signal density
var.to.plot <- "log2SignalDensity"
spe.set.1$polygons[,var.to.plot] <- spe.set.1[,var.to.plot]
colLabel <- "Signal density"
p2 <- ggplot() +
  geom_sf(data=spe.set.1$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="mako") +
  scale_color_viridis_c(option="mako", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.set.1.orig$polygons)[1],
           xmax = st_bbox(spe.set.1.orig$polygons)[3],
           ymin = st_bbox(spe.set.1.orig$polygons)[2],
           ymax = st_bbox(spe.set.1.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

### Proportion of control probe counts
var.to.plot <- "log2Ctrl_total_ratio"
spe.set.1$polygons[,var.to.plot] <- spe.set.1[,var.to.plot]
colLabel <- "Proportion of control probe counts"
p3 <- ggplot() +
  geom_sf(data=spe.set.1$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  colorspace::scale_fill_continuous_sequential(palette = "BuPu", 
                                               begin = 0.15, 
                                               breaks = scales::breaks_pretty(n = 5),
                                               limits = c(-8, -2),
                                               oob = scales::squish) + 
  colorspace::scale_color_continuous_sequential(palette = "BuPu", 
                                                begin = 0.15, 
                                                guide = "none",
                                                limits = c(-8, -2),
                                                oob = scales::squish) + 
  labs(title = colLabel,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.set.1.orig$polygons)[1],
           xmax = st_bbox(spe.set.1.orig$polygons)[3],
           ymin = st_bbox(spe.set.1.orig$polygons)[2],
           ymax = st_bbox(spe.set.1.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

### Border effect
var.to.plot <- "abs_log2AspectRatio_dist"
spe.set.1$polygons[,var.to.plot] <- ifelse(spe.set.1$dist_border<50, abs(spe.set.1$log2AspectRatio), NA)
colLabel <- "Aspect ratio"
p4 <- ggplot() +
  geom_sf(data=spe.set.1$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_gradientn(colors = heat.colors(50, rev = TRUE), breaks = scales::breaks_pretty(n = 5), na.value = "grey90") +
  scale_color_gradientn(colors = heat.colors(50, rev = TRUE), guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.set.1.orig$polygons)[1],
           xmax = st_bbox(spe.set.1.orig$polygons)[3],
           ymin = st_bbox(spe.set.1.orig$polygons)[2],
           ymax = st_bbox(spe.set.1.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

### Quality score
var.to.plot <- "QC_score"
spe.set.1$polygons[,var.to.plot] <- spe.set.1[,var.to.plot]
colLabel <- "Quality score"
p5 <- ggplot() +
  geom_sf(data=spe.set.1$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.set.1.orig$polygons)[1],
           xmax = st_bbox(spe.set.1.orig$polygons)[3],
           ymin = st_bbox(spe.set.1.orig$polygons)[2],
           ymax = st_bbox(spe.set.1.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

## add the scale legend
p5 <- plotScaleBar(p5, spe.set.1.orig)

### QS per cell type
var.cell.types <- "cell_types"
var.to.plot <- "QC_score"
types.order <- c("Not assigned",
                 "Alpha cells",
                 "Beta cells",
                 "Gamma cells",
                 "Delta cells",
                 "Epsilon cells",
                 "Acinar.1",
                 "Acinar.2",
                 "Ductal",
                 "Active stellate",
                 "Quiescent stellate",
                 "Macrophage")
spe.set.1$polygons[,var.cell.types] <- factor(spe.set.1[,var.cell.types],
                                        levels = types.order)
colLabel <- "Cell types"
p.viol <- ggplot(spe.set.1$polygons, aes(x = !!sym(var.cell.types), 
                                   y = !!sym(var.to.plot), 
                                   fill = cell_types)) +
  geom_violin(trim = T, color = "black", 
              scale="width",
              alpha = 0.9,
              linewidth = 0.2) +
  geom_boxplot(width = 0.25, 
               outlier.shape = NA, 
               alpha = 0.7, 
               color = "black", 
               linewidth = 0.2) +
  scale_fill_manual(values = celltype_palette) +
  scale_x_discrete(expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     expand = c(0.01,0.01)) +
  theme_minimal() +
  xlab(colLabel) +
  ylab("Quality score") +
  theme_minimal() +
  theme(legend.position="NONE",
        panel.grid.minor = element_blank(),
        panel.border=element_rect(color = "black", fill = "transparent"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### create the final plot
row1 <- p5 | p.viol
row2 <- p1 | p2
row3 <- p3 | p4

# Combine rows with equal heights
final_plot <- row1 / row2 / row3 + 
  plot_layout(
    widths = c(1, 1),          # equal width for columns
    heights = c(1, 1, 1)       # equal height for rows
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

# Save PDF
pdf(file.path("FigureSupp14.pdf"), 
    width = 10, 
    height = 16, 
    bg = "transparent")
print(final_plot)
dev.off()

# Supplementary figure 15

### Area/Volume
var.to.plot <- "Area_um"
spe.set.2$polygons[,var.to.plot] <- spe.set.2[,var.to.plot]
colLabel <- "Area"
p1 <- ggplot() +
  geom_sf(data=spe.set.2$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="rocket", direction = -1) +
  scale_color_viridis_c(option="rocket", direction = -1, guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.set.2.orig$polygons)[1],
           xmax = st_bbox(spe.set.2.orig$polygons)[3],
           ymin = st_bbox(spe.set.2.orig$polygons)[2],
           ymax = st_bbox(spe.set.2.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

### Signal density
var.to.plot <- "log2SignalDensity"
spe.set.2$polygons[,var.to.plot] <- spe.set.2[,var.to.plot]
colLabel <- "Signal density"
p2 <- ggplot() +
  geom_sf(data=spe.set.2$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="mako") +
  scale_color_viridis_c(option="mako", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.set.2.orig$polygons)[1],
           xmax = st_bbox(spe.set.2.orig$polygons)[3],
           ymin = st_bbox(spe.set.2.orig$polygons)[2],
           ymax = st_bbox(spe.set.2.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

### Proportion of control probe counts
var.to.plot <- "log2Ctrl_total_ratio"
spe.set.2$polygons[,var.to.plot] <- spe.set.2[,var.to.plot]
colLabel <- "Proportion of control probe counts"
p3 <- ggplot() +
  geom_sf(data=spe.set.2$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  colorspace::scale_fill_continuous_sequential(palette = "BuPu", 
                                               begin = 0.15, 
                                               breaks = scales::breaks_pretty(n = 5),
                                               limits = c(-8, -2),
                                               oob = scales::squish) + 
  colorspace::scale_color_continuous_sequential(palette = "BuPu", 
                                                begin = 0.15, 
                                                guide = "none",
                                                limits = c(-8, -2),
                                                oob = scales::squish) + 
  labs(title = colLabel,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.set.2.orig$polygons)[1],
           xmax = st_bbox(spe.set.2.orig$polygons)[3],
           ymin = st_bbox(spe.set.2.orig$polygons)[2],
           ymax = st_bbox(spe.set.2.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

### Border effect
var.to.plot <- "abs_log2AspectRatio_dist"
spe.set.2$polygons[,var.to.plot] <- ifelse(spe.set.2$dist_border<50, abs(spe.set.2$log2AspectRatio), NA)
colLabel <- "Aspect ratio"
p4 <- ggplot() +
  geom_sf(data=spe.set.2$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_gradientn(colors = heat.colors(50, rev = TRUE), breaks = scales::breaks_pretty(n = 5), na.value = "grey90") +
  scale_color_gradientn(colors = heat.colors(50, rev = TRUE), guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.set.2.orig$polygons)[1],
           xmax = st_bbox(spe.set.2.orig$polygons)[3],
           ymin = st_bbox(spe.set.2.orig$polygons)[2],
           ymax = st_bbox(spe.set.2.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

### Quality score   
var.to.plot <- "QC_score"
spe.set.2$polygons[,var.to.plot] <- spe.set.2[,var.to.plot]
colLabel <- "Quality score"
p5 <- ggplot() +
  geom_sf(data=spe.set.2$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.set.2.orig$polygons)[1],
           xmax = st_bbox(spe.set.2.orig$polygons)[3],
           ymin = st_bbox(spe.set.2.orig$polygons)[2],
           ymax = st_bbox(spe.set.2.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

## add the scale legend
p5 <- plotScaleBar(p5, spe.set.2.orig)

### QS per cell type
var.cell.types <- "cell_types"
var.to.plot <- "QC_score"
types.order <- c("Not assigned",
                 "Alpha cells",
                 "Beta cells",
                 "Gamma cells",
                 "Delta cells",
                 "Epsilon cells",
                 "Acinar.1",
                 "Acinar.2",
                 "Ductal",
                 "Active stellate",
                 "Quiescent stellate",
                 "Macrophage")
spe.set.2$polygons[,var.cell.types] <- factor(spe.set.2[,var.cell.types],
                                        levels = types.order)
colLabel <- "Cell types"
p.viol <- ggplot(spe.set.2$polygons, aes(x = !!sym(var.cell.types), 
                                   y = !!sym(var.to.plot), 
                                   fill = cell_types)) +
  geom_violin(trim = T, color = "black", 
              scale="width",
              alpha = 0.9,
              linewidth = 0.2) +
  geom_boxplot(width = 0.25, 
               outlier.shape = NA, 
               alpha = 0.7, 
               color = "black", 
               linewidth = 0.2) +
  scale_fill_manual(values = celltype_palette) +
  scale_x_discrete(expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     expand = c(0.01,0.01)) +
  theme_minimal() +
  xlab(colLabel) +
  ylab("Quality score") +
  theme_minimal() +
  theme(legend.position="NONE",
        panel.grid.minor = element_blank(),
        panel.border=element_rect(color = "black", fill = "transparent"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### create the final plot
row1 <- p5 | p.viol
row2 <- p1 | p2
row3 <- p3 | p4

# Combine rows with equal heights
final_plot <- row1 / row2 / row3 + 
  plot_layout(
    widths = c(1, 1),          # equal width for columns
    heights = c(1, 1, 1)       # equal height for rows
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

## Save PDF
pdf(file.path("FigureSupp15.pdf"), 
    width = 10, height = 16, bg = "transparent")
print(final_plot)
dev.off()

sessionInfo()
# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
# LAPACK version 3.12.1
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: Europe/Rome
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.3.2             sf_1.0-24                   dplyr_1.1.4                 ggplot2_4.0.1              
# [5] SpaceTrooper_1.1.3          SpatialExperiment_1.20.0    SingleCellExperiment_1.32.0 SummarizedExperiment_1.40.0
# [9] Biobase_2.70.0              GenomicRanges_1.62.1        Seqinfo_1.0.0               IRanges_2.44.0             
# [13] S4Vectors_0.48.0            BiocGenerics_0.56.0         generics_0.1.4              MatrixGenerics_1.22.0      
# [17] matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.3                 gridExtra_2.3             rlang_1.1.7               magrittr_2.0.4           
# [5] scater_1.38.0             e1071_1.7-17              compiler_4.5.1            DelayedMatrixStats_1.32.0
# [9] SpatialExperimentIO_1.2.0 sfheaders_0.4.5           vctrs_0.7.1               pkgconfig_2.0.3          
# [13] shape_1.4.6.1             backports_1.5.0           magick_2.9.0              XVector_0.50.0           
# [17] scuttle_1.20.0            ggbeeswarm_0.7.3          purrr_1.2.1               bit_4.6.0                
# [21] glmnet_4.1-10             beachmat_2.26.0           rhdf5filters_1.22.0       DelayedArray_0.36.0      
# [25] Rhdf5lib_1.32.0           BiocParallel_1.44.0       broom_1.0.11              irlba_2.3.5.1            
# [29] parallel_4.5.1            R6_2.6.1                  RColorBrewer_1.1-3        limma_3.66.0             
# [33] car_3.1-3                 Rcpp_1.1.1                assertthat_0.2.1          iterators_1.0.14         
# [37] R.utils_2.13.0            Matrix_1.7-4              splines_4.5.1             tidyselect_1.2.1         
# [41] viridis_0.6.5             rstudioapi_0.17.1         dichromat_2.0-0.1         abind_1.4-8              
# [45] codetools_0.2-20          lattice_0.22-7            tibble_3.3.1              withr_3.0.2              
# [49] S7_0.2.1                  survival_3.8-3            units_1.0-0               proxy_0.4-29             
# [53] pillar_1.11.1             ggpubr_0.6.2              carData_3.0-5             KernSmooth_2.23-26       
# [57] foreach_1.5.2             sparseMatrixStats_1.22.0  scales_1.4.0              class_7.3-23             
# [61] glue_1.8.0                tools_4.5.1               BiocNeighbors_2.4.0       robustbase_0.99-6        
# [65] data.table_1.18.0         ScaledMatrix_1.18.0       locfit_1.5-9.12           ggsignif_0.6.4           
# [69] cowplot_1.2.0             rhdf5_2.54.1              grid_4.5.1                tidyr_1.3.2              
# [73] DropletUtils_1.30.0       edgeR_4.8.2               colorspace_2.1-2          beeswarm_0.4.0           
# [77] BiocSingular_1.26.1       HDF5Array_1.38.0          vipor_0.4.7               Formula_1.2-5            
# [81] cli_3.6.5                 rsvd_1.0.5                viridisLite_0.4.2         S4Arrays_1.10.1          
# [85] arrow_23.0.0              gtable_0.3.6              DEoptimR_1.1-4            R.methodsS3_1.8.2        
# [89] rstatix_0.7.3             classInt_0.4-11           ggrepel_0.9.6             SparseArray_1.10.8       
# [93] dqrng_0.4.1               rjson_0.2.23              farver_2.1.2              R.oo_1.27.1              
# [97] lifecycle_1.0.5           h5mread_1.2.1             statmod_1.5.1             bit64_4.6.0-1