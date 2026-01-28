# Clean, minimal, and reproducible script to reproduce main figure 4 and related supplementary 18
# To reproduce the analysis be sure to have the required packages installed

library(SpaceTrooper)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(sf)
library(patchwork)

# R = 4.5.1

# Load CosMx Protein dataset

dataset <- "CosMx_protein_tonsil"
fov.insert.ery <- c(24)
fov.insert.poly <- c(201)

source("service.R")

### Define ad hoc plot theme
my.theme <- theme(panel.border=element_blank(),
                  panel.background=element_rect(fill="transparent", color=NA),
                  plot.background=element_rect(fill="transparent", color=NA),
                  plot.title =element_text(face = "bold", size = 12,  hjust = 0.5),
                  axis.line=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.text.y=element_blank(),
                  axis.text.x=element_blank(),
                  legend.background=element_rect(fill="transparent", color=NA),
                  legend.title= element_text(color = "black", size = 10),
                  legend.text=element_text(color="black", size = 6),
                  panel.grid=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.grid.major=element_blank())

# Create SpatialExperiment as shown in 
# bicciatolab_data/SpaceTrooper_data_and_scripts/Scripts/Scripts_for_spe/Xenium_lung_cancer_analysis_script.R

spe.orig <- readRDS("CosMx_protein_tonsil_spe.rds")
sub.spe.ery.orig <- spe.orig[,spe.orig$fov%in%fov.insert.ery]
sub.spe.poly.orig <- spe.orig[,spe.orig$fov%in%fov.insert.poly]
spe <- spe.orig@colData
sub.spe.ery <- spe.orig[,spe.orig$fov%in%fov.insert.ery]@colData
sub.spe.poly <- spe.orig[,spe.orig$fov%in%fov.insert.poly]@colData

### Coordinates of the erythrocytes insert
bbox <- st_bbox(sub.spe.ery$polygons)
insert.ery.x <- c(bbox[1],bbox[3])
insert.ery.y <- c(bbox[2],bbox[4])

### Coordinates of the poly insert
bbox <- st_bbox(sub.spe.poly$polygons)
insert.poly.x <- c(bbox[1],bbox[3])
insert.poly.y <- c(bbox[2],bbox[4])

# celltype palette

celltype_palette <- c(
  "Activated.B.Cell" = "deeppink",
  "B.Cell" = "blueviolet",
  "CD4.T.cell" = "limegreen",
  "dropped" = "grey30",
  "Epithelial"= "#6665DD",
  "ITGAX Slan-like"= "darkturquoise", # type of CD16+ monocytes expressing slan, a specific sugar structure on their surface
  "Macrophages.1" = "blue",
  "Macrophages.2" ="dodgerblue4",
  "Mesenchymal" = "lightgoldenrod1",
  "NK" = "darkolivegreen1",
  "Proliferating.B.Cell"= "orchid", 
  "T.Cell" = "darkseagreen1", 
  "Tfh"= "#134502",
  "unknown" = "#c0c8cf"
)

# Figure 4A

### Quality score
var.to.plot <- "QC_score"
spe$polygons[,var.to.plot] <- spe[,var.to.plot]
colLabel <- "Quality score"
p1 <- ggplot(data=spe$polygons) + 
  geom_sf(data=spe$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", guide = "none") +
  labs(title = NULL,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.orig$polygons)[1],
           xmax = st_bbox(spe.orig$polygons)[3],
           ymin = st_bbox(spe.orig$polygons)[2],
           ymax = st_bbox(spe.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05) +
  annotate("rect",
           xmin = insert.ery.x[1],
           xmax = insert.ery.x[2],
           ymin = insert.ery.y[1],
           ymax = insert.ery.y[2],
           fill = NA,
           color = "dodgerblue",
           linewidth = 1) 
#+ annotate("rect",
#           xmin = insert.poly.x[1],
#           xmax = insert.poly.x[2],
#           ymin = insert.poly.y[1],
#           ymax = insert.poly.y[2],
#           fill = NA,
#           color = "black",
#           linewidth = 0.5)
## add the scale legend
p1 <- plotScaleBar(p1, spe.orig)

### save the final figure ###
pdf(file.path("Figure4_A.pdf"), 
    width = 10, 
    height = 6,
    bg = "transparent")
p1
dev.off()

# Figure 4B, 4C, 4F, 4G, 4H, 4E

### Quality score in erythrocytes insert
var.to.plot <- "QC_score"
sub.spe.ery$polygons[,var.to.plot] <- sub.spe.ery[,var.to.plot]
colLabel <- "Quality score"
p1 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_viridis_c(option="plasma", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", 
                        guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  theme(legend.position = "none") +
  annotate("rect",
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p1 <- plotScaleBar(p1, sub.spe.ery.orig)

### Proportion of control probe counts
var.to.plot <- "log2Ctrl_total_ratio"
sub.spe.ery$polygons[,var.to.plot] <- sub.spe.ery[,var.to.plot]
colLabel <- "Proportion of control probe counts"
p2 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  colorspace::scale_fill_continuous_sequential(palette = "BuPu", begin = 0.15, 
                                               breaks = scales::breaks_pretty(n = 5) ) + 
  colorspace::scale_color_continuous_sequential(palette = "BuPu", 
                                                begin = 0.15, 
                                                guide = "none") + 
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p2 <- plotScaleBar(p2, sub.spe.ery.orig)

### Area equal to zero
var.to.filt <- "NucArea"
var.to.plot <- "nuc_area_zero"
sub.spe.ery$polygons[[var.to.plot]] <- sub.spe.ery[[var.to.filt]] == 0  
colLabel <- "Is nuclear area zero?"
p3 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, aes(fill=!!sym(var.to.plot), 
                                         color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_manual(values = c("TRUE"  = "violetred4",
                               "FALSE" = "grey95"),
                    breaks = c("TRUE", "FALSE"),
                    na.value = "transparent") +
  scale_color_manual(values = c("FALSE"  = "#B3C2F2",
                                "TRUE" = "#B3C2F2"),
                     breaks = c("FALSE", "TRUE"),
                     na.value = "transparent",
                     guide = "none") +
  labs(title = NULL,
       fill  = colLabel) +    
  my.theme + 
  theme(legend.position = "bottom",
        legend.justification = c(0.95, 0.5),
        legend.margin = margin(t = -15)) +
  annotate("rect",
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p3 <- plotScaleBar(p3, sub.spe.ery.orig)

### Pan-CK expression level
var.to.plot <- "Channel-PanCK"
sub.spe.ery$polygons[,var.to.plot] <- log2(sub.spe.ery.orig@assays@data$counts[var.to.plot,])
colLabel <- "Pan-CK"
p4 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1) +
  scale_fill_gradient(low="white", high="blue",
                      breaks = scales::breaks_pretty(n = 5)) +
  scale_color_gradient(low="white", high="blue",
                       breaks = scales::breaks_pretty(n = 5),
                       guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p4 <- plotScaleBar(p4, sub.spe.ery.orig)

### CD68 expression level
var.to.plot <- "CD68"
sub.spe.ery$polygons[,var.to.plot] <- log2(sub.spe.ery.orig@assays@data$counts[var.to.plot,])
colLabel <- "CD68"
p5 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1) +
  scale_fill_gradient(low="white", high="blue",
                      breaks = scales::breaks_pretty(n = 5)) +
  scale_color_gradient(low="white", high="blue",
                       breaks = scales::breaks_pretty(n = 5),
                       guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p5 <- plotScaleBar(p5, sub.spe.ery.orig)

### SMA expression level
var.to.plot <- "SMA"
sub.spe.ery$polygons[,var.to.plot] <- log2(sub.spe.ery.orig@assays@data$counts[var.to.plot,])
colLabel <- "SMA"
p6 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1) +
  scale_fill_gradient(low="white", high="blue",
                      breaks = scales::breaks_pretty(n = 5)) +
  scale_color_gradient(low="white", high="blue",
                       breaks = scales::breaks_pretty(n = 5),
                       guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p6 <- plotScaleBar(p6, sub.spe.ery.orig)

### create the final plot
row1 <- p1 | p2 | p3
row2 <- p4 | p5 | p6

# Combine rows with equal heights
final_plot <- row1 / row2 + 
  plot_layout(
    widths = c(1, 1),          # equal width for columns
    heights = c(1, 1, 1)       # equal height for rows
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

# Save PDF
pdf(file.path("Figure4_BCEFGH.pdf"), 
    width = 10, height = 12, bg = "transparent")
print(final_plot)
dev.off()

# Supplementary figure 18

### QS per cell type
var.cell.types <- "cell_type"
var.to.plot <- "QC_score"
types.order <- c("dropped",
                 "unknown",
                 "B.Cell",
                 "Proliferating.B.Cell",
                 "Activated.B.Cell",
                 "T.Cell",
                 "CD4.T.cell",
                 "Tfh",
                 "NK",
                 "ITGAX Slan-like",
                 "Macrophages.1",
                 "Macrophages.2",
                 "Mesenchymal",
                 "Epithelial")
spe$polygons[,var.cell.types] <- factor(spe[,var.cell.types],
                                            levels = types.order)
spe$polygons[,var.to.plot] <- spe[,var.to.plot]
spe.filt <- spe$polygons %>% 
  filter(!is.na(!!sym(var.to.plot))) %>%
  group_by(!!sym(var.cell.types)) %>%
  filter(n() >= 2) %>%
  ungroup()
colLabel <- "Cell types"
p.viol <- ggplot(spe.filt, aes(x = !!sym(var.cell.types), 
                               y = !!sym(var.to.plot), 
                               fill = !!sym(var.cell.types))) +
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

### Area/Volume in erythrocytes insert
var.to.plot <- "Area_um"
sub.spe.ery$polygons[,var.to.plot] <- sub.spe.ery[,var.to.plot]
colLabel <- "Area"
p1 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="rocket", 
                       direction = -1) +
  scale_color_viridis_c(option="rocket", 
                        direction = -1, 
                        guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p1 <- plotScaleBar(p1, sub.spe.ery.orig)

### Signal density in erythrocytes insert
var.to.plot <- "log2SignalDensity"
sub.spe.ery$polygons[,var.to.plot] <- sub.spe.ery[,var.to.plot]
colLabel <- "Signal density"
p2 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_viridis_c(option="mako") +
  scale_color_viridis_c(option="mako", 
                        guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p2 <- plotScaleBar(p2, sub.spe.ery.orig)

### Border effect in erythrocytes insert
var.to.plot <- "abs_log2AspectRatio_dist"
sub.spe.ery$polygons[,var.to.plot] <- ifelse(sub.spe.ery$dist_border<50, abs(sub.spe.ery$log2AspectRatio), NA)
colLabel <- "Aspect ratio"
p3 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, aes(fill=!!sym(var.to.plot), 
                                         color = !!sym(var.to.plot)), 
          lwd = 0.01)+
  scale_fill_gradientn(colors = heat.colors(50, rev = TRUE), 
                       breaks = scales::breaks_pretty(n = 5), 
                       na.value = "grey90") +
  scale_color_gradientn(colors = heat.colors(50, rev = TRUE), 
                        guide = "none") +
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
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p3 <- plotScaleBar(p3, sub.spe.ery.orig)

### GITR/FOXP3 expression level in erythrocytes insert
#var.to.plot <- "GITR"
var.to.plot <- "FOXP3"
sub.spe.ery$polygons[,var.to.plot] <- log2(sub.spe.ery.orig@assays@data$counts[var.to.plot,])
colLabel <- "GITR"
p4 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_gradient(low="white", high="blue",
                      breaks = scales::breaks_pretty(n = 5)) +
  scale_color_gradient(low="white", high="blue",
                       breaks = scales::breaks_pretty(n = 5),
                       guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p4 <- plotScaleBar(p4, sub.spe.ery.orig)

### EpCAM expression level in erythrocytes insert
var.to.plot <- "EpCAM"
sub.spe.ery$polygons[,var.to.plot] <- log2(sub.spe.ery.orig@assays@data$counts[var.to.plot,])
colLabel <- "EpCAM"
p5 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_gradient(low="white", high="blue",
                      breaks = scales::breaks_pretty(n = 5)) +
  scale_color_gradient(low="white", high="blue",
                       breaks = scales::breaks_pretty(n = 5),
                       guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(sub.spe.ery$polygons)[1],
           xmax = st_bbox(sub.spe.ery$polygons)[3],
           ymin = st_bbox(sub.spe.ery$polygons)[2],
           ymax = st_bbox(sub.spe.ery$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p5 <- plotScaleBar(p5, sub.spe.ery.orig)

### Cell types in erythrocytes insert
var.to.plot <- "cell_type"
sub.spe.ery$polygons[,var.to.plot] <- factor(sub.spe.ery[,var.to.plot],
                                             levels = types.order)
sub.spe.ery$polygons[,var.to.plot] <- sub.spe.ery[,var.to.plot]
p6 <- ggplot() +
  geom_sf(data=sub.spe.ery$polygons, 
          aes(fill=!!sym(var.to.plot),
              color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_manual(values = celltype_palette,
                    name = colLabel) +
  scale_color_manual(values = celltype_palette,
                     guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme
p6 <- plotScaleBar(p6, sub.spe.ery.orig)

### create the final plot
row1 <- p.viol 
row2 <- p1 | p2 | p3
row3 <- p4 | p5 | p6

# Combine rows with equal heights
final_plot <- row1 / row2 / row3 + 
  plot_layout(
    widths = c(1, 1, 1),          # equal width for columns
    heights = c(1, 1, 1)       # equal height for rows
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

# Save PDF
pdf(file.path("SuppFigure18.pdf"), 
    width = 10, 
    height = 16, 
    bg = "transparent")
print(final_plot)
dev.off()






