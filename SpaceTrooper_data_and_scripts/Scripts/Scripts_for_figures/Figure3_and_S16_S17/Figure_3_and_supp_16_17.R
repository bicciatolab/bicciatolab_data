# Clean, minimal, and reproducible script to reproduce main figure 3 and related supplementaries 16 and 17
# To reproduce the analysis be sure to have the required packages installed

library(SpaceTrooper)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(sf)
library(patchwork)
library(arrow)

# R = 4.5.1

# Load Xenium dataset

dataset <- "Xenium_lung_cancer"
fov.insert.1 <- "D18"
fov.insert.2 <- "B18"
fov.insert.3 <- "E9"
fov.insert.4 <- "D12"

source("service/service.R")

### Define ad hoc plot theme
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
# bicciatolab_data/SpaceTrooper_data_and_scripts/Scripts/Scripts_for_spe/Xenium_lung_cancer_analysis_script.R

spe.orig <- readRDS("Xenium_lung_cancer_spe.rds")

sub.spe.orig.1 <- spe.orig[,spe.orig$fov%in%fov.insert.1]
sub.spe.orig.2 <- spe.orig[,spe.orig$fov%in%fov.insert.2]
sub.spe.orig.3 <- spe.orig[,spe.orig$fov%in%fov.insert.3]
sub.spe.orig.4 <- spe.orig[,spe.orig$fov%in%fov.insert.4]
spe <- spe.orig@colData
sub.spe.1 <- sub.spe.orig.1@colData
sub.spe.2 <- sub.spe.orig.2@colData
sub.spe.3 <- sub.spe.orig.3@colData
sub.spe.4 <- sub.spe.orig.4@colData

### Coordinates of insert 1
bbox <- st_bbox(sub.spe.1$polygons)
insert.1.x <- c(bbox[1],bbox[3])
insert.1.y <- c(bbox[2],bbox[4])

### Coordinates of insert 2
bbox <- st_bbox(sub.spe.2$polygons)
insert.2.x <- c(bbox[1],bbox[3])
insert.2.y <- c(bbox[2],bbox[4])

### Coordinates of insert 3
bbox <- st_bbox(sub.spe.3$polygons)
insert.3.x <- c(bbox[1],bbox[3])
insert.3.y <- c(bbox[2],bbox[4])

### Coordinates of insert 4
bbox <- st_bbox(sub.spe.4$polygons)
insert.4.x <- c(bbox[1],bbox[3])
insert.4.y <- c(bbox[2],bbox[4])

# load transcript file downloaded from 
# https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard

tx <- read_parquet("Xenium_lung_cancer/transcripts.parquet")
tx <- st_as_sf(tx, coords = c("x_location", "y_location"))

tx.insert.1 <- tx[tx$fov_name==fov.insert.1,]
tx.insert.2 <- tx[tx$fov_name==fov.insert.2,]

### Shift the insert 1
attr(tx.insert.1$geometry, "bbox") <- st_bbox(sub.spe.1$polygons)

### Shift the insert 2
attr(tx.insert.2$geometry, "bbox") <- st_bbox(sub.spe.2$polygons)

# celltype_simple_palette

celltype_palette_simple <- c(
  "Alveolar" = "#6665DD",
  "B" = "blueviolet",
  "Basal" = "chocolate",
  "Ciliated" = "burlywood",
  "Club" = "darkgoldenrod",
  "DC" = "dodgerblue4",
  "Fibroblast" = "goldenrod1",
  "Goblet" = "lightsalmon",
  "Lymphatic" = "violetred4",
  "Macrophage" = "darkturquoise",
  "Mast" = "slategray1", 
  "Mesothelial" = "deeppink",
  "Monocyte" = "blue",
  "Pericyte" = "lightgoldenrod1",
  "PNEC" = "orange", #pulmonary neuroendocrine cells
  "SMC" = "lightgoldenrod4",
  "T/NK" = "limegreen",
  "undefined" = "grey80",
  "Vasc.Endo" = "red")

# Figure 3A

### Quality score
var.to.plot <- "QC_score"
spe$polygons[,var.to.plot] <- spe[,var.to.plot]
colLabel <- "Quality score"
p1 <- ggplot() +
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
           xmin = insert.1.x[1],
           xmax = insert.1.x[2],
           ymin = insert.1.y[1],
           ymax = insert.1.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  annotate("rect",
           xmin = insert.2.x[1],
           xmax = insert.2.x[2],
           ymin = insert.2.y[1],
           ymax = insert.2.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2)

## add the scale legend
p1 <- plotScaleBar(p1, spe.orig)

### save the final figure ###
pdf(file.path("Figure_3A.pdf"), 
    width = 10, 
    height = 6,
    bg = "transparent")
p1
dev.off()

# Figure 3B and 3C

### Quality score insert 1
var.to.plot <- "QC_score"
colLabel <- "Quality score"
sub.spe.1$polygons[,var.to.plot] <- sub.spe.1[,var.to.plot]
p1 <- ggplot() + 
  geom_sf(data=sub.spe.1$polygons, 
          aes(fill=!!sym(var.to.plot),
              color = !!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="plasma", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", 
                        guide = "none") +
  annotate("rect",
           xmin = insert.1.x[1],
           xmax = insert.1.x[2],
           ymin = insert.1.y[1],
           ymax = insert.1.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  labs(title = NULL,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2)))

### Quality score insert 2
sub.spe.2$polygons[,var.to.plot] <- sub.spe.2[,var.to.plot]
p2 <- ggplot() + 
  geom_sf(data=sub.spe.2$polygons, 
          aes(fill=!!sym(var.to.plot),
              color = !!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="plasma", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", 
                        guide = "none") +
  annotate("rect",
           xmin = insert.2.x[1],
           xmax = insert.2.x[2],
           ymin = insert.2.y[1],
           ymax = insert.2.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  labs(title = NULL,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2)))
### Signal density insert 1
var.to.plot <- "log2SignalDensity"
colLabel <- "Signal density"
sub.spe.1$polygons[,var.to.plot] <- sub.spe.1[,var.to.plot]
p3 <- ggplot() + 
  geom_sf(data=sub.spe.1$polygons, 
          aes(fill=!!sym(var.to.plot),
              color = !!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="mako", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako", 
                        guide = "none") +
  annotate("rect",
           xmin = insert.1.x[1],
           xmax = insert.1.x[2],
           ymin = insert.1.y[1],
           ymax = insert.1.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  labs(title = NULL,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2)))

### Quality score insert 2
sub.spe.2$polygons[,var.to.plot] <- sub.spe.2[,var.to.plot]
p4 <- ggplot() + 
  geom_sf(data=sub.spe.2$polygons, 
          aes(fill=!!sym(var.to.plot),
              color = !!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="mako", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako", 
                        guide = "none") +
  annotate("rect",
           xmin = insert.2.x[1],
           xmax = insert.2.x[2],
           ymin = insert.2.y[1],
           ymax = insert.2.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  labs(title = NULL,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2)))
p4 <- plotScaleBar(p4,  sub.spe.orig.2)

### Final plot and save ###  
final_plot <- ((p1 | p3) /
                 (p2 | p4)) +
  plot_layout(heights = c(1, 1, 1), widths = c(1, 1))
# --- Build rows ---
row1 <- p1 | p3
row2 <- p2 | p4
# --- Add equal spacing between rows ---
final_plot <- (row1 / plot_spacer() /
                 row2) +
  plot_layout(
    heights = c(1, 0.05, 1),
    widths = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

### save the final figure ###
pdf(file.path("Figure3B_C.pdf"), 
    width = 10, 
    height = 6,
    bg = "transparent")
final_plot
dev.off()

# Supplementary figure 16A and 16D

### QS per cell type
var.cell.types <- "InSituType_Simple"
var.to.plot <- "QC_score"
types.order <- c("undefined",
                 "PNEC",
                 "Mesothelial",
                 "Alveolar",
                 "Basal",
                 "Ciliated",
                 "Club",
                 "Goblet",
                 "B",
                 "DC",
                 "T/NK",
                 "Macrophage",
                 "Monocyte",
                 "Mast",
                 "Fibroblast",
                 "SMC",
                 "Vasc.Endo",
                 "Pericyte",
                 "Lymphatic")
spe$polygons[,var.cell.types] <- factor(spe[,var.cell.types],
                                            levels = types.order)
spe$polygons[,var.to.plot] <- spe[,var.to.plot]
spe.filt <- spe$polygons %>% 
  filter(!is.na(!!sym(var.to.plot)))
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
  scale_fill_manual(values = celltype_palette_simple) +
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

# Supplementary figure 16D

### Quality score in border insert
var.to.plot <- "QC_score"
sub.spe.3$polygons[,var.to.plot] <- sub.spe.3[,var.to.plot]
colLabel <- "Quality score"
p1 <- ggplot() +
  geom_sf(data=sub.spe.3$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_viridis_c(option="plasma", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", 
                        guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  annotate("rect",
           xmin = insert.3.x[1],
           xmax = insert.3.x[2],
           ymin = insert.3.y[1],
           ymax = insert.3.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) 
p1 <- plotScaleBar(p1, sub.spe.orig.3)

### Signal density in border insert
var.to.plot <- "log2SignalDensity"
colLabel <- "Signal density"
sub.spe.3$polygons[,var.to.plot] <- sub.spe.3[,var.to.plot]
p2 <- ggplot() +
  geom_sf(data=sub.spe.3$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_viridis_c(option="mako", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako", 
                        guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = insert.3.x[1],
           xmax = insert.3.x[2],
           ymin = insert.3.y[1],
           ymax = insert.3.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2)))

### create the final plot
row1 <- p.viol 
row2 <- p1 | p2

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
pdf(file.path("SuppFigure16A_D.pdf"), 
    width = 10, 
    height = 16, 
    bg = "transparent")
print(final_plot)
dev.off()


# Supplementary figure 16C

### Signal density insert 1
var.to.plot <- "log2SignalDensity"
colLabel <- "Signal density"
sub.spe.1$polygons[,var.to.plot] <- sub.spe.1[,var.to.plot]
p1 <- ggplot() + 
  geom_sf(data=sub.spe.1$polygons, 
          aes(fill=!!sym(var.to.plot),
              color = !!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="mako", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako", 
                        guide = "none") +
  geom_sf(data = tx.insert.1,
          color = "firebrick",
          size = 0.005,
          alpha = 0.4) +
  coord_sf(xlim = c(insert.1.x[1], insert.1.x[2]-10),
           ylim = c(insert.1.y[1], insert.1.y[2]-15),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.1.x[1],
           xmax = insert.1.x[2]-10,
           ymin = insert.1.y[1],
           ymax = insert.1.y[2]-15,
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  labs(title = NULL,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2)))

### Signal density insert 2
sub.spe.2$polygons[,var.to.plot] <- sub.spe.2[,var.to.plot]
p2 <- ggplot() + 
  geom_sf(data=sub.spe.2$polygons, 
          aes(fill=!!sym(var.to.plot),
              color = !!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="mako", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako", 
                        guide = "none") +
  geom_sf(data = tx.insert.2,
          color = "firebrick",
          size = 0.005,
          alpha = 0.4) +
  coord_sf(xlim = c(insert.2.x[1], insert.2.x[2]-23),
           ylim = c(insert.2.y[1], insert.2.y[2]-8),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.2.x[1],
           xmax = insert.2.x[2]-23,
           ymin = insert.2.y[1],
           ymax = insert.2.y[2]-8,
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  labs(title = NULL,fill = NULL) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2)))

### Final plot and save ###  
final_plot <- ((p1 | p2)) +
  plot_layout(heights = c(1, 1), widths = c(1, 1))
# --- Build rows ---
row1 <- p1 | p2
# --- Add equal spacing between rows ---
final_plot <- (row1) +
  plot_layout(
    heights = c(1, 0.05, 1),
    widths = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

### save the final figure ###
pdf(file.path("SuppFigure16C.pdf"), 
    width = 10, 
    height = 6,
    bg = "transparent")
final_plot
dev.off()

# Supplementary figure 16E

xmin <- 6456.5
xmax <- 6546
ymin <- 2226    
ymax <- 2314

### Quality score in border insert
var.to.plot <- "QC_score"
sub.spe.4$polygons[,var.to.plot] <- sub.spe.4[,var.to.plot]
colLabel <- "Quality score"
p1 <- ggplot() +
  geom_sf(data=sub.spe.4$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_viridis_c(option="plasma", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", 
                        guide = "none") +
  coord_sf(xlim = c(xmin, xmax),
           ylim = c(ymin, ymax),
           expand = FALSE) +
  labs(title = NULL,
       fill = NULL) +
  annotate("rect",
           xmin = xmin,
           xmax = xmax,
           ymin = ymin,
           ymax = ymax,
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2)))

# Save PDF
pdf(file.path("SuppFigure16E.pdf"), 
    width = 10, 
    height = 16, 
    bg = "transparent")
print(p1)
dev.off()

# Load MERFISH dataset

dataset <- "MERFISH_mouse_liver"
fov.insert.bord <- c(1699:1706, 1736:1749, 1773:1780)
fov.insert.triad <- c(593:600,653:660)

# Create SpatialExperiment as shown in 
# bicciatolab_data/SpaceTrooper_data_and_scripts/Scripts/Scripts_for_spe/Xenium_lung_cancer_analysis_script.R

spe.orig <- readRDS("MERFISH_mouse_liver_spe.rds")

sub.spe.orig.bord <- spe.orig[,spe.orig$fov%in%fov.insert.bord]
sub.spe.orig.triad <- spe.orig[,spe.orig$fov%in%fov.insert.triad]
spe <- spe.orig@colData
sub.spe.bord <- sub.spe.orig.bord@colData
sub.spe.triad <- sub.spe.orig.triad@colData

### Coordinates of the board insert
bbox <- st_bbox(sub.spe.bord$polygons)
insert.bord.x <- c(bbox[1],bbox[3])
insert.bord.y <- c(bbox[2],bbox[4])

### Coordinates of the triad insert
bbox <- st_bbox(sub.spe.triad$polygons)
insert.triad.x <- c(bbox[1],bbox[3])
insert.triad.y <- c(bbox[2],bbox[4])

celltype_palette_simple <- c(
  "B cells" = "blueviolet",
  "Basophils" = "slategray1",
  "DC" = "navyblue",
  "Cholangiocytes" = "burlywood",
  "Endothelial cells" = "red", 
  "Fibroblast" = "yellow",
  "Hepatocytes" = "chocolate", 
  "HsPCs" = "firebrick",
  "ILC1s" = "violetred4",
  "Kupffer cells" = "#6665DD",
  "Mesothelial cells" = "#eb37b2",
  "Myeloid" = "royalblue1",
  "Neutrophils" = "#9B9ECE", 
  "T/NK cells" = "darkolivegreen1",
  "pDCs" = "mediumpurple",
  "Stellate cells" = "tomato",
  "undefined" = "grey80",
  "VSMCs" = "lightgoldenrod4"
)

# Figure 3D 

### Quality score
var.to.plot <- "QC_score"
spe$polygons[,var.to.plot] <- spe[,var.to.plot]
colLabel <- "Quality score"
p1 <- ggplot() +
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
           xmin = insert.bord.x[1],
           xmax = insert.bord.x[2],
           ymin = insert.bord.y[1],
           ymax = insert.bord.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.5) +
  annotate("rect",
           xmin = insert.triad.x[1],
           xmax = insert.triad.x[2],
           ymin = insert.triad.y[1],
           ymax = insert.triad.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.5)
## add the scale legend
p1 <- plotScaleBar(p1, spe.orig)

### save the final figure ###
pdf(file.path("Figure3_D.pdf"), 
    width = 10, 
    height = 6,
    bg = "transparent")
p1
dev.off()

# Figure 3E and 3F

### Quality score in border insert
var.to.plot <- "QC_score"
sub.spe.bord$polygons[,var.to.plot] <- sub.spe.bord[,var.to.plot]
colLabel <- "Quality score"
p1 <- ggplot() +
  geom_sf(data=sub.spe.bord$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_viridis_c(option="plasma", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", 
                        guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  annotate("rect",
           xmin = insert.bord.x[1],
           xmax = insert.bord.x[2],
           ymin = insert.bord.y[1],
           ymax = insert.bord.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  my.theme +
  theme(legend.position = "none")
p1 <- plotScaleBar(p1, sub.spe.orig.bord)

### Quality score in triad insert
sub.spe.triad$polygons[,var.to.plot] <- sub.spe.triad[,var.to.plot]
colLabel <- "Quality score"
p2 <- ggplot() +
  geom_sf(data=sub.spe.triad$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_viridis_c(option="plasma", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma",
                        guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  annotate("rect",
           xmin = insert.triad.x[1],
           xmax = insert.triad.x[2],
           ymin = insert.triad.y[1],
           ymax = insert.triad.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  my.theme +
  theme(legend.position = "none")
p2 <- plotScaleBar(p2, sub.spe.orig.triad)

### Signal density in border insert
var.to.plot <- "log2SignalDensity"
colLabel <- "Signal density"
sub.spe.bord$polygons[,var.to.plot] <- sub.spe.bord[,var.to.plot]
p3 <- ggplot() +
  geom_sf(data=sub.spe.bord$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_viridis_c(option="mako", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako", 
                        guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = insert.bord.x[1],
           xmax = insert.bord.x[2],
           ymin = insert.bord.y[1],
           ymax = insert.bord.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2)))
p3 <- plotScaleBar(p3, sub.spe.orig.bord)

### Density signal in triad insert
sub.spe.triad$polygons[,var.to.plot] <- sub.spe.triad[,var.to.plot]
p4 <- ggplot() +
  geom_sf(data=sub.spe.triad$polygons, 
          aes(fill=!!sym(var.to.plot), 
              color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_viridis_c(option="mako", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako",
                        guide = "none") +
  labs(title = NULL,
       fill = NULL) +
  annotate("rect",
           xmin = insert.triad.x[1],
           xmax = insert.triad.x[2],
           ymin = insert.triad.y[1],
           ymax = insert.triad.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  my.theme +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2)))
p4 <- plotScaleBar(p4, sub.spe.orig.triad)

### create the final plot
row1 <- p1 | p3
row2 <- p2 | p4

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
pdf(file.path("Figure3_E-F.pdf"), 
    width = 10, height = 12, bg = "transparent")
print(final_plot)
dev.off()

# Supplementary figure 17A, 17B and 17C

### QS per cell type
var.cell.types <- "InSituType_Simple"
var.to.plot <- "QC_score"
types.order <- c("undefined",
                 "HsPCs",
                 "Fibroblast",
                 "VSMCs",
                 "Endothelial cells",
                 "Stellate cells",
                 "ILC1s",
                 "Mesothelial cells",
                 "Cholangiocytes",
                 "Hepatocytes",
                 "B cells",
                 "T/NK",
                 "DC",
                 "Kupffer cells",
                 "Myeloid",
                 "Basophils",
                 "Neutrophils")
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
  scale_fill_manual(values = celltype_palette_simple) +
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

pdf(file.path("SuppFigure17A.pdf"), 
    width = 8, 
    height = 6, 
    bg = "transparent")
p.viol
dev.off()

### Cell type in border insert
var.to.plot <- "InSituType_Simple"
sub.spe.bord$polygons[,var.to.plot] <- factor(sub.spe.bord[,var.to.plot],
                                              levels = types.order)
sub.spe.bord$polygons[,var.to.plot] <- sub.spe.bord[,var.to.plot]
spe.filt <- sub.spe.bord$polygons %>% 
  filter(!is.na(!!sym(var.to.plot))) %>%
  group_by(!!sym(var.to.plot)) %>%
  filter(n() >= 2) %>%
  ungroup()
colLabel <- "Cell types"
p1 <- ggplot() +
  geom_sf(data=spe.filt, 
          aes(fill=!!sym(var.to.plot),
              color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_manual(values = celltype_palette_simple,
                    name = colLabel) +
  scale_color_manual(values = celltype_palette_simple,
                     guide = "none")+
  labs(title = NULL,
       fill = NULL) +
  annotate("rect",
           xmin = insert.bord.x[1],
           xmax = insert.bord.x[2],
           ymin = insert.bord.y[1],
           ymax = insert.bord.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  my.theme 
p1 <- plotScaleBar(p1, sub.spe.orig.bord)

# Save PDF
pdf(file.path("SuppFigure17B.pdf"), 
    width = 10, 
    height = 16, 
    bg = "transparent")
p1
dev.off()

### Cell type in triad insert
var.to.plot <- "InSituType_Simple"
sub.spe.triad$polygons[,var.to.plot] <- factor(sub.spe.triad[,var.to.plot],
                                               levels = types.order)
sub.spe.triad$polygons[,var.to.plot] <- sub.spe.triad[,var.to.plot]
spe.filt <- sub.spe.triad$polygons %>% 
  filter(!is.na(!!sym(var.to.plot))) %>%
  group_by(!!sym(var.to.plot)) %>%
  filter(n() >= 2) %>%
  ungroup()
colLabel <- "Cell types"
p2 <- ggplot() +
  geom_sf(data=spe.filt, 
          aes(fill=!!sym(var.to.plot),
          color =!!sym(var.to.plot)), lwd = 0.2)+
  scale_fill_manual(values = celltype_palette_simple,
                    name = colLabel) +
  scale_color_manual(values = celltype_palette_simple,
                     guide = "none")+
  labs(title = NULL,
       fill = NULL) +
  annotate("rect",
           xmin = insert.triad.x[1],
           xmax = insert.triad.x[2],
           ymin = insert.triad.y[1],
           ymax = insert.triad.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  my.theme 

# Save PDF
pdf(file.path("SuppFigure17C.pdf"), 
    width = 10, 
    height = 16, 
    bg = "transparent")
p2
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
#   [1] arrow_23.0.0                    data.table_1.18.0               SpatialFeatureExperiment_1.12.1 SFEData_1.12.0                 
# [5] patchwork_1.3.2                 sf_1.0-24                       dplyr_1.1.4                     ggplot2_4.0.1                  
# [9] SpaceTrooper_1.1.3              SpatialExperiment_1.20.0        SingleCellExperiment_1.32.0     SummarizedExperiment_1.40.0    
# [13] Biobase_2.70.0                  GenomicRanges_1.62.1            Seqinfo_1.0.0                   IRanges_2.44.0                 
# [17] S4Vectors_0.48.0                BiocGenerics_0.56.0             generics_0.1.4                  MatrixGenerics_1.22.0          
# [21] matrixStats_1.5.0              
# 
# loaded via a namespace (and not attached):
#   [1] splines_4.5.1             bitops_1.0-9              filelock_1.0.3            tibble_3.3.1              R.oo_1.27.1              
# [6] lifecycle_1.0.5           httr2_1.2.2               rstatix_0.7.3             edgeR_4.8.2               MASS_7.3-65              
# [11] lattice_0.22-7            backports_1.5.0           magrittr_2.0.4            limma_3.66.0              yaml_2.3.12              
# [16] otel_0.2.0                sp_2.2-0                  cowplot_1.2.0             DBI_1.2.3                 RColorBrewer_1.1-3       
# [21] multcomp_1.4-29           abind_1.4-8               spatialreg_1.4-2          purrr_1.2.1               R.utils_2.13.0           
# [26] RCurl_1.98-1.17           TH.data_1.1-5             sandwich_3.1-1            rappdirs_0.3.4            ggrepel_0.9.6            
# [31] irlba_2.3.5.1             terra_1.8-93              units_1.0-0               dqrng_0.4.1               DelayedMatrixStats_1.32.0
# [36] codetools_0.2-20          DropletUtils_1.30.0       DelayedArray_0.36.0       scuttle_1.20.0            tidyselect_1.2.1         
# [41] shape_1.4.6.1             farver_2.1.2              ScaledMatrix_1.18.0       viridis_0.6.5             BiocFileCache_3.0.0      
# [46] BiocNeighbors_2.4.0       e1071_1.7-17              Formula_1.2-5             survival_3.8-3            scater_1.38.0            
# [51] iterators_1.0.14          systemfonts_1.3.1         foreach_1.5.2             tools_4.5.1               ragg_1.5.0               
# [56] Rcpp_1.1.1                glue_1.8.0                gridExtra_2.3             SparseArray_1.10.8        EBImage_4.52.0           
# [61] HDF5Array_1.38.0          withr_3.0.2               BiocManager_1.30.27       fastmap_1.2.0             boot_1.3-32              
# [66] rhdf5filters_1.22.0       spData_2.3.4              digest_0.6.39             rsvd_1.0.5                R6_2.6.1                 
# [71] textshaping_1.0.4         colorspace_2.1-2          wk_0.9.5                  LearnBayes_2.15.1         jpeg_0.1-11              
# [76] dichromat_2.0-0.1         RSQLite_2.4.5             R.methodsS3_1.8.2         h5mread_1.2.1             tidyr_1.3.2              
# [81] robustbase_0.99-6         class_7.3-23              httr_1.4.7                htmlwidgets_1.6.4         S4Arrays_1.10.1          
# [86] spdep_1.4-1               pkgconfig_2.0.3           gtable_0.3.6              blob_1.3.0                S7_0.2.1                 
# [91] XVector_0.50.0            htmltools_0.5.9           carData_3.0-5             fftwtools_0.9-11          scales_1.4.0             
# [96] png_0.1-8                 rstudioapi_0.17.1         rjson_0.2.23              nlme_3.1-168              coda_0.19-4.1            
# [101] curl_7.0.0                zoo_1.8-15                proxy_0.4-29              cachem_1.1.0              rhdf5_2.54.1             
# [106] BiocVersion_3.22.0        KernSmooth_2.23-26        parallel_4.5.1            vipor_0.4.7               AnnotationDbi_1.72.0     
# [111] s2_1.1.9                  pillar_1.11.1             grid_4.5.1                vctrs_0.7.1               ggpubr_0.6.2             
# [116] car_3.1-3                 BiocSingular_1.26.1       dbplyr_2.5.1              beachmat_2.26.0           sfheaders_0.4.5          
# [121] beeswarm_0.4.0            SpatialExperimentIO_1.2.0 zeallot_0.2.0             magick_2.9.0              mvtnorm_1.3-3            
# [126] cli_3.6.5                 locfit_1.5-9.12           compiler_4.5.1            rlang_1.1.7               crayon_1.5.3             
# [131] ggsignif_0.6.4            labeling_0.4.3            classInt_0.4-11           ggbeeswarm_0.7.3          viridisLite_0.4.2        
# [136] deldir_2.0-4              BiocParallel_1.44.0       assertthat_0.2.1          Biostrings_2.78.0         tiff_0.1-12              
# [141] glmnet_4.1-10             Matrix_1.7-4              ExperimentHub_3.0.0       sparseMatrixStats_1.22.0  bit64_4.6.0-1            
# [146] Rhdf5lib_1.32.0           KEGGREST_1.50.0           statmod_1.5.1             AnnotationHub_3.99.6      broom_1.0.11             
# [151] memoise_2.0.1             DEoptimR_1.1-4            bit_4.6.0
