# Clean, minimal, and reproducible script to reproduce figure 1 and related supplementaries 8-11
# To reproduce the analysis be sure to have the required packages installed

library(SpaceTrooper)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(sf)
library(patchwork)
library(data.table)
library(colorspace)
library(ggrastr)

# R = 4.5.1

dataset <- "CosMx_rna_DBKERO"
fov.13 <- c(13)
fov.11 <- c(11)
fov.11.12 <- c(11,12)
fov.31 <- c(31)
fov.38.39 <- c(38:39)

# load plot theme

source("service/service.R")

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

# Create SpatialExperiment as shown in 
# bicciatolab_data/SpaceTrooper_data_and_scripts/Scripts/Scripts_for_spe/CosMx_DBKERO_1k_BC_analysis_script.R

spe.orig <- readRDS("CosMx_DBKERO_1k_BC_spe.rds")

# subset for each FoV/pair of FoV to show
spe.13.orig <- spe.orig[,rownames(spe.orig@colData[spe.orig@colData$fov%in%fov.13,])]
spe.11.orig <- spe.orig[,rownames(spe.orig@colData[spe.orig@colData$fov%in%fov.11,])]
spe.11.12.orig <- spe.orig[,rownames(spe.orig@colData[spe.orig@colData$fov%in%fov.11.12,])]
spe.31.orig <- spe.orig[,rownames(spe.orig@colData[spe.orig@colData$fov%in%fov.31,])]
spe.38.39.orig <- spe.orig[,rownames(spe.orig@colData[spe.orig@colData$fov%in%fov.38.39,])]
spe <- spe.orig@colData
spe.13 <- spe.13.orig@colData
spe.11 <- spe.11.orig@colData
spe.11.12 <- spe.11.12.orig@colData
spe.31 <- spe.31.orig@colData
spe.38.39 <- spe.38.39.orig@colData

# Load transcript file available at https://kero.hgc.jp/Breast_Cancer_Spatial.html
tx <- fread("Run5810_Case2_tx_file.csv")

# Create subset of transcript file with only nuclear transcripts
sub.tx <- tx[tx$target%in%c("NEAT1", "MALAT1")&
                             tx$fov%in%fov.11,]

# selection of only control probes
ctrl_list <- c("NegPrb")

idxlist <- lapply(ctrl_list, function(ng){
  grep(paste0("^", ng), tx$target)
})

ctrl.tx <- tx[unlist(idxlist),]

### select all non-nuclear transcripts
no.nucl.tx <- tx[tx$target!="MALAT1"& tx$target!="NEAT1",]
  
# Load cells flagged with FastReseg, file available at bicciatolab_data/SpaceTrooper_data_and_scripts/Scripts/Scripts_for_figures/data/Figure1
FastRes <- read.csv(file.path("data/Figure1/CosMx_rna_DBKERO_FastReseg_flags.csv"))

# celltype palette
celltype_palette <- c(
  "TAMs" =	"blue",
  "DCs" =	"dodgerblue4",
  "Mast cells" =	"lightblue1",
  "T cells" =	"darkseagreen1",
  "NK cells" = "darkolivegreen1",
  "Mural cells" =	"goldenrod4",
  "Myoepithelial cells" =	"goldenrod1",
  "Blood ECs" =	"lightgoldenrod1",
  "CAFs" =	"lightgoldenrod4",
  "Plasma cells" =	"deeppink1",
  "Mix BC cells TAMs" = "purple",
  "BC cells" =	"orangered1"
)

# Figure 1B

var.to.plot <- "QC_score"
xlabel <- "Quality score"
  
metadata <- data.frame(spe)
fov.positions <- metadata(spe.orig)$fov_positions

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
  geom_vline(xintercept=median(metadata[,var.to.plot], na.rm = TRUE) -3*mad(metadata[,var.to.plot], na.rm = TRUE),
             col="black",
             linewidth = 0.4,
             linetype = "dashed") + 
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
                     expand = expansion(mult = c(0, 0.05))) +
  ylab("Number of cells") + 
  xlab(xlabel) + 
  theme_minimal() + 
  light_theme +
  theme(plot.margin = unit(c(0,20,0,0), "pt"),
        aspect.ratio = 1)

### save the final figure ###
pdf(file.path("Figure1B.pdf"), 
    width = 10, 
    height = 10,
    bg = "transparent")
pp
dev.off()

# Figure 1C

### Quality score
var.to.plot <- "QC_score"
spe.11.12$polygons[,var.to.plot] <- spe.11.12[,var.to.plot]
colLabel <- "Quality score"
p1 <- ggplot() +
  geom_sf(data=spe.11.12$polygons, aes(fill=!!sym(var.to.plot), 
                                 color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_viridis_c(option="plasma", 
                       breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", 
                        guide = "none") +
  labs(title = NULL,fill = NULL) +
  theme(panel.border=element_blank(),
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
        panel.grid.major=element_blank()) +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.11.12.orig$polygons)[1],
           xmax = st_bbox(spe.11.12.orig$polygons)[3],
           ymin = st_bbox(spe.11.12.orig$polygons)[2],
           ymax = st_bbox(spe.11.12.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
## add the scale legend
p1 <- plotScaleBar(p1, spe.11.12.orig)

### save the QS figure ###
pdf(file.path("Figure1C.pdf"), 
    width = 10, 
    height = 10,
    bg = "transparent")
p1
dev.off()

# Figure 1D

### Coordinates of the insert
xmin <- 3456
xmax <- 4256
ymin <- 750
ymax <- 1750
fov.size <- 4256
if (length(fov.11) == 1){
  x.global <- metadata(spe.orig)$fov_positions[metadata(spe.orig)$fov_positions$fov%in%fov.11,"x_global_px"]
  y.global <- metadata(spe.orig)$fov_positions[metadata(spe.orig)$fov_positions$fov%in%fov.11,"y_global_px"]
  insert.x <- c(xmin,xmax)+x.global
  insert.y <- y.global + fov.size -c(ymin,ymax)
}

### Area/Volume
var.to.plot <- "Area_um"
spe.11$polygons[,var.to.plot] <- spe.11[,var.to.plot]
colLabel <- "Area"
p1 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
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
           xmin = st_bbox(spe.11.orig$polygons)[1],
           xmax = st_bbox(spe.11.orig$polygons)[3],
           ymin = st_bbox(spe.11.orig$polygons)[2],
           ymax = st_bbox(spe.11.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 1)

### insert area/volume  
p2 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot)),
          color = "grey90",
          lwd = 0.2)+
  scale_fill_viridis_c(option="rocket", direction = -1) +
  #scale_color_viridis_c(option="rocket", direction = -1, guide = "none") +
  labs(title = colLabel,fill = NULL) +
  coord_sf(xlim = c(insert.x[1], insert.x[2]),
           ylim = c(insert.y[1], insert.y[2]),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  theme(panel.border=element_blank(),
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

### insert QS
var.to.plot <- "QC_score"
spe.11$polygons[,var.to.plot] <- spe.11[,var.to.plot]
colLabel <- "Quality score"
p3 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot)), 
          color = "grey90",
          lwd = 0.2)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  coord_sf(xlim = c(insert.x[1], insert.x[2]),
           ylim = c(insert.y[1], insert.y[2]),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  theme(panel.border=element_blank(),
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

### create the final plot with the legends
p1 <- plotScaleBar(p1, spe.11.orig)

# --- Build rows ---
row1 <- p1 
row2 <- p2 | p3 
# --- Add equal spacing between rows ---
final_plot <- (row1 / plot_spacer() /
                 row2) +
  plot_layout(
    heights = c(1.5, 0.05, 0.7),
    widths = c(1, 0.01)
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

### save the final figure ###
pdf(file.path("Figure1D.pdf"), 
    width = 10, 
    height = 6,
    bg = "transparent")
final_plot
dev.off()

# Figure 1E

### Coordinates of the insert
xmin <- 3456
xmax <- 4256
ymin <- 750
ymax <- 1750
fov.size <- 4256
if (length(fov.11) == 1){
  x.global <- metadata(spe.orig)$fov_positions[metadata(spe.orig)$fov_positions$fov%in%fov.11,"x_global_px"]
  y.global <- metadata(spe.orig)$fov_positions[metadata(spe.orig)$fov_positions$fov%in%fov.11,"y_global_px"]
  insert.x <- c(xmin,xmax)+x.global
  insert.y <- y.global + fov.size -c(ymin,ymax)
}

### Fast Reseg flag
var.to.plot <- "fastreseg_flag"
spe.11$polygons[,var.to.plot] <- NA
sub <- FastRes[grep(paste0("f",fov.11),FastRes$cell_ID),]
spe.11$polygons[sub$cell_ID,var.to.plot] <- sub$fastreseg_flag
colLabel <- "FastReseg"
p1 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot), 
                                 color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_manual(values = c("TRUE"  = "red",
                               "FALSE" = "grey95"),
                    breaks = c("TRUE", "FALSE"),
                    na.value = "transparent") +
  scale_color_manual(values = c("FALSE"  = "#B3C2F2",
                                "TRUE" = "#B3C2F2"),
                     breaks = c("FALSE", "TRUE"),
                     na.value = "transparent",
                     guide = "none") +
  labs(title = colLabel,fill  = NULL) +    
  my.theme + 
  theme(legend.position = "bottom",
        legend.justification = c(0.95, 0.5),
        legend.margin = margin(t = -15)) +
  annotate("rect",
           xmin = st_bbox(spe.11.orig$polygons)[1],
           xmax = st_bbox(spe.11.orig$polygons)[3],
           ymin = st_bbox(spe.11.orig$polygons)[2],
           ymax = st_bbox(spe.11.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

## Fast Reseg flag insert
p2 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot), 
                                 color = !!sym(var.to.plot)), 
          lwd = 0.1)+
  scale_fill_manual(values = c("TRUE"  = "red",
                               "FALSE" = "grey95"),
                    breaks = c("TRUE", "FALSE"),
                    na.value = "transparent") +
  scale_color_manual(values = c("FALSE"  = "#B3C2F2",
                                "TRUE" = "#B3C2F2"),
                     breaks = c("FALSE", "TRUE"),
                     na.value = "transparent",
                     guide = "none") +
  geom_point(data = sub.tx, # plotting nuclear transcripts
             aes(x = x_global_px, y = y_global_px),
             size = 0.5,
             color = "dodgerblue") +
  coord_sf(xlim = c(insert.x[1], insert.x[2]),
           ylim = c(insert.y[1], insert.y[2]),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  theme(panel.border=element_blank(),
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
                    legend.position = "none",
                    legend.justification = c(0.95, 0.5),
                    legend.box.just = "right",
                    legend.margin = margin(t = -15)  ,
                    panel.grid=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_blank())

### QS
var.to.plot <- "QC_score"
spe.11$polygons[,var.to.plot] <- spe.11[,var.to.plot]
colLabel <- "Quality score"
p3 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(spe.11.orig$polygons)[1],
           xmax = st_bbox(spe.11.orig$polygons)[3],
           ymin = st_bbox(spe.11.orig$polygons)[2],
           ymax = st_bbox(spe.11.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

## QS insert
p4 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", guide = "none") +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot)), 
          color = "#B3C2F2",
          lwd = 0.2)+ 
  geom_point(data = sub.tx, # plotting nuclear transcripts
             aes(x = x_global_px, y = y_global_px),
             size = 0.4,
             color = "dodgerblue") +
  coord_sf(xlim = c(insert.x[1], insert.x[2]),
           ylim = c(insert.y[1], insert.y[2]),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  theme(panel.border=element_blank(),
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
        legend.position = "none",
        legend.justification = c(0.95, 0.5),
        legend.box.just = "right",
        legend.margin = margin(t = -15)  ,
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

### create the final plot with the legends
p1 <- plotScaleBar(p1, spe.11.orig)
p3 <- plotScaleBar(p3, spe.11.orig)

# --- Build rows ---
row1 <- p1 | p3
row2 <- p2 | p4
# --- Add equal spacing between rows ---
final_plot <- (row1 / plot_spacer() /
                 row2) +
  plot_layout(
    heights = c(1.5, 0.05, 1),
    widths = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

### save the final figure ###
pdf(file.path("Figure1E.pdf"), 
    width = 10, 
    height = 6,
    bg = "white")
final_plot
dev.off()

# Supplementary figure 8

label_coord <- metadata(spe.orig)$fov_positions
label_coord$y_global_px <- label_coord$y_global_px+4256

### Quality score
fov.xdim <- 4256
fov.ydim <- 4256
offset <- 700
var.to.plot <- "QC_score"
spe$polygons[,var.to.plot] <- spe[,var.to.plot]
colLabel <- "Quality score"
p1 <- ggplot() +
  geom_sf(data=spe$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", guide = "none") +
  geom_text(aes(x=label_coord$x_global_px+offset,
                y=label_coord$y_global_px-offset,
                label = label_coord$fov), color="black", fontface = "bold",
            size= 2.5) +
  labs(title = NULL,fill = NULL) +
  theme(panel.border=element_blank(),
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
        panel.grid.major=element_blank()) +
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
           linewidth = 0.05)
## add the scale legend
p1 <- plotScaleBar(p1, spe.orig)

### save the QS figure ###
pdf(file.path("SuppFigure8.pdf"), 
    width = 10, 
    height = 10,
    bg = "transparent")
p1
dev.off()

# Supplementary figure 9

### Area/Volume
var.to.plot <- "Area_um"
spe.11.12$polygons[,var.to.plot] <- spe.11.12[,var.to.plot]
colLabel <- "Area"
p1 <- ggplot() +
  geom_sf(data=spe.11.12$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="rocket", direction = -1) +
  scale_color_viridis_c(option="rocket", direction = -1, guide = "none") +
  labs(title = colLabel,fill = NULL) +
  theme(panel.border=element_blank(),
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
                    panel.grid.major=element_blank()) +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.11.12.orig$polygons)[1],
           xmax = st_bbox(spe.11.12.orig$polygons)[3],
           ymin = st_bbox(spe.11.12.orig$polygons)[2],
           ymax = st_bbox(spe.11.12.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
## add the scale legend
p1 <- plotScaleBar(p1, spe.11.12.orig)

### Signal density
var.to.plot <- "log2SignalDensity"
spe.11.12$polygons[,var.to.plot] <- spe.11.12[,var.to.plot]
colLabel <- "Signal density"
p2 <- ggplot() +
  geom_sf(data=spe.11.12$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_viridis_c(option="mako") +
  scale_color_viridis_c(option="mako", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  theme(panel.border=element_blank(),
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
        panel.grid.major=element_blank()) +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.11.12.orig$polygons)[1],
           xmax = st_bbox(spe.11.12.orig$polygons)[3],
           ymin = st_bbox(spe.11.12.orig$polygons)[2],
           ymax = st_bbox(spe.11.12.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
## add the scale legend
p2 <- plotScaleBar(p2, spe.11.12.orig)

### Proportion of control probe counts
var.to.plot <- "log2Ctrl_total_ratio"
spe.11.12$polygons[,var.to.plot] <- spe.11.12[,var.to.plot]
colLabel <- "Proportion of control probe counts"
p3 <- ggplot() +
  geom_sf(data=spe.11.12$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_continuous_sequential(palette = "BuPu", begin = 0.15, breaks = scales::breaks_pretty(n = 5) ) + 
  scale_color_continuous_sequential(palette = "BuPu", begin = 0.15, guide = "none") + 
  labs(title = colLabel,fill = NULL) +
  theme(panel.border=element_blank(),
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
        panel.grid.major=element_blank()) +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.11.12.orig$polygons)[1],
           xmax = st_bbox(spe.11.12.orig$polygons)[3],
           ymin = st_bbox(spe.11.12.orig$polygons)[2],
           ymax = st_bbox(spe.11.12.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
## add the scale legend
p3 <- plotScaleBar(p3, spe.11.12.orig)

### Border effect
var.to.plot <- "abs_log2AspectRatio_dist"
spe.11.12$polygons[,var.to.plot] <- ifelse(spe.11.12$dist_border<50, abs(spe.11.12$log2AspectRatio), NA)
colLabel <- "Aspect ratio"
p4 <- ggplot() +
  geom_sf(data=spe.11.12$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_gradientn(colors = heat.colors(50, rev = TRUE), breaks = scales::breaks_pretty(n = 5), na.value = "grey90") +
  scale_color_gradientn(colors = heat.colors(50, rev = TRUE), guide = "none") +
  labs(title = colLabel,fill = NULL) +
  theme(panel.border=element_blank(),
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
        panel.grid.major=element_blank()) +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.11.12.orig$polygons)[1],
           xmax = st_bbox(spe.11.12.orig$polygons)[3],
           ymin = st_bbox(spe.11.12.orig$polygons)[2],
           ymax = st_bbox(spe.11.12.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
## add the scale legend
p4 <- plotScaleBar(p4, spe.11.12.orig)

# --- Build rows ---
row1 <- p1 | p2
row2 <- p3 | p4

# --- Add equal spacing between rows ---
final_plot <- (row1 / plot_spacer() /
                 row2) +
  plot_layout(
    heights = c(1.5, 0.05, 1),
    widths = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

### save the final figure ###
pdf(file.path("SuppFigure9.pdf"), 
    width = 16, 
    height = 10,
    bg = "transparent")
final_plot
dev.off()

# Supplementary figure 10A

### Coordinates of the insert
xmin <- 136550
xmax <- 138300
ymin <- 43500
ymax <- 44900
insert.x <- c(xmin,xmax)
insert.y <- c(ymin,ymax)

### Density
var.to.plot <- "log2SignalDensity"
spe.11$polygons[,var.to.plot] <- spe.11[,var.to.plot]
colLabel <- "Signal density"

p1 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="mako", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(spe.11.orig$polygons)[1],
           xmax = st_bbox(spe.11.orig$polygons)[3],
           ymin = st_bbox(spe.11.orig$polygons)[2],
           ymax = st_bbox(spe.11.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05) 
p1 <- plotScaleBar(p1, spe.11.orig)

p2 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="mako", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(spe.11.orig$polygons)[1],
           xmax = st_bbox(spe.11.orig$polygons)[3],
           ymin = st_bbox(spe.11.orig$polygons)[2],
           ymax = st_bbox(spe.11.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05) +
  annotate("rect",
           xmin = xmin,
           xmax = xmax,
           ymin = ymin,
           ymax = ymax,
           fill = NA,
           color = "black",
           linewidth = 1)

p2 <- plotScaleBar(p2, spe.11.orig)

## Density insert with non-nuclear transcripts
p3 <- ggplot() + 
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2) +
  scale_fill_viridis_c(option="mako", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako", guide = "none") +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot)), 
          color = "#B3C2F2",
          lwd = 0.2)+
  coord_sf(xlim = c(xmin, xmax),
           ylim = c(ymin, ymax),
           expand = FALSE) +
  geom_point_rast(data = no.nucl.tx,
             aes(x = x_global_px, y = y_global_px),
             size = 0.01,
             color = "firebrick",
             raster.dpi=1000) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  theme(panel.border=element_blank(),
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
        legend.position = "none",
        legend.justification = c(0.95, 0.5),
        legend.box.just = "right",
        legend.margin = margin(t = -15),
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

## QS insert 
var.to.plot <- "QC_score"
spe$polygons[,var.to.plot] <- spe[,var.to.plot]
colLabel <- "Quality score"
p4 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", guide = "none") +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot)), 
          color = "#B3C2F2",
          lwd = 0.2) +
  coord_sf(xlim = c(xmin, xmax),
           ylim = c(ymin, ymax),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  theme(panel.border=element_blank(),
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
        legend.justification = c(0.95, 0.5),
        legend.box.just = "right",
        legend.margin = margin(t = -15),
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

### Proportion of control probe counts with only control probes
var.to.plot <- "log2Ctrl_total_ratio"
spe.11$polygons[,var.to.plot] <- spe.11[,var.to.plot]
colLabel <- "Proportion of control probe counts"
p5 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_continuous_sequential(palette = "BuPu", begin = 0.15, breaks = scales::breaks_pretty(n = 5) ) + 
  scale_color_continuous_sequential(palette = "BuPu", begin = 0.15, guide = "none") + 
  geom_point_rast(data = ctrl.tx[ctrl.tx$fov==11,],
             aes(x = x_global_px, y = y_global_px),
             size = 0.05,
             color = "hotpink", raster.dpi=1000) +
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
           xmin = st_bbox(spe.11.orig$polygons)[1],
           xmax = st_bbox(spe.11.orig$polygons)[3],
           ymin = st_bbox(spe.11.orig$polygons)[2],
           ymax = st_bbox(spe.11.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05) + 
  annotate("rect",
           xmin = xmin,
           xmax = xmax,
           ymin = ymin,
           ymax = ymax,
           fill = NA,
           color = "black",
           linewidth = 1)

p5 <- plotScaleBar(p5, spe.11.orig)

## Proportion of control probe counts in insert with only control probes
p6 <- ggplot() +
  geom_sf(data=spe.11$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_continuous_sequential(palette = "BuPu", begin = 0.15, breaks = scales::breaks_pretty(n = 5) ) + 
  scale_color_continuous_sequential(palette = "BuPu", begin = 0.15, guide = "none") + 
  geom_point_rast(data = ctrl.tx[ctrl.tx$fov==11,],
             aes(x = x_global_px, y = y_global_px),
             size = 0.05,
             color = "hotpink",
             raster.dpi=1000) +
  labs(title = colLabel,fill = NULL) +
  coord_sf(xlim = c(xmin, xmax),
           ylim = c(ymin, ymax),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  theme(panel.border=element_blank(),
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
        legend.position = "none",
        legend.justification = c(0.95, 0.5),
        legend.box.just = "right",
        legend.margin = margin(t = -15),
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

### create the final plot with the legends
row1 <- p1 | p3
row2 <- p2 | p4
row3 <- p5 | p6

# --- Assemble with spacers for equal vertical spacing ---
final_plot <- (row1 / plot_spacer() /
                 row2 / plot_spacer() /
                 row3) +
  plot_layout(
    heights = c(1, 0.05, 1, 0.05, 1),
    widths  = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "."
  ) & 
  theme(
    plot.tag = element_text(size = 14, face = "bold")
  )

### Save the final figure ----
pdf(file.path("FigureSupp10A.pdf"),
  width  = 10,
  height = 10,
  bg     = "white"
)
final_plot
dev.off()

# Supplementary figure 10B

### Coordinates of the insert
xmin <- 540
xmax <- 2104
ymin <- 803
ymax <- 2367
fov.size <- 4256
if (length(fov.31) == 1){
  x.global <- metadata(spe.orig)$fov_positions[metadata(spe.orig)$fov_positions$fov%in%fov.31,"x_global_px"]
  y.global <- metadata(spe.orig)$fov_positions[metadata(spe.orig)$fov_positions$fov%in%fov.31,"y_global_px"]
  insert.x <- c(xmin,xmax)+x.global
  insert.y <- y.global + fov.size -c(ymin,ymax)
}

### Density
var.to.plot <- "log2SignalDensity"
spe.31$polygons[,var.to.plot] <- spe.31[,var.to.plot]
colLabel <- "Signal density"
p1 <- ggplot() +
  geom_sf(data=spe.31$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="mako", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="mako", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(spe.31.orig$polygons)[1],
           xmax = st_bbox(spe.31.orig$polygons)[3],
           ymin = st_bbox(spe.31.orig$polygons)[2],
           ymax = st_bbox(spe.31.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 1)
  
## Density insert 
p2 <- ggplot() + 
  geom_sf(data=spe.31$polygons, aes(fill=!!sym(var.to.plot)), 
          color = "#B3C2F2",
          lwd = 0.2)+
  coord_sf(xlim = c(insert.x[1], insert.x[2]),
           ylim = c(insert.y[1], insert.y[2]),
           expand = FALSE) +
  scale_fill_viridis_c(option="mako", breaks = scales::breaks_pretty(n = 5)) +
  coord_sf(xlim = c(insert.x[1], insert.x[2]),
           ylim = c(insert.y[1], insert.y[2]),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  theme(panel.border=element_blank(),
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
        legend.title = element_blank(),
        legend.text=element_text(color="black", size = 6),
        legend.position = "right",
        legend.justification = c(0.95, 0.5),
        legend.box.just = "right",
        legend.margin = margin(t = -15)  ,
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) 

### Quality score
var.to.plot <- "QC_score"
spe.31$polygons[,var.to.plot] <- spe.31[,var.to.plot]
colLabel <- "Quality score"
p3 <- ggplot() +
  geom_sf(data=spe.31$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(spe.31.orig$polygons)[1],
           xmax = st_bbox(spe.31.orig$polygons)[3],
           ymin = st_bbox(spe.31.orig$polygons)[2],
           ymax = st_bbox(spe.31.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

## QS insert 
p4 <- p3 + 
  geom_sf(data=spe.31$polygons, aes(fill=!!sym(var.to.plot)), 
          color = "#B3C2F2",
          lwd = 0.2)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  #scale_color_viridis_c(option="plasma", guide = "none") +
  coord_sf(xlim = c(insert.x[1], insert.x[2]),
           ylim = c(insert.y[1], insert.y[2]),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  theme(panel.border=element_blank(),
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

### Cell types
var.to.plot <- "InSituType_Simple"
spe.31$polygons[,var.to.plot] <- as.factor(spe.31[,var.to.plot])
colLabel <- "Cell types"
pp.ct <- ggplot() +
  geom_sf(data=spe.31$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_manual(values = celltype_palette,
                    name = colLabel) +
  scale_color_manual(values = celltype_palette,
                     guide = "none") +
  labs(title = colLabel, fill = NULL) +
  my.theme +
  annotate("rect",
           xmin = st_bbox(spe.31.orig$polygons)[1],
           xmax = st_bbox(spe.31.orig$polygons)[3],
           ymin = st_bbox(spe.31.orig$polygons)[2],
           ymax = st_bbox(spe.31.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

## Cell types insert
p5 <- pp.ct + 
  geom_sf(data=spe.31$polygons, aes(fill=!!sym(var.to.plot)), 
          color = "#B3C2F2",
          lwd = 0.2)+
  scale_fill_manual(values = celltype_palette,
                    name = colLabel) +
  coord_sf(xlim = c(insert.x[1], insert.x[2]),
           ylim = c(insert.y[1], insert.y[2]),
           expand = FALSE) +
  annotate("rect",
           xmin = insert.x[1],
           xmax = insert.x[2],
           ymin = insert.y[1],
           ymax = insert.y[2],
           fill = NA,
           color = "black",
           linewidth = 0.2) +
  theme(panel.border=element_blank(),
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

### create the final plot with the legends
row1 <- p1 | p3
row2 <- p2 | p4
row3 <- p5

# --- Assemble with spacers for equal vertical spacing ---
final_plot <- (row1 / plot_spacer() /
                 row2 / plot_spacer() /
                 row3) +
  plot_layout(
    heights = c(1, 0.05, 1, 0.05, 1),
    widths  = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "."
  ) & 
  theme(
    plot.tag = element_text(size = 14, face = "bold")
  )

### Save the final figure ----
pdf(file.path("FigureSupp10B.pdf"),
  width  = 10,
  height = 10,
  bg     = "white"
)
final_plot
dev.off()

# Supplementary figure 11

### Insert
xmin <- 164500 
xmax <- 165300
ymin <- -38350
ymax <- -37800

### Border effect
var.to.plot <- "abs_log2AspectRatio_dist"
spe.38.39$polygons[,var.to.plot] <- ifelse(spe.38.39$dist_border<50, abs(spe.38.39$log2AspectRatio), NA)
colLabel <- "Aspect ratio"
p1 <- ggplot() +
  geom_sf(data=spe.38.39$polygons, aes(fill=!!sym(var.to.plot), color = !!sym(var.to.plot)), lwd = 0.1)+
  scale_fill_gradientn(colors = heat.colors(50, rev = TRUE), breaks = scales::breaks_pretty(n = 5), na.value = "grey90") +
  scale_color_gradientn(colors = heat.colors(50, rev = TRUE), guide = "none") +
  labs(title = colLabel,fill = NULL) +
  theme(panel.border=element_blank(),
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
        panel.grid.major=element_blank())  +
  guides(fill = guide_colorbar(title.position = "top", 
                               title.hjust = 0.5,
                               barwidth = 0.5,
                               barheight = 4,
                               frame.colour = "black",
                               frame.linewidth = 0.02,
                               ticks = element_line(linewidth = 0.2))) +
  annotate("rect",
           xmin = st_bbox(spe.38.39.orig$polygons)[1],
           xmax = st_bbox(spe.38.39.orig$polygons)[3],
           ymin = st_bbox(spe.38.39.orig$polygons)[2],
           ymax = st_bbox(spe.38.39.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)
p2 <- p1 +
  annotate("rect",
           xmin = xmin,
           xmax = xmax,
           ymin = ymin,
           ymax = ymax,
           fill = NA,
           color = "black",
           linewidth = 1)

p2 <- plotScaleBar(p2, spe.38.39.orig)

### Border effect on insert 
p3 <- p1 + 
  geom_sf(data=spe.38.39$polygons, aes(fill=!!sym(var.to.plot)), 
          color = "#B3C2F2",
          lwd = 0.2)+
  coord_sf(xlim = c(xmin, xmax),
           ylim = c(ymin, ymax),
           expand = FALSE) +
  theme(panel.border=element_blank(),
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
        legend.position = "none",
        legend.justification = c(0.95, 0.5),
        legend.box.just = "right",
        legend.margin = margin(t = -15)  ,
        panel.grid=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) 

### Quality score
var.to.plot <- "QC_score"
spe.38.39$polygons[,var.to.plot] <- spe.38.39[,var.to.plot]
colLabel <- "Quality score"
pp <- ggplot() +
  geom_sf(data=spe.38.39$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  scale_color_viridis_c(option="plasma", guide = "none") +
  labs(title = colLabel,fill = NULL) +
  theme(panel.border=element_blank(),
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
        panel.grid.major=element_blank())  +
  annotate("rect",
           xmin = st_bbox(spe.38.39.orig$polygons)[1],
           xmax = st_bbox(spe.38.39.orig$polygons)[3],
           ymin = st_bbox(spe.38.39.orig$polygons)[2],
           ymax = st_bbox(spe.38.39.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

pp <- plotScaleBar(pp, spe.38.39.orig)

## QS insert 
p4 <- pp + 
  geom_sf(data=spe.38.39$polygons, aes(fill=!!sym(var.to.plot)), 
          color = "#B3C2F2",
          lwd = 0.2)+
  scale_fill_viridis_c(option="plasma", breaks = scales::breaks_pretty(n = 5)) +
  #scale_color_viridis_c(option="plasma", guide = "none") +
  coord_sf(xlim = c(xmin, xmax),
           ylim = c(ymin, ymax),
           expand = FALSE) +
  theme(panel.border=element_blank(),
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

### Cell types
var.to.plot <- "InSituType_Simple"
spe.38.39$polygons[,var.to.plot] <- as.factor(spe.38.39[,var.to.plot])
colLabel <- "Cell types"
pp.2 <- ggplot() +
  geom_sf(data=spe.38.39$polygons, aes(fill=!!sym(var.to.plot),
                                 color =!!sym(var.to.plot)), 
          lwd = 0.2)+
  scale_fill_manual(values = celltype_palette,
                    name = colLabel) +
  scale_color_manual(values = celltype_palette,
                     guide = "none") +
  labs(title = colLabel, fill = NULL) +
  theme(panel.border=element_blank(),
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
        panel.grid.major=element_blank()) +
  annotate("rect",
           xmin = st_bbox(spe.38.39.orig$polygons)[1],
           xmax = st_bbox(spe.38.39.orig$polygons)[3],
           ymin = st_bbox(spe.38.39.orig$polygons)[2],
           ymax = st_bbox(spe.38.39.orig$polygons)[4],
           fill = NA,
           color = "black",
           linewidth = 0.05)

pp.2 <- plotScaleBar(pp.2, spe.38.39.orig)

## Cell types insert
p6 <- pp.2 + 
  geom_sf(data=spe.38.39$polygons, aes(fill=!!sym(var.to.plot)), 
          color = "#B3C2F2",
          lwd = 0.2)+
  scale_fill_manual(values = celltype_palette,
                    name = colLabel) +
  coord_sf(xlim = c(xmin, xmax),
           ylim = c(ymin, ymax),
           expand = FALSE) +
  theme(panel.border=element_blank(),
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

# --- Build rows ---
row1 <- p2 | pp | pp.2
row2 <- p3 | p4 | p6

# --- Add equal spacing between rows ---
final_plot <- (row1 / plot_spacer() /
                 row2) +
  plot_layout(
    heights = c(2, 0.05, 0.7),
    widths = c(1, 1, 1)
  ) 

### save the final figure
pdf(file.path("SuppFigure11.pdf"), 
    width = 18, 
    height = 10,
    bg = "transparent")
final_plot
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
#   [1] ggrastr_1.0.2               colorspace_2.1-2            data.table_1.18.0           patchwork_1.3.2            
# [5] sf_1.0-24                   dplyr_1.1.4                 ggplot2_4.0.1               SpaceTrooper_1.1.3         
# [9] SpatialExperiment_1.20.0    SingleCellExperiment_1.32.0 SummarizedExperiment_1.40.0 Biobase_2.70.0             
# [13] GenomicRanges_1.62.1        Seqinfo_1.0.0               IRanges_2.44.0              S4Vectors_0.48.0           
# [17] BiocGenerics_0.56.0         generics_0.1.4              MatrixGenerics_1.22.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.3                 gridExtra_2.3             rlang_1.1.7               magrittr_2.0.4            scater_1.38.0            
# [6] e1071_1.7-17              compiler_4.5.1            DelayedMatrixStats_1.32.0 SpatialExperimentIO_1.2.0 sfheaders_0.4.5          
# [11] vctrs_0.7.1               pkgconfig_2.0.3           shape_1.4.6.1             backports_1.5.0           magick_2.9.0             
# [16] XVector_0.50.0            scuttle_1.20.0            ggbeeswarm_0.7.3          purrr_1.2.1               bit_4.6.0                
# [21] glmnet_4.1-10             beachmat_2.26.0           rhdf5filters_1.22.0       DelayedArray_0.36.0       Rhdf5lib_1.32.0          
# [26] BiocParallel_1.44.0       broom_1.0.11              irlba_2.3.5.1             parallel_4.5.1            R6_2.6.1                 
# [31] RColorBrewer_1.1-3        limma_3.66.0              car_3.1-3                 Rcpp_1.1.1                assertthat_0.2.1         
# [36] iterators_1.0.14          R.utils_2.13.0            Matrix_1.7-4              splines_4.5.1             tidyselect_1.2.1         
# [41] viridis_0.6.5             rstudioapi_0.17.1         dichromat_2.0-0.1         abind_1.4-8               codetools_0.2-20         
# [46] lattice_0.22-7            tibble_3.3.1              withr_3.0.2               S7_0.2.1                  survival_3.8-3           
# [51] units_1.0-0               proxy_0.4-29              pillar_1.11.1             BiocManager_1.30.27       ggpubr_0.6.2             
# [56] carData_3.0-5             KernSmooth_2.23-26        foreach_1.5.2             sparseMatrixStats_1.22.0  scales_1.4.0             
# [61] class_7.3-23              glue_1.8.0                tools_4.5.1               BiocNeighbors_2.4.0       robustbase_0.99-6        
# [66] ScaledMatrix_1.18.0       locfit_1.5-9.12           ggsignif_0.6.4            Cairo_1.7-0               cowplot_1.2.0            
# [71] rhdf5_2.54.1              grid_4.5.1                tidyr_1.3.2               DropletUtils_1.30.0       edgeR_4.8.2              
# [76] beeswarm_0.4.0            BiocSingular_1.26.1       HDF5Array_1.38.0          vipor_0.4.7               Formula_1.2-5            
# [81] cli_3.6.5                 rsvd_1.0.5                viridisLite_0.4.2         S4Arrays_1.10.1           arrow_23.0.0             
# [86] gtable_0.3.6              DEoptimR_1.1-4            R.methodsS3_1.8.2         rstatix_0.7.3             classInt_0.4-11          
# [91] ggrepel_0.9.6             SparseArray_1.10.8        dqrng_0.4.1               rjson_0.2.23              farver_2.1.2             
# [96] R.oo_1.27.1               lifecycle_1.0.5           h5mread_1.2.1             statmod_1.5.1             bit64_4.6.0-1
