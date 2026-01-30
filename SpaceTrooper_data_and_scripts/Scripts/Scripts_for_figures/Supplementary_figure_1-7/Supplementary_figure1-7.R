# Clean, minimal, and reproducible script to reproduce supplementary figures 1-7
# To reproduce the analysis be sure to have the required packages installed

library(ggplot2)
library(patchwork)
library(dplyr)
library(robustbase)

# R = 4.5.1

# Specify spatial transcriptomics datasets to load metadata files available at: 
# bicciatolab_data/SpaceTrooper_data_and_scripts/Spe_metadata"

metadata_dir <- "bicciatolab_data/SpaceTrooper_data_and_scripts/Spe_metadata"

datasets <- c("CosMx_rna_DBKERO",
              "Xenium_lung_cancer",
              "MERFISH_mouse_liver")

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

# Supplementary figure 1

### read metadata table and create the plots ###
pair_matrix <- matrix(1:6, ncol = 2, byrow = TRUE)

for (k in 1:length(datasets)){
  metadata <- readRDS(file.path(metadata_dir, paste0(datasets[k],"_metadata.rds")))
  ### Area/Volume
  var.to.plot <- "Area_um"
  xlabel <- expression("Area ("*µm^2*")")
  x.bins <- 100
  file.name <- "Area"
  if (datasets[k] == "MERFISH_mouse_liver") {
    var.to.plot <- "volume"
    xlabel <- expression("Volume ("*µm^3*")")
    x.bins <- 1000
    file.name <- "Volume"
  }
  ## histogram non logged axes
  pp <- ggplot(metadata, aes(!!sym(var.to.plot))) + 
    geom_histogram(bins = 50, 
                   fill = "#B3C2F2",
                   color = "grey40",
                   linewidth = 0.1) + 
    geom_vline(xintercept=median(metadata[,var.to.plot]), 
               col="black",
               linewidth = 0.8) + 
    geom_vline(xintercept=median(metadata[,var.to.plot]) +3*mad(metadata[,var.to.plot]),
               col="tomato",
               linewidth = 0.8) + 
    geom_vline(xintercept=adjbox(metadata[,var.to.plot], plot = FALSE)$fence[2,1],
               col="firebrick",
               linewidth = 0.8) + 
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                       expand = c(0,0)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
                       expand = expansion(mult = c(0, 0.05))) +
    ylab("Number of cells") + 
    xlab(xlabel) + 
    theme_minimal() + 
    light_theme +
    theme(plot.margin = unit(c(0,20,0,0), "pt"))
  assign(paste0("p",pair_matrix[k,1]),pp)
  ## scatter against total
  metadata <- metadata %>%
    mutate(area_mad_color = case_when(
      !!sym(var.to.plot) > adjbox(!!sym(var.to.plot), plot = FALSE)$fence[2] ~ "HIGH_MC",
      !!sym(var.to.plot) > median(!!sym(var.to.plot)) + 3 * mad(!!sym(var.to.plot)) ~ "HIGH_MAD",
      TRUE ~ "NO"
    ))
  pp <- ggplot(metadata, aes(x = !!sym(var.to.plot), y = total)) + 
    ggrastr::geom_point_rast(aes(color = area_mad_color, 
                                 size = area_mad_color), 
                             shape = 16, 
                             raster.dpi = 1000) + 
    scale_color_manual(values = c("HIGH_MAD" = "tomato", 
                                  "HIGH_MC" = "firebrick",
                                  "NO" = "#c0c8cf")) + 
    scale_size_manual(values = c("HIGH_MAD" = 0.6, 
                                 "HIGH_MC" = 0.6, 
                                 "NO" = 0)) +
    geom_smooth(col = "blue", 
                se = FALSE,
                linewidth = 0.8) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 8),
      expand = expansion(mult = c(0, 0.05))) + 
    scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
                       expand = expansion(mult = c(0, 0.05))) +
    xlab(xlabel) + 
    ylab("Total probe counts") +
    theme_minimal() + 
    theme(aspect.ratio = 1) +
    light_theme +
    theme(plot.margin = unit(c(0,0,0,20), "pt"))
  assign(paste0("p",pair_matrix[k,2]),pp)
}

# --- Build rows ---
row1 <- p1 | p2
row2 <- p3 | p4
row3 <- p5 | p6
# --- Add equal spacing between rows ---
final_plot <- (row1 / plot_spacer() /
                 row2 / plot_spacer() /
                 row3) +
  plot_layout(
    heights = c(1, 0.05, 1, 0.05, 1),
    widths = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

### save the final figure ###
pdf(file.path("SuppFigure_1.pdf"), 
    width = 10, 
    height = 14,
    bg = "transparent")
final_plot
dev.off()

# Supplementary figure 2

### read metadata table and create the plots ###
pair_matrix <- matrix(1:6, ncol = 2, byrow = TRUE)

for (k in 1:length(datasets)){
  metadata <- readRDS(file.path(metadata_dir, paste0(datasets[k],"_metadata.rds")))
  ### Total counts and signal density
  var.to.plot <- "total"
  ## histogram logged axes
  pp <- ggplot(metadata, aes(!!sym(var.to.plot))) + 
    geom_histogram(bins = 50, 
                   fill = "#B3C2F2",
                   color = "grey40",
                   linewidth = 0.1) +
    geom_vline(xintercept=median(metadata[,var.to.plot]), 
               col="black",
               linewidth = 0.8) + 
    scale_x_continuous(trans = "log10",
                       breaks = scales::breaks_log(n = 6)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                       expand = expansion(mult = c(0, 0.05))) +
    ylab("Number of cells") + 
    xlab("Total probe counts") + 
    theme_minimal() + 
    light_theme +
    theme(plot.margin = unit(c(0,20,0,0), "pt"))
  assign(paste0("p",pair_matrix[k,1]),pp)
  
  ## histogram on signal density
  var.to.plot <- "log2CountArea"
  if (datasets[k] == "MERFISH_mouse_liver") {
    var.to.plot <- "log2CountVolume"
  }
  pp <- ggplot(metadata, aes(!!sym(var.to.plot))) + 
    geom_histogram(bins = 50, 
                   fill = "#B3C2F2",
                   color = "grey40",
                   linewidth = 0.1) + 
    geom_vline(xintercept = median(metadata[,var.to.plot]), 
               col="black",
               linewidth = 0.8) + 
    geom_vline(xintercept = median(metadata[,var.to.plot]) - 3*mad(metadata[,var.to.plot]),
               col = "magenta",
               linewidth = 0.8) + 
    geom_vline(xintercept = adjbox(metadata[,var.to.plot], plot = FALSE)$fence[1],
               col = "purple3",
               linewidth = 0.8) + 
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                       expand = expansion(mult = c(0, 0.05))) +
    ylab("Number of cells") + 
    xlab("Signal density") + 
    theme_minimal() + 
    light_theme +
    theme(plot.margin = unit(c(0,0,0,20), "pt"))
  assign(paste0("p",pair_matrix[k,2]),pp)
}

# --- Build rows ---
row1 <- p1 | p2
row2 <- p3 | p4
row3 <- p5 | p6
# --- Add equal spacing between rows ---
final_plot <- (row1 / plot_spacer() /
                 row2 / plot_spacer() /
                 row3) +
  plot_layout(
    heights = c(1, 0.05, 1, 0.05, 1),
    widths = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

### save the final figure ###
pdf(file.path("SuppFigure_2.pdf"), 
    width = 10, 
    height = 14,
    bg = "transparent")
final_plot
dev.off()

# Supplementary figure 3

### read metadata table and create the plots ###

pair_matrix <- matrix(1:6, ncol = 2, byrow = TRUE)

for (k in 1:length(datasets)){
  metadata <- readRDS(file.path(metadata_dir, paste0(datasets[k],"_metadata.rds")))
  ### Control probe counts
  var.to.filter <- "ctrl_total_ratio"
  var.to.plot <- "log2Ctrl_total_ratio"
  filter.metadata <- metadata |> dplyr::filter(!!sym(var.to.filter)!=0)
  text.y <- 100
  if (datasets[k] == "MERFISH_mouse_liver") {
    text.y <- 5000
  }
  ## histogram
  pp <- ggplot(filter.metadata, aes(!!sym(var.to.plot))) + 
    geom_histogram(bins = 50, 
                   fill = "#B3C2F2",
                   color = "grey40",
                   linewidth = 0.1) + 
    geom_vline(xintercept = median(filter.metadata[,var.to.plot]), 
               col = "black",
               linewidth = 0.8) + 
    geom_vline(xintercept = median(filter.metadata[,var.to.plot]) +3*mad(filter.metadata[,var.to.plot]),
               col = "tomato",
               linewidth = 0.8) + 
    geom_vline(xintercept = adjbox(filter.metadata[,var.to.plot],plot = FALSE)$fence[2],
               col = "firebrick",
               linewidth = 0.8) + 
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                       expand = expansion(mult = c(0, 0.05))) +
    ylab("Number of cells") + 
    xlab("Proportion of control probe counts") + 
    theme_minimal() + 
    light_theme +
    theme(plot.margin = unit(c(0,20,0,0), "pt"))
  assign(paste0("p",pair_matrix[k,1]),pp)
  
  ## scatter plot
  if (adjbox(filter.metadata[,var.to.plot],plot = FALSE)$fence[2] > median(filter.metadata[,var.to.plot]) + 3 * mad(filter.metadata[,var.to.plot])){
    filter.metadata <- filter.metadata %>%
      mutate(log2Ctrl_mad_color = case_when(
        !!sym(var.to.plot) > adjbox(!!sym(var.to.plot), plot = FALSE)$fence[2] ~ "HIGH_MC",
        !!sym(var.to.plot) > median(!!sym(var.to.plot)) + 3 * mad(!!sym(var.to.plot)) ~ "HIGH_MAD",
        TRUE ~ "NO"
      ))
  } else {
    filter.metadata <- filter.metadata %>%
      mutate(log2Ctrl_mad_color = case_when(
        !!sym(var.to.plot) > median(!!sym(var.to.plot)) + 3 * mad(!!sym(var.to.plot)) ~ "HIGH_MAD",
        !!sym(var.to.plot) > adjbox(!!sym(var.to.plot), plot = FALSE)$fence[2] ~ "HIGH_MC",
        TRUE ~ "NO"
      ))
  }
  pp <- ggplot(filter.metadata, aes(x = total, y = !!sym(var.to.plot))) + 
    ggrastr::geom_point_rast(aes(color = log2Ctrl_mad_color, 
                                 size = log2Ctrl_mad_color), 
                             shape = 16, 
                             raster.dpi = 1000) + 
    scale_color_manual(values = c("HIGH_MAD" = "tomato", 
                                  "HIGH_MC" = "firebrick",
                                  "NO" = "#c0c8cf")) + 
    scale_size_manual(values = c("HIGH_MAD" = 0.6, 
                                 "HIGH_MC" = 0.6, 
                                 "NO" = 0)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                       expand = expansion(mult = c(0, 0.05))) + 
    scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                       expand = expansion(mult = c(0, 0.05))) +
    xlab("Total probe counts") + 
    ylab("Proportion of control probe counts") +
    theme_minimal() + 
    theme(aspect.ratio = 1) +
    light_theme +
    theme(plot.margin = unit(c(0,0,0,20), "pt"))
  assign(paste0("p",pair_matrix[k,2]),pp)
}

# --- Build rows ---
row1 <- p1 | p2
row2 <- p3 | p4
row3 <- p5 | p6
# --- Add equal spacing between rows ---
final_plot <- (row1 / plot_spacer() /
                 row2 / plot_spacer() /
                 row3) +
  plot_layout(
    heights = c(1, 0.05, 1, 0.05, 1),
    widths = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

### save the final figure ###
pdf(file.path("SuppFigure_3.pdf"), 
    width = 10, 
    height = 14,
    bg = "transparent")
final_plot
dev.off()

# Supplementary figure 4

### switching to CosMx Protein assay dataset
dataset <- "CosMx_protein_tonsil"

### read metadata table and create the plots ###
metadata <- readRDS(file.path(metadata_dir, paste0(dataset,"_metadata.rds")))
# metadata$SignalDensity corresponds to metadata$total
# metadata$log2SignalDensity corresponds to log2(metadata$SignalDensity)
# metadata$total_intensity corresponds to metadata$total*metadata$Area_um

### Area/Volume
var.to.plot <- "Area_um"
xlabel <- expression("Area ("*µm^2*")")
## histogram non logged axes
p1 <- ggplot(metadata, aes(!!sym(var.to.plot))) + 
  geom_histogram(binwidth = 30, 
                 fill = "#B3C2F2",
                 color = "grey40",
                 linewidth = 0.1) + 
  geom_vline(xintercept=median(metadata[,var.to.plot]), 
             col="black",
             linewidth = 0.8) + 
  geom_vline(xintercept=median(metadata[,var.to.plot]) +3*mad(metadata[,var.to.plot]),
             col="tomato",
             linewidth = 0.8) + 
  geom_vline(xintercept=adjbox(metadata[,var.to.plot], plot = FALSE)$fence[2,1],
             col="firebrick",
             linewidth = 0.8) + 
  scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = expansion(mult = c(0, 0.05)),
                     labels = function(x) {
                       sapply(x, function(xi) {
                         if (!is.na(xi) && xi == 0) {
                           "0"
                         } else if (!is.na(xi)) {
                           parse(text = gsub("e\\+?", " %*% 10^", scales::scientific_format()(xi)))
                         } else {
                           NA
                         }
                       })
                     }) +
  ylab("Number of cells") + 
  xlab(xlabel) + 
  theme_minimal() + 
  light_theme +
  theme(plot.margin = unit(c(0,20,0,0), "pt"))

## scatter against total intensity
metadata <- metadata %>%
  mutate(area_mad_color = case_when(
    !!sym(var.to.plot) > adjbox(!!sym(var.to.plot), plot = FALSE)$fence[2] ~ "HIGH_MC",
    !!sym(var.to.plot) > median(!!sym(var.to.plot)) + 3 * mad(!!sym(var.to.plot)) ~ "HIGH_MAD",
    TRUE ~ "NO"
  ))
p2 <- ggplot(metadata, aes(x = !!sym(var.to.plot), y = total_intensity)) + 
  ggrastr::geom_point_rast(aes(color = area_mad_color, 
                               size = area_mad_color), 
                           shape = 16, 
                           raster.dpi = 1000) + 
  scale_color_manual(values = c("HIGH_MAD" = "tomato", 
                                "HIGH_MC" = "firebrick",
                                "NO" = "#c0c8cf")) + 
  scale_size_manual(values = c("HIGH_MAD" = 0.6, 
                               "HIGH_MC" = 0.6, 
                               "NO" = 0)) +
  geom_smooth(col = "blue", 
              se = FALSE,
              linewidth = 0.8) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = expansion(mult = c(0, 0.05))) + 
  scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = expansion(mult = c(0, 0.05)),
                     labels = function(x) {
                       sapply(x, function(xi) {
                         if (!is.na(xi) && xi == 0) {
                           "0"
                         } else if (!is.na(xi)) {
                           parse(text = gsub("e\\+?", " %*% 10^", scales::scientific_format()(xi)))
                         } else {
                           NA
                         }
                       })
                     }) +
  xlab(xlabel) + 
  ylab("Total intensity") +
  theme_minimal() + 
  theme(aspect.ratio = 1) +
  light_theme +
  theme(plot.margin = unit(c(0,0,0,20), "pt"))

### Signal density
## histogram on signal density
var.to.plot <- "log2SignalDensity"
p3 <- ggplot(metadata, aes(!!sym(var.to.plot))) + 
  geom_histogram(bins=50,
                 fill = "#B3C2F2",
                 color = "grey40",
                 linewidth = 0.1) + 
  geom_vline(xintercept = median(metadata[,var.to.plot]), 
             col="black",
             linewidth = 0.8) + 
  geom_vline(xintercept = median(metadata[,var.to.plot]) - 3*mad(metadata[,var.to.plot]),
             col = "magenta",
             linewidth = 0.8) + 
  geom_vline(xintercept = adjbox(metadata[,var.to.plot], plot = FALSE)$fence[1],
             col = "purple3",
             linewidth = 0.8) + 
  scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = expansion(mult = c(0, 0.05))) +
  ylab("Number of cells") + 
  xlab("Signal density") + 
  theme_minimal() + 
  light_theme +
  theme(plot.margin = unit(c(0,20,0,0), "pt"))

### Control probe intensity
var.to.filter <- "ctrl_total_ratio"
var.to.plot <- "log2Ctrl_total_ratio"
filter.metadata <- metadata |> dplyr::filter(!!sym(var.to.filter)!=0)
text.y <- 100
## histogram
p5 <- ggplot(filter.metadata, aes(!!sym(var.to.plot))) + 
  geom_histogram(binwidth = 0.125, 
                 fill = "#B3C2F2",
                 color = "grey40",
                 linewidth = 0.1) + 
  geom_vline(xintercept = median(filter.metadata[,var.to.plot]), 
             col = "black",
             linewidth = 0.8) + 
  geom_vline(xintercept = median(filter.metadata[,var.to.plot]) +3*mad(filter.metadata[,var.to.plot]),
             col = "tomato",
             linewidth = 0.8) + 
  geom_vline(xintercept = adjbox(filter.metadata[,var.to.plot],plot = FALSE)$fence[2],
             col = "firebrick",
             linewidth = 0.8) + 
  scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = expansion(mult = c(0, 0.05))) +
  ylab("Number of cells") + 
  xlab("Proportion of control probe counts") + 
  theme_minimal() + 
  light_theme +
  theme(plot.margin = unit(c(0,20,0,0), "pt"))

## scatter plot
if (adjbox(filter.metadata[,var.to.plot],plot = FALSE)$fence[2] > median(filter.metadata[,var.to.plot]) + 3 * mad(filter.metadata[,var.to.plot])){
  filter.metadata <- filter.metadata %>%
    mutate(log2Ctrl_mad_color = case_when(
      !!sym(var.to.plot) > adjbox(!!sym(var.to.plot), plot = FALSE)$fence[2] ~ "HIGH_MC",
      !!sym(var.to.plot) > median(!!sym(var.to.plot)) + 3 * mad(!!sym(var.to.plot)) ~ "HIGH_MAD",
      TRUE ~ "NO"
    ))
} else {
  filter.metadata <- filter.metadata %>%
    mutate(log2Ctrl_mad_color = case_when(
      !!sym(var.to.plot) > median(!!sym(var.to.plot)) + 3 * mad(!!sym(var.to.plot)) ~ "HIGH_MAD",
      !!sym(var.to.plot) > adjbox(!!sym(var.to.plot), plot = FALSE)$fence[2] ~ "HIGH_MC",
      TRUE ~ "NO"
    ))
}
p6 <- ggplot(filter.metadata, aes(x = total_intensity, y = !!sym(var.to.plot))) + 
  ggrastr::geom_point_rast(aes(color = log2Ctrl_mad_color, 
                               size = log2Ctrl_mad_color), 
                           shape = 16, 
                           raster.dpi = 1000) + 
  scale_color_manual(values = c("HIGH_MAD" = "tomato", 
                                "HIGH_MC" = "firebrick",
                                "NO" = "#c0c8cf")) + 
  scale_size_manual(values = c("HIGH_MAD" = 0.6, 
                               "HIGH_MC" = 0.6, 
                               "NO" = 0)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = expansion(mult = c(0, 0.05)),
                     labels = function(x) {
                       sapply(x, function(xi) {
                         if (!is.na(xi) && xi == 0) {
                           "0"
                         } else if (!is.na(xi)) {
                           parse(text = gsub("e\\+?", " %*% 10^", scales::scientific_format()(xi)))
                         } else {
                           NA
                         }
                       })
                     }) + 
  scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = expansion(mult = c(0, 0.05))) +
  xlab("Total intensity") + 
  ylab("Proportion of control probe counts") +
  theme_minimal() + 
  theme(aspect.ratio = 1) +
  light_theme +
  theme(plot.margin = unit(c(0,0,0,20), "pt"))

### create the final plot
p4 <- NULL
final_plot <- ((p1 | p2) /
                 (p3 | p4) /
                 (p5 | p6)) +
  plot_layout(heights = c(1, 1, 1), widths = c(1, 1))
# --- Build rows ---
row1 <- p1 | p2
row2 <- p3 | p4
row3 <- p5 | p6
# --- Add equal spacing between rows ---
final_plot <- (row1 / plot_spacer() /
                 row2 / plot_spacer() /
                 row3) +
  plot_layout(
    heights = c(1, 0.05, 1, 0.05, 1),
    widths = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A", 
                  tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 14, face = "bold"))


### save the final figure ###
pdf(file.path("SuppFigure_4.pdf"), 
    width = 10, 
    height = 14,
    bg = "transparent")
final_plot
dev.off()

# Supplementary figure 5 and 6

### plotting function   
p <- function(df, x, y) {
  # df <- df %>% sample_n(min(5000, nrow(df)))
  ggplot(df, aes(.data[[x]], .data[[y]])) +
    geom_point(color="grey", size=0.1) +
    geom_smooth() + 
    theme(panel.background = element_rect(fill="white", color=NA),
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
} 

### read the metadata 
dataset <- "CosMx_rna_DBKERO"

### To obtain supplementary figure 6, switch to protein dataset and then run
### all the following functions until the code chunk to produce the pdf of 
### Supplementary figure 6:
dataset <- "CosMx_protein_tonsil"

metadata <- readRDS(file.path(metadata_dir, paste0(dataset,"_metadata.rds")))

### create the plot
vars <- c("total",
          "Area_um",
          "log2AspectRatio",
          "log2AspectRatio")
bord.dist <- c("dist_border",
               "dist_border",
               "dist_border_hor",
               "dist_border_vert")
y.labels <- c("Total probe counts",
              expression("Area ("*µm^2*")"),
              "Aspect ratio",
              "Aspect ratio")
x.labels <- c("Distance from border (px)",
              "Distance from border (px)",
              "Distance from horizontal border (px)",
              "Distance from vertical border (px)")
if (dataset == "CosMx_protein_tonsil"){
  vars <- c("total_intensity",
            "Area_um",
            "log2AspectRatio",
            "log2AspectRatio")
  y.labels <- c("Total intensity",
                expression("Area ("*µm^2*")"),
                "Aspect ratio",
                "Aspect ratio")
}

out.deriv <- matrix(data = NA, nrow = 4, ncol = 3)
plotList <- list()
## plot the various metrics against the distance from the border
for (k in 1:2){
  plotList[[k]] <- p(metadata, bord.dist[k], vars[k]) +
    xlim(c(0,500))
  smt <- ggplot_build(plotList[[k]])$data[[2]]
  smm <- data.frame(x = smt$x,
                    y = smt$y)
  min.y <- min(metadata[metadata[bord.dist[k]]<=500,vars[k]])
  max.y <- max(metadata[metadata[bord.dist[k]]<=500,vars[k]])
  plotList[[k]] <- plotList[[k]] + 
    scale_x_continuous(
      limits = c(0, 500),
      breaks = scales::breaks_pretty(n = 8),
      expand = expansion(mult = c(0, 0), add = c(0, 2.7))) +
    scale_y_continuous(
      limits = c(min.y, max.y),
      breaks = scales::breaks_pretty(n = 6),
      expand = c(0, 0.05)) +
    labs(x = x.labels[k],
         y = y.labels[k])
  ## get the derivative 
  der <- c(NA, diff(smm$y)/diff(smm$x))
  chg.idx <- -2
  if (bord.dist[k] == "dist_border_hor") {chg.idx <- 2}
  der0 <- smm$x[which(diff(sign(der)) == chg.idx)[1]]
  der0_p1 <- smm$x[which(diff(sign(der)) == chg.idx)[1]+1]
  extr.x <- smm[which(smm$y==max(smm$y, na.rm = T)),"x"]
  if (bord.dist[k] == "dist_border_hor") {
    extr.x <- smm[which(smm$y==min(smm$y, na.rm = T)),"x"]
  }
  out.deriv[k,] <- c(extr.x, der0, der0_p1)
}
colnames(out.deriv) <- c("extrm","der0","der0_p1")
# limit.2get <- "extrm" # "der0_p1"
# border.limit <- round(max(out.deriv[,limit.2get],na.rm = T),0)

## add the outlier detection thresholds
k=3
if (adjbox(metadata[,vars[k]],plot = FALSE)$fence[2] > median(metadata[,vars[k]]) + 3 * mad(metadata[,vars[k]])){
  metadata <- metadata %>%
    mutate(aspRatcolor = case_when(
      (!!sym(vars[k]) > adjbox(!!sym(vars[k]), plot = FALSE)$fence[2])&(!!sym(bord.dist[k])<50) ~ "HIGH_MC",
      (!!sym(vars[k]) > median(!!sym(vars[k])) + 3 * mad(!!sym(vars[k])))&(!!sym(bord.dist[k])<50) ~ "HIGH_MAD",
      TRUE ~ "NO"
    ))
} else {
  metadata <- metadata %>%
    mutate(aspRatcolor = case_when(
      (!!sym(vars[k]) > median(!!sym(vars[k])) + 3 * mad(!!sym(vars[k])))&(!!sym(bord.dist[k])<50) ~ "HIGH_MAD",
      (!!sym(vars[k]) > adjbox(!!sym(vars[k]), plot = FALSE)$fence[2])&(!!sym(bord.dist[k])<50) ~ "HIGH_MC",
      TRUE ~ "NO"
    ))
}
plotList[[k]] <- ggplot(metadata, aes(!!sym(bord.dist[k]), !!sym(vars[k]))) + 
  geom_point(aes(color = aspRatcolor, 
                 size = aspRatcolor), size=0.5) +
  geom_smooth() + 
  xlim(c(0,500)) +
  scale_color_manual(values = c("HIGH_MAD" = "tomato", 
                                "HIGH_MC" = "firebrick",
                                "NO" = "#c0c8cf")) + 
  scale_size_manual(values = c("HIGH_MAD" = 0.6, 
                               "HIGH_MC" = 0.6, 
                               "NO" = 0)) +
  geom_hline(yintercept = median(metadata[,vars[k]]) + 3*mad(metadata[,vars[k]]), 
             linewidth = 0.5,
             color = "tomato") +
  geom_hline(yintercept = adjbox(metadata[,vars[k]], plot = FALSE)$fence[2,1], 
             linewidth = 0.5,
             color = "firebrick") +
  theme(panel.background = element_rect(fill="white", color=NA),
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

smt <- ggplot_build(plotList[[k]])$data[[2]]
smm <- data.frame(x = smt$x,
                  y = smt$y)
min.y <- min(metadata[metadata[bord.dist[k]]<=500,vars[k]])
max.y <- max(metadata[metadata[bord.dist[k]]<=500,vars[k]])
plotList[[k]] <- plotList[[k]] + 
  scale_x_continuous(
    limits = c(0, 500),
    breaks = scales::breaks_pretty(n = 8),
    expand = expansion(mult = c(0, 0), add = c(0, 2.7))) +
  scale_y_continuous(
    limits = c(min.y, max.y),
    breaks = scales::breaks_pretty(n = 6),
    expand = c(0, 0.05)) +
  labs(x = x.labels[k],
       y = y.labels[k])

## get the derivative 
der <- c(NA, diff(smm$y)/diff(smm$x))
chg.idx <- -2
if (bord.dist[k] == "dist_border_hor") {chg.idx <- 2}
der0 <- smm$x[which(diff(sign(der)) == chg.idx)[1]]
der0_p1 <- smm$x[which(diff(sign(der)) == chg.idx)[1]+1]
extr.x <- smm[which(smm$y==max(smm$y, na.rm = T)),"x"]
if (bord.dist[k] == "dist_border_hor") {
  extr.x <- smm[which(smm$y==min(smm$y, na.rm = T)),"x"]
}
out.deriv[k,] <- c(extr.x, der0, der0_p1)
assign(paste0("p",k),plotList[[k]])

## add the outlier detection thresholds
k=4
if (adjbox(metadata[,vars[k]],plot = FALSE)$fence[1] > median(metadata[,vars[k]]) - 3 * mad(metadata[,vars[k]])){
  metadata <- metadata %>%
    mutate(aspRatcolor = case_when(
      (!!sym(vars[k]) < median(!!sym(vars[k])) - 3 * mad(!!sym(vars[k])))&(!!sym(bord.dist[k])<50) ~ "HIGH_MAD",
      (!!sym(vars[k]) < adjbox(!!sym(vars[k]), plot = FALSE)$fence[1])&(!!sym(bord.dist[k])<50) ~ "HIGH_MC",
      TRUE ~ "NO"
    ))
} else {
  metadata <- metadata %>%
    mutate(aspRatcolor = case_when(
      (!!sym(vars[k]) < adjbox(!!sym(vars[k]), plot = FALSE)$fence[1])&(!!sym(bord.dist[k])<50) ~ "HIGH_MC",
      (!!sym(vars[k]) < median(!!sym(vars[k])) - 3 * mad(!!sym(vars[k])))&(!!sym(bord.dist[k])<50) ~ "HIGH_MAD",
      TRUE ~ "NO"
    ))
}
plotList[[k]] <- ggplot(metadata, aes(!!sym(bord.dist[k]), !!sym(vars[k]))) + 
  geom_point(aes(color = aspRatcolor, 
                 size = aspRatcolor), size=0.5) +
  geom_smooth() + 
  xlim(c(0,500)) +
  scale_color_manual(values = c("HIGH_MAD" = "tomato", 
                                "HIGH_MC" = "firebrick",
                                "NO" = "#c0c8cf")) + 
  scale_size_manual(values = c("HIGH_MAD" = 0.6, 
                               "HIGH_MC" = 0.6, 
                               "NO" = 0)) +
  geom_hline(yintercept = median(metadata[,vars[k]]) - 3*mad(metadata[,vars[k]]), 
             linewidth = 0.5,
             color = "tomato") +
  geom_hline(yintercept = adjbox(metadata[,vars[k]], plot = FALSE)$fence[1,1], 
             linewidth = 0.5,
             color = "firebrick") +
  theme(panel.background = element_rect(fill="white", color=NA),
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

smt <- ggplot_build(plotList[[k]])$data[[2]]
smm <- data.frame(x = smt$x,
                  y = smt$y)
min.y <- min(metadata[metadata[bord.dist[k]]<=500,vars[k]])
max.y <- max(metadata[metadata[bord.dist[k]]<=500,vars[k]])
plotList[[k]] <- plotList[[k]] + 
  scale_x_continuous(
    limits = c(0, 500),
    breaks = scales::breaks_pretty(n = 8),
    expand = expansion(mult = c(0, 0), add = c(0, 2.7))) +
  scale_y_continuous(
    limits = c(min.y, max.y),
    breaks = scales::breaks_pretty(n = 6),
    expand = c(0, 0.05)) +
  labs(x = x.labels[k],
       y = y.labels[k])

## get the derivative 
der <- c(NA, diff(smm$y)/diff(smm$x))
chg.idx <- -2
if (bord.dist[k] == "dist_border_hor") {chg.idx <- 2}
der0 <- smm$x[which(diff(sign(der)) == chg.idx)[1]]
der0_p1 <- smm$x[which(diff(sign(der)) == chg.idx)[1]+1]
extr.x <- smm[which(smm$y==max(smm$y, na.rm = T)),"x"]
if (bord.dist[k] == "dist_border_hor") {
  extr.x <- smm[which(smm$y==min(smm$y, na.rm = T)),"x"]
}
out.deriv[k,] <- c(extr.x, der0, der0_p1)
assign(paste0("p",k),plotList[[k]])

## add the vertical lines
colnames(out.deriv) <- c("extrm","der0","der0_p1")
plotList_add <- NULL
for (i in 1:length(plotList)){
  plotList_add[[i]] <- plotList[[i]] +
    # vertical line at 50 px
    geom_vline(xintercept = 50, 
               linewidth = 0.5,
               color = "grey40") +
    # vertical line at the smooth der=0
    geom_vline(xintercept = as.numeric(out.deriv[i,"der0"]),
               linewidth = 0.5,
               linetype = "dashed",
               color = "black") +
    # vertical line at average cell radius
    geom_vline(xintercept = sqrt(mean(metadata$Area)/pi), 
               linewidth = 0.5,
               linetype = "dotdash",
               color = "grey40")
  assign(paste0("p",i),plotList_add[[i]])
}

if(dataset == "CosMx_protein_tonsil" ){
  p1 <- p1 +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 6),
      expand = c(0, 0.05),
      labels = function(x) {
        sapply(x, function(xi) {
          if (!is.na(xi) && xi == 0) {
            "0"
          } else if (!is.na(xi)) {
            parse(text = gsub("e\\+?", " %*% 10^", scales::scientific_format()(xi)))
          } else {
            NA
          }
        })
      })
}

### create the final plot
final_plot <- ((p1 | p2) /
                 (p3 | p4)) +
  plot_layout(heights = c(1, 1, 1), widths = c(1, 1))
# --- Build rows ---
row1 <- p1 | p2
row2 <- p3 | p4
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
pdf(file.path("SuppFigure_5.pdf"), 
    width = 10, 
    height = 6,
    bg = "transparent")
final_plot
dev.off()

### If you used Protein dataset: 

### save the final figure ###
pdf(file.path("SuppFigure_6.pdf"), 
    width = 10, 
    height = 6,
    bg = "transparent")
final_plot
dev.off()


# Supplementary figure 7

datasets <- c("CosMx_rna_pancreas",
              "Xenium_lung_cancer",
              "MERFISH_mouse_liver",
              "CosMx_protein_tonsil")

var.to.plot <- "QC_score"
xlabel <- "Quality score"
for (k in 1:length(datasets)){
  metadata <- readRDS(file.path(metadata_dir,(paste0(datasets[k],"_metadata.rds"))))
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
    ylab("Number of cells") + 
    xlab(xlabel) + 
    theme_minimal() + 
    light_theme +
    theme(plot.margin = unit(c(0,20,0,0), "pt"),
          aspect.ratio = 1)
  assign(paste0("p",k),pp)
} 

### create the final plot
final_plot <- ((p1 | p2) /
                 (p3 | p4)) +
  plot_layout(heights = c(1, 1, 1), widths = c(1, 1))
# --- Build rows ---
row1 <- p1 | p2
row2 <- p3 | p4
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
pdf(file.path("SuppFigure7.pdf"), 
    width = 14, 
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] robustbase_0.99-6 dplyr_1.1.4       patchwork_1.3.2   ggplot2_4.0.1    
# 
# loaded via a namespace (and not attached):
#   [1] vctrs_0.7.1        ggbeeswarm_0.7.3   nlme_3.1-168       cli_3.6.5          rlang_1.1.7        DEoptimR_1.1-4     generics_0.1.4    
# [8] Cairo_1.7-0        S7_0.2.1           labeling_0.4.3     glue_1.8.0         scales_1.4.0       grid_4.5.1         tibble_3.3.1      
# [15] lifecycle_1.0.5    compiler_4.5.1     vipor_0.4.7        ggrastr_1.0.2      RColorBrewer_1.1-3 pkgconfig_2.0.3    mgcv_1.9-4        
# [22] rstudioapi_0.17.1  beeswarm_0.4.0     lattice_0.22-7     farver_2.1.2       viridisLite_0.4.2  R6_2.6.1           dichromat_2.0-0.1 
# [29] tidyselect_1.2.1   splines_4.5.1      pillar_1.11.1      magrittr_2.0.4     Matrix_1.7-4       tools_4.5.1        withr_3.0.2       
# [36] gtable_0.3.6      
