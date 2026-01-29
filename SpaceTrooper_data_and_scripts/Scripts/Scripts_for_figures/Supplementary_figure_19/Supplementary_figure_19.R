# Clean, minimal, and reproducible script to reproduce supplementary figure 19
# To reproduce the analysis be sure to have the required packages installed

library(ggplot2)
library(patchwork)

# R = 4.5.1

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

### read data table for full runtime ###
dataset <- "full_run_data"
full.run.data <- readRDS(file.path("data","scalability",(paste0(dataset,".rds"))))

### plot the full runtime vs. cells
y_range_time <- diff(range(full.run.data$time_sec, na.rm = TRUE))
point_jitter_time <- ifelse(is.na(y_range_time) || y_range_time == 0, 0.5, 0.01 * y_range_time)
p1 <- ggplot(full.run.data, aes(n, time_sec)) +
  geom_point(shape = 21,
             fill = "red",
             colour = "red",
             alpha = 0.4,
             size = 3,
             stroke = 0.6,
             position = position_jitter(width = 5000, 
                                        height = point_jitter_time)) +
  geom_smooth(color = "blue",
              method = "lm",
              se = FALSE,
              linetype = "dashed",
              linewidth = 0.6,
              alpha = 0.9) +
  scale_x_continuous(breaks = c(10000,seq(100000, 1500000, by = 100000)),
                     labels = scales::comma,
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
                     expand = expansion(mult = c(0.02, 0.05))) +
  xlab("Number of cells") +
  ylab("Runtime (sec)") +
  light_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

### read data table for full peak memory ###
dataset <- "full_peak_data"
full.peak.data <- readRDS(file.path("data","scalability",(paste0(dataset,".rds"))))

### plot the RAM vs. cells
y_range_time <- diff(range(full.peak.data$peak_RAM_GB, na.rm = TRUE))
point_jitter_y <- ifelse(is.na(y_range_time) || y_range_time == 0, 0.5, 0.01 * y_range_time)
p2 <-  ggplot(full.peak.data, aes(x = n, y = peak_RAM_GB)) +
  geom_point(shape = 21,
             fill = "red",
             color = "red",
             alpha = 0.4,
             size = 3,
             stroke = 0.6,
             position = position_jitter(width = 0.009, height = point_jitter_y)) +
  geom_smooth(color = "blue",
              method = "lm",
              se = FALSE,
              linetype = "dashed",
              linewidth = 0.6,
              alpha = 0.9) +
  #scale_x_log10(expand = expansion(mult = c(0.025, 0.025)),
  #              labels = scales::label_number(),
  #              breaks = c(10000, seq(100000, 1400000, by = 100000))) +
  scale_x_continuous(breaks = c(10000,seq(100000, 1500000, by = 100000)),
                     labels = scales::comma,
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
                     expand = expansion(mult = c(0.02, 0.05))) +
  xlab("Number of cells") +
  ylab("Peak memory usage (GB)") +
  light_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

### read data table for core runtime ###
dataset <- "core_run_data"
core.run.data <- readRDS(file.path(main.dir,"data","scalability",(paste0(dataset,".rds"))))

### plot the core runtime vs. cells
y_range_time <- diff(range(core.run.data$time_sec, na.rm = TRUE))
point_jitter_time <- ifelse(is.na(y_range_time) || y_range_time == 0, 0.5, 0.01 * y_range_time)
p3 <- ggplot(core.run.data, aes(n, time_sec)) +
  geom_point(shape = 21,
             fill = "red",
             colour = "red",
             alpha = 0.4,
             size = 3,
             stroke = 0.6,
             position = position_jitter(width = 5000, 
                                        height = point_jitter_time)) +
  geom_smooth(color = "blue",
              method = "lm",
              se = FALSE,
              linetype = "dashed",
              linewidth = 0.6,
              alpha = 0.9) +
  #scale_x_log10(expand = expansion(mult = c(0.025, 0.025)),
  #              labels = scales::label_number(),
  #                breaks = c(10000, seq(100000, 1400000, by = 100000))) +
  scale_x_continuous(breaks = c(10000,seq(100000, 1500000, by = 100000)),
                     labels = scales::comma,
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
                     expand = expansion(mult = c(0.02, 0.05))) +
  xlab("Number of cells") +
  ylab("Runtime (sec)") +
  light_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

### read data table for core peak memory ###
dataset <- "core_peak_data"
core.peak.data <- readRDS(file.path(main.dir,"data","scalability",(paste0(dataset,".rds"))))

### plot the RAM vs. cells
y_range_time <- diff(range(core.peak.data$peak_RAM_GB, na.rm = TRUE))
point_jitter_y <- ifelse(is.na(y_range_time) || y_range_time == 0, 0.5, 0.01 * y_range_time)
p4 <-  ggplot(core.peak.data, aes(x = n, y = peak_RAM_GB)) +
  geom_point(shape = 21,
             fill = "red",
             color = "red",
             alpha = 0.4,
             size = 3,
             stroke = 0.6,
             position = position_jitter(width = 0.009, height = point_jitter_y)) +
  geom_smooth(color = "blue",
              method = "lm",
              se = FALSE,
              linetype = "dashed",
              linewidth = 0.6,
              alpha = 0.9) +
  #scale_x_log10(expand = expansion(mult = c(0.025, 0.025)),
  #              labels = scales::label_number(),
  #              breaks = c(10000, seq(100000, 1400000, by = 100000))) +
  scale_x_continuous(breaks = c(10000,seq(100000, 1500000, by = 100000)),
                     labels = scales::comma,
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 6),
                     expand = expansion(mult = c(0.02, 0.05))) +
  xlab("Number of cells") +
  ylab("Peak memory usage (GB)") +
  light_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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
pdf(file.path("SuppFigure19.pdf"), 
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
#   [1] patchwork_1.3.2 ggplot2_4.0.1  
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.4     gtable_0.3.6      
# [7] glue_1.8.0         tibble_3.3.1       dichromat_2.0-0.1  pkgconfig_2.0.3    generics_0.1.4     dplyr_1.1.4       
# [13] lifecycle_1.0.5    cli_3.6.5          S7_0.2.1           scales_1.4.0       grid_4.5.1         vctrs_0.7.1       
# [19] withr_3.0.2        compiler_4.5.1     rstudioapi_0.17.1  tools_4.5.1        pillar_1.11.1      colorspace_2.1-2  
# [25] rlang_1.1.7 