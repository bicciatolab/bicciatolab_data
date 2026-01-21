# Clean, minimal, and reproducible script to reproduce figures 2 and related supplementaries
# To reproduce the analysis be sure to have the required packages installed
# We suggest using a conda environment created using the environment file provided in https://github.com/bicciatolab/bicciatolab_data/blob/main/SpaceTrooper_data_and_scripts/Scripts/Seurat_v5.yml

# R = 4.4.3
# Usage: set wdir, data_dir, and DBkero_metadata, then run

library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(pbapply)
library(future)
options(
  future.globals.maxSize = 20 * 1024^3,   # 20 Gb
  future.seed = TRUE
  )
future::plan("multisession", workers = 16)

# install ggvenn if not already installed
if (!requireNamespace("ggvenn", quietly = TRUE)) {
  install.packages("ggvenn")
}
library(ggvenn)

sessionInfo()
# R version 4.4.3 (2025-02-28)
# Platform: x86_64-conda-linux-gnu
# Running under: Ubuntu 24.04.1 LTS
# 
# Matrix products: default
# BLAS/LAPACK: R version 4.4.3 (2025-02-28)
# Platform: x86_64-conda-linux-gnu
# Running under: Ubuntu 24.04.1 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /data01/Programs/anaconda3/envs/test/lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0
# 
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
# 
# time zone: Etc/UTC
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
#  [1] future_1.68.0      pbapply_1.7-4      cowplot_1.2.0      viridis_0.6.5
#  [5] viridisLite_0.4.2  ggplot2_4.0.1      tidyr_1.3.1        data.table_1.17.8
#  [9] dplyr_1.1.4        Seurat_5.3.1       SeuratObject_5.2.0 sp_2.2-0
# 
# loaded via a namespace (and not attached):
#  [1] deldir_2.0-4           gridExtra_2.3          rlang_1.1.6
#  [4] magrittr_2.0.4         RcppAnnoy_0.0.22       otel_0.2.0
#  [7] matrixStats_1.5.0      ggridges_0.5.7         compiler_4.4.3
# [10] spatstat.geom_3.6-1    png_0.1-8              vctrs_0.6.5
# [13] reshape2_1.4.5         stringr_1.6.0          pkgconfig_2.0.3
# [16] fastmap_1.2.0          promises_1.5.0         purrr_1.2.0
# [19] jsonlite_2.0.0         goftest_1.2-3          later_1.4.4
# [22] spatstat.utils_3.2-0   irlba_2.3.5.1          parallel_4.4.3
# [25] cluster_2.1.8.1        R6_2.6.1               ica_1.0-3
# [28] stringi_1.8.7          RColorBrewer_1.1-3     spatstat.data_3.1-9
# [31] reticulate_1.44.1      parallelly_1.45.1      spatstat.univar_3.1-5
# [34] lmtest_0.9-40          scattermore_1.2        Rcpp_1.1.0
# [37] tensor_1.5.1           future.apply_1.20.0    zoo_1.8-14
# [40] sctransform_0.4.2      httpuv_1.6.16          Matrix_1.7-4
# [43] splines_4.4.3          igraph_2.1.4           tidyselect_1.2.1
# [46] abind_1.4-8            codetools_0.2-20       spatstat.random_3.4-3
# [49] miniUI_0.1.2           spatstat.explore_3.6-0 listenv_0.10.0
# [52] lattice_0.22-7         tibble_3.3.0           plyr_1.8.9
# [55] withr_3.0.2            shiny_1.11.1           S7_0.2.1
# [58] ROCR_1.0-11            Rtsne_0.17             fastDummies_1.7.5
# [61] survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.2-4
# [64] pillar_1.11.1          KernSmooth_2.23-26     plotly_4.11.0
# [67] generics_0.1.4         RcppHNSW_0.6.0         scales_1.4.0
# [70] globals_0.18.0         xtable_1.8-4           glue_1.8.0
# [73] lazyeval_0.2.2         tools_4.4.3            RSpectra_0.16-2
# [76] ggvenn_0.1.19          RANN_2.6.2             dotCall64_1.2
# [79] grid_4.4.3             nlme_3.1-168           patchwork_1.3.2
# [82] cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1
# [85] uwot_0.2.4             gtable_0.3.6           digest_0.6.39
# [88] progressr_0.18.0       ggrepel_0.9.6          htmlwidgets_1.6.4
# [91] farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.4
# [94] httr_1.4.7             mime_0.13              MASS_7.3-65


# --- User settings: set these before running ----
wdir <- ".../SpaceTrooper"                              # working directory (project root)
data_dir <- ".../CosMx_Breast/CosMX_data_Case2"         # Nanostring data directory (for LoadNanostring)
DBkero_metadata <- ".../CosMx_rna_DBKERO_metadata.rds"  # file path to SpatialExperiment metadata (data from https://github.com/bicciatolab/bicciatolab_data/blob/main/SpaceTrooper_data_and_scripts/Spe_metadata/CosMx_rna_DBKERO_metadata.rds)
# ------------------------------------------------

wdir <- "/NAS06/work/matteo/SpaceTrooper"
data_dir <- "/NAS06/work/Spatial_omics/DBKero/CosMx_Breast/CosMX_data_Case2"
DBkero_metadata <- "/NAS06/work/matteo/SpaceTrooper/CosMx_rna_DBKERO_metadata.rds"


if (wdir == "" || data_dir == "") stop("Please set 'wdir' and 'data_dir' at the top of the script.")
setwd(wdir)

outdir <- file.path(wdir, "Figure_S13")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

metadata <- readRDS(DBkero_metadata)  # load pre-existing metadata (from https://github.com/bicciatolab/bicciatolab_data/blob/main/SpaceTrooper_data_and_scripts/Spe_metadata/CosMx_rna_DBKERO_metadata.rds)
rownames(metadata) <- paste0(metadata$cellID, "_", metadata$fov)


# Cell type palette
celltype_palette <- c(
  "TAMs" =  "blue",
  "DCs" =   "dodgerblue4",
  "Mast cells" =    "lightblue1",
  "T cells" =   "darkseagreen1",
  "NK cells" = "darkolivegreen1",
  "Mural cells" =   "goldenrod4",
  "Myoepithelial cells" =   "goldenrod1",
  "Blood ECs" = "lightgoldenrod1",
  "CAFs" =  "lightgoldenrod4",
  "Plasma cells" =  "deeppink1",
  "Mix BC cells TAMs" = "purple",
  "BC cells" =  "orangered1"
)

cl_levels <- c(
  "Mixed", 
  "T cells",
  "NK cells",
  "Plasma cells",
  "Mast cells",
  "DCs",
  "TAMs",
  "Blood ECs",
  "Mural cells",
  "CAFs",
  "Myoepithelial cells",
  "Mix BC cells TAMs",
  "BC cells"
)


# --------------------------
# FIGURES (publication-ready)
# - 1) Bar plots: proportion of cell types filtered by different QC metrics
# - 2) Bar plots: proportion of cell types filtered by QS thresholds
# --------------------------

# Define themes
theme_barplot <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, hjust = 1, margin = margin(r = 2)),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    aspect.ratio = 0.5,
    panel.border = element_rect(colour = "#444444", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.margin = margin(0, 0, 5.5, 5.5)  # Adjust plot margins
  )


### 1) Bar plots: proportion of cell types filtered by different QC metrics ----
# define QC metric thresholds
area_thr <- median(metadata$Area_um) + 3 * mad(metadata$Area_um) # 216.1425

metadata$Area <- NULL

qc_thresholds <- list(
  Total_probe_counts = list("total" = 20, "direction" = "above"),
  Area = list("Area_um" = area_thr, "direction" = "below"),
  Control_probe_ratio = list("ctrl_total_ratio" = 0.01, "direction" = "below")
)

# Function to classify cells based on QC metrics
classify_qc <- function(metadata, qc_thresholds) {
  qc_results <- lapply(names(qc_thresholds), function(metric) {
    metric_def <- qc_thresholds[[metric]]
    # first element should be the named threshold (e.g. "nCount_Nanostring" = 20)
    colname <- names(metric_def)[1]
    if (is.null(colname) || colname == "" || colname == "direction") {
      stop("Invalid qc_thresholds format for metric: ", metric)
    }
    threshold_value <- metric_def[[1]]
    direction <- metric_def[["direction"]]
    if (is.null(direction)) {
      stop("Missing 'direction' for metric: ", metric)
    }
    if (direction == "above") {
      pass <- metadata[[colname]] >= threshold_value
    } else if (direction == "below") {
      pass <- metadata[[colname]] <= threshold_value
    } else {
      stop("Invalid direction in qc_thresholds")
    }
    return(pass)
  })
  qc_df <- as.data.table(do.call(cbind, qc_results))
  setnames(qc_df, names(qc_thresholds))
  qc_df[, Combined := rowSums(.SD) == length(qc_thresholds)]
  return(qc_df)
}

# Classify cells based on QC metrics
qc_classification <- classify_qc(metadata, qc_thresholds)
metadata <- cbind(metadata, qc_classification)


# Plot together QC classification proportion in a single plot, per celltype
# Combined plot: proportion of filtered cells per cell type for each QC metric
metric_palette <- c(
      "Total probe counts" = "#b7d3fcff",
      "Area" = "#f3a8daff",
      "Control probe ratio" = "#5468a0ff",
      "Combined" = "#bd0000ff"
    )

# For each QC classification plot the proportion of cell types passing/failing
qc_metrics <- names(qc_thresholds)
qc_metrics <- c(qc_metrics, "Combined")

metadata$InSituType_Simple <- factor(metadata$InSituType_Simple, levels = cl_levels)
metadata <- metadata[!is.na(metadata$InSituType_Simple), ]  # remove empty cells (with NA as InSituType_Simple)

metadata_grouped <- metadata %>%
  group_by(InSituType_Simple) %>%
  mutate(InSituType_Simple_plusCount = paste0(InSituType_Simple, " - ", n(), " cells"))

plot_data_list <- lapply(qc_metrics, function(metric) {
  df <- metadata_grouped %>%
    group_by(InSituType_Simple_plusCount) %>%
    summarize(
      Total = n(),
      Pass = sum(get(metric)),
      .groups = "drop"
    ) %>% 
    mutate(Filtered_ratio = (Total - Pass) / Total)
  df$QC_Metric <- metric
  return(df)
})
plot_data <- bind_rows(plot_data_list)
plot_data$QC_Metric <- factor(gsub("_", " ", plot_data$QC_Metric), levels = c("Combined", "Area", "Control probe ratio", "Total probe counts"))

max_y <- 0.2

p_combined <- ggplot(plot_data, aes(x = InSituType_Simple_plusCount, y = Filtered_ratio, fill = QC_Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(
    aes(label = Total - Pass),
    position = position_dodge(width = 0.9),
    vjust = 0.5,
    hjust = -1,
    size = 3
  ) +
  labs(
    title = "Proportion of Cells Filtered by QC Metrics (by Cell Type)",
    y = "Proportion of Cells Filtered",
    x = "Cell Type"
  ) +
  scale_fill_manual(
    values = metric_palette,
    name = "QC Metric",
    breaks = rev(levels(plot_data$QC_Metric))
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max_y), breaks = seq(0, max_y, by = 0.02), labels = scales::percent_format(accuracy = 2)) +
  theme_barplot +
  theme(
    panel.grid.major.x = element_line(color = "#dddddd", linewidth = 0.25),
    aspect.ratio = 0.75)

pdf(file = file.path(outdir, "QC_metrics_celltype_filtered_proportions_combined.pdf"), width = 13, height = 8)
print(p_combined)
dev.off()


# Create a venn plot using ggplot2 to show overlaps between rejected cells for different QC metrics, excluding "Combined"
venn_data <- lapply(names(qc_thresholds), function(metric) {
  which(!metadata[[metric]])  # Use `!metadata[[metric]]` to get indices of rejected cells
})
names(venn_data) <- gsub("_", " ", qc_metrics[qc_metrics != "Combined"])
pdf(file = file.path(outdir, "QC_metrics_venn_plot_rejected.pdf"), width = 8, height = 6)
ggvenn(
  venn_data,
  fill_color = metric_palette[gsub("_", " ", qc_metrics[qc_metrics != "Combined"])],
  stroke_size = 0.5,
  set_name_size = 4,
  auto_scale = FALSE
)
ggvenn(
  venn_data[c("Total probe counts", "Control probe ratio")],
  fill_color = metric_palette[c("Total probe counts", "Control probe ratio")],
  stroke_size = 0.5,
  set_name_size = 4,
  auto_scale = TRUE
)
ggvenn(
  venn_data[c("Control probe ratio", "Area")],
  fill_color = metric_palette[c("Control probe ratio", "Area")],
  stroke_size = 0.5,
  set_name_size = 4,
  auto_scale = TRUE
)
dev.off()


# Plot cell type proportion of rejected cells for each QC metric (x axis) using a stacked bar plot
plot_data_list <- lapply(qc_metrics, function(metric) {
  df <- metadata %>%
    group_by(InSituType_Simple) %>%
    summarize(
      Total = n(),
      Filtered = Total - sum(get(metric)),
      .groups = "drop"
    )
  df$QC_Metric <- metric
  return(df)
})
plot_data <- bind_rows(plot_data_list) %>%
  group_by(QC_Metric) %>%
  mutate(Proportion = Filtered / sum(Filtered)) %>%
  ungroup() %>%
  mutate(
    QC_Metric = factor(gsub("_", " ", QC_Metric), levels = c("Combined", "Area", "Control probe ratio", "Total probe counts")),
    InSituType_Simple = factor(InSituType_Simple, levels = rev(cl_levels))
  )

overall_plot_data <- metadata %>%
  group_by(InSituType_Simple) %>%
  summarize(Proportion = n() / nrow(metadata)) %>%
  mutate(InSituType_Simple = factor(InSituType_Simple, levels = rev(cl_levels)))


pdf(file = file.path(outdir, "QC_metrics_celltype_proportions_filtered_stackedbar.pdf"), width = 10, height = 6)

p <- ggplot(plot_data, aes(x = QC_Metric, y = Proportion, fill = InSituType_Simple)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.5) +
  labs(
    title = paste("Cell Type Proportions of Filtered Cells - QC Metric:", metric),
    y = "Proportion",
    x = ""
  ) +
  scale_fill_manual(values = celltype_palette) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = scales::percent_format(accuracy = 2)) +
  theme_barplot
print(p)

p_overall <- ggplot(overall_plot_data, aes(x = "Control probe ratio", y = Proportion, fill = InSituType_Simple)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.5) +
  labs(
    title = "Overall Cell Type Proportions in the Sample",
    y = "Proportion",
    x = ""
  ) +
  scale_fill_manual(values = celltype_palette) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = scales::percent_format(accuracy = 2)) +
  theme_barplot +
  theme(
    aspect.ratio = 0.125,
    axis.text.y = element_text(color = "white"),
  )
print(p_overall)

dev.off()





### 2) Bar plots: proportion of cell types filtered by QS thresholds ----

qs_eq <- quantile(metadata$QC_score, probs = 1 - mean(metadata$Combined)) # 0.3402
qs_mad <- median(metadata$QC_score) - 3 * mad(metadata$QC_score) # 0.6340

qc_thresholds <- list(
  Combined = NULL,
  QS_equivalent = list("QC_score" = qs_eq, "direction" = "above"),
  QS_suggested = list("QC_score" = qs_mad, "direction" = "above")
)

# Function to classify cells based on QC metrics
classify_qc <- function(metadata, qc_thresholds) {
  if ("Combined" %in% names(qc_thresholds)) {
    qc_results <- lapply(names(qc_thresholds), function(metric) {
      if (metric == "Combined") {
        # Use the custom filter previously defined
        return(metadata$Combined)
      } else {
        metric_def <- qc_thresholds[[metric]]
        # first element should be the named threshold (e.g. "nCount_Nanostring" = 20)
        colname <- names(metric_def)[1]
        if (is.null(colname) || colname == "" || colname == "direction") {
          stop("Invalid qc_thresholds format for metric: ", metric)
        }
        threshold_value <- metric_def[[1]]
        direction <- metric_def[["direction"]]
        if (is.null(direction)) {
          stop("Missing 'direction' for metric: ", metric)
        }
        if (direction == "above") {
          pass <- metadata[[colname]] >= threshold_value
        } else if (direction == "below") {
          pass <- metadata[[colname]] <= threshold_value
        } else {
          stop("Invalid direction in qc_thresholds")
        }
        return(pass)
      }
    })
  } else {
    qc_results <- lapply(names(qc_thresholds), function(metric) {
      metric_def <- qc_thresholds[[metric]]
      colname <- names(metric_def)[1]
      if (is.null(colname) || colname == "" || colname == "direction") {
        stop("Invalid qc_thresholds format for metric: ", metric)
      }
      threshold_value <- metric_def[[1]]
      direction <- metric_def[["direction"]]
      if (is.null(direction)) {
        stop("Missing 'direction' for metric: ", metric)
      }
      if (direction == "above") {
        pass <- metadata[[colname]] >= threshold_value
      } else if (direction == "below") {
        pass <- metadata[[colname]] <= threshold_value
      } else {
        stop("Invalid direction in qc_thresholds")
      }
      return(pass)
    })
  }
  qc_df <- as.data.table(do.call(cbind, qc_results))
  setnames(qc_df, names(qc_thresholds))
  if (!"Combined" %in% names(qc_thresholds)) {
    qc_df[, Combined := rowSums(.SD) == length(qc_thresholds)]
  }
  return(qc_df)
}

# Classify cells based on QC metrics
qc_classification <- classify_qc(metadata, qc_thresholds)
metadata <- cbind(metadata, qc_classification[, !"Combined", with = FALSE])

# For each QC classification plot the proportion of cell types passing/failing
qc_metrics <- names(qc_thresholds)

# Plot together QC classification proportion in a single plot, per celltype
# Combined plot: proportion of filtered cells per cell type for each QC metric
metric_palette <- c(
      "QS equivalent" = "#8a3bbeff",
      "QS suggested" = "#5d009bff",
      "Combined" = "#bd0000ff"
    )

metadata_grouped <- metadata %>%
  group_by(InSituType_Simple) %>%
  mutate(InSituType_Simple_plusCount = paste0(InSituType_Simple, " - ", n(), " cells"))

plot_data_list <- lapply(qc_metrics, function(metric) {
  df <- metadata_grouped %>%
    group_by(InSituType_Simple_plusCount) %>%
    summarize(
      Total = n(),
      Pass = sum(get(metric)),
      .groups = "drop"
    ) %>% 
    mutate(Filtered_ratio = (Total - Pass) / Total)
  df$QC_Metric <- metric
  return(df)
})
plot_data <- bind_rows(plot_data_list)
plot_data$QC_Metric <- factor(gsub("_", " ", plot_data$QC_Metric), levels = rev(c("Combined", "QS equivalent", "QS suggested")))

max_y <- 0.5

p_combined <- ggplot(plot_data, aes(x = InSituType_Simple_plusCount, y = Filtered_ratio, fill = QC_Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(
    aes(label = Total - Pass),
    position = position_dodge(width = 0.9),
    vjust = 0.5,
    hjust = -1,
    size = 3
  ) +
  labs(
    title = "Proportion of Cells Filtered by QS threshold (by Cell Type)",
    y = "Proportion of Cells Filtered",
    x = "Cell Type"
  ) +
  scale_fill_manual(
    values = metric_palette,
    name = "QS Threshold",
    breaks = rev(levels(plot_data$QC_Metric))
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max_y), breaks = seq(0, max_y, by = 0.1), labels = scales::percent_format(accuracy = 2)) +
  theme_barplot +
  theme(
    panel.grid.major.x = element_line(color = "#dddddd", linewidth = 0.25),
    aspect.ratio = 0.75)

pdf(file = file.path(outdir, "QS_thr_celltype_filtered_proportions_combined.pdf"), width = 13, height = 8)
print(p_combined)
dev.off()


# Create a venn plot using ggplot2 to show overlaps between rejected cells for different QC metrics, excluding "Combined"
venn_data <- lapply(names(qc_thresholds), function(metric) {
  which(!metadata[[metric]])  # Use `!metadata[[metric]]` to get indices of rejected cells
})
names(venn_data) <- gsub("_", " ", qc_metrics)
pdf(file = file.path(outdir, "QS_thr_venn_plot_rejected.pdf"), width = 8, height = 6)
ggvenn(
  venn_data[names(venn_data) != "QS suggested"],
  fill_color = metric_palette[gsub("_", " ", qc_metrics[qc_metrics != "QS suggested"])],
  stroke_size = 0.5,
  set_name_size = 4,
  auto_scale = TRUE
)
ggvenn(
  venn_data[names(venn_data) != "QS equivalent"],
  fill_color = metric_palette[gsub("_", " ", qc_metrics[qc_metrics != "QS equivalent"])],
  stroke_size = 0.5,
  set_name_size = 4,
  auto_scale = TRUE
)
dev.off()


# Plot cell type proportion of rejected cells for each QC metric (x axis) using a stacked bar plot
plot_data_list <- lapply(qc_metrics, function(metric) {
  df <- metadata %>%
    group_by(InSituType_Simple) %>%
    summarize(
      Total = n(),
      Filtered = Total - sum(get(metric)),
      .groups = "drop"
    )
  df$QC_Metric <- metric
  return(df)
})
plot_data <- bind_rows(plot_data_list) %>%
  group_by(QC_Metric) %>%
  mutate(Proportion = Filtered / sum(Filtered)) %>%
  ungroup() %>%
  mutate(
    QC_Metric = factor(gsub("_", " ", QC_Metric), levels = c("Combined", "QS equivalent", "QS suggested")),
    InSituType_Simple = factor(InSituType_Simple, levels = rev(cl_levels))
  )

overall_plot_data <- metadata %>%
  group_by(InSituType_Simple) %>%
  summarize(Proportion = n() / nrow(metadata)) %>%
  mutate(InSituType_Simple = factor(InSituType_Simple, levels = rev(cl_levels)))


pdf(file = file.path(outdir, "QS_thr_celltype_proportions_filtered_stackedbar.pdf"), width = 10, height = 6)

p <- ggplot(plot_data, aes(x = QC_Metric, y = Proportion, fill = InSituType_Simple)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.5) +
  labs(
    title = paste("Cell Type Proportions of Filtered Cells - QC Metric:", metric),
    y = "Proportion",
    x = ""
  ) +
  scale_fill_manual(values = celltype_palette) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = scales::percent_format(accuracy = 2)) +
  theme_barplot
print(p)

p_overall <- ggplot(overall_plot_data, aes(x = "Control probe ratio", y = Proportion, fill = InSituType_Simple)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.5) +
  labs(
    title = "Overall Cell Type Proportions in the Sample",
    y = "Proportion",
    x = ""
  ) +
  scale_fill_manual(values = celltype_palette) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), labels = scales::percent_format(accuracy = 2)) +
  theme_barplot +
  theme(
    aspect.ratio = 0.167,
    axis.text.y = element_text(color = "white"),
  )
print(p_overall)

dev.off()