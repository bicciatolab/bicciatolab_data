# Clean, minimal, and reproducible script to reproduce figures 2 and related supplementaries
# To reproduce the analysis be sure to have the required packages installed
# We suggest using a conda environment created using the environment file provided in https://github.com/bicciatolab/bicciatolab_data/blob/main/SpaceTrooper_data_and_scripts/Scripts/Seurat_v5.yml

# R = 4.4.3
# Usage: set wdir, data_dir, DBkero_metadata, and downstream_results, then run

library(Seurat)
library(dplyr)
library(data.table)
library(tidyr)
library(igraph)
library(ggplot2)
library(viridis)
library(cowplot)
library(pbapply)
library(future)
options(
  future.globals.maxSize = 100 * 1024^3,   # 100 Gb
  future.seed = TRUE
  )
future::plan("multisession", workers = 8)
set.seed(42)

sessionInfo()
# R version 4.4.3 (2025-02-28)
# Platform: x86_64-conda-linux-gnu
# Running under: Ubuntu 24.04.1 LTS
# 
# Matrix products: default
# BLAS/LAPACK: libopenblasp-r0.3.30.so;  LAPACK version 3.12.0
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
#  [5] viridisLite_0.4.2  ggplot2_4.0.1      igraph_2.1.4       tidyr_1.3.1
#  [9] data.table_1.17.8  dplyr_1.1.4        Seurat_5.3.1       SeuratObject_5.2.0
# [13] sp_2.2-0
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
# [43] splines_4.4.3          tidyselect_1.2.1       abind_1.4-8
# [46] codetools_0.2-20       spatstat.random_3.4-3  miniUI_0.1.2
# [49] spatstat.explore_3.6-0 listenv_0.10.0         lattice_0.22-7
# [52] tibble_3.3.0           plyr_1.8.9             withr_3.0.2
# [55] shiny_1.11.1           S7_0.2.1               ROCR_1.0-11
# [58] Rtsne_0.17             fastDummies_1.7.5      survival_3.8-3
# [61] polyclip_1.10-7        fitdistrplus_1.2-4     pillar_1.11.1
# [64] KernSmooth_2.23-26     plotly_4.11.0          generics_0.1.4
# [67] RcppHNSW_0.6.0         scales_1.4.0           globals_0.18.0
# [70] xtable_1.8-4           glue_1.8.0             lazyeval_0.2.2
# [73] tools_4.4.3            RSpectra_0.16-2        RANN_2.6.2
# [76] dotCall64_1.2          grid_4.4.3             nlme_3.1-168
# [79] patchwork_1.3.2        cli_3.6.5              spatstat.sparse_3.1-0
# [82] spam_2.11-1            uwot_0.2.4             gtable_0.3.6
# [85] digest_0.6.39          progressr_0.18.0       ggrepel_0.9.6
# [88] htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1
# [91] lifecycle_1.0.4        httr_1.4.7             mime_0.13
# [94] MASS_7.3-65



# --- User settings: set these before running ----
wdir <- ".../SpaceTrooper"                              # working directory (project root)
data_dir <- ".../CosMx_Breast/CosMX_data_Case2"         # Nanostring data directory (for LoadNanostring)
DBkero_metadata <- ".../CosMx_rna_DBKERO_metadata.rds"  # file path to SpatialExperiment metadata (data from https://github.com/bicciatolab/bicciatolab_data/blob/main/SpaceTrooper_data_and_scripts/Spe_metadata/CosMx_rna_DBKERO_metadata.rds)
downstream_results <- ".../CosMx_rna_DBKERO_downstream_effects_results.rds" # file path to downstream effects results (data from https://github.com/bicciatolab/bicciatolab_data/blob/main/SpaceTrooper_data_and_scripts/Spe_metadata/CosMx_rna_DBKERO_downstream_effects_results.rds)
# ------------------------------------------------

if (wdir == "" || data_dir == "") stop("Please set 'wdir' and 'data_dir' at the top of the script.")
setwd(wdir)

outdir <- file.path(wdir, "Figure_2")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


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

# Simplified palette for some plots (colored by immune, stromal, endothelial, and epithelial)
celltype_palette_simplified <- c(
  "TAMs" =  "#6868f6ff", # Immune
  "DCs" =   "#6868f6ff", # Immune
  "Mast cells" =    "#6868f6ff", # Immune
  "T cells" =   "#6868f6ff", # Immune
  "NK cells" = "#6868f6ff", # Immune
  "Mural cells" =   "#ffa467ff", # Stromal
  "Myoepithelial cells" =   "#ffa467ff", # Stromal
  "Blood ECs" = "#fed678ff", # Endothelial
  "CAFs" =  "#ffa467ff", # Stromal
  "Plasma cells" =  "#6868f6ff", # Immune
  "Mix BC cells TAMs" = "#c56ffaff", # Mixed
  "BC cells" =  "#fb5454ff" # Epithelial (tumor)
)

# Define homophily function
homophily <- function(graph, attr_name, ego = 1) {
  n <- vcount(graph)
  outs <- numeric(n)
  v_names <- V(graph)$name
  if (is.null(v_names)) v_names <- as.character(seq_len(n))
  for (i in seq_len(n)) {
    nn_vs <- ego(graph, order = ego, nodes = V(graph)[i], mode = "all")[[1]]
    nn_idx <- as.integer(nn_vs)
    attrs <- vertex_attr(graph, attr_name, index = nn_idx)
    v_attr <- vertex_attr(graph, attr_name, index = i)
    outs[i] <- sum(attrs == v_attr) / length(nn_idx)
  }
  names(outs) <- v_names
  return(outs)
}

# --- Load nanostring / Seurat object ----
message("Loading Nanostring data (Seurat)...")
seurat_obj <- LoadNanostring(data.dir = data_dir, fov = "DBkero")
if (is.null(seurat_obj)) stop("Failed to load Seurat object with LoadNanostring()")

# Load new metadata from existing object
metadata <- readRDS(DBkero_metadata)
rownames(metadata) <- paste0(metadata$cellID, "_", metadata$fov)
metadata <- metadata[colnames(seurat_obj), ]  # ensure same order

# Add new metadata with QC info
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata %>% select(QC_score, InSituType_Simple, InSituType_Simple.prob, Area_um, ctrl_total_ratio))

# Remove negative controls (probe names containing "Neg")
message("Removing negative control probes from assay features (if present)...")
feat_keep <- setdiff(rownames(seurat_obj), grep("Neg", rownames(seurat_obj), value = TRUE))
seurat_obj <- seurat_obj[feat_keep, ]

seurat_obj
# An object of class Seurat 
# 1000 features across 59234 samples within 1 assay 
# Active assay: Nanostring (1000 features, 0 variable features)
#  1 layer present: counts
#  1 spatial field of view present: DBkero


# Compute MAD-3 threshold for QC_score and define thresholds to test
mad_th <- median(seurat_obj$QC_score, na.rm = TRUE) - 3 * mad(seurat_obj$QC_score, na.rm = TRUE)
threshold_values <- list("0" = 0, "mad3" = mad_th, "custom" = "custom")
labels <- c("0" = "All cells", "mad3" = paste0("QS >= 3MAD (", round(mad_th, 3), ")"), "custom" = "Combined")

results_list <- readRDS(downstream_results)



# --------------------------
# FIGURES (publication-ready)
# - 1) Faceted UMAP: clusters and cell types
# - 2) Faceted UMAP: homophily, QC score, cell type probability
# - 3) Violin plots of homophily, QC score, cell type probability by cluster
# - 4) Violin plots of homophily, QC score, cell type probability by cell types
# - 5) Cell type composition by cluster
# - 6) Jitter plots of mean homophily by cluster
# --------------------------

# Build combined dataframe across thresholds
df_list <- list()
for (thr in names(results_list)) {
  res <- results_list[[thr]]
  cl <- res$clustering
  umap <- res$umap
  hs <- res$homophily_scores
  cells <- intersect(names(cl), rownames(umap))
  if (length(cells) == 0) next
  d <- data.frame(
    cell_id = cells,
    threshold = thr,
    cluster = as.character(cl[cells]),
    homophily = as.numeric(hs[cells]),
    UMAP_1 = umap[cells, 1],
    UMAP_2 = umap[cells, 2],
    stringsAsFactors = FALSE
  )
  # attach cell type from seurat_obj metadata (if present)
  d$cell_type <- seurat_obj$InSituType_Simple[match(d$cell_id, colnames(seurat_obj))]
  # attach QC and prob if present
  d$QC_score <- seurat_obj$QC_score[match(d$cell_id, colnames(seurat_obj))]
  d$celltype_prob <- seurat_obj$InSituType_Simple.prob[match(d$cell_id, colnames(seurat_obj))]
  df_list[[thr]] <- d
}
bigdf <- bind_rows(df_list)

new_labels  <- c("0" = "All cells", "custom" = "Combined", "mad3" = "QS >= 3MAD")
bigdf$threshold <- factor(
  mutate(bigdf, threshold = recode(threshold, !!!new_labels)) %>% pull(threshold),
  levels = new_labels)
bigdf <- data.table::as.data.table(bigdf)

# Define themes
theme_umap <- theme_minimal() +
    theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1
    )

scatter_theme <- theme(
    panel.background = element_rect(fill="white", color=NA),
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)
  )




### 1) Faceted UMAP: clusters and cell types

# plot clusters on UMAP for each threshold in a faceted plot with cluster labels at centroids
umap_plots <- list()

for (thr_name in unique(bigdf$threshold)) {
  df <- bigdf[threshold == thr_name, .(UMAP_1, UMAP_2, cluster, cell_id)]
  df$cluster <- factor(df$cluster, levels = sort(as.numeric(unique(df$cluster))))
  centroids <- df[, .(
    UMAP_1 = mean(UMAP_1, na.rm = TRUE),
    UMAP_2 = mean(UMAP_2, na.rm = TRUE)
  ), by = cluster]
  original_cluster_colors <- scales::hue_pal()(length(levels(df$cluster)))
  names(original_cluster_colors) <- levels(df$cluster)
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point(size = 0.1, alpha = 0.65) +
    geom_text(
      data = centroids,
      aes(x = UMAP_1, y = UMAP_2, label = as.character(cluster)),
      inherit.aes = FALSE,
      size = 3,
      color = "black",
      fontface = "bold",
      alpha = 0.8
    ) +
    scale_color_manual(values = original_cluster_colors) +
    labs(title = paste("UMAP Clusters -", as.character(thr_name)),
          x = "UMAP 1",
          y = "UMAP 2") +
    theme_umap
  umap_plots[[thr_name]] <- p
}

umap_plots[["legend"]] <- NULL
pdf(file.path(outdir, "UMAP_clusters_all_thresholds_with_centroids.pdf"), width = 12, height = 12, useDingbats = FALSE)
print(cowplot::plot_grid(plotlist = umap_plots, ncol = 2))
dev.off()


# plot InSituType_Simple on UMAP for each threshold in a faceted plot
label_to_plot <- "cell_type" # InSituType simplified labels
cell_palette <- celltype_palette
umap_plots <- list()

for (thr_name in unique(bigdf$threshold)) {
  df <- bigdf[threshold == thr_name, .(UMAP_1, UMAP_2, color_value = get(label_to_plot), cell_id)]
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_value)) +
    geom_point(size = 0.1, alpha = 0.65) +
    labs(title = paste("UMAP -", as.character(thr_name)),
         x = "UMAP 1",
         y = "UMAP 2",
         color = label_to_plot) +
    theme_umap 
  if(!is.null(cell_palette)){
    p <- p + scale_color_manual(values = cell_palette)
  }
  umap_plots[[thr_name]] <- p
}
umap_plots[["legend"]] <- cowplot::get_legend(umap_plots[[1]] + geom_point(size = 5) + theme(legend.position = "right"))

pdf(file.path(outdir, paste0("UMAP_", label_to_plot, "_all_thresholds.pdf")), width = 12, height = 12, useDingbats = FALSE)
print(cowplot::plot_grid(plotlist = umap_plots, ncol = 2))
dev.off()



### 2) Faceted UMAP: homophily, QC score, cell type probability

# plot homophily on UMAP for each threshold in a faceted plot
label_to_plot <- "homophily" # homophily scores
cell_palette <- c("#040404", "#214878ff", "#6394d3ff", "#a6b8dc", "#dedede")
umap_plots <- list()

for (thr_name in unique(bigdf$threshold)) {
  df <- bigdf[threshold == thr_name, .(UMAP_1, UMAP_2, color_value = get(label_to_plot), cell_id)]
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_value)) +
    geom_point(size = 0.1, alpha = 0.65) +
    labs(title = paste("UMAP -", as.character(thr_name)),
         x = "UMAP 1",
         y = "UMAP 2",
         color = label_to_plot) +
    theme_umap 
  if(!is.null(cell_palette)){
    p <- p + scale_color_gradientn(colors = cell_palette, limits = c(0,1), breaks = seq(0, 1, by = 0.2))
  }
  umap_plots[[thr_name]] <- p
}
umap_plots[["legend"]] <- cowplot::get_legend(umap_plots[[1]] + geom_point(size = 5) + theme(legend.position = "right"))

pdf(file.path(outdir, paste0("UMAP_", label_to_plot, "_all_thresholds.pdf")), width = 12, height = 12, useDingbats = FALSE)
print(cowplot::plot_grid(plotlist = umap_plots, ncol = 2))
dev.off()


# plot QC score on UMAP for each threshold in a faceted plot
label_to_plot <- "QC_score" # QC scores
cell_palette <- viridis::viridis(256, option = "plasma")
umap_plots <- list()

for (thr_name in unique(bigdf$threshold)) {
  df <- bigdf[threshold == thr_name, .(UMAP_1, UMAP_2, color_value = get(label_to_plot), cell_id)]
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_value)) +
    geom_point(size = 0.1, alpha = 0.65) +
    labs(title = paste("UMAP -", as.character(thr_name)),
         x = "UMAP 1",
         y = "UMAP 2",
         color = label_to_plot) +
    theme_umap 
  if(!is.null(cell_palette)){
    p <- p + scale_color_gradientn(colors = cell_palette, limits = c(0,1), breaks = seq(0, 1, by = 0.2))
  }
  umap_plots[[thr_name]] <- p
}
umap_plots[["legend"]] <- cowplot::get_legend(umap_plots[[1]] + geom_point(size = 5) + theme(legend.position = "right"))

pdf(file.path(outdir, paste0("UMAP_", label_to_plot, "_all_thresholds.pdf")), width = 12, height = 12, useDingbats = FALSE)
print(cowplot::plot_grid(plotlist = umap_plots, ncol = 2))
dev.off()


# plot InSituType probability on UMAP for each threshold in a faceted plot
label_to_plot <- "celltype_prob" # InSituType probabilities
cell_palette <- viridis::viridis(256, option = "viridis")
umap_plots <- list()

for (thr_name in unique(bigdf$threshold)) {
  df <- bigdf[threshold == thr_name, .(UMAP_1, UMAP_2, color_value = get(label_to_plot), cell_id)]
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = color_value)) +
    geom_point(size = 0.1, alpha = 0.65) +
    labs(title = paste("UMAP -", as.character(thr_name)),
         x = "UMAP 1",
         y = "UMAP 2",
         color = label_to_plot) +
    theme_umap 
  if(!is.null(cell_palette)){
    p <- p + scale_color_gradientn(colors = cell_palette, limits = c(0,1), breaks = seq(0, 1, by = 0.2))
  }
  umap_plots[[thr_name]] <- p
}
umap_plots[["legend"]] <- cowplot::get_legend(umap_plots[[1]] + geom_point(size = 5) + theme(legend.position = "right"))

pdf(file.path(outdir, paste0("UMAP_", label_to_plot, "_all_thresholds.pdf")), width = 12, height = 12, useDingbats = FALSE)
print(cowplot::plot_grid(plotlist = umap_plots, ncol = 2))
dev.off()



### 3) Violin plots of homophily, QC score, cell type probability by cluster

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

# For every threshold, find the most represented cell type if it is above 30%, otherwise it is "Mixed"
cl_order <- list()
for (thr in unique(bigdf$threshold)) {
  df <- bigdf[threshold == thr, .(cluster, cell_type)]
  composition_df <- as.data.frame.matrix(table(df$cluster, df$cell_type))
  composition_df_prop <- composition_df / rowSums(composition_df)
  dominant_ct <- apply(composition_df_prop, 1, function(x) {
    if (max(x) >= 0.3) {
      return(names(which.max(x)))
    } else {
      return("Mixed")
    }
  })
  abund <- sapply(names(dominant_ct), function(cl) {
    ct <- dominant_ct[cl]
    if (ct == "Mixed") {
      # Use the maximum proportion in that cluster when labeled Mixed
      return(max(as.numeric(composition_df_prop[cl, ]), na.rm = TRUE))
    } else if (ct %in% colnames(composition_df_prop)) {
      return(as.numeric(composition_df_prop[cl, ct]))
    } else {
      return(0)
    }
  })

  ct_rank <- sapply(dominant_ct, function(ct) {
    idx <- match(as.character(ct), as.character(cl_levels))
    if (is.na(idx)) idx <- length(cl_levels) + 1
    idx
  })

  ord <- order(ct_rank, abund, decreasing = FALSE)

  cluster_ids <- names(dominant_ct)
  cl_order[[as.character(thr)]] <- cluster_ids[ord]
}

# Create a levels list with clusters order for every threshold
levels_list <- lapply(cl_order, function(order) factor(order, levels = order))

# plot violin plots of homophily by cluster for each threshold in a faceted plot
create_violin_plot <- function(metric, y_label, file_name) {
  violin_plots <- list()
  for (thr_name in levels(bigdf$threshold)) {
    metric_values <- bigdf[threshold == thr_name, get(metric)]
    clusters <- bigdf[threshold == thr_name, cluster]
    df <- data.frame(
      value = metric_values,
      cluster = factor(clusters, levels = levels_list[[thr_name]])
    )
    original_cluster_colors <- scales::hue_pal()(length(levels(df$cluster)))
    original_cluster_colors <- original_cluster_colors[as.numeric(levels(df$cluster)) + 1]
    names(original_cluster_colors) <- levels(df$cluster)
    p <- ggplot(df, aes(x = cluster, y = value, fill = cluster)) +
      geom_violin(trim = TRUE, color = "black", scale = "width", alpha = 0.9, linewidth = 0.2) +
      geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.7, color = "black", linewidth = 0.2) +
      labs(title = paste("Violin Plot -", as.character(thr_name)),
           x = "Cluster",
           y = y_label) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      scale_fill_manual(values = original_cluster_colors)
    violin_plots[[thr_name]] <- p
  }
  pdf(file.path(outdir, file_name), width = 21, height = 4, useDingbats = FALSE)
  print(cowplot::plot_grid(plotlist = violin_plots, ncol = 3))
  dev.off()
}

# Create violin plots for homophily, QC score, and cell type probability
create_violin_plot("homophily", "Homophily Score", "Violin_Homophily_by_cluster_all_thresholds.pdf")
create_violin_plot("QC_score", "QC Score", "Violin_QC_score_by_cluster_all_thresholds.pdf")
create_violin_plot("celltype_prob", "Cell Type Probability", "Violin_CellType_Prob_by_cluster_all_thresholds.pdf")



### 4) Violin plots of homophily, QC score, cell type probability by cell types
create_violin_plot_by_celltype <- function(metric, y_label, file_name) {
  violin_plots <- list()
  for (thr_name in levels(bigdf$threshold)) {
    metric_values <- bigdf[threshold == thr_name, get(metric)]
    cell_types <- bigdf[threshold == thr_name, cell_type]
    df <- data.frame(
      value = metric_values,
      cell_type = factor(cell_types, levels = cl_levels)
    )
    p <- ggplot(df, aes(x = cell_type, y = value, fill = cell_type)) +
      geom_violin(trim = TRUE, color = "black", scale = "width", alpha = 0.9, linewidth = 0.2) +
      geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.7, color = "black", linewidth = 0.2) +
      labs(title = paste("Violin Plot -", as.character(thr_name)),
           x = "Cell Type",
           y = y_label) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      scale_fill_manual(values = celltype_palette)
    violin_plots[[thr_name]] <- p
  }
  pdf(file.path(outdir, file_name), width = 21, height = 4, useDingbats = FALSE)
  print(cowplot::plot_grid(plotlist = violin_plots, ncol = 3))
  dev.off()
}

# Create violin plots for homophily, QC score, and cell type probability
create_violin_plot_by_celltype("homophily", "Homophily Score", "Violin_Homophily_by_celltype_all_thresholds.pdf")
create_violin_plot_by_celltype("QC_score", "QC Score", "Violin_QC_score_by_celltype_all_thresholds.pdf")
create_violin_plot_by_celltype("celltype_prob", "Cell Type Probability", "Violin_CellType_Prob_by_celltype_all_thresholds.pdf")



### 5) Cell type composition by cluster
create_composition_plots <- function(palette = celltype_palette, label = "Cell Type Composition by Cluster") {
  for (thr_name in levels(bigdf$threshold)) {
    df <- bigdf[threshold == thr_name, .(cluster, cell_type)]
    composition_df <- as.data.frame.matrix(table(df$cluster, df$cell_type))
    composition_df <- composition_df / rowSums(composition_df)
    composition_df_melt <- reshape2::melt(as.matrix(composition_df))
    colnames(composition_df_melt) <- c("Cluster", "CellType", "Proportion")
    p <- ggplot(composition_df_melt, aes(x = factor(Cluster, levels = levels_list[[thr_name]]), y = Proportion, fill = factor(CellType, levels = cl_levels))) +
      geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.1) +
      labs(title = paste(label, "-", as.character(thr_name)),
           x = "Cluster",
           y = "Proportion",
           fill = "Cell Type") +
      scale_fill_manual(values = palette) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
  }
}

pdf(file.path(outdir, paste0("CellType_composition_by_cluster_QS_thresholds.pdf")), width = 8, height = 6, useDingbats = FALSE)
create_composition_plots(celltype_palette, "Cell Type Composition by Cluster")
dev.off()

pdf(file.path(outdir, paste0("CellType_composition_by_cluster_QS_thresholds_simplifiedColors.pdf")), width = 8, height = 6, useDingbats = FALSE)
create_composition_plots(celltype_palette_simplified, "Lineage Composition by Cluster")
dev.off()



### 6) Jitter plots of mean homophily by cluster
df_list <- list()
for (thr_name in levels(bigdf$threshold)) {
  df <- bigdf[threshold == thr_name, .(homophily, cluster)]
  df <- df[, .(mean_homophily = mean(homophily, na.rm = TRUE)), by = cluster]
  df$threshold <- thr_name
  df_list[[thr_name]] <- df
}
df <- data.table::rbindlist(df_list)

create_jitter_plot <- function(df) {
  p <- ggplot(df, aes(x = threshold, y = mean_homophily, fill = mean_homophily)) +
    geom_jitter(color = "#272964ff", shape = 21, alpha = 0.8, size = 3, stroke = 0.6, position = position_jitter(width = 0.1, height = 0)) +
    labs(x = "Threshold used for QC filtering",
         y = "Mean Homophily Score") +
    scatter_theme +
    scale_fill_gradientn(limits = c(0, 1), colors = c("#040404", "#214878ff", "#6394d3ff", "#a6b8dc", "#f5f5f5ff"), breaks = seq(0, 1, by = 0.2)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
  print(p)
}

pdf(file.path(outdir, "Mean_homophily_by_cluster_all_thresholds.pdf"), width = 5.5, height = 4.5, useDingbats = FALSE)
create_jitter_plot(df)
dev.off()

# make only for "All cells" and "QS >= 3MAD"
df_list <- list()
for (thr_name in c("All cells", "QS >= 3MAD")) {
  df <- bigdf[threshold == thr_name, .(homophily, cluster)]
  df <- df[, .(mean_homophily = mean(homophily, na.rm = TRUE)), by = cluster]
  df$threshold <- thr_name
  df_list[[thr_name]] <- df
}
df <- data.table::rbindlist(df_list)

pdf(file.path(outdir, "Mean_homophily_by_cluster_noCustom.pdf"), width = 4.5, height = 4.5, useDingbats = FALSE)
create_jitter_plot(df)
dev.off()





### Significance testing of QS and homophily differences between thresholds, pairwise between clusters
stats_results <- list()
for (thr_name in levels(bigdf$threshold)) {
  message(crayon::green(paste0("\nComputing statistics for threshold: ", thr_name)))
  df <- bigdf[threshold == thr_name, .(homophily, QC_score, cluster, celltype_prob)]
  homophily_values <- df$homophily
  qc_values <- df$QC_score
  celltype_prob_values <- df$celltype_prob
  clusters <- df$cluster
  cluster_ids <- unique(clusters)
  # Prepare dataframes to store results
  homophily_stats <- data.frame(
    Cluster1 = character(),
    Cluster2 = character(),
    p_value = numeric(),
    p_value_adj = numeric(),
    stringsAsFactors = FALSE
  )
  qc_stats <- data.frame(
    Cluster1 = character(),
    Cluster2 = character(),
    p_value = numeric(),
    p_value_adj = numeric(),
    stringsAsFactors = FALSE
  )
  celltype_prob_stats <- data.frame(
    Cluster1 = character(),
    Cluster2 = character(),
    p_value = numeric(),
    p_value_adj = numeric(),
    stringsAsFactors = FALSE
  )
  # Pairwise comparisons between clusters
  for (i in 0:(length(cluster_ids) - 1)) {
    for (j in 0:(length(cluster_ids) - 1)) {
      if (i == j) next
      cl1 <- cluster_ids[i]
      cl2 <- cluster_ids[j]
      homophily_cl1 <- homophily_values[clusters == cl1]
      homophily_cl2 <- homophily_values[clusters == cl2]
      qc_cl1 <- qc_values[clusters == cl1]
      qc_cl2 <- qc_values[clusters == cl2]
      celltype_prob_cl1 <- celltype_prob_values[clusters == cl1]
      celltype_prob_cl2 <- celltype_prob_values[clusters == cl2]
      if (length(homophily_cl1) > 0 && length(homophily_cl2) > 0) {
        hom_p <- wilcox.test(homophily_cl1, homophily_cl2, alternative = "less")$p.value
        homophily_stats <- rbind(homophily_stats, data.frame(Cluster1 = cl1, Cluster2 = cl2, p_value = hom_p, p_value_adj = NA, stringsAsFactors = FALSE))
      }
      if (length(qc_cl1) > 0 && length(qc_cl2) > 0) {
        qc_p <- wilcox.test(qc_cl1, qc_cl2, alternative = "less")$p.value
        qc_stats <- rbind(qc_stats, data.frame(Cluster1 = cl1, Cluster2 = cl2, p_value = qc_p, p_value_adj = NA, stringsAsFactors = FALSE))
      }
      if (length(celltype_prob_cl1) > 0 && length(celltype_prob_cl2) > 0) {
        prob_p <- wilcox.test(celltype_prob_cl1, celltype_prob_cl2, alternative = "less")$p.value
        celltype_prob_stats <- rbind(celltype_prob_stats, data.frame(Cluster1 = cl1, Cluster2 = cl2, p_value = prob_p, p_value_adj = NA, stringsAsFactors = FALSE))
      }
    }
  }
  # adjust p values with BH method
  homophily_stats$p_value_adj <- p.adjust(homophily_stats$p_value, method = "BH")
  qc_stats$p_value_adj <- p.adjust(qc_stats$p_value, method = "BH")
  celltype_prob_stats$p_value_adj <- p.adjust(celltype_prob_stats$p_value, method = "BH")
  stats_results[[thr_name]] <- list(
    homophily_stats = homophily_stats,
    qc_stats = qc_stats,
    celltype_prob_stats = celltype_prob_stats
  )
}

stats_results[["All cells"]]$qc_stats %>% filter(Cluster1 == "12") %>% arrange(p_value_adj)
#    Cluster1 Cluster2       p_value   p_value_adj
# 1        12        7  0.000000e+00  0.000000e+00
# 2        12        8  0.000000e+00  0.000000e+00
# 3        12        4  0.000000e+00  0.000000e+00
# 4        12        0  0.000000e+00  0.000000e+00
# 5        12        9  0.000000e+00  0.000000e+00
# 6        12        1  0.000000e+00  0.000000e+00
# 7        12        5  0.000000e+00  0.000000e+00
# 8        12        3  0.000000e+00  0.000000e+00
# 9        12       14 2.553928e-281 1.595223e-279
# 10       12       17 1.802845e-270 9.759402e-269
# 11       12       11 2.695343e-266 1.367887e-264
# 12       12        6 2.881054e-264 1.376127e-262
# 13       12       19 2.412458e-254 1.031008e-252
# 14       12        2 3.479345e-251 1.345347e-249
# 15       12       15 1.291991e-227 4.371237e-226
# 16       12       18 2.170966e-201 6.528979e-200
# 17       12       21 2.468259e-186 6.911126e-185
# 18       12       22 4.381989e-156 1.046522e-154
# 19       12       16 5.769738e-153 1.338579e-151
# 20       12       20 1.892027e-129 3.840815e-128
# 21       12       10 1.903199e-125 3.679518e-124
# 22       12       13 6.792919e-125 1.282756e-123
# 23       12       23  6.162548e-92  8.778928e-91
# 24       12       25  1.419524e-86  1.746445e-85
# 25       12       28  3.009569e-71  2.944302e-70
# 26       12       26  9.203345e-35  5.225955e-34
# 27       12       27  1.404333e-32  7.757265e-32
# 28       12       24  2.243188e-24  1.023297e-23

stats_results[["All cells"]]$homophily_stats %>% filter(Cluster1 == "12") %>% arrange(p_value_adj)
#    Cluster1 Cluster2       p_value   p_value_adj
# 1        12        6  0.000000e+00  0.000000e+00
# 2        12        2  0.000000e+00  0.000000e+00
# 3        12        7  0.000000e+00  0.000000e+00
# 4        12        8  0.000000e+00  0.000000e+00
# 5        12        4  0.000000e+00  0.000000e+00
# 6        12        0  0.000000e+00  0.000000e+00
# 7        12       14  0.000000e+00  0.000000e+00
# 8        12        9  0.000000e+00  0.000000e+00
# 9        12       17  0.000000e+00  0.000000e+00
# 10       12        1  0.000000e+00  0.000000e+00
# 11       12       18  0.000000e+00  0.000000e+00
# 12       12        5  0.000000e+00  0.000000e+00
# 13       12        3  0.000000e+00  0.000000e+00
# 14       12       20  0.000000e+00  0.000000e+00
# 15       12       15  0.000000e+00  0.000000e+00
# 16       12       19 7.581342e-274 7.892371e-273
# 17       12       22 4.433915e-255 4.235693e-254
# 18       12       11 7.279491e-212 5.738783e-211
# 19       12       13 4.578548e-195 3.379801e-194
# 20       12       16 2.809502e-181 1.983753e-180
# 21       12       23 2.971831e-154 1.856251e-153
# 22       12       21 1.215843e-148 7.423039e-148
# 23       12       10 1.100027e-146 6.616457e-146
# 24       12       25 1.424692e-139 8.033678e-139
# 25       12       26 3.181908e-129 1.722473e-128
# 26       12       27  2.305448e-85  1.040013e-84
# 27       12       24  3.678600e-77  1.588842e-76
# 28       12       28  2.335224e-49  8.698173e-49

stats_results[["All cells"]]$celltype_prob_stats %>% filter(Cluster1 == "12") %>% arrange(p_value_adj)
#    Cluster1 Cluster2       p_value   p_value_adj
# 1        12        6  0.000000e+00  0.000000e+00
# 2        12       13  0.000000e+00  0.000000e+00
# 3        12        2  0.000000e+00  0.000000e+00
# 4        12       11  0.000000e+00  0.000000e+00
# 5        12        7  0.000000e+00  0.000000e+00
# 6        12        8  0.000000e+00  0.000000e+00
# 7        12        4  0.000000e+00  0.000000e+00
# 8        12        0  0.000000e+00  0.000000e+00
# 9        12       14  0.000000e+00  0.000000e+00
# 10       12        9  0.000000e+00  0.000000e+00
# 11       12       17  0.000000e+00  0.000000e+00
# 12       12        1  0.000000e+00  0.000000e+00
# 13       12       18  0.000000e+00  0.000000e+00
# 14       12        5  0.000000e+00  0.000000e+00
# 15       12        3  0.000000e+00  0.000000e+00
# 16       12       20  0.000000e+00  0.000000e+00
# 17       12       19  0.000000e+00  0.000000e+00
# 18       12       15  0.000000e+00  0.000000e+00
# 19       12       10 8.512289e-307 1.382396e-305
# 20       12       22 3.869527e-261 5.610814e-260
# 21       12       16 2.077423e-239 2.811446e-238
# 22       12       21 1.167959e-215 1.374468e-214
# 23       12       23 1.102330e-181 1.133028e-180
# 24       12       25 3.515025e-148 3.102392e-147
# 25       12       26 3.624459e-131 2.857341e-130
# 26       12       24 2.535359e-118 1.871556e-117
# 27       12       27  7.198751e-97  4.602666e-96
# 28       12       28  1.019918e-85  5.915524e-85