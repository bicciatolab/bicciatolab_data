# Clean, minimal, and reproducible script to reproduce supplementary figure 24 
# To reproduce the analysis be sure to have the required packages installed

library(SpaceTrooper)
library(ggplot2)
library(dplyr)
library(pROC)

# Supplementary figure 24A

# CosMx DBKERO - Evaluate QS model performance across different data splits

# 5-fold single random split (50/50)
# procedure:
# 1) QS computation on full DBKERO dataset
# 2) threshold definition for low-quality cell classification based on full dataset QS using median - 3*mad
# 3) split DBKERO dataset into training and test set (50/50) using a random seed
# 4) model transfer from training to test set -> pseudoQS
# 5) repeat 5 times using each time a different random seed to generate the split
# 6) ROC curve computation using the pseudoQS as the predictor and the low-quality cell classification based on full dataset QS

# Modes:
# - `k_folds = 1`: single random split using `train_fraction`
# - `k_folds > 1`: k-fold cross-validation
#
# Important notes:
# - `train_fraction` is ignored when `k_folds > 1`
# - `threshold` is used only for classification summaries, not for model fitting

# We were interested in showing the roc curve for the 5-fold cross-validation

# Create train/test indices for a single random split.
#
# `train_fraction` is only used when `k_folds = 1`. It controls the proportion
# of cells assigned to the training set, allowing e.g. 50/50, 70/30, 80/20.

make_random_split <- function(n_cells, train_fraction, seed) {
    stopifnot(length(n_cells) == 1L, n_cells > 1L)
    stopifnot(length(train_fraction) == 1L, !is.na(train_fraction))

    if (train_fraction <= 0 || train_fraction >= 1) {
        stop("train_fraction must be strictly between 0 and 1")
    }

    set.seed(seed)

    # Sample the training cells at random.
    n_train <- floor(n_cells * train_fraction)
    if (n_train < 1L || n_train >= n_cells) {
        stop("train_fraction produces an empty training or test set")
    }

    train_idx <- sample(seq_len(n_cells), size=n_train, replace=FALSE)
    test_idx <- setdiff(seq_len(n_cells), train_idx)

    list(train_idx=train_idx, test_idx=test_idx)
}

safe_classification_labels <- function(spe_subset, verbose=FALSE) {
    tryCatch(
        list(
            data={
                df <- as.data.frame(colData(spe_subset))
                df
            },
            error=NULL
        ),
        error=function(e) list(data=NULL, error=conditionMessage(e))
    )
}

# Train the QS model on one training split/fold and evaluate it on both train
# and test cells.
#
# This is the core routine used by both single-split validation and k-fold CV.
run_one_fold <- function(spe, test_idx, verbose=TRUE) {
    # Partition the input SpatialExperiment into training and test subsets.
    spe_train <- spe[, -test_idx]
    spe_test <- spe[, test_idx]

    # Compute QS on training cells only
    spe_train <- computeQScore(spe_train, verbose=verbose)
    train_model <- metadata(spe_train)$QScore_model

    # apply the training model with its coefficients to the test set
    spe_test <- .applyQScoreModel(
        spe=spe_test,
        qsModel=train_model,
        scoreName="QScore_train_model"
    )

    train_eval <- data.frame(colData(spe_train)) |> select(cell_id, QScore)

    # safe function to return error
    test_ref <- safe_classification_labels(spe_test, verbose=verbose)

    test_eval <- if (is.null(test_ref$data)) {
        NULL
    } else {
        test_eval <- test_ref$data |>
            select(cell_id, QScore_train_model, is_low_full)   
    } 

    test_roc <- if (!is.null(test_eval) && nrow(test_eval) > 0L) {
        roc(test_eval$is_low_full, test_eval$QScore_train_model)
    } else {
        NULL
    }

    # Return both summaries and raw objects/tables for further inspection.
    list(
        model_formula=train_model,
        train_class_balance=table(train_eval$is_low_full),
        train_label_error=NA_character_,
        test_label_error=test_ref$error,
        spe_train=spe_train,
        spe_test=spe_test,
        train_eval=train_eval,
        test_eval=test_eval,
        test_roc=test_roc,
        test_auc=if (!is.null(test_roc)) pROC::auc(test_roc) else NA_real_
    )
}

# k-fold cross-validation

# Assign each cell to one of the k folds used in cross-validation.
#
# Each fold acts as test set once, while the remaining folds are used for
# training. This function is only used when `k_folds > 1`.

make_fold_ids <- function(n_cells, k_folds, seed) {
    stopifnot(length(n_cells) == 1L, n_cells > 1L)
    stopifnot(length(k_folds) == 1L, !is.na(k_folds), k_folds >= 1)

    k_folds <- as.integer(k_folds)
    if (k_folds < 1L) {
        stop("k_folds must be >= 1")
    }

    set.seed(seed)

    if (k_folds > n_cells) {
        stop("k_folds cannot be greater than the number of cells")
    }

    # Shuffle cells and then distribute them approximately evenly across folds.
    shuffled <- sample(seq_len(n_cells))
    fold_id <- integer(n_cells)
    fold_id[shuffled] <- rep(seq_len(k_folds), length.out=n_cells)
    fold_id
}

# Safely extract an optional scalar character field from a list.
#
# This is used when collecting fold-level error messages at the end of k-fold
# evaluation. Missing, NULL, or NA entries are converted to `NA_character_`.
extract_optional_chr1 <- function(x, name) {
    value <- x[[name]]

    if (is.null(value) || length(value) == 0L || all(is.na(value))) {
        return(NA_character_)
    }

    as.character(value[[1]])
}

evaluate_qs_split <- function(
    data_dir="D:\\CosMx_DBKERO_input_data",
    sample_name="CosMx_DBKERO",
    seed=1998,
    k_folds=5,
    train_fraction=0.5,
    verbose=TRUE
) {
    # Validate the high-level input parameters.
    stopifnot(length(k_folds) == 1L, !is.na(k_folds), k_folds >= 1)
    k_folds <- as.integer(k_folds)
    stopifnot(length(train_fraction) == 1L, !is.na(train_fraction))

    if (train_fraction <= 0 || train_fraction >= 1) {
        stop("train_fraction must be strictly between 0 and 1")
    }

    # Read the CosMx CosMx dataset
    spe <- readCosmxSPE(data_dir, sampleName=sample_name)

    # Compute the per-cell QC metrics required by the QS model.
    spe <- spatialPerCellQC(spe)

    # Compute full dataset QS

    set.seed(1998) # this is fixed
    spe <- computeQScore(spe, modelFormula=NULL, verbose=verbose)

    threshold <- median(spe$QScore, na.rm = TRUE) - 3*mad(spe$QScore, na.rm = TRUE)
    message("Using threshold on QS computed on full dataset= ", round(threshold, 4), "\n for low-quality cell classification")
    spe$is_low_full <- spe$QScore < threshold

    # Count the number of cells remaining after QC preprocessing.
    n_cells <- ncol(spe)

    if (k_folds == 1L) {
        # Single random split: build train/test indices according to the
        # requested training fraction.
        split_idx <- make_random_split(
            n_cells=n_cells,
            train_fraction=train_fraction,
            seed=seed
        )

        res <- run_one_fold(
            spe=spe,
            test_idx=split_idx$test_idx,
            verbose=verbose
        )

        return(c(
            list(
                split_type="random_split",
                k_folds=k_folds,
                seed=seed,
                train_fraction=train_fraction,
                test_fraction=1 - train_fraction,
                n_cells=n_cells,
                n_train=length(split_idx$train_idx),
                n_test=length(split_idx$test_idx)
            ),
            res
        ))
    }

    # K-fold mode: assign each cell to a fold once.
    fold_id <- make_fold_ids(n_cells=n_cells, k_folds=k_folds, seed=seed)

    fold_results <- lapply(seq_len(k_folds), function(i) {
        if (verbose) {
            message("Running fold ", i, " of ", k_folds)
        }

        # Fold i is the test set; all other folds are the training set.
        run_one_fold(
            spe=spe,
            test_idx=which(fold_id == i),
            verbose=verbose
        )
    })

    test_mean_auc <- mean(unlist(lapply(fold_results, function(x) x$test_auc)), na.rm=TRUE)
    test_label_errors <- vapply(
        fold_results,
        extract_optional_chr1,
        character(1),
        name = "test_label_error"
    )

    train_label_errors <- vapply(
        fold_results,
        extract_optional_chr1,
        character(1),
        name = "train_label_error"
    )

    # if you only want real errors later, drop missing values explicitly
    test_label_errors <- stats::na.omit(test_label_errors)
    train_label_errors <- stats::na.omit(train_label_errors)

    list(
        split_type="k_fold_cv",
        k_folds=k_folds,
        seed=seed,
        n_cells=n_cells,
        fold_assignment=fold_id,
        test_mean_auc=test_mean_auc,
        test_label_errors=test_label_errors,
        train_label_errors=train_label_errors,
        folds=fold_results
    )
}

seeds <- c(1998, 1995, 1993, 2000, 1713)

res_list <- list()
for(i in 1:length(seeds)){
    res_list[[i]] <- evaluate_qs_split(
    data_dir="D:\\CosMx_DBKERO_input_data",
    sample_name="CosMx_DBKERO",
    seed=seeds[i],
    k_folds=1,
    train_fraction=0.5,
    verbose=TRUE)
}

# res_kfold_5$folds

# Plot ROC curves for a single split or across folds.
#
# In k-fold mode the function overlays one ROC curve per fold and reports the
# average AUC in the subtitle.
plot_roc_eval <- function(result) {

    roc_list <- lapply(seq_along(result), function(i) {
        df_eval <- result[[i]]$test_roc
        roc_df <- data.frame(
                fpr = 1 - df_eval$specificities,
                tpr = df_eval$sensitivities)
        if (nrow(roc_df) == 0L) {
            return(NULL)
        }
        roc_df$fold <- i
        roc_df
    })

    roc_list <- Filter(Negate(is.null), roc_list)

    if (length(roc_list) == 0L) {
        stop("No ROC curves available for the selected dataset")
    }

    roc_df <- do.call(rbind, roc_list)
    auc_mean <- mean(unlist(lapply(result, function(x) x$test_auc)), na.rm=TRUE)

    return(
        ggplot(roc_df,
            aes(x=fpr, y=tpr, colour=factor(fold), group=fold)) +
            geom_abline(intercept=0, slope=1, linetype=2,
                colour="red4") +
            geom_path(linewidth=0.7, alpha=0.8) +
            coord_fixed() +
            theme_bw() +
            labs(
                title=paste0("ROC curves across folds (test set)"),
                subtitle=paste0("Mean AUC = ", round(auc_mean, 4)),
                x="False positive rate",
                y="True positive rate",
                colour="Fold"
            )
    )
}

pdf(file.path("SuppFigure_24A.pdf"), 
    width = 8, 
    height = 8,
    bg = "transparent")
plot_roc_eval(res_list)
dev.off()

# Supplementary figure 24B

# correlation and R-squared between full dataset QS and pseudoQS computed on test set

# r-squared function
rsq <- function(x, y) summary(lm(y~x))$r.squared

rsq_cor_list <- lapply(res_list, function(x) {
    ref_score <- x$spe_test$QScore
    score <- x$spe_test$QScore_train_model
    rsquared  <-  rsq(ref_score, score)
    spearman <- cor(ref_score,
                score,
                use="complete.obs",
                method="spearman")

    c(rsquared=rsquared, spearman_cor=spearman)
})

# coefficient boxplot across folds
coef_list <- lapply(seq_along(res_list), function(i) {
    coef_value <- res_list[[i]]$model_formula$coefficients[,1][-2]
    coef_value <- as.data.frame(coef_value)
    coef_value$coef_name <- case_when(
        rownames(coef_value) == "(Intercept)" ~ "Intercept",
        rownames(coef_value) == "log2SignalDensity" ~ "Signal density",
        rownames(coef_value) == "Area_um" ~ "Cell size",
        rownames(coef_value) == "I(abs(log2AspectRatio) * as.numeric(dist_border < 50))" ~ "Border effect",
        rownames(coef_value) == "log2Ctrl_total_ratio" ~ "Background signal",
        rownames(coef_value) == "log2SignalDensity:Area_um" ~ "Signal density:Cell size",
        rownames(coef_value) == "log2SignalDensity:I(abs(log2AspectRatio) * as.numeric(dist_border < 50))" ~ "Signal density:border effect",
        rownames(coef_value) == "log2SignalDensity:log2Ctrl_total_ratio" ~ "Signal density:Background signal",
        rownames(coef_value) == "Area_um:I(abs(log2AspectRatio) * as.numeric(dist_border < 50))" ~ "Cell size:border effect",
        rownames(coef_value) == "Area_um:log2Ctrl_total_ratio" ~ "Cell size:Background signal",
        rownames(coef_value) == "I(abs(log2AspectRatio) * as.numeric(dist_border < 50)):log2Ctrl_total_ratio" ~ "Border effect:Background signal"
    )
    coef_value$fold <- i
    coef_value
})

coef_df <- do.call(rbind, coef_list)
coef_df$coef_name <- factor(coef_df$coef_name, levels=c("Intercept", "Cell size", "Signal density", "Background signal", "Border effect", "Signal density:Cell size", "Signal density:border effect", "Signal density:Background signal", "Cell size:border effect", "Cell size:Background signal", "Border effect:Background signal"))

pdf(file.path("SuppFigure_24B.pdf"), 
    width = 10, 
    height = 6,
    bg = "transparent")
ggplot()+
   geom_boxplot(data=coef_df, aes(y=coef_value, x=coef_name), fill="#B3C2F2", outlier.shape=NA, alpha=0.5) +
   geom_jitter(data=coef_df, aes(y=coef_value, x=coef_name), color = "black", width=0.2) +
   geom_hline(yintercept=0, color = "black") +
   guides(fill="none", color=guide_legend(title="Coefficient")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
        title=paste0("Model coefficients across folds"),
        x="Folds",
        y="Coefficient Value"
)
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
#  [1] pROC_1.19.0.1               dplyr_1.1.4                
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
#  [39] labeling_0.4.3            SpatialExperimentIO_1.2.0
#  [41] abind_1.4-8               compiler_4.5.1           
#  [43] proxy_0.4-29              remotes_2.5.0            
#  [45] bit64_4.6.0-1             withr_3.0.2              
#  [47] S7_0.2.1                  backports_1.5.0          
#  [49] BiocParallel_1.44.0       carData_3.0-5            
#  [51] viridis_0.6.5             DBI_1.2.3                
#  [53] pkgbuild_1.4.8            HDF5Array_1.38.0         
#  [55] R.utils_2.13.0            ggsignif_0.6.4           
#  [57] DelayedArray_0.36.0       sessioninfo_1.2.3        
#  [59] rjson_0.2.23              classInt_0.4-11          
#  [61] tools_4.5.1               units_1.0-0              
#  [63] vipor_0.4.7               otel_0.2.0               
#  [65] beeswarm_0.4.0            R.oo_1.27.1              
#  [67] glue_1.8.0                h5mread_1.2.1            
#  [69] rhdf5filters_1.22.0       grid_4.5.1               
#  [71] sf_1.0-24                 gtable_0.3.6             
#  [73] R.methodsS3_1.8.2         class_7.3-23             
#  [75] tidyr_1.3.2               data.table_1.18.0        
#  [77] BiocSingular_1.26.1       ScaledMatrix_1.18.0      
#  [79] car_3.1-3                 XVector_0.50.0           
#  [81] ggrepel_0.9.6             foreach_1.5.2            
#  [83] pillar_1.11.1             limma_3.66.0             
#  [85] robustbase_0.99-6         splines_4.5.1            
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
# [115] purrr_1.2.1               rlang_1.1.7              
# [117] cowplot_1.2.0            
