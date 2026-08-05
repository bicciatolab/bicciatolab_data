# Clean, minimal, and reproducible script to reproduce supplementary figure 23 
# To reproduce the analysis be sure to have the required packages installed

library(SpaceTrooper)
library(ggplot2)
library(patchwork)

# R = 4.5.1

# define ablation formulas

formulas <- list(
    full="~(log2SignalDensity + Area_um + I(abs(log2AspectRatio) * as.numeric(dist_border < 50)) + log2Ctrl_total_ratio)^2",

    no_Area_um=
        "~(log2SignalDensity + I(abs(log2AspectRatio) * as.numeric(dist_border < 50)) + log2Ctrl_total_ratio)^2",

    no_log2SignalDensity=
    "~(Area_um + I(abs(log2AspectRatio) * as.numeric(dist_border < 50)) + log2Ctrl_total_ratio)^2",   

    no_log2Ctrl_total_ratio=
    "~(log2SignalDensity + Area_um + I(abs(log2AspectRatio) * as.numeric(dist_border < 50)))^2", 

    no_border_effect=
        "~(log2SignalDensity + Area_um + log2Ctrl_total_ratio)^2",

    no_interaction_terms=
    "~log2SignalDensity + Area_um + I(abs(log2AspectRatio) * as.numeric(dist_border < 50)) + log2Ctrl_total_ratio"    
)

# load CosMx DBKERO spe

# dirname pointing to directory containing count matrix, fov position,
# polygons and metadata files as downloaded from https://kero.hgc.jp/Breast_Cancer_Spatial.html

dirname <- "CosMx_DBKERO_input_data"
sample_name <- "CosMx_DBKERO"

dbkero_spe <- readCosmxSPE(dirname, sampleName=sample_name)
dbkero_spe <- spatialPerCellQC(dbkero_spe)

# for ablation study, we need to compute the model also in absence of log2SignalDensity
# which would generally result in an error.
# It is necessary to change a function to avoid this error

.qscoreSupportedPredictors <- function() {
    c(
        "log2SignalDensity",
        "Area_um",
        "log2AspectRatio",
        "log2Ctrl_total_ratio"
    )
}

.validateQScoreFormula <- function(modelFormula, dataNames, technology) {
    model_formula <- tryCatch(
        stats::as.formula(modelFormula),
        error=function(e) {
            stop("Unable to parse 'modelFormula': ", conditionMessage(e))
        }
    )
    if (length(model_formula) != 2L) {
        stop("'modelFormula' must be a one-sided formula.")
    }

    required_variables <- base::all.vars(model_formula)
    supported <- .qscoreSupportedPredictors()
    allowed_variables <- c(supported, "dist_border")
    unsupported_variables <- setdiff(
        required_variables,
        allowed_variables
    )
    if (length(unsupported_variables) > 0L) {
        stop(
            "Unsupported Quality Score predictor(s): ",
            paste(unsupported_variables, collapse=", "),
            ". Formula structure may be customised, but predictors and ",
            "transformations are restricted to: ",
            paste(supported, collapse=", "),
            "."
        )
    }

    term_labels <- attr(stats::terms(model_formula), "term.labels")
    if (length(term_labels) == 0L) {
        stop("'modelFormula' must contain at least one supported predictor.")
    }
    normalise <- function(x) gsub("[[:space:]]+", "", x)
    border_term <- .qscoreBorderTerm()
    unsupported_terms <- character()
    for (term_label in term_labels) {
        components <- strsplit(
            normalise(term_label),
            ":",
            fixed=TRUE
        )[[1]]
        valid_components <- components %in% supported |
            components == border_term
        if (!all(valid_components)) {
            unsupported_terms <- c(unsupported_terms, term_label)
        }
    }
    unsupported_terms <- unique(unsupported_terms)
    if (length(unsupported_terms) > 0L) {
        stop(
            "Unsupported Quality Score formula term(s) or transformation(s): ",
            paste(unsupported_terms, collapse=", "),
            ". Formula structure and interactions may be customised, but ",
            "predictors and transformations may not be extended. Supported ",
            "predictors are: ",
            paste(supported, collapse=", "),
            "; the supported border-effect expression is ",
            "'I(abs(log2AspectRatio) * as.numeric(dist_border < 50))'."
        )
    }

    missing_variables <- setdiff(required_variables, dataNames)
    if (length(missing_variables) > 0L) {
        stop(
            "Missing variables required by 'modelFormula': ",
            paste(missing_variables, collapse=", "),
            ". Run 'spatialPerCellQC()' before computing Quality Score."
        )
    }

    is_cosmx <- technology %in%
        c("Nanostring_CosMx", "Nanostring_CosMx_Protein")
    uses_border_predictor <- any(
        c("log2AspectRatio", "dist_border") %in% required_variables
    )
    if (uses_border_predictor && !is_cosmx) {
        stop(
            "'log2AspectRatio' and the border-effect expression are ",
            "supported only for Nanostring CosMx datasets; detected ",
            "technology: '", technology, "'."
        )
    }

    formula_text <- if (is.character(modelFormula)) {
        paste(modelFormula, collapse=" ")
    } else {
        paste(deparse(model_formula, width.cutoff=500L), collapse=" ")
    }
    training_metrics <- required_variables[
        required_variables %in% supported
    ]
    list(
        formula=model_formula,
        text=formula_text,
        variables=required_variables,
        terms=term_labels,
        training_metrics=training_metrics
    )
}

checkOutliers <- function(spe, verbose=FALSE) {
    warnstopmsg <- function(var, warnstop=c("w","s")) {
        warnstop <- match.arg(warnstop)
        m1 <- paste0("Not enough outlier cells for ", var, ".\n")
        m2 <- switch(warnstop,
                s="In this version of SpaceTrooper, presence of outliers
                for log2SignalDensity is required to compute Quality Score.",
                w="This variable will not be used in the final formula")
        return(paste0(m1, m2))
    }
    out_var <- metadata(spe)$formula_variables
    cd <- colData(spe)
    if (verbose) {
        for (i in names(out_var)) {
            message("Outliers found for ", i, ":")
            outlier_counts <- table(cd[[out_var[i]]])
            for (k in seq_along(outlier_counts)) {
                message(names(outlier_counts)[k], ": ", outlier_counts[k])
            }
        }
    }
    # stopifnot(
    #     "log2SignalDensity is not included in the Quality Score formula.\n
    #     Quality Score cannot be computed"=
    #         "log2SignalDensity" %in% names(out_var)
    # )
    cfg <- list(
        log2SignalDensity=list(pattern="log2SignalDensity_outlier",
            remove="log2SignalDensity_outlier_train", label="LOW", act=stop,
            code="s"),
        Area_um=list(pattern="Area_um_outlier", remove="Area_um_outlier",
            label="HIGH", act=warning, code="w"),
        log2Ctrl_total_ratio=list(pattern="log2Ctrl_total_ratio_outlier",
            remove="log2Ctrl_total_ratio_outlier_train", label="HIGH",
            act=warning, code="w")
    )
    for (v in intersect(names(cfg), names(out_var))) {
        r <- cfg[[v]]
        gi <- grep(r$pattern, out_var)
        if (length(gi)==0L) next
        col <- out_var[gi]
        labs <- factor(cd[[col]], levels=c("LOW","HIGH","NO"))
        tab <- table(labs)
        cnt <- tab[r$label]
        if (is.null(cnt) || is.na(cnt)) cnt <- 0L
        if (cnt < ncol(spe)*0.001) {
            r$act(warnstopmsg(v, r$code))
            out_var <- out_var[-grep(r$remove, out_var)]
        }
    }
    var <- "log2AspectRatio"
    pat <- paste0(var, "_outlier")
    is_cosmx <- metadata(spe)$technology %in%
                c("Nanostring_CosMx", "Nanostring_CosMx_Protein")
    idx <- grep(pat, out_var)            
    if (is_cosmx && (var %in% names(out_var))) {
        if (length(idx)) {
            col <- out_var[idx]
            labs <- factor(cd[[col]], levels=c("LOW","HIGH","NO"))
            tab <- table(labs)
            nmin <- ncol(spe) * 0.001

            low  <- tab[["LOW"]];  if (is.na(low))  low  <- 0L
            high <- tab[["HIGH"]]; if (is.na(high)) high <- 0L

            if (low < nmin && high < nmin) {
                warning(warnstopmsg(var, "w"))
                out_var <- out_var[-idx]
            }
        }
    } else {
        if (length(idx)) {
            out_var <- out_var[-idx]
        }
    }
    metadata(spe)$formula_variables <- out_var
    return(spe)
}

.prepQCContext <- function(spe, metricList=c("log2SignalDensity", "Area_um",
    "log2AspectRatio", "log2Ctrl_total_ratio"), verbose=FALSE) {

    stopifnot(is(spe,"SpatialExperiment"))

    # remove zero-count cells once
    zerocells <- spe$total==0
    if (sum(zerocells) > 0) {
        warning(sum(zerocells),
            " cells with 0 counts were found. These cells will be removed.")
        spe <- spe[, !zerocells]
    }
    if("log2CountArea" %in% names(colData(spe)))
    {
        stop("log2CountArea found in colData.\n",
        "Please updated the object with the latest version of spatialPerCellQC.")
    }
    spe1 <- computeOutliersQScore(spe, metricList)
    spe1 <- checkOutliers(spe1, verbose)

    out_var <- metadata(spe1)$formula_variables
    df <- as.data.frame(colData(spe1))
    tech <- metadata(spe1)$technology

    return(list(df=df, out_var=out_var, tech=tech))
}

.filterCompleteModelCases <- function(df, modelFormula, response=NULL,
    context="cells") {

    vars <- base::all.vars(stats::as.formula(modelFormula))

    missing_vars <- setdiff(vars, colnames(df))

    if (length(missing_vars) > 0L) {
        stop(
            "Missing variables required by model formula: ",
            paste(missing_vars, collapse=", ")
        )
    }
    ## Missing values may come from variables used in `modelFormula`, for
    ## example `log2AspectRatio` when aspect ratio cannot be computed from
    ## polygons.
    ok <- complete.cases(df[, vars, drop=FALSE])

    if (!is.null(response)) {
        ok <- ok & !is.na(response)
    }

    if (!all(ok)) {
        warning(
            sum(!ok),
            " ", context,
            " contain missing model variables."
        )
    }

    return(ok)
}

# I was wrong, also this function checks for log2SignalDensity absence

computeTrainDF <- function(colData, formulaVars, tech, verbose=FALSE) {
    out_var <- formulaVars
    df <- as.data.frame(colData)
    train_bad_var <- character()
    train_good_var <- character()

    # stopifnot("log2SignalDensity is not included in the Quality Score formula.\nQuality Score cannot be computed"="log2SignalDensity" %in% names(out_var))

    cfg <- list(
        log2SignalDensity=list(bad="LOW", good=c(0.90,0.99)),
        Area_um=list(bad="HIGH", good=c(0.25,0.75)),
        log2Ctrl_total_ratio=list(bad="HIGH", good=NULL)
    )

    vars <- intersect(names(cfg), names(out_var))

    for (v in vars) {
        r <- cfg[[v]]
        out_col <- out_var[names(out_var) == v]
        if (length(out_col) == 0L) next

        # BAD ids
        bad_ids <- df |>
            dplyr::filter(.data[[out_col]] == r$bad) |>
            dplyr::pull(cell_id)
        train_bad_var <- unique(c(train_bad_var, bad_ids))

        # GOOD ids (if a band is defined)
        if (!is.null(r$good)) {
            qs <- stats::quantile(df[[v]], probs=r$good)
            good_ids <- df |>
                dplyr::filter(.data[[v]] > qs[1] & .data[[v]] < qs[2]) |>
                dplyr::pull(cell_id)
            train_good_var <- unique(c(train_good_var, good_ids))
        }
    }

    # CosMx-specific handling for log2AspectRatio
    is_cosmx <- tech %in% c("Nanostring_CosMx", "Nanostring_CosMx_Protein")

    if (is_cosmx && "log2AspectRatio" %in% names(out_var)) {

        if (!("dist_border" %in% colnames(df))) {
            stop("dist_border column is required for CosMx handling")
        }

        out_col <- out_var[names(out_var) == "log2AspectRatio"]

        bad_ids <- df |>
            dplyr::filter(.data[[out_col]] %in% c("HIGH","LOW") &
                            dist_border < 50) |>
            dplyr::pull(cell_id)
        train_bad_var <- unique(c(train_bad_var, bad_ids))

        qs <- stats::quantile(df$log2AspectRatio, probs=c(0.25,0.75))
        good_ids <- df |>
            dplyr::filter(log2AspectRatio > qs[1] &
                            log2AspectRatio < qs[2] &
                            dist_border > 50) |>
            dplyr::pull(cell_id)
        train_good_var <- unique(c(train_good_var, good_ids))

        idx <- grep(pattern="log2AspectRatio_outlier", x=out_var)
        if (length(idx)) {
            names(out_var)[idx] <-
                "I(abs(log2AspectRatio) * as.numeric(dist_border < 50))"
        }
    }

    train_bad <- df |>
        dplyr::filter(cell_id %in% train_bad_var) |>
        dplyr::mutate(QScore_train=0)

    train_good <- df |>
        dplyr::filter(cell_id %in% train_good_var) |>
        dplyr::mutate(
            QScore_train=1,
            is_a_bad_boy=cell_id %in% train_bad$cell_id
        )

    train_bad <- train_bad |>
        dplyr::distinct(cell_id, .keep_all=TRUE)

    if (verbose) message(paste0("Chosen low qual examples: ", nrow(train_bad)))

    # good example duplicates removal without any warning to the user
    train_good <- train_good |>
        dplyr::distinct(cell_id, .keep_all=TRUE)

    train_good <- train_good[!train_good$is_a_bad_boy, ]

    n_bad <- nrow(train_bad)
    stopifnot("Not enough good examples to match bad examples"=
                nrow(train_good) > n_bad)
    idx <- sample(seq_len(nrow(train_good)), n_bad, replace=FALSE)
    train_good <- train_good[idx, ]

    if (verbose) {
        message(paste0(
            "Chosen good quality examples (should match bad): ",
            nrow(train_good)
        ))
    }

    train_good$is_a_bad_boy <- NULL
    train_df <- rbind(train_bad, train_good)
    train_df <- train_df |>
        dplyr::distinct(cell_id, .keep_all=TRUE)
    ## `qcscore_train` remains available because `computeTrainDF()` is a
    ## published API. Canonical internals use `QScore_train`.
    train_df$qcscore_train <- train_df$QScore_train
    return(train_df)
}

computeQScore <- function(spe, bestLambda=NULL, modelFormula=NULL, verbose=FALSE) {
    stopifnot(is(spe, "SpatialExperiment"))
    if (dim(spe[,spe$total == 0])[2] != 0) {
        warning(paste0(dim(spe[,spe$total == 0])[2],
            " cells with 0 counts were found. These cells will be removed."))
        spe <- spe[,spe$total > 0]
    }
    formula_info <- NULL
    if (is.null(modelFormula)) {
        metric_list <- intersect(
            .qscoreSupportedPredictors(),
            names(colData(spe))
        )
    } else {
        formula_info <- .validateQScoreFormula(
            modelFormula=modelFormula,
            dataNames=names(colData(spe)),
            technology=metadata(spe)$technology
        )
        metric_list <- formula_info$training_metrics
    }

    # if (!"log2SignalDensity" %in% metric_list) {
    #     stop(
    #         "'log2SignalDensity' is required to construct Quality Score ",
    #         "training labels."
    #     )
    # }

    ctx <- .prepQCContext(spe, metric_list, verbose)
    df <- ctx$df; out_var <- ctx$out_var; tech <- ctx$tech

    if (is.null(modelFormula)) {
        model_formula <- getModelFormula(names(out_var))
        model_formula_object <- stats::as.formula(model_formula)
    } else {
        model_formula <- formula_info$text
        model_formula_object <- formula_info$formula
    }

    if (verbose) {
        message("Using final model formula:")
        message(model_formula)
    }

    train_df <- computeTrainDF(df, out_var, tech, verbose)

    train_ok <- .filterCompleteModelCases(
        df=train_df,
        modelFormula=model_formula_object,
        response=train_df$QScore_train,
        context="training cells"
    )

    train_df <- train_df[train_ok, , drop=FALSE]

    model_matrix <- stats::model.matrix(
        model_formula_object,
        data=train_df
    )
    model <- trainModel(model_matrix, train_df)
    
    if(is.null(bestLambda)) {
        bestLambda <- .computeLambda(
            modelMatrix=model_matrix,
            response=train_df$QScore_train
        )
    }

    coefs <- coef(model, s=bestLambda)

    if (verbose) {
        message("Model coefficients for every term used in the formula:")
        message( paste( rownames(coefs), round(as.numeric(coefs), 2), sep="=",
            collapse=" "))
    }

    full_ok <- .filterCompleteModelCases(
        df=df,
        modelFormula=model_formula_object,
        context="cells")

    full_matrix <- stats::model.matrix(
        model_formula_object,
        data=df[full_ok, , drop=FALSE]
    )

    ## NAs may come from model variables used in `model_formula`, e.g. cells with
    ## missing `log2AspectRatio` when aspect ratio cannot be computed from polygons.
    ## Predict only complete cases; incomplete cells keep QScore = NA.
    qscore <- rep(NA_real_, nrow(df))
    qscore[full_ok] <- as.vector(
        predict(model, s=bestLambda, newx=full_matrix, type="response")
    )

    spe$QScore <- qscore

    metadata(spe)$QScore_model <- list(
        model=model,
        bestLambda=bestLambda,
        model_formula=model_formula,
        formula_variables=out_var,
        model_matrix_colnames=colnames(model_matrix),
        coefficients=coefs,
        coefficients_table=data.frame(
            term=rownames(coefs),
            coefficient=as.numeric(coefs),
            row.names=NULL
        )
    )
    return(spe)
}

# define ablation function

run_qs_ablation <- function(spe, formulas, verbose=TRUE) {
    out <- list()

    for (nm in names(formulas)) {
        message("Running model: ", nm)

        spe_i <- computeQScore(
            spe=spe,
            modelFormula=formulas[[nm]],
            verbose=verbose
        )

        out[[nm]] <- list(
            spe=spe_i,
            score=spe_i$QScore,
            model=metadata(spe_i)$QScore_model
        )
    }

    return(out)
}

# run ablation on CosMx DBKERO

set.seed(1998)

abl_spe_list <- run_qs_ablation(
    spe=dbkero_spe,
    formulas=formulas,
    verbose=TRUE
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
        geom_point(alpha=0.35, size=0.6, color="#c0c8cf") +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
        labs(x=x_label, y=y_label, title=title) +
        theme_minimal() +
        annotate(
            "text",
            x=quantile(df$x, 0.99, na.rm=TRUE),
            y=quantile(df$y, 0.01, na.rm=TRUE),
            label=paste0("R² = ",
                formatC(rsq, digits=3, format="f")),
            hjust=1,
            vjust=0,
            color="black"
        ) +
        annotate(
            "text",
            x=quantile(df$x, 0.99, na.rm=TRUE),
            y=quantile(df$y, 0.005, na.rm=TRUE),
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

ablation_names <- setdiff(names(abl_spe_list), "full")

plots_cosmx <- lapply(ablation_names, function(nm) {
    plot_qscore_comparison(
        x=abl_spe_list$full$score,
        y=abl_spe_list[[nm]]$score,
        x_label="Full QS",
        y_label=paste0("QS: ", nm),
        title=paste0("CosMx: full QS vs ", nm),
        cor_method="spearman"
    )
})

names(plots_cosmx) <- ablation_names

plot_list <- lapply(plots_cosmx, function(x) x$plot)

# --- Build rows ---
row1 <- plot_list[[1]] | plot_list[[2]] 
row2 <- plot_list[[3]] | plot_list[[4]] 
row3 <- plot_list[[5]] | plot_spacer()

final_plot <- (row1 / plot_spacer() /
                 row2 / plot_spacer() /
                 row3) +
  plot_layout(
    heights = c(1, 0.01, 1, 0.01, 1),  # smaller spacing between rows
    widths = c(1, 1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "."
  ) &
  theme(plot.tag = element_text(size = 14, face = "bold"))

## Print plots
pdf(file.path("SuppFigure_23.pdf"), 
    width = 10, 
    height = 15,
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
#  [1] patchwork_1.3.2             ggplot2_4.0.1              
#  [3] SpaceTrooper_1.1.9          testthat_3.3.2             
#  [5] SpatialExperiment_1.20.0    SingleCellExperiment_1.32.0
#  [7] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#  [9] GenomicRanges_1.62.1        Seqinfo_1.0.0              
# [11] IRanges_2.44.0              S4Vectors_0.48.0           
# [13] BiocGenerics_0.56.0         generics_0.1.4             
# [15] MatrixGenerics_1.22.0       matrixStats_1.5.0          

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
 


