# ═══════════════════════════════════════════════════════════════════════════════
# run_lodo_validation
#
# Leave-One-Dataset-Out (LODO) cross-dataset validation.
#
# For each held-out dataset:
#   - Train glmnet, glm, and RF on all *other* datasets combined
#   - Predict on the held-out dataset
#   - Compute AUC, BalAcc(opt), sensitivity, specificity per model
#
# This tests whether biomarker signals generalise across cohorts, which is a
# stronger validity claim than within-cohort LOOCV. Age and sex are included
# as covariates in all models — they are risk factors, not nuisance variables,
# so they remain in the model and their contribution is part of the signal.
#
# Arguments:
#   dataset_list  — named list, each element is list(X, y, meta_use)
#                   X: sample × gene matrix (same genes across all datasets)
#                   y: binary outcome vector (1=AD, 0=CTRL)
#                   meta_use: data.frame with Age, Sex columns
#   ensembl_ids   — character vector of all gene Ensembl IDs (cols of X)
#   hgnc_symbols  — matching HGNC symbols
#   output_dir    — output directory
#   today         — date string for filenames
#   alpha_val     — glmnet alpha (used if not tuning; tuning is always on)
#   B_boot        — bootstrap iterations for AUC CI
#   gene_chr_map  — named character vector: names=ensembl_id, values=chromosome
#   num.trees     — ranger num.trees
#   node_grid     — RF min.node.size search grid
#   alpha_grid    — glmnet alpha search grid
#   nfolds_inner  — inner CV folds for tuning on combined training data
#
# Returns a list with:
#   lodo_summary     — per held-out dataset × model performance table
#   lodo_gene_scores — per gene × held-out dataset AUC (from single-gene models)
#   p_lodo_summary   — ggplot comparison across datasets
# ═══════════════════════════════════════════════════════════════════════════════

run_lodo_validation <- function(dataset_list,
                                ensembl_ids,
                                hgnc_symbols,
                                output_dir,
                                today,
                                alpha_val    = 0.5,
                                B_boot       = 100,
                                gene_chr_map = NULL,
                                num.trees    = 500,
                                node_grid    = c(1, 3, 5, 10),
                                alpha_grid   = c(0, 0.25, 0.5, 0.75, 1.0),
                                nfolds_inner = 5) {

  library(ranger)
  dataset_names <- names(dataset_list)
  if (is.null(dataset_names) || any(dataset_names == ""))
    stop("dataset_list must be a fully named list (names = dataset labels)")
  if (length(dataset_list) < 3)
    warning("LODO with fewer than 3 datasets gives limited generalisation evidence")

  # Helper: align gene columns across datasets — intersect only
  align_genes <- function(mats) {
    common <- Reduce(intersect, lapply(mats, colnames))
    lapply(mats, function(m) m[, common, drop = FALSE])
  }

  # ── Per-gene LODO: single-gene models trained on n-1 datasets ─────────────
  # For each gene, for each held-out dataset:
  #   fit glmnet + glm on pooled training data → predict on held-out
  cat("\n── LODO: per-gene single-gene models ──\n")

  # Build full gene list from intersection
  all_X <- lapply(dataset_list, `[[`, "X")
  common_genes <- Reduce(intersect, lapply(all_X, colnames))
  cat("  Genes in common across all datasets:", length(common_genes), "\n")

  # Results: gene × held-out-dataset matrix for AUC (glmnet) and AUC (glm)
  lodo_gene_rows <- list()

  for (g in common_genes) {
    sym <- hgnc_symbols[match(g, ensembl_ids)]
    sym <- if (!is.na(sym)) sym else g

    chr_annot  <- if (!is.null(gene_chr_map)) gene_chr_map[g] else NA_character_
    on_sex_chr <- !is.na(chr_annot) & chr_annot %in% c("X", "chrX", "Y", "chrY")

    row_vals <- list(ensembl_id = g, gene_symbol = sym,
                     chromosome = ifelse(is.na(chr_annot), "unknown", chr_annot),
                     on_sex_chr = on_sex_chr)

    for (held_out in dataset_names) {
      train_names <- setdiff(dataset_names, held_out)

      # Pool training data
      train_X_list   <- lapply(train_names, function(d) dataset_list[[d]]$X[, g, drop = FALSE])
      train_y_list   <- lapply(train_names, function(d) dataset_list[[d]]$y)
      train_meta_list <- lapply(train_names, function(d) dataset_list[[d]]$meta_use)

      X_tr_gene  <- do.call(rbind, train_X_list)[, 1]
      y_tr       <- do.call(c, train_y_list)
      age_tr     <- as.numeric(do.call(c, lapply(train_meta_list, `[[`, "Age")))
      sex_tr     <- ifelse(toupper(do.call(c, lapply(train_meta_list, `[[`, "Sex"))) == "M", 1, 0)

      # Test data
      X_te_gene  <- dataset_list[[held_out]]$X[, g]
      y_te       <- dataset_list[[held_out]]$y
      age_te     <- as.numeric(dataset_list[[held_out]]$meta_use$Age)
      sex_te     <- ifelse(toupper(dataset_list[[held_out]]$meta_use$Sex) == "M", 1, 0)

      if (var(X_tr_gene) == 0 || length(unique(y_tr)) < 2) {
        row_vals[[paste0("auc_enet_", held_out)]] <- NA_real_
        row_vals[[paste0("auc_glm_",  held_out)]] <- NA_real_
        next
      }

      # glmnet — tune alpha on training pool
      X_tr_mat <- cbind(gene = X_tr_gene, age = age_tr, sex = sex_tr)
      X_te_mat <- matrix(c(X_te_gene, age_te, sex_te), ncol = 3,
                         dimnames = list(NULL, c("gene", "age", "sex")))
      pf       <- c(gene = 1, age = 0, sex = 0)

      tune_res <- tune_glmnet_alpha(X_tr_mat, y_tr,
                                    alpha_grid     = alpha_grid,
                                    penalty.factor = pf,
                                    nfolds_inner   = min(nfolds_inner, length(y_tr)))
      fit_enet <- tryCatch(
        cv.glmnet(x = X_tr_mat, y = y_tr, family = "binomial",
                  alpha = tune_res$alpha, penalty.factor = pf,
                  nfolds = min(nfolds_inner, length(y_tr)),
                  type.measure = "class", standardize = TRUE),
        error = function(e) NULL
      )
      pred_enet_te <- if (!is.null(fit_enet)) {
        gc_min <- tryCatch(as.numeric(coef(fit_enet, s = "lambda.min")["gene", ]), error = function(e) 0)
        gc_1se <- tryCatch(as.numeric(coef(fit_enet, s = "lambda.1se")["gene", ]), error = function(e) 0)
        s_use  <- if (!is.na(gc_min) && gc_min != 0) "lambda.min" else
                  if (!is.na(gc_1se) && gc_1se != 0) "lambda.1se" else NA
        if (is.na(s_use)) rep(0.5, length(y_te)) else
          tryCatch(as.numeric(predict(fit_enet, newx = X_te_mat, s = s_use, type = "response")),
                   error = function(e) rep(0.5, length(y_te)))
      } else rep(0.5, length(y_te))

      # glm
      has_sex_var <- length(unique(sex_tr)) > 1
      df_tr_glm   <- data.frame(y = y_tr, gene = X_tr_gene, age = age_tr, sex = sex_tr)
      df_te_glm   <- data.frame(gene = X_te_gene, age = age_te, sex = sex_te)
      formula_glm <- if (has_sex_var) y ~ gene + age + sex else y ~ gene + age
      fit_glm     <- tryCatch(
        suppressWarnings(glm(formula_glm, data = df_tr_glm, family = binomial())),
        error = function(e) NULL
      )
      pred_glm_te <- if (!is.null(fit_glm)) {
        tryCatch(as.numeric(predict(fit_glm, newdata = df_te_glm, type = "response")),
                 error = function(e) rep(0.5, length(y_te)))
      } else rep(0.5, length(y_te))

      # AUC and BalAcc on held-out
      auc_enet_ho <- tryCatch(
        as.numeric(auc(roc(response = y_te, predictor = pred_enet_te, quiet = TRUE))),
        error = function(e) NA_real_
      )
      auc_glm_ho <- tryCatch(
        as.numeric(auc(roc(response = y_te, predictor = pred_glm_te, quiet = TRUE))),
        error = function(e) NA_real_
      )
      balacc_enet_ho <- find_opt_threshold(pred_enet_te, y_te)$balacc
      balacc_glm_ho  <- find_opt_threshold(pred_glm_te,  y_te)$balacc

      row_vals[[paste0("auc_enet_",    held_out)]] <- round(auc_enet_ho,    3)
      row_vals[[paste0("auc_glm_",     held_out)]] <- round(auc_glm_ho,     3)
      row_vals[[paste0("balacc_enet_", held_out)]] <- round(balacc_enet_ho, 3)
      row_vals[[paste0("balacc_glm_",  held_out)]] <- round(balacc_glm_ho,  3)
    }
    lodo_gene_rows[[g]] <- as.data.frame(row_vals, stringsAsFactors = FALSE)
  }

  lodo_gene_scores <- do.call(rbind, lodo_gene_rows)

  # Mean AUC columns across held-out datasets (per model)
  enet_cols <- paste0("auc_enet_", dataset_names)
  glm_cols  <- paste0("auc_glm_",  dataset_names)
  # Mean BalAcc and AUC across held-out datasets (per model)
  enet_balacc_cols <- paste0("balacc_enet_", dataset_names)
  glm_balacc_cols  <- paste0("balacc_glm_",  dataset_names)

  lodo_gene_scores$mean_auc_enet_lodo    <- rowMeans(lodo_gene_scores[, enet_cols], na.rm = TRUE)
  lodo_gene_scores$mean_auc_glm_lodo     <- rowMeans(lodo_gene_scores[, glm_cols],  na.rm = TRUE)
  lodo_gene_scores$mean_auc_lodo         <- rowMeans(
    lodo_gene_scores[, c(enet_cols, glm_cols)], na.rm = TRUE
  )
  # BalAcc columns (added in loop below — back-fill if present)
  if (all(enet_balacc_cols %in% names(lodo_gene_scores))) {
    lodo_gene_scores$mean_balacc_enet_lodo <- rowMeans(lodo_gene_scores[, enet_balacc_cols], na.rm = TRUE)
    lodo_gene_scores$mean_balacc_glm_lodo  <- rowMeans(lodo_gene_scores[, glm_balacc_cols],  na.rm = TRUE)
    lodo_gene_scores$mean_balacc_lodo      <- rowMeans(
      lodo_gene_scores[, c(enet_balacc_cols, glm_balacc_cols)], na.rm = TRUE
    )
    # Primary sort: mean BalAcc; tiebreaker: mean AUC
    lodo_gene_scores <- lodo_gene_scores[order(-lodo_gene_scores$mean_balacc_lodo,
                                               -lodo_gene_scores$mean_auc_lodo), ]
  } else {
    lodo_gene_scores <- lodo_gene_scores[order(-lodo_gene_scores$mean_auc_lodo), ]
  }

  write.csv(lodo_gene_scores,
            file.path(output_dir, paste0("lodo_per_gene_auc_", today, ".csv")),
            row.names = FALSE)

  # ── Multi-gene LODO: train on pooled n-1 datasets, test on held-out ───────
  cat("\n── LODO: multi-gene models ──\n")

  lodo_summary_rows <- list()

  for (held_out in dataset_names) {
    cat("  Held-out:", held_out, "\n")
    train_names <- setdiff(dataset_names, held_out)

    # Pool and align training matrices
    train_X_aligned <- align_genes(lapply(train_names, function(d) dataset_list[[d]]$X))
    common_g <- colnames(train_X_aligned[[1]])

    X_tr       <- do.call(rbind, train_X_aligned)
    y_tr       <- do.call(c, lapply(train_names, function(d) dataset_list[[d]]$y))
    age_tr     <- as.numeric(do.call(c, lapply(train_names,
                                               function(d) dataset_list[[d]]$meta_use$Age)))
    sex_tr     <- ifelse(toupper(do.call(c, lapply(train_names,
                                                   function(d) dataset_list[[d]]$meta_use$Sex))) == "M", 1, 0)

    # Test
    X_te   <- dataset_list[[held_out]]$X[, common_g, drop = FALSE]
    y_te   <- dataset_list[[held_out]]$y
    age_te <- as.numeric(dataset_list[[held_out]]$meta_use$Age)
    sex_te <- ifelse(toupper(dataset_list[[held_out]]$meta_use$Sex) == "M", 1, 0)
    n_te   <- length(y_te)

    cov_tr     <- cbind(Age = age_tr, Sex = sex_tr)
    cov_te     <- cbind(Age = age_te, Sex = sex_te)
    X_tr_full  <- cbind(X_tr, cov_tr)
    X_te_full  <- cbind(X_te, cov_te)
    pf_multi   <- c(rep(1, ncol(X_tr)), rep(0, ncol(cov_tr)))

    # ── Multi-gene glmnet ──────────────────────────────────────────────────
    tune_multi_enet <- tune_glmnet_alpha(X_tr_full, y_tr,
                                         alpha_grid     = alpha_grid,
                                         penalty.factor = pf_multi,
                                         nfolds_inner   = min(nfolds_inner, length(y_tr)))
    fit_multi_enet <- tryCatch(
      cv.glmnet(x = X_tr_full, y = y_tr, family = "binomial",
                alpha = tune_multi_enet$alpha, penalty.factor = pf_multi,
                nfolds = min(nfolds_inner, length(y_tr)),
                type.measure = "class", standardize = TRUE),
      error = function(e) NULL
    )
    pred_multi_enet <- if (!is.null(fit_multi_enet)) {
      tryCatch(as.numeric(predict(fit_multi_enet, newx = X_te_full,
                                  s = "lambda.min", type = "response")),
               error = function(e) rep(0.5, n_te))
    } else rep(0.5, n_te)

    # ── Multi-gene RF ──────────────────────────────────────────────────────
    df_tr_rf      <- as.data.frame(X_tr_full)
    df_te_rf      <- as.data.frame(X_te_full)
    tune_multi_rf <- tune_rf(df_tr_rf, y_tr,
                             node_grid    = node_grid,
                             nfolds_inner = min(nfolds_inner, length(y_tr)),
                             num.trees    = num.trees)
    df_tr_rf$y <- factor(y_tr, levels = c(0, 1))
    fit_multi_rf <- tryCatch(
      ranger::ranger(
        y             ~ .,
        data          = df_tr_rf,
        num.trees     = num.trees,
        mtry          = max(1, floor(sqrt(ncol(df_tr_rf) - 1))),
        min.node.size = tune_multi_rf$min.node.size,
        probability   = TRUE,
        importance    = "permutation",
        num.threads   = 1
      ),
      error = function(e) NULL
    )
    pred_multi_rf <- if (!is.null(fit_multi_rf)) {
      tryCatch(predict(fit_multi_rf, data = df_te_rf)$predictions[, "1"],
               error = function(e) rep(0.5, n_te))
    } else rep(0.5, n_te)

    # ── Null model on training data → predict held-out ────────────────────
    df_null_tr <- data.frame(y = y_tr, age = age_tr, sex = sex_tr)
    df_null_te <- data.frame(age = age_te, sex = sex_te)
    has_sex_v  <- length(unique(sex_tr)) > 1
    fit_null   <- tryCatch(
      suppressWarnings(glm(if (has_sex_v) y ~ age + sex else y ~ age,
                           data = df_null_tr, family = binomial())),
      error = function(e) NULL
    )
    pred_null_te <- if (!is.null(fit_null)) {
      tryCatch(as.numeric(predict(fit_null, newdata = df_null_te, type = "response")),
               error = function(e) rep(0.5, n_te))
    } else rep(0.5, n_te)

    # ── Compute metrics for all three + null ──────────────────────────────
    compute_metrics <- function(pred, y, label_str) {
      roc_obj  <- tryCatch(roc(response = y, predictor = pred, quiet = TRUE),
                           error = function(e) NULL)
      auc_val  <- if (!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA_real_
      opt_res  <- find_opt_threshold(pred, y)
      data.frame(
        model           = label_str,
        auc             = round(auc_val,            3),
        balacc_opt      = round(opt_res$balacc,     3),
        sensitivity_opt = round(opt_res$sensitivity, 3),
        specificity_opt = round(opt_res$specificity, 3),
        opt_threshold   = round(opt_res$threshold,  3),
        n_test          = length(y),
        n_AD_test       = sum(y == 1),
        n_CTRL_test     = sum(y == 0),
        stringsAsFactors = FALSE
      )
    }

    held_rows <- rbind(
      compute_metrics(pred_multi_enet, y_te, "glmnet_multi"),
      compute_metrics(pred_multi_rf,   y_te, "RF_multi"),
      compute_metrics(pred_null_te,    y_te, "null_age_sex")
    )
    held_rows$held_out_dataset   <- held_out
    held_rows$train_datasets     <- paste(train_names, collapse = "+")
    held_rows$best_alpha_glmnet  <- tune_multi_enet$alpha
    held_rows$best_node_rf       <- tune_multi_rf$min.node.size
    lodo_summary_rows[[held_out]] <- held_rows

    cat("    glmnet AUC:", round(held_rows$auc[held_rows$model == "glmnet_multi"], 3),
        " RF AUC:", round(held_rows$auc[held_rows$model == "RF_multi"], 3),
        " null AUC:", round(held_rows$auc[held_rows$model == "null_age_sex"], 3), "\n")
  }

  lodo_summary <- do.call(rbind, lodo_summary_rows)
  rownames(lodo_summary) <- NULL

  write.csv(lodo_summary,
            file.path(output_dir, paste0("lodo_multigen_summary_", today, ".csv")),
            row.names = FALSE)

  # ── LODO summary plot ─────────────────────────────────────────────────────
  # One panel per held-out dataset; models as coloured dots; null as dashed line
  null_df <- lodo_summary[lodo_summary$model == "null_age_sex", ]
  main_df <- lodo_summary[lodo_summary$model != "null_age_sex", ]
  main_df$model <- factor(main_df$model, levels = c("glmnet_multi", "RF_multi"))

  model_colors_lodo <- c("glmnet_multi" = "#377eb8", "RF_multi" = "#e41a1c")
  model_labels_lodo <- c("glmnet_multi" = "Multi-gene glmnet", "RF_multi" = "Multi-gene RF")

  p_lodo <- ggplot(main_df, aes(x = model, color = model)) +
    # null BalAcc reference line per facet (dashed = primary metric)
    geom_hline(data = null_df,
               aes(yintercept = balacc_opt),
               linetype = "dashed", color = "grey50", linewidth = 0.8,
               inherit.aes = FALSE) +
    # null AUC reference line (dotted = secondary metric)
    geom_hline(data = null_df,
               aes(yintercept = auc),
               linetype = "dotted", color = "grey70", linewidth = 0.5,
               inherit.aes = FALSE) +
    # BalAcc(opt) — primary metric — large filled circle
    geom_point(aes(y = balacc_opt), size = 5, shape = 16) +
    # AUC — secondary — open circle, slightly smaller
    geom_point(aes(y = auc), size = 3.5, shape = 1, stroke = 1.1) +
    # Segment connecting BalAcc to AUC for same model
    geom_segment(aes(xend = model, y = balacc_opt, yend = auc),
                 linewidth = 0.5, alpha = 0.4) +
    # Sens/Spec as small triangles
    geom_point(aes(y = sensitivity_opt), size = 2.5, shape = 24, fill = "white") +
    geom_point(aes(y = specificity_opt), size = 2.5, shape = 25, fill = "white") +
    facet_wrap(~ held_out_dataset, nrow = 1) +
    scale_color_manual(values = model_colors_lodo, labels = model_labels_lodo,
                       name = "Model") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    geom_hline(yintercept = 0.5, color = "grey85", linewidth = 0.4) +
    geom_hline(yintercept = 0.7, linetype = "dotted", color = "orange",  linewidth = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "red",     linewidth = 0.5) +
    annotate("text", x = 0.6, y = 0.03,
             label = paste0("● BalAcc(opt)   ○ AUC   △ Sens   ▽ Spec",
                            "   dashed = null BalAcc   dotted = null AUC"),
             hjust = 0, size = 2.5, color = "grey35") +
    labs(
      title    = "Leave-One-Dataset-Out (LODO) Validation — Multi-gene models",
      subtitle = paste0("Primary metric: BalAcc(opt) ●   Secondary: AUC ○",
                        "  |  Each panel: trained on all other datasets, tested on held-out"),
      x = NULL, y = "Metric value"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 10),
      plot.subtitle   = element_text(size = 8, color = "grey30"),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      legend.position = "bottom",
      strip.text      = element_text(face = "bold")
    )

  ggsave(file.path(output_dir, paste0("lodo_summary_plot_", today, ".pdf")),
         plot   = p_lodo,
         width  = max(6, length(dataset_names) * 2.5),
         height = 5,
         units  = "in",
         device = cairo_pdf)

  # ── LODO per-gene heatmap: mean LODO AUC across held-out datasets ─────────
  top30_lodo <- head(lodo_gene_scores, 30)

  # Primary heatmap: per-gene BalAcc (glmnet) per held-out dataset
  balacc_enet_cols_avail <- paste0("balacc_enet_", dataset_names)
  heatmap_mat <- as.matrix(top30_lodo[, balacc_enet_cols_avail])
  rownames(heatmap_mat) <- top30_lodo$gene_symbol
  colnames(heatmap_mat) <- dataset_names

  # Column annotation: mean AUC per dataset (for context)
  col_annot <- data.frame(
    Mean_AUC = colMeans(as.matrix(top30_lodo[, enet_cols]), na.rm = TRUE),
    row.names = dataset_names
  )

  # Row annotation: sex-chromosome and mean BalAcc
  row_annot <- data.frame(
    Sex_chr      = ifelse(top30_lodo$on_sex_chr, "Yes", "No"),
    Mean_BalAcc  = round(top30_lodo$mean_balacc_lodo, 2),
    row.names    = top30_lodo$gene_symbol
  )

  pheatmap::pheatmap(
    heatmap_mat,
    color           = colorRampPalette(c("#d73027", "white", "#4575b4"))(100),
    breaks          = seq(0, 1, length.out = 101),
    display_numbers = TRUE, number_format = "%.2f", fontsize_number = 7,
    cluster_rows    = TRUE, cluster_cols  = FALSE,
    annotation_row  = row_annot,
    annotation_col  = col_annot,
    annotation_colors = list(
      Sex_chr     = c("Yes" = "purple", "No" = "grey90"),
      Mean_BalAcc = colorRampPalette(c("white", "#4575b4"))(50),
      Mean_AUC    = colorRampPalette(c("white", "tomato"))(50)
    ),
    main     = "LODO per-gene BalAcc(opt) — glmnet — top 30 by mean BalAcc",
    filename = file.path(output_dir,
                         paste0("lodo_per_gene_heatmap_", today, ".pdf")),
    width  = max(5, length(dataset_names) * 1.2 + 3),
    height = max(6, 30 * 0.28 + 2)
  )

  cat("\nLODO validation complete.\n")
  cat("  Genes evaluated:", nrow(lodo_gene_scores), "\n")
  cat("  Datasets used as held-out:", paste(dataset_names, collapse = ", "), "\n")

  invisible(list(
    lodo_summary     = lodo_summary,
    lodo_gene_scores = lodo_gene_scores,
    p_lodo           = p_lodo
  ))
}
