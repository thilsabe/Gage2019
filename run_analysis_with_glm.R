run_analysis <- function(X, y, meta_use, dname, label, ensembl_ids, hgnc_symbols,
                         output_dir, today, alpha_val, B_boot,
                         gene_chr_map = NULL) {
  
  gene_names   <- colnames(X)
  gene_symbols <- hgnc_symbols[match(gene_names, ensembl_ids)]
  n            <- length(y)
  auc_vec      <- numeric(length(gene_names))
  ci_low_vec   <- numeric(length(gene_names))
  ci_high_vec  <- numeric(length(gene_names))
  names(auc_vec) <- names(ci_low_vec) <- names(ci_high_vec) <- gene_names
  
  pred_store     <- list()
  pred_glm_store <- list()   # GLM predictions stored in parallel
  boxplot_list   <- list()
  gene_n         <- 0
  
  # ── Null model (age + sex only) LOOCV — computed once, reused per gene ──
  sex_numeric_all  <- ifelse(toupper(meta_use$Sex) == "M", 1, 0)
  pred_null        <- numeric(n)
  has_sex_var_null <- length(unique(sex_numeric_all)) > 1
  
  for (i in 1:n) {
    train_idx <- setdiff(1:n, i)
    df_null   <- data.frame(
      y   = y[train_idx],
      age = as.numeric(meta_use$Age[train_idx]),
      sex = sex_numeric_all[train_idx]
    )
    df_null_test <- data.frame(
      age = as.numeric(meta_use$Age[i]),
      sex = sex_numeric_all[i]
    )
    formula_null <- if (has_sex_var_null) y ~ age + sex else y ~ age
    fit_null <- tryCatch(
      suppressWarnings(glm(formula_null, data = df_null, family = binomial())),
      error = function(e) NULL
    )
    pred_null[i] <- if (!is.null(fit_null)) {
      tryCatch(predict(fit_null, newdata = df_null_test, type = "response"),
               error = function(e) 0.5)
    } else 0.5
  }
  null_res    <- find_opt_threshold(pred_null, y)
  null_balacc <- null_res$balacc
  cat("  Null model (age+sex) BalAcc:", round(null_balacc * 100, 1), "%\n")
  
  # ── Per-gene LOOCV loop ───────────────────────────────────────────────────
  for (g in gene_names) {
    gene_n   <- gene_n + 1
    pred_all <- numeric(n)
    pred_glm <- numeric(n)
    
    if (var(X[, g]) == 0) {
      auc_vec[g] <- NA; ci_low_vec[g] <- NA; ci_high_vec[g] <- NA
      next
    }
    
    for (i in 1:n) {
      train_idx <- setdiff(1:n, i)
      
      sex_train <- ifelse(toupper(meta_use$Sex[train_idx]) == "M", 1, 0)
      sex_test  <- ifelse(toupper(meta_use$Sex[i])         == "M", 1, 0)
      
      X_train <- cbind(
        gene = X[train_idx, g],
        age  = as.numeric(meta_use$Age[train_idx]),
        sex  = sex_train
      )
      X_test <- matrix(c(
        X[i, g],
        as.numeric(meta_use$Age[i]),
        sex_test
      ), nrow = 1, dimnames = list(NULL, c("gene", "age", "sex")))
      
      y_train <- y[train_idx]
      
      if (var(X_train[, "gene"]) == 0 || length(unique(y_train)) < 2) {
        pred_all[i] <- 0.5
        pred_glm[i] <- 0.5
        next
      }
      
      # ── glmnet prediction ──────────────────────────────────────────────
      pf <- c(gene = 1, age = 0, sex = 0)
      
      fit_cv <- tryCatch(
        cv.glmnet(
          x              = X_train,
          y              = y_train,
          family         = "binomial",
          alpha          = alpha_val,
          penalty.factor = pf,
          nfolds         = min(5, length(y_train)),
          type.measure   = "class",
          standardize    = TRUE
        ),
        error = function(e) NULL
      )
      
      if (is.null(fit_cv)) {
        pred_all[i] <- 0.5
      } else {
        gene_coef_min <- tryCatch(
          as.numeric(coef(fit_cv, s = "lambda.min")["gene", ]),
          error = function(e) 0
        )
        gene_coef_1se <- tryCatch(
          as.numeric(coef(fit_cv, s = "lambda.1se")["gene", ]),
          error = function(e) 0
        )
        s_use <- if (!is.na(gene_coef_min) && gene_coef_min != 0) {
          "lambda.min"
        } else if (!is.na(gene_coef_1se) && gene_coef_1se != 0) {
          "lambda.1se"
        } else {
          NA
        }
        pred_all[i] <- if (is.na(s_use)) 0.5 else tryCatch(
          as.numeric(predict(fit_cv, newx = X_test, s = s_use, type = "response")),
          error = function(e) 0.5
        )
      }
      
      # ── GLM prediction (unpenalized, for comparison) ───────────────────
      has_sex_var <- length(unique(sex_train)) > 1
      df_tr <- data.frame(
        y    = y_train,
        gene = X_train[, "gene"],
        age  = X_train[, "age"],
        sex  = sex_train
      )
      df_te <- data.frame(
        gene = X_test[, "gene"],
        age  = X_test[, "age"],
        sex  = sex_test
      )
      formula_glm <- if (has_sex_var) y ~ gene + age + sex else y ~ gene + age
      fit_glm <- tryCatch(
        suppressWarnings(glm(formula_glm, data = df_tr, family = binomial())),
        error = function(e) NULL
      )
      pred_glm[i] <- if (!is.null(fit_glm)) {
        tryCatch(
          as.numeric(predict(fit_glm, newdata = df_te, type = "response")),
          error = function(e) 0.5
        )
      } else 0.5
    }
    
    # Store completed predictions after inner loop
    pred_store[[g]]     <- pred_all
    pred_glm_store[[g]] <- pred_glm
    
    roc_g <- tryCatch(
      roc(response = y, predictor = pred_all, quiet = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(roc_g)) {
      auc_vec[g] <- NA; ci_low_vec[g] <- NA; ci_high_vec[g] <- NA; next
    }
    
    auc_vec[g]     <- as.numeric(auc(roc_g))
    ci_g           <- bootstrap_auc_ci(y = y, p_hat = pred_all, B = B_boot, conf = 0.95)
    ci_low_vec[g]  <- ci_g$ci_low
    ci_high_vec[g] <- ci_g$ci_high
    
    sym <- if (!is.na(gene_symbols[gene_n])) gene_symbols[gene_n] else g
    
    # ── Optimal threshold (glmnet) ────────────────────────────────────────
    opt_res          <- find_opt_threshold(pred_all, y)
    opt_thresh       <- opt_res$threshold
    sensitivity_opt  <- opt_res$sensitivity
    specificity_opt  <- opt_res$specificity
    balanced_acc_opt <- opt_res$balacc
    thresh_inverted  <- isTRUE(opt_res$inverted)
    
    ppv_opt <- if (!thresh_inverted) {
      sum((pred_all >= opt_thresh) & (y == 1)) / max(1, sum(pred_all >= opt_thresh))
    } else {
      sum((pred_all <  opt_thresh) & (y == 1)) / max(1, sum(pred_all <  opt_thresh))
    }
    npv_opt <- if (!thresh_inverted) {
      sum((pred_all <  opt_thresh) & (y == 0)) / max(1, sum(pred_all <  opt_thresh))
    } else {
      sum((pred_all >= opt_thresh) & (y == 0)) / max(1, sum(pred_all >= opt_thresh))
    }
    predicted_class_opt <- if (!thresh_inverted) {
      ifelse(pred_all >= opt_thresh, "AD", "CTRL")
    } else {
      ifelse(pred_all <  opt_thresh, "AD", "CTRL")
    }
    correct_opt <- predicted_class_opt == ifelse(y == 1, "AD", "CTRL")
    
    # ── Stats at fixed 0.5 threshold ──────────────────────────────────────
    sensitivity_05  <- sum((pred_all >= 0.5) & (y == 1)) / max(1, sum(y == 1))
    specificity_05  <- sum((pred_all <  0.5) & (y == 0)) / max(1, sum(y == 0))
    balanced_acc_05 <- (sensitivity_05 + specificity_05) / 2
    correct_05      <- ifelse(pred_all >= 0.5, "AD", "CTRL") == ifelse(y == 1, "AD", "CTRL")
    
    # ── GLM stats ─────────────────────────────────────────────────────────
    glm_opt_res      <- find_opt_threshold(pred_glm, y)
    glm_balacc_opt   <- glm_opt_res$balacc
    glm_sens_opt     <- glm_opt_res$sensitivity
    glm_spec_opt     <- glm_opt_res$specificity
    
    # ── Incremental gain over null model (descriptive only — not a filter) ──
    balacc_gain <- balanced_acc_opt - null_balacc

    # ── Adjusted expression residuals ─────────────────────────────────────
    df_full <- data.frame(
      x   = as.numeric(X[, g]),
      sex = as.factor(meta_use$Sex),
      age = as.numeric(meta_use$Age)
    )
    has_sex_variance_full <- length(unique(df_full$sex)) > 1
    lm_cov <- tryCatch(
      if (has_sex_variance_full) lm(x ~ sex + age, data = df_full)
      else                       lm(x ~ age,       data = df_full),
      error = function(e) NULL
    )
    adj_expression <- if (!is.null(lm_cov)) residuals(lm_cov) else df_full$x - mean(df_full$x)
    
    # ── Build plot data frame ─────────────────────────────────────────────
    plot_df <- data.frame(
      raw_expression = as.numeric(X[, g]),
      adj_expression = as.numeric(adj_expression),
      predicted_prob = pred_all,
      predicted_glm  = pred_glm,
      actual_diag    = ifelse(y == 1, "AD", "CTRL"),
      group          = as.factor(meta_use[["Diag"]]),
      sex            = as.factor(meta_use$Sex),
      age            = as.numeric(meta_use$Age),
      sample         = rownames(meta_use),
      stringsAsFactors = FALSE
    )
    plot_df$predicted_opt <- predicted_class_opt
    plot_df$predicted_05  <- ifelse(plot_df$predicted_prob >= 0.5, "AD", "CTRL")
    plot_df$correct_opt   <- plot_df$predicted_opt == plot_df$actual_diag
    plot_df$correct_05    <- plot_df$predicted_05  == plot_df$actual_diag
    
    # ── Panel A: Raw expression ───────────────────────────────────────────
    p_raw <- ggplot(plot_df, aes(x = group, y = raw_expression, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(aes(shape = sex, color = age), width = 0.2, size = 2.5, alpha = 0.9) +
      scale_shape_manual(values = c("M" = 16, "F" = 17), na.value = 15) +
      scale_color_gradient(low = "lightblue", high = "darkred", name = "Age") +
      labs(title = "A: Raw expression",
           x = "Diagnosis", y = "Expression",
           fill = "Diagnosis", shape = "Sex") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # ── Panel B: Adjusted expression ─────────────────────────────────────
    p_adj <- ggplot(plot_df, aes(x = group, y = adj_expression, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(aes(shape = sex, color = age), width = 0.2, size = 2.5, alpha = 0.9) +
      scale_shape_manual(values = c("M" = 16, "F" = 17), na.value = 15) +
      scale_color_gradient(low = "lightblue", high = "darkred", name = "Age") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      labs(title = "B: Adjusted expression (sex + age residuals)",
           x = "Diagnosis", y = "Adjusted Expression",
           fill = "Diagnosis", shape = "Sex") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # ── Panel C: Predicted probability per sample ─────────────────────────
    p_pred <- ggplot(plot_df,
                     aes(x     = reorder(sample, predicted_prob),
                         y     = predicted_prob,
                         color = actual_diag,
                         shape = correct_opt)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = opt_thresh, ymax = 1.0,
               fill = "lightgreen", alpha = 0.08) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.0, ymax = opt_thresh,
               fill = "lightyellow", alpha = 0.08) +
      geom_hline(yintercept = opt_thresh, linetype = "dashed",
                 color = "darkgreen", linewidth = 0.8) +
      geom_hline(yintercept = 0.5, linetype = "dotted",
                 color = "grey40", linewidth = 0.5) +
      geom_point(size = 3, alpha = 0.9) +
      # GLM predictions overlaid as grey crosses
      geom_point(aes(y = predicted_glm), shape = 3, size = 2,
                 color = "grey50", alpha = 0.7) +
      annotate("text", x = 1, y = opt_thresh + 0.03,
               label = paste0("Optimal threshold: ", round(opt_thresh, 2)),
               color = "darkgreen", size = 3, hjust = 0) +
      scale_color_manual(values = c("AD" = "tomato", "CTRL" = "steelblue"),
                         name = "Actual Diagnosis") +
      scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                         name   = "Prediction (optimal threshold)",
                         labels = c("TRUE" = "Correct", "FALSE" = "Wrong")) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(title = "C: LOOCV predicted P(AD) per sample  [+ = glm]",
           x     = "Sample (ordered by glmnet predicted probability)",
           y     = "Predicted P(AD)") +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x     = element_text(angle = 90, hjust = 1, size = 7),
        plot.title      = element_text(face = "bold", size = 10),
        legend.position = "bottom"
      )
    
    # ── Panel D: Dual confusion matrix (optimal vs 0.5) ───────────────────
    make_conf_df <- function(predicted, thresh_label) {
      df <- data.frame(
        actual_diag    = ifelse(y == 1, "AD", "CTRL"),
        predicted_diag = predicted,
        stringsAsFactors = FALSE
      ) %>%
        group_by(actual_diag, predicted_diag) %>%
        summarize(n = n(), .groups = "drop") %>%
        mutate(correct   = actual_diag == predicted_diag,
               threshold = thresh_label)
      all_combos <- expand.grid(
        actual_diag    = c("AD", "CTRL"),
        predicted_diag = c("AD", "CTRL"),
        threshold      = thresh_label,
        stringsAsFactors = FALSE
      )
      df <- merge(all_combos, df, all.x = TRUE)
      df$n[is.na(df$n)]             <- 0
      df$correct[is.na(df$correct)] <- FALSE
      df
    }
    
    conf_both <- bind_rows(
      make_conf_df(plot_df$predicted_opt,
                   paste0("glmnet Optimal (", round(opt_thresh, 2), ")\n",
                          "BalAcc=", round(balanced_acc_opt * 100, 1), "%")),
      make_conf_df(plot_df$predicted_05,
                   paste0("glmnet Fixed (0.5)\n",
                          "BalAcc=", round(balanced_acc_05 * 100, 1), "%"))
    )
    
    p_conf <- ggplot(conf_both,
                     aes(x = predicted_diag, y = actual_diag,
                         fill = correct, label = n)) +
      geom_tile(alpha = 0.7) +
      geom_text(size = 7, fontface = "bold") +
      facet_wrap(~ threshold) +
      scale_fill_manual(values = c("TRUE" = "lightgreen", "FALSE" = "lightsalmon"),
                        guide = "none") +
      labs(title = "D: Confusion matrix", x = "Predicted", y = "Actual") +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold", size = 8))
    
    # ── Combine all 4 panels ──────────────────────────────────────────────
    p_combined <- cowplot::plot_grid(
      p_raw, p_adj, p_pred, p_conf,
      ncol = 2, rel_heights = c(1, 1.3)
    )
    
    title_grob <- cowplot::ggdraw() +
      cowplot::draw_label(
        paste0(sym, " (", g, ")",
               "  |  glmnet AUC: ",   round(auc_vec[g], 3),
               " [", round(ci_low_vec[g], 3), ", ", round(ci_high_vec[g], 3), "]",
               "  |  glm BalAcc: ",   round(glm_balacc_opt * 100, 1), "%",
               "  |  glmnet BalAcc(opt): ", round(balanced_acc_opt * 100, 1), "%",
               "  |  BalAcc(0.5): ",  round(balanced_acc_05  * 100, 1), "%",
               "  |  Sens: ",         round(sensitivity_opt  * 100, 1), "%",
               "  |  Spec: ",         round(specificity_opt  * 100, 1), "%",
               "  |  Gain: ",         round(balacc_gain      * 100, 1), "%",
               "  |  ", label),
        fontface = "bold", size = 8, x = 0.01, hjust = 0
      )
    
    boxplot_list[[g]] <- cowplot::plot_grid(
      title_grob, p_combined,
      ncol = 1, rel_heights = c(0.04, 1)
    )
    
    cat("  Gene", gene_n, "/", length(gene_names), ":", sym,
        "| glmnet AUC:", round(auc_vec[g], 3),
        "| glmnet BalAcc:", round(balanced_acc_opt * 100, 1), "%",
        "| glm BalAcc:", round(glm_balacc_opt * 100, 1), "%",
        "| Gain over null:", round(balacc_gain * 100, 1), "%\n")
  }
  
  # ── Save all gene plots to one combined PDF ───────────────────────────────
  combined_pdf_path <- file.path(output_dir,
                                 paste0(label, "_all_gene_plots_", today, ".pdf"))
  message("Saving combined gene plot PDF: ", combined_pdf_path)
  
  balacc_tmp <- sapply(names(boxplot_list), function(g) {
    pred_all <- pred_store[[g]]
    if (is.null(pred_all)) return(NA_real_)
    find_opt_threshold(pred_all, y)$balacc
  })
  sort_df <- data.frame(
    g      = names(boxplot_list),
    balacc = balacc_tmp,
    auc    = auc_vec[names(boxplot_list)],   # retained for tiebreaking only
    stringsAsFactors = FALSE
  )
  sort_df          <- sort_df[order(-sort_df$balacc, -sort_df$auc, na.last = TRUE), ]
  plot_order_valid <- sort_df$g
  
  pdf(combined_pdf_path, width = 14, height = 10)
  for (g in plot_order_valid) print(boxplot_list[[g]])
  dev.off()
  message("Saved ", length(boxplot_list), " gene plots to ", combined_pdf_path)
  
  # ── Build gene summary ────────────────────────────────────────────────────
  gene_summary_rows <- lapply(gene_names, function(g) {
    pred_all <- pred_store[[g]]
    pred_glm <- pred_glm_store[[g]]
    if (is.null(pred_all)) return(NULL)
    opt_check <- find_opt_threshold(pred_all, y)
    if (is.na(opt_check$balacc)) return(NULL)
    
    roc_g <- tryCatch(
      roc(response = y, predictor = pred_all, quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_g)) return(NULL)
    
    # glmnet stats
    opt_res          <- find_opt_threshold(pred_all, y)
    opt_thresh       <- opt_res$threshold
    sensitivity_opt  <- opt_res$sensitivity
    specificity_opt  <- opt_res$specificity
    balanced_acc_opt <- opt_res$balacc
    thresh_inverted  <- isTRUE(opt_res$inverted)
    
    ppv_opt <- if (!thresh_inverted) {
      sum((pred_all >= opt_thresh) & (y == 1)) / max(1, sum(pred_all >= opt_thresh))
    } else {
      sum((pred_all <  opt_thresh) & (y == 1)) / max(1, sum(pred_all <  opt_thresh))
    }
    npv_opt <- if (!thresh_inverted) {
      sum((pred_all <  opt_thresh) & (y == 0)) / max(1, sum(pred_all <  opt_thresh))
    } else {
      sum((pred_all >= opt_thresh) & (y == 0)) / max(1, sum(pred_all >= opt_thresh))
    }
    predicted_class_opt <- if (!thresh_inverted) {
      ifelse(pred_all >= opt_thresh, "AD", "CTRL")
    } else {
      ifelse(pred_all <  opt_thresh, "AD", "CTRL")
    }
    correct_opt <- predicted_class_opt == ifelse(y == 1, "AD", "CTRL")
    
    sensitivity_05  <- sum((pred_all >= 0.5) & (y == 1)) / max(1, sum(y == 1))
    specificity_05  <- sum((pred_all <  0.5) & (y == 0)) / max(1, sum(y == 0))
    balanced_acc_05 <- (sensitivity_05 + specificity_05) / 2
    correct_05      <- ifelse(pred_all >= 0.5, "AD", "CTRL") == ifelse(y == 1, "AD", "CTRL")
    
    # GLM stats
    glm_opt_res    <- if (!is.null(pred_glm)) find_opt_threshold(pred_glm, y) else
                      list(threshold = NA, sensitivity = NA, specificity = NA, balacc = NA)
    roc_glm        <- tryCatch(
      roc(response = y, predictor = pred_glm, quiet = TRUE),
      error = function(e) NULL
    )
    auc_glm        <- if (!is.null(roc_glm)) as.numeric(auc(roc_glm)) else NA_real_
    ci_glm         <- if (!is.null(pred_glm) && !is.na(auc_glm)) {
      bootstrap_auc_ci(y = y, p_hat = pred_glm, B = B_boot, conf = 0.95)
    } else list(ci_low = NA_real_, ci_high = NA_real_)
    
    sym            <- hgnc_symbols[match(g, ensembl_ids)]
    sym            <- if (!is.na(sym)) sym else g
    balacc_gain    <- balanced_acc_opt - null_balacc

    # Chromosome annotation — looked up from gene_chr_map passed into function
    chr_annot  <- if (!is.null(gene_chr_map)) gene_chr_map[g] else NA_character_
    on_sex_chr <- !is.na(chr_annot) & chr_annot %in% c("X", "chrX", "Y", "chrY")

    data.frame(
      ensembl_id        = g,
      gene_symbol       = sym,
      chromosome        = ifelse(is.na(chr_annot), "unknown", chr_annot),
      on_sex_chr        = on_sex_chr,
      # glmnet metrics
      auc               = round(auc_vec[g],          3),
      ci_low            = round(ci_low_vec[g],        3),
      ci_high           = round(ci_high_vec[g],       3),
      opt_threshold     = round(opt_thresh,           3),
      sensitivity_opt   = round(sensitivity_opt,      3),
      specificity_opt   = round(specificity_opt,      3),
      ppv_opt           = round(ppv_opt,              3),
      npv_opt           = round(npv_opt,              3),
      balanced_acc_opt  = round(balanced_acc_opt,     3),
      accuracy_opt      = round(mean(correct_opt),    3),
      n_correct_opt     = sum(correct_opt),
      sensitivity_05    = round(sensitivity_05,       3),
      specificity_05    = round(specificity_05,       3),
      balanced_acc_05   = round(balanced_acc_05,      3),
      accuracy_05       = round(mean(correct_05),     3),
      n_correct_05      = sum(correct_05),
      # GLM metrics
      auc_glm           = round(auc_glm,              3),
      ci_low_glm        = round(ci_glm$ci_low,        3),
      ci_high_glm       = round(ci_glm$ci_high,       3),
      balacc_opt_glm    = round(glm_opt_res$balacc,   3),
      sensitivity_glm   = round(glm_opt_res$sensitivity, 3),
      specificity_glm   = round(glm_opt_res$specificity, 3),
      # Null comparison (descriptive)
      n_total           = length(y),
      n_AD              = sum(y == 1),
      n_CTRL            = sum(y == 0),
      null_balacc       = round(null_balacc,          3),
      balacc_gain       = round(balacc_gain,          3),
      stringsAsFactors  = FALSE
    )
  })
  
  gene_summary_df <- do.call(rbind, Filter(Negate(is.null), gene_summary_rows))
  gene_summary_df <- gene_summary_df[order(-gene_summary_df$balanced_acc_opt), ]
  
  write.csv(gene_summary_df,
            file.path(output_dir, paste0("per_gene_summary_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  # ── Build per_gene_df ─────────────────────────────────────────────────────
  per_gene_df <- gene_summary_df[, c("ensembl_id", "gene_symbol",
                                     "chromosome", "on_sex_chr",
                                     "auc",     "ci_low",     "ci_high",
                                     "auc_glm", "ci_low_glm", "ci_high_glm",
                                     "balanced_acc_opt", "balanced_acc_05", "balacc_opt_glm",
                                     "opt_threshold", "sensitivity_opt", "specificity_opt",
                                     "sensitivity_glm", "specificity_glm",
                                     "null_balacc", "balacc_gain")]
  
  per_gene_df_sorted <- per_gene_df[!is.na(per_gene_df$balanced_acc_opt), ]
  per_gene_df_sorted <- per_gene_df_sorted[order(-per_gene_df_sorted$balanced_acc_opt,
                                                  -per_gene_df_sorted$auc,
                                                  na.last = TRUE), ]
  # Strength bands based on BalAcc(opt) — primary metric
  per_gene_df_sorted$strength <- cut(
    per_gene_df_sorted$balanced_acc_opt,
    breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels = c("Poor (<0.6)", "Fair (0.6-0.7)", "Good (0.7-0.8)",
               "Strong (0.8-0.9)", "Excellent (>0.9)"),
    include.lowest = TRUE
  )
  
  # ── Top biomarker dot plot ────────────────────────────────────────────────
  top20            <- head(per_gene_df_sorted, 20)
  top20_sex_chr    <- top20[!is.na(top20$on_sex_chr) & top20$on_sex_chr, ]

  p_biomarker <- ggplot(top20, aes(y = reorder(gene_symbol, balanced_acc_opt),
                                   color = strength)) +
    # glmnet AUC + CI (solid circle)
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.25, linewidth = 0.7) +
    geom_point(aes(x = auc), size = 4, shape = 16) +
    # glm AUC + CI (dashed, grey triangle)
    geom_errorbarh(aes(xmin = ci_low_glm, xmax = ci_high_glm),
                   height = 0.15, linewidth = 0.5, linetype = "dashed",
                   color = "grey55") +
    geom_point(aes(x = auc_glm), size = 3, shape = 17, color = "grey45") +
    # BalAcc: glmnet optimal (filled diamond), glmnet @0.5 (filled square), glm optimal (open diamond)
    geom_point(aes(x = balanced_acc_opt, fill = strength),
               size = 4, shape = 23, color = "grey30") +
    geom_point(aes(x = balanced_acc_05, fill = strength),
               size = 3, shape = 22, color = "grey30") +
    geom_point(aes(x = balacc_opt_glm),
               size = 3, shape = 5, color = "grey45") +
    scale_fill_manual(
      values = c("Poor (<0.6)"      = "#d9d9d9",
                 "Fair (0.6-0.7)"   = "#fc8d59",
                 "Good (0.7-0.8)"   = "#fee090",
                 "Strong (0.8-0.9)" = "#91bfdb",
                 "Excellent (>0.9)" = "#4575b4"),
      guide = "none"
    ) +
    # Sex-chromosome genes annotated by chromosome label, not penalised
    {if (nrow(top20_sex_chr) > 0)
      geom_label(data = top20_sex_chr,
                 aes(x = 0.415, y = reorder(gene_symbol, balanced_acc_opt),
                     label = chromosome),
                 size = 2.5, color = "purple", fill = "white",
                 label.padding = unit(0.12, "lines"),
                 label.size = 0.3, inherit.aes = FALSE)
    } +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0.7, linetype = "dotted", color = "orange") +
    geom_vline(xintercept = 0.8, linetype = "dotted", color = "red") +
    annotate("text", x = 0.46, y = -Inf,
             label = paste0("● glmnet AUC+CI   ▲ glm AUC+CI   ",
                            "◆ glmnet BalAcc(opt)   ■ glmnet BalAcc(0.5)   ",
                            "◇ glm BalAcc(opt)   purple = sex chromosome gene"),
             hjust = 0, vjust = -0.5, size = 2.8, color = "grey30") +
    scale_color_discrete(name = "AUC strength") +
    scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
    labs(
      title    = paste0("Top AD Biomarkers — ", label),
      subtitle = paste0("Null model (age+sex) BalAcc: ", round(null_balacc * 100, 1),
                        "%  |  balacc_gain = gene increment over demographic baseline",
                        "  |  Solid/colour = glmnet   Grey = glm"),
      x        = "Metric value",
      y        = "Gene (sorted by glmnet balanced accuracy)",
      color    = "AUC strength"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 10),
      plot.subtitle   = element_text(size = 8, color = "grey30"),
      legend.position = "bottom"
    )
  
  ggsave(file.path(output_dir, paste0("top_biomarkers_", label, "_", today, ".pdf")),
         plot = p_biomarker, width = 9, height = 7, units = "in", device = cairo_pdf)
  
  write.csv(head(per_gene_df_sorted, 20),
            file.path(output_dir, paste0("top20_biomarkers_", label, "_", today, ".csv")),
            row.names = FALSE)
  write.csv(per_gene_df,
            file.path(output_dir, paste0("per_gene_LOOCV_auc_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  # ── Elastic net ───────────────────────────────────────────────────────────
  pred_multi <- run_loocv_elastic_net(X, y, meta_use = meta_use, alpha_val = alpha_val)
  roc_multi  <- roc(response = y, predictor = pred_multi, quiet = TRUE)
  auc_multi  <- as.numeric(auc(roc_multi))
  boot_res   <- bootstrap_auc_ci(y, pred_multi, B = B_boot, conf = 0.95)
  
  boot_df <- data.frame(iteration = 1:B_boot, auc = boot_res$auc_boot)
  write.csv(boot_df,
            file.path(output_dir, paste0("multigene_bootstrap_auc_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  # ── Elastic net stats ─────────────────────────────────────────────────────
  enet_opt_res    <- find_opt_threshold(pred_multi, y)
  enet_opt_thresh <- enet_opt_res$threshold
  enet_sens_opt   <- enet_opt_res$sensitivity
  enet_spec_opt   <- enet_opt_res$specificity
  enet_balacc_opt <- enet_opt_res$balacc
  enet_sens_05    <- sum((pred_multi >= 0.5) & (y == 1)) / max(1, sum(y == 1))
  enet_spec_05    <- sum((pred_multi <  0.5) & (y == 0)) / max(1, sum(y == 0))
  enet_balacc_05  <- (enet_sens_05 + enet_spec_05) / 2
  enet_boot_res   <- bootstrap_auc_ci(y = y, p_hat = pred_multi, B = B_boot, conf = 0.95)
  enet_balacc_strength <- cut(enet_balacc_opt, breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
                              labels = c("Poor", "Fair", "Good", "Strong", "Excellent"),
                              include.lowest = TRUE)
  
  sex_numeric     <- ifelse(toupper(meta_use$Sex) == "M", 1, 0)
  cov_mat         <- cbind(Age = as.numeric(meta_use$Age), Sex = sex_numeric)
  X_full          <- cbind(X, cov_mat)
  penalty_factors <- c(rep(1, ncol(X)), rep(0, ncol(cov_mat)))
  
  cvfit_enet <- cv.glmnet(
    x = X_full, y = y, family = "binomial", alpha = alpha_val,
    nfolds = min(5, length(y)), type.measure = "class",
    penalty.factor = penalty_factors
  )
  cvfit_lasso <- cv.glmnet(
    x = X_full, y = y, family = "binomial", alpha = 1,
    nfolds = min(5, length(y)), type.measure = "class",
    penalty.factor = penalty_factors
  )
  
  cat("\n", label, "- Elastic net (alpha=", alpha_val, "):\n")
  cat("  lambda.min:", sum(coef(cvfit_enet,  s = "lambda.min")[-1] != 0), "genes\n")
  cat("  lambda.1se:", sum(coef(cvfit_enet,  s = "lambda.1se")[-1] != 0), "genes\n")
  cat("Lasso:\n")
  cat("  lambda.min:", sum(coef(cvfit_lasso, s = "lambda.min")[-1] != 0), "genes\n")
  cat("  lambda.1se:", sum(coef(cvfit_lasso, s = "lambda.1se")[-1] != 0), "genes\n")
  
  pdf(file.path(output_dir, paste0("glmnet_CV_curves_", label, "_", today, ".pdf")),
      width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot(cvfit_enet,  main = paste("Elastic net -", label))
  plot(cvfit_lasso, main = paste("Lasso -", label))
  dev.off()
  
  best_fit  <- if (sum(coef(cvfit_enet, s = "lambda.min")[-1] != 0) > 0) cvfit_enet else cvfit_lasso
  coef_full <- coef(best_fit, s = "lambda.min")
  coef_df   <- data.frame(
    feature = rownames(coef_full)[-1],
    coef    = as.numeric(coef_full)[-1]
  )
  coef_df <- coef_df[coef_df$coef != 0, ]
  coef_df <- coef_df[!coef_df$feature %in% c("Age", "Sex"), ]
  coef_df$gene_symbol <- stable_genes$hgnc_symbol[match(coef_df$feature, stable_genes$ensembl_id)]
  coef_df$gene_symbol[is.na(coef_df$gene_symbol)] <- coef_df$feature[is.na(coef_df$gene_symbol)]
  coef_df <- coef_df[order(-abs(coef_df$coef)), ]
  
  cat("Final selected genes:", nrow(coef_df), "\n")
  cat("Selected:", paste(coef_df$gene_symbol, collapse = ", "), "\n")
  write.csv(coef_df,
            file.path(output_dir, paste0("elastic_net_coefficients_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  # ── Collinearity plots ────────────────────────────────────────────────────
  top20_genes <- head(per_gene_df_sorted$ensembl_id[
    per_gene_df_sorted$ensembl_id %in% colnames(X)], 20)
  top20_sym <- stable_genes$hgnc_symbol[match(top20_genes, stable_genes$ensembl_id)]
  top20_sym[is.na(top20_sym)] <- top20_genes[is.na(top20_sym)]
  
  if (length(top20_genes) > 1) {
    X_top20 <- X[, top20_genes, drop = FALSE]
    colnames(X_top20) <- top20_sym
    cor_mat <- cor(X_top20, use = "pairwise.complete.obs")
    pheatmap::pheatmap(
      cor_mat,
      color           = colorRampPalette(c("steelblue", "white", "tomato"))(100),
      breaks          = seq(-1, 1, length.out = 101),
      display_numbers = TRUE, number_format = "%.2f", fontsize_number = 7,
      main            = paste0("Top gene correlation - ", label),
      filename        = file.path(output_dir,
                                  paste0("collinearity_clustered_", label, "_", today, ".pdf")),
      width  = max(6, length(top20_genes) * 0.4),
      height = max(5, length(top20_genes) * 0.4)
    )
    if (nrow(coef_df) > 1) {
      sel_genes <- coef_df$feature[coef_df$feature %in% colnames(X)]
      sel_sym   <- coef_df$gene_symbol[coef_df$feature %in% colnames(X)]
      X_sel     <- X[, sel_genes, drop = FALSE]
      colnames(X_sel) <- sel_sym
      cor_sel   <- cor(X_sel, use = "pairwise.complete.obs")
      pheatmap::pheatmap(
        cor_sel,
        color           = colorRampPalette(c("steelblue", "white", "tomato"))(100),
        breaks          = seq(-1, 1, length.out = 101),
        display_numbers = TRUE, number_format = "%.2f", fontsize_number = 8,
        main            = paste0("Selected gene correlation - ", label),
        filename        = file.path(output_dir,
                                    paste0("selected_gene_correlation_", label, "_", today, ".pdf")),
        width  = max(4, nrow(coef_df) * 0.5),
        height = max(4, nrow(coef_df) * 0.5)
      )
    }
  }
  
  # ── Summary plots ─────────────────────────────────────────────────────────
  p_gene_bar <- ggplot(per_gene_df_sorted,
                       aes(x = reorder(gene_symbol, balanced_acc_opt), y = balanced_acc_opt)) +
    geom_col(fill = "tomato", alpha = 0.8) +
    geom_point(aes(y = auc),          color = "steelblue", size = 2) +
    geom_point(aes(y = balacc_opt_glm), color = "grey45",  size = 2, shape = 17) +
    geom_hline(yintercept = 0.7, linetype = "dotted", color = "orange") +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "red") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(title    = paste("Per-gene balanced accuracy (bars) + AUC (dots) -", label),
         subtitle = "Bars=glmnet BalAcc(opt)  Blue dot=glmnet AUC  Grey triangle=glm BalAcc(opt)",
         x = "Gene", y = "Metric value")
  
  ggsave(file.path(output_dir, paste0("per_gene_LOOCV_auc_bar_", label, "_", today, ".png")),
         p_gene_bar, width = 6, height = 8, dpi = 300)
  
  balacc_strength <- cut(enet_balacc_opt, breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
                         labels = c("Poor", "Fair", "Good", "Strong", "Excellent"),
                         include.lowest = TRUE)

  p_boot <- ggplot(boot_df, aes(x = auc)) +
    geom_histogram(color = "black", fill = "steelblue", bins = 30) +
    geom_vline(xintercept = auc_multi,        color = "red",     linewidth = 1.1) +
    geom_vline(xintercept = boot_res$ci_low,  color = "red",     linetype = "dashed") +
    geom_vline(xintercept = boot_res$ci_high, color = "red",     linetype = "dashed") +
    geom_vline(xintercept = enet_balacc_opt,  color = "tomato",  linetype = "dashed", linewidth = 1.0) +
    geom_vline(xintercept = enet_balacc_05,   color = "salmon",  linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = 0.7,              color = "orange",  linetype = "dotted") +
    geom_vline(xintercept = 0.8,              color = "darkred", linetype = "dotted") +
    annotate("text", x = enet_balacc_opt, y = Inf,
             label = paste0("BalAcc(opt)=", round(enet_balacc_opt, 3), " (", balacc_strength, ")"),
             vjust = 2, hjust = -0.1, color = "tomato", size = 3.5) +
    annotate("text", x = auc_multi, y = Inf,
             label = paste0("AUC=", round(auc_multi, 3)),
             vjust = 4, hjust = -0.1, color = "red", size = 3.2) +
    annotate("text", x = enet_balacc_05, y = Inf,
             label = paste0("BalAcc(@0.5)=", round(enet_balacc_05, 3)),
             vjust = 6, hjust = -0.05, color = "salmon", size = 3.2) +
    theme_bw(base_size = 10) +
    labs(title = paste0("Multigene LOOCV — ", label,
                        "\nBalAcc(opt): ", round(enet_balacc_opt * 100, 1), "%",
                        " (", balacc_strength, ")",
                        "  |  Sens: ",     round(enet_sens_opt  * 100, 1), "%",
                        "  |  Spec: ",     round(enet_spec_opt  * 100, 1), "%",
                        "  |  BalAcc(0.5): ", round(enet_balacc_05 * 100, 1), "%",
                        "\nAUC: ", round(boot_res$mean, 3),
                        " CI [", round(boot_res$ci_low, 3), ", ", round(boot_res$ci_high, 3), "]"),
         x = "Bootstrap AUC", y = "Count")
  ggsave(file.path(output_dir, paste0("multigene_bootstrap_auc_hist_", label, "_", today, ".png")),
         p_boot, width = 6, height = 4, dpi = 300)
  
  # ── Structured comparison: glmnet + glm per-gene vs elastic net ───────────
  needed_cols  <- c("gene_symbol", "auc", "ci_low", "ci_high",
                    "balanced_acc_opt", "balanced_acc_05",
                    "sensitivity_opt", "specificity_opt")
  missing_cols <- setdiff(needed_cols, names(gene_summary_df))
  if (length(missing_cols) > 0)
    stop("gene_summary_df missing columns: ", paste(missing_cols, collapse = ", "))
  
  # Top 10 glmnet per-gene
  top10_glmnet <- head(gene_summary_df[order(-gene_summary_df$balanced_acc_opt), ], 10)[, needed_cols]
  top10_glmnet$model_type <- "Single gene (glmnet)"
  
  # Top 10 GLM per-gene (sorted by glm balacc)
  top10_glm <- head(gene_summary_df[order(-gene_summary_df$balacc_opt_glm), ], 10)
  top10_glm_compare <- data.frame(
    gene_symbol      = top10_glm$gene_symbol,
    auc              = top10_glm$auc_glm,
    ci_low           = top10_glm$ci_low_glm,
    ci_high          = top10_glm$ci_high_glm,
    balanced_acc_opt = top10_glm$balacc_opt_glm,
    balanced_acc_05  = NA_real_,
    sensitivity_opt  = top10_glm$sensitivity_glm,
    specificity_opt  = top10_glm$specificity_glm,
    model_type       = "Single gene (glm)",
    stringsAsFactors = FALSE
  )
  
  # Elastic net row
  enet_row <- data.frame(
    gene_symbol      = paste0("Elastic net (", nrow(coef_df), " genes)"),
    auc              = round(auc_multi,             3),
    ci_low           = round(enet_boot_res$ci_low,  3),
    ci_high          = round(enet_boot_res$ci_high, 3),
    balanced_acc_opt = round(enet_balacc_opt,       3),
    balanced_acc_05  = round(enet_balacc_05,        3),
    sensitivity_opt  = round(enet_sens_opt,         3),
    specificity_opt  = round(enet_spec_opt,         3),
    model_type       = "Elastic net (LOOCV)",
    stringsAsFactors = FALSE
  )
  
  comparison_df         <- rbind(top10_glmnet, top10_glm_compare, enet_row)
  comparison_df$is_enet <- comparison_df$model_type == "Elastic net (LOOCV)"
  
  write.csv(comparison_df,
            file.path(output_dir, paste0("comparison_singlegene_vs_enet_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  p_compare <- ggplot(comparison_df,
                      aes(y = reorder(gene_symbol, balanced_acc_opt),
                          color = model_type)) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.3, linewidth = 0.7) +
    geom_point(aes(x = auc), size = 4, shape = 16) +
    geom_point(aes(x = balanced_acc_opt, fill = model_type),
               size = 4, shape = 23, color = "grey30") +
    geom_point(aes(x = balanced_acc_05, fill = model_type),
               size = 3, shape = 22, color = "grey30") +
    scale_color_manual(
      name   = "Model type",
      values = c("Single gene (glmnet)" = "steelblue",
                 "Single gene (glm)"    = "grey45",
                 "Elastic net (LOOCV)"  = "tomato")
    ) +
    scale_fill_manual(
      name   = "Model type",
      values = c("Single gene (glmnet)" = "steelblue",
                 "Single gene (glm)"    = "grey45",
                 "Elastic net (LOOCV)"  = "tomato"),
      guide  = "none"
    ) +
    annotate("text", x = 0.41, y = -Inf,
             label = "● AUC+CI   ◆ BalAcc(optimal)   ■ BalAcc(0.5)",
             hjust = 0, vjust = -0.5, size = 3, color = "grey30") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0.7, linetype = "dotted", color = "orange") +
    geom_vline(xintercept = 0.8, linetype = "dotted", color = "red") +
    scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
    labs(
      title    = paste0("Single-gene (glmnet + glm) vs elastic net — ", label),
      subtitle = paste0(
        "Enet: AUC=",     round(auc_multi, 3),
        " [",             round(enet_boot_res$ci_low, 3), ", ",
                          round(enet_boot_res$ci_high, 3), "]",
        "  BalAcc(opt)=", round(enet_balacc_opt * 100, 1), "%",
        "  Sens=",        round(enet_sens_opt   * 100, 1), "%",
        "  Spec=",        round(enet_spec_opt   * 100, 1), "%"
      ),
      x = "Metric value", y = "Model"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 10),
      plot.subtitle   = element_text(size = 8, color = "grey30"),
      legend.position = "bottom"
    )
  
  ggsave(file.path(output_dir, paste0("comparison_singlegene_vs_enet_", label, "_", today, ".pdf")),
         plot   = p_compare,
         width  = 9,
         height = max(5, nrow(comparison_df) * 0.35),
         units  = "in",
         device = cairo_pdf)
  
  # ── Bootstrap CIs for elastic net coefficients ────────────────────────────
  set.seed(42)
  coef_boot_list <- lapply(1:B_boot, function(b) {
    idx_b    <- sample(1:length(y), replace = TRUE)
    y_b      <- y[idx_b]
    X_full_b <- X_full[idx_b, , drop = FALSE]
    if (length(unique(y_b)) < 2) return(NULL)
    fit_b <- tryCatch(
      cv.glmnet(x = X_full_b, y = y_b, family = "binomial",
                alpha = alpha_val, nfolds = min(5, length(y_b)),
                type.measure = "class", penalty.factor = penalty_factors,
                standardize = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit_b)) return(NULL)
    cf <- coef(fit_b, s = "lambda.min")
    data.frame(feature = rownames(cf)[-1], coef = as.numeric(cf)[-1],
               stringsAsFactors = FALSE)
  })
  
  coef_boot_df <- do.call(rbind, Filter(Negate(is.null), coef_boot_list))
  
  coef_boot_summary <- coef_boot_df %>%
    group_by(feature) %>%
    summarise(
      coef_mean      = mean(coef,            na.rm = TRUE),
      coef_ci_low    = quantile(coef, 0.025, na.rm = TRUE),
      coef_ci_high   = quantile(coef, 0.975, na.rm = TRUE),
      selection_freq = mean(coef != 0,       na.rm = TRUE),
      .groups = "drop"
    )
  
  coef_df <- merge(coef_df, coef_boot_summary, by = "feature", all.x = TRUE)
  coef_df <- coef_df[order(-abs(coef_df$coef)), ]
  
  write.csv(coef_df,
            file.path(output_dir, paste0("elastic_net_coefficients_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  if (nrow(coef_df) > 0) {
    coef_df$direction <- ifelse(coef_df$coef > 0, "Higher in AD", "Lower in AD")
    
    p_coef <- ggplot(coef_df,
                     aes(x = coef, y = reorder(gene_symbol, abs(coef)),
                         fill = direction)) +
      geom_col(alpha = 0.6, width = 0.6) +
      geom_errorbarh(aes(xmin = coef_ci_low, xmax = coef_ci_high),
                     height = 0.3, linewidth = 0.7, color = "grey30") +
      geom_point(aes(x = coef, color = direction), size = 2.5) +
      geom_text(aes(x    = ifelse(coef >= 0, coef_ci_high, coef_ci_low),
                    label = paste0(round(selection_freq * 100), "%")),
                hjust = ifelse(coef_df$coef >= 0, -0.15, 1.15),
                size = 3, color = "grey30") +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
      scale_fill_manual(values  = c("Higher in AD" = "tomato", "Lower in AD" = "steelblue")) +
      scale_color_manual(values = c("Higher in AD" = "tomato", "Lower in AD" = "steelblue"),
                         guide  = "none") +
      labs(
        title    = paste0("Elastic Net Selected Genes — ", label),
        subtitle = paste0(
          "AUC: ",          round(auc_multi, 3),
          " [",             round(boot_res$ci_low, 3), ", ", round(boot_res$ci_high, 3), "]",
          "  BalAcc(opt)=", round(enet_balacc_opt * 100, 1), "%",
          "  Sens=",        round(enet_sens_opt   * 100, 1), "%",
          "  Spec=",        round(enet_spec_opt   * 100, 1), "%\n",
          "Bars=bootstrap 95% CI | % label=selection frequency across ", B_boot, " bootstraps"
        ),
        x = "Coefficient (full-data fit)", y = "Gene", fill = "Direction"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title      = element_text(face = "bold", size = 10),
        plot.subtitle   = element_text(size = 8, color = "grey30"),
        legend.position = "bottom"
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))
    
    ggsave(file.path(output_dir, paste0("elastic_net_coef_", label, "_", today, ".pdf")),
           plot   = p_coef,
           width  = 7,
           height = max(4, nrow(coef_df) * 0.45),
           units  = "in",
           device = cairo_pdf)
  }
  
  cat("Dataset", label, "done.\n")
  cat("  Multigene LOOCV AUC:", round(auc_multi, 3), "\n")
  cat("  Bootstrap mean:", round(boot_res$mean, 3),
      "CI [", round(boot_res$ci_low, 3), ",", round(boot_res$ci_high, 3), "]\n")
  cat("  Elastic net selected", nrow(coef_df), "genes\n")
  
  list(
    per_gene_df       = per_gene_df,
    gene_summary_df   = gene_summary_df,
    auc_multi         = auc_multi,
    boot_mean         = boot_res$mean,
    boot_ci_low       = boot_res$ci_low,
    boot_ci_high      = boot_res$ci_high,
    elastic_net_genes = coef_df
  )
}
