# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: tune_glmnet_alpha
# Selects best alpha from a grid via inner 5-fold CV on the training set.
# Returns list(alpha, lambda, balacc) — best balanced accuracy on inner CV.
# ═══════════════════════════════════════════════════════════════════════════════
tune_glmnet_alpha <- function(X_train, y_train, alpha_grid = c(0, 0.25, 0.5, 0.75, 1),
                              penalty.factor = NULL, nfolds_inner = 5) {
  if (is.null(penalty.factor))
    penalty.factor <- rep(1, ncol(X_train))
  
  nfolds_inner <- min(nfolds_inner, length(y_train))
  
  best_alpha  <- alpha_grid[1]
  best_lambda <- NULL
  best_balacc <- -Inf
  
  for (a in alpha_grid) {
    fit <- tryCatch(
      cv.glmnet(
        x              = X_train,
        y              = y_train,
        family         = "binomial",
        alpha          = a,
        penalty.factor = penalty.factor,
        nfolds         = nfolds_inner,
        type.measure   = "class",
        standardize    = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    # Evaluate balanced accuracy at lambda.min on inner CV predictions
    # cv.glmnet stores cvm = mean CV error; approximate balacc from 1 - cvm
    # (class error → accuracy → use as proxy; exact balacc needs fold predictions)
    balacc_proxy <- 1 - min(fit$cvm, na.rm = TRUE)
    if (balacc_proxy > best_balacc) {
      best_balacc <- balacc_proxy
      best_alpha  <- a
      best_lambda <- fit$lambda.min
    }
  }
  list(alpha = best_alpha, lambda = best_lambda, balacc = best_balacc)
}

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: tune_rf
# Selects best min.node.size via inner 5-fold CV on the training set.
# Returns list(min.node.size, balacc).
# ═══════════════════════════════════════════════════════════════════════════════
tune_rf <- function(df_train, y_train,
                    node_grid   = c(1, 3, 5, 10),
                    nfolds_inner = 5,
                    num.trees    = 500) {
  nfolds_inner <- min(nfolds_inner, nrow(df_train))
  fold_ids     <- sample(rep(1:nfolds_inner, length.out = nrow(df_train)))
  
  best_node   <- node_grid[1]
  best_balacc <- -Inf
  
  for (node in node_grid) {
    fold_preds <- numeric(nrow(df_train))
    
    for (f in 1:nfolds_inner) {
      tr <- df_train[fold_ids != f, , drop = FALSE]
      te <- df_train[fold_ids == f, , drop = FALSE]
      y_tr <- y_train[fold_ids != f]
      y_te <- y_train[fold_ids == f]
      
      if (length(unique(y_tr)) < 2) { fold_preds[fold_ids == f] <- 0.5; next }
      
      tr$y <- factor(y_tr, levels = c(0, 1))
      fit  <- tryCatch(
        ranger::ranger(
          y             ~ .,
          data          = tr,
          num.trees     = num.trees,
          mtry          = max(1, floor(sqrt(ncol(tr) - 1))),
          min.node.size = node,
          probability   = TRUE,
          num.threads   = 1
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) { fold_preds[fold_ids == f] <- 0.5; next }
      
      preds <- tryCatch(
        predict(fit, data = te)$predictions[, "1"],
        error = function(e) rep(0.5, nrow(te))
      )
      fold_preds[fold_ids == f] <- preds
    }
    
    res     <- find_opt_threshold(fold_preds, y_train)
    balacc  <- res$balacc
    if (balacc > best_balacc) {
      best_balacc <- balacc
      best_node   <- node
    }
  }
  list(min.node.size = best_node, balacc = best_balacc)
}

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN: run_analysis_rf
# Mirrors run_analysis but uses Random Forest (ranger) as the primary model.
# Includes:
#   - Per-gene LOOCV with RF (tuned min.node.size per fold)
#   - Parallel glmnet (tuned alpha+lambda) and glm per fold for comparison
#   - Multi-gene RF with permutation importance
#   - Null model comparison and sex confounding flag
#   - Comparison plot: glm vs glmnet vs RF per gene + multi-gene RF vs elastic net
# ═══════════════════════════════════════════════════════════════════════════════
options(ranger.num.threads = 10)
run_analysis_rf <- function(X, y, meta_use, dname, label, ensembl_ids, hgnc_symbols,
                            output_dir, today, alpha_val, B_boot, B_perm = 1000, 
                            num.trees    = 500,
                            node_grid    = c(1, 3, 5, 10),
                            alpha_grid   = c(0, 0.25, 0.5, 0.75, 1.0),
                            nfolds_inner = 5,
                            gene_chr_map = NULL) {
  
  library(ranger)
  
  gene_names   <- colnames(X)
  gene_symbols <- hgnc_symbols[match(gene_names, ensembl_ids)]
  n            <- length(y)
  
  # Per-gene metric stores
  auc_vec_rf      <- numeric(length(gene_names));  names(auc_vec_rf)      <- gene_names
  auc_vec_enet    <- numeric(length(gene_names));  names(auc_vec_enet)    <- gene_names
  auc_vec_glm     <- numeric(length(gene_names));  names(auc_vec_glm)     <- gene_names
  ci_low_rf       <- numeric(length(gene_names));  names(ci_low_rf)       <- gene_names
  ci_high_rf      <- numeric(length(gene_names));  names(ci_high_rf)      <- gene_names
  ci_low_enet     <- numeric(length(gene_names));  names(ci_low_enet)     <- gene_names
  ci_high_enet    <- numeric(length(gene_names));  names(ci_high_enet)    <- gene_names
  ci_low_glm      <- numeric(length(gene_names));  names(ci_low_glm)      <- gene_names
  ci_high_glm     <- numeric(length(gene_names));  names(ci_high_glm)     <- gene_names
  best_alpha_vec  <- numeric(length(gene_names));  names(best_alpha_vec)  <- gene_names
  best_node_vec   <- numeric(length(gene_names));  names(best_node_vec)   <- gene_names
  
  pred_store_rf   <- list()
  pred_store_enet <- list()
  pred_store_glm  <- list()
  boxplot_list    <- list()
  gene_n          <- 0
  
  # ── Null model LOOCV (age + sex only) ─────────────────────────────────────
  sex_numeric_all  <- ifelse(toupper(meta_use$Sex) == "M", 1, 0)
  pred_null        <- numeric(n)
  has_sex_var_null <- length(unique(sex_numeric_all)) > 1
  
  for (i in 1:n) {
    train_idx <- setdiff(1:n, i)
    df_null   <- data.frame(y   = y[train_idx],
                            age = as.numeric(meta_use$Age[train_idx]),
                            sex = sex_numeric_all[train_idx])
    df_null_test <- data.frame(age = as.numeric(meta_use$Age[i]),
                               sex = sex_numeric_all[i])
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
    gene_n      <- gene_n + 1
    pred_rf     <- numeric(n)
    pred_enet   <- numeric(n)
    pred_glm_g  <- numeric(n)
    
    if (var(X[, g]) == 0) {
      auc_vec_rf[g] <- NA; auc_vec_enet[g] <- NA; auc_vec_glm[g] <- NA; next
    }
    
    fold_best_alpha <- numeric(n)
    fold_best_node  <- numeric(n)
    
    for (i in 1:n) {
      train_idx <- setdiff(1:n, i)
      sex_train <- ifelse(toupper(meta_use$Sex[train_idx]) == "M", 1, 0)
      sex_test  <- ifelse(toupper(meta_use$Sex[i])         == "M", 1, 0)
      age_train <- as.numeric(meta_use$Age[train_idx])
      age_test  <- as.numeric(meta_use$Age[i])
      gene_train <- X[train_idx, g]
      gene_test  <- X[i, g]
      y_train    <- y[train_idx]
      
      if (var(gene_train) == 0 || length(unique(y_train)) < 2) {
        pred_rf[i] <- 0.5; pred_enet[i] <- 0.5; pred_glm_g[i] <- 0.5; next
      }
      
      # Shared design objects
      X_train_mat <- cbind(gene = gene_train, age = age_train, sex = sex_train)
      X_test_mat  <- matrix(c(gene_test, age_test, sex_test), nrow = 1,
                            dimnames = list(NULL, c("gene", "age", "sex")))
      df_train    <- data.frame(gene = gene_train, age = age_train,
                                sex = sex_train, stringsAsFactors = FALSE)
      df_test     <- data.frame(gene = gene_test,  age = age_test,
                                sex = sex_test,  stringsAsFactors = FALSE)
      has_sex_var <- length(unique(sex_train)) > 1
      
      # ── Tune + fit glmnet ────────────────────────────────────────────────
      pf         <- c(gene = 1, age = 0, sex = 0)
      tune_enet  <- tune_glmnet_alpha(X_train_mat, y_train,
                                      alpha_grid     = alpha_grid,
                                      penalty.factor = pf,
                                      nfolds_inner   = min(nfolds_inner, length(y_train)))
      fold_best_alpha[i] <- tune_enet$alpha
      
      fit_enet <- tryCatch(
        cv.glmnet(x = X_train_mat, y = y_train, family = "binomial",
                  alpha = tune_enet$alpha, penalty.factor = pf,
                  nfolds = min(nfolds_inner, length(y_train)),
                  type.measure = "class", standardize = TRUE),
        error = function(e) NULL
      )
      if (!is.null(fit_enet)) {
        gc_min <- tryCatch(as.numeric(coef(fit_enet, s = "lambda.min")["gene", ]), error = function(e) 0)
        gc_1se <- tryCatch(as.numeric(coef(fit_enet, s = "lambda.1se")["gene", ]), error = function(e) 0)
        s_use  <- if (!is.na(gc_min) && gc_min != 0) "lambda.min" else
          if (!is.na(gc_1se) && gc_1se != 0) "lambda.1se" else NA
        pred_enet[i] <- if (is.na(s_use)) 0.5 else
          tryCatch(as.numeric(predict(fit_enet, newx = X_test_mat,
                                      s = s_use, type = "response")),
                   error = function(e) 0.5)
      } else pred_enet[i] <- 0.5
      
      # ── Fit glm ──────────────────────────────────────────────────────────
      df_tr_glm   <- data.frame(y = y_train, gene = gene_train,
                                age = age_train, sex = sex_train)
      formula_glm <- if (has_sex_var) y ~ gene + age + sex else y ~ gene + age
      fit_glm     <- tryCatch(
        suppressWarnings(glm(formula_glm, data = df_tr_glm, family = binomial())),
        error = function(e) NULL
      )
      pred_glm_g[i] <- if (!is.null(fit_glm)) {
        tryCatch(as.numeric(predict(fit_glm, newdata = df_test, type = "response")),
                 error = function(e) 0.5)
      } else 0.5
      
      # ── Tune + fit RF ────────────────────────────────────────────────────
      tune_res       <- tune_rf(df_train, y_train,
                                node_grid    = node_grid,
                                nfolds_inner = min(nfolds_inner, length(y_train)),
                                num.trees    = num.trees)
      fold_best_node[i] <- tune_res$min.node.size
      
      df_train_rf        <- df_train
      df_train_rf$y      <- factor(y_train, levels = c(0, 1))
      
      fit_rf <- tryCatch(
        ranger::ranger(
          y             ~ .,
          data          = df_train_rf,
          num.trees     = num.trees,
          mtry          = max(1, floor(sqrt(ncol(df_train_rf) - 1))),
          min.node.size = tune_res$min.node.size,
          probability   = TRUE,
          num.threads   = 1
        ),
        error = function(e) NULL
      )
      pred_rf[i] <- if (!is.null(fit_rf)) {
        tryCatch(predict(fit_rf, data = df_test)$predictions[, "1"],
                 error = function(e) 0.5)
      } else 0.5
    }
    
    # Store
    pred_store_rf[[g]]   <- pred_rf
    pred_store_enet[[g]] <- pred_enet
    pred_store_glm[[g]]  <- pred_glm_g
    best_alpha_vec[g]    <- round(mean(fold_best_alpha, na.rm = TRUE), 3)
    best_node_vec[g]     <- round(mean(fold_best_node,  na.rm = TRUE), 1)
    
    # ── AUC + CI for all three models ────────────────────────────────────
    compute_auc_ci <- function(pred, y, B) {
      roc_obj <- tryCatch(roc(response = y, predictor = pred, quiet = TRUE),
                          error = function(e) NULL)
      if (is.null(roc_obj)) return(list(auc = NA, ci_low = NA, ci_high = NA))
      ci <- bootstrap_auc_ci(y = y, p_hat = pred, B = B, conf = 0.95)
      list(auc = as.numeric(auc(roc_obj)), ci_low = ci$ci_low, ci_high = ci$ci_high)
    }
    
    rf_ci   <- compute_auc_ci(pred_rf,    y, B_boot)
    enet_ci <- compute_auc_ci(pred_enet,  y, B_boot)
    glm_ci  <- compute_auc_ci(pred_glm_g, y, B_boot)
    
    auc_vec_rf[g]   <- rf_ci$auc;    ci_low_rf[g]   <- rf_ci$ci_low;   ci_high_rf[g]   <- rf_ci$ci_high
    auc_vec_enet[g] <- enet_ci$auc;  ci_low_enet[g] <- enet_ci$ci_low; ci_high_enet[g] <- enet_ci$ci_high
    auc_vec_glm[g]  <- glm_ci$auc;   ci_low_glm[g]  <- glm_ci$ci_low;  ci_high_glm[g]  <- glm_ci$ci_high
    
    sym <- if (!is.na(gene_symbols[gene_n])) gene_symbols[gene_n] else g
    
    # ── Threshold stats for all three models ──────────────────────────────
    rf_opt   <- find_opt_threshold(pred_rf,    y)
    enet_opt <- find_opt_threshold(pred_enet,  y)
    glm_opt  <- find_opt_threshold(pred_glm_g, y)
    
    balacc_gain    <- rf_opt$balacc - null_balacc
    
    cat("  Gene", gene_n, "/", length(gene_names), ":", sym,
        "| RF AUC:", round(rf_ci$auc, 3),
        "| RF BalAcc:", round(rf_opt$balacc * 100, 1), "%",
        "| glmnet BalAcc:", round(enet_opt$balacc * 100, 1), "%",
        "| glm BalAcc:", round(glm_opt$balacc * 100, 1), "%",
        "| Gain:", round(balacc_gain * 100, 1), "%",
        "| best_alpha:", best_alpha_vec[g],
        "| best_node:", best_node_vec[g], "\n")
    
    # ── Per-gene 4-panel plot (RF as primary) ─────────────────────────────
    rf_thresh   <- rf_opt$threshold
    rf_inverted <- isTRUE(rf_opt$inverted)
    rf_class    <- if (!rf_inverted) ifelse(pred_rf >= rf_thresh, "AD", "CTRL") else
      ifelse(pred_rf <  rf_thresh, "AD", "CTRL")
    rf_correct  <- rf_class == ifelse(y == 1, "AD", "CTRL")
    
    adj_expression <- {
      df_full <- data.frame(x = as.numeric(X[, g]), sex = as.factor(meta_use$Sex),
                            age = as.numeric(meta_use$Age))
      has_sv  <- length(unique(df_full$sex)) > 1
      lm_cov  <- tryCatch(if (has_sv) lm(x ~ sex + age, data = df_full) else
        lm(x ~ age, data = df_full), error = function(e) NULL)
      if (!is.null(lm_cov)) residuals(lm_cov) else df_full$x - mean(df_full$x)
    }
    
    plot_df <- data.frame(
      raw_expression = as.numeric(X[, g]),
      adj_expression = as.numeric(adj_expression),
      pred_rf        = pred_rf,
      pred_enet      = pred_enet,
      pred_glm       = pred_glm_g,
      actual_diag    = ifelse(y == 1, "AD", "CTRL"),
      group          = as.factor(meta_use[["Diag"]]),
      sex            = as.factor(meta_use$Sex),
      age            = as.numeric(meta_use$Age),
      sample         = rownames(meta_use),
      rf_class       = rf_class,
      rf_correct     = rf_correct,
      stringsAsFactors = FALSE
    )
    
    p_raw <- ggplot(plot_df, aes(x = group, y = raw_expression, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(aes(shape = sex, color = age), width = 0.2, size = 2.5, alpha = 0.9) +
      scale_shape_manual(values = c("M" = 16, "F" = 17), na.value = 15) +
      scale_color_gradient(low = "lightblue", high = "darkred", name = "Age") +
      labs(title = "A: Raw expression", x = "Diagnosis", y = "Expression",
           fill = "Diagnosis", shape = "Sex") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    p_adj <- ggplot(plot_df, aes(x = group, y = adj_expression, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(aes(shape = sex, color = age), width = 0.2, size = 2.5, alpha = 0.9) +
      scale_shape_manual(values = c("M" = 16, "F" = 17), na.value = 15) +
      scale_color_gradient(low = "lightblue", high = "darkred", name = "Age") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      labs(title = "B: Adjusted expression (sex + age residuals)",
           x = "Diagnosis", y = "Adjusted Expression", fill = "Diagnosis", shape = "Sex") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Panel C: all three model predictions overlaid
    plot_df_long <- rbind(
      data.frame(sample = plot_df$sample, pred = pred_rf,    model = "RF",
                 actual_diag = plot_df$actual_diag, correct = rf_correct),
      data.frame(sample = plot_df$sample, pred = pred_enet,  model = "glmnet",
                 actual_diag = plot_df$actual_diag,
                 correct = (ifelse(pred_enet >= enet_opt$threshold, "AD", "CTRL") ==
                              plot_df$actual_diag)),
      data.frame(sample = plot_df$sample, pred = pred_glm_g, model = "glm",
                 actual_diag = plot_df$actual_diag,
                 correct = (ifelse(pred_glm_g >= glm_opt$threshold, "AD", "CTRL") ==
                              plot_df$actual_diag))
    )
    
    p_pred <- ggplot(plot_df_long,
                     aes(x = reorder(sample, pred), y = pred,
                         color = actual_diag, shape = model)) +
      geom_hline(yintercept = rf_thresh, linetype = "dashed",
                 color = "darkgreen", linewidth = 0.7) +
      geom_hline(yintercept = 0.5, linetype = "dotted",
                 color = "grey40", linewidth = 0.5) +
      geom_point(size = 2.5, alpha = 0.85) +
      annotate("text", x = 1, y = rf_thresh + 0.03,
               label = paste0("RF threshold: ", round(rf_thresh, 2)),
               color = "darkgreen", size = 2.8, hjust = 0) +
      scale_color_manual(values = c("AD" = "tomato", "CTRL" = "steelblue"),
                         name = "Actual Diagnosis") +
      scale_shape_manual(values = c("RF" = 16, "glmnet" = 17, "glm" = 15),
                         name = "Model") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(title = "C: LOOCV P(AD) — RF (●), glmnet (▲), glm (■)",
           x = "Sample (ordered by RF probability)", y = "Predicted P(AD)") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
            legend.position = "bottom")
    
    # Panel D: confusion matrix (RF only)
    make_conf_df <- function(predicted, thresh_label) {
      df <- data.frame(actual_diag = ifelse(y == 1, "AD", "CTRL"),
                       predicted_diag = predicted, stringsAsFactors = FALSE) %>%
        group_by(actual_diag, predicted_diag) %>%
        summarize(n = n(), .groups = "drop") %>%
        mutate(correct = actual_diag == predicted_diag, threshold = thresh_label)
      all_combos <- expand.grid(actual_diag = c("AD", "CTRL"),
                                predicted_diag = c("AD", "CTRL"),
                                threshold = thresh_label, stringsAsFactors = FALSE)
      df <- merge(all_combos, df, all.x = TRUE)
      df$n[is.na(df$n)] <- 0
      df$correct[is.na(df$correct)] <- FALSE
      df
    }
    
    conf_rf   <- make_conf_df(rf_class,
                              paste0("RF (", round(rf_thresh, 2), ")\n",
                                     "BalAcc=", round(rf_opt$balacc * 100, 1), "%"))
    conf_enet <- make_conf_df(
      ifelse(pred_enet >= enet_opt$threshold, "AD", "CTRL"),
      paste0("glmnet (", round(enet_opt$threshold, 2), ")\n",
             "BalAcc=", round(enet_opt$balacc * 100, 1), "%")
    )
    conf_both <- bind_rows(conf_rf, conf_enet)
    
    p_conf <- ggplot(conf_both,
                     aes(x = predicted_diag, y = actual_diag, fill = correct, label = n)) +
      geom_tile(alpha = 0.7) +
      geom_text(size = 6, fontface = "bold") +
      facet_wrap(~ threshold) +
      scale_fill_manual(values = c("TRUE" = "lightgreen", "FALSE" = "lightsalmon"),
                        guide = "none") +
      labs(title = "D: Confusion matrix (RF vs glmnet)", x = "Predicted", y = "Actual") +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold", size = 8))
    
    p_combined <- cowplot::plot_grid(p_raw, p_adj, p_pred, p_conf,
                                     ncol = 2, rel_heights = c(1, 1.3))
    title_grob <- cowplot::ggdraw() +
      cowplot::draw_label(
        paste0(sym, " (", g, ")",
               "  |  RF AUC:", round(rf_ci$auc, 3),
               " [", round(rf_ci$ci_low, 3), ", ", round(rf_ci$ci_high, 3), "]",
               "  |  RF BalAcc:", round(rf_opt$balacc * 100, 1), "%",
               "  |  glmnet BalAcc:", round(enet_opt$balacc * 100, 1), "%",
               "  |  glm BalAcc:", round(glm_opt$balacc * 100, 1), "%",
               "  |  Gain:", round(balacc_gain * 100, 1), "%",
               "  |  alpha:", best_alpha_vec[g], "| node:", best_node_vec[g],
               "  |  ", label),
        fontface = "bold", size = 7.5, x = 0.01, hjust = 0
      )
    boxplot_list[[g]] <- cowplot::plot_grid(title_grob, p_combined,
                                            ncol = 1, rel_heights = c(0.04, 1))
  }
  
  # ── Save per-gene PDF ─────────────────────────────────────────────────────
  combined_pdf_path <- file.path(output_dir,
                                 paste0(label, "_rf_all_gene_plots_", today, ".pdf"))
  message("Saving RF gene plot PDF: ", combined_pdf_path)
  balacc_tmp <- sapply(names(boxplot_list), function(g) {
    p_rf   <- pred_store_rf[[g]]
    p_enet <- pred_store_enet[[g]]
    p_glm  <- pred_store_glm[[g]]
    if (is.null(p_rf)) return(NA_real_)
    mean(c(
      find_opt_threshold(p_rf,   y)$balacc,
      find_opt_threshold(p_enet, y)$balacc,
      find_opt_threshold(p_glm,  y)$balacc
    ), na.rm = TRUE)
  })
  sort_df <- data.frame(
    g           = names(boxplot_list),
    mean_balacc = balacc_tmp,
    mean_auc    = rowMeans(cbind(auc_vec_rf[names(boxplot_list)],
                                 auc_vec_enet[names(boxplot_list)],
                                 auc_vec_glm[names(boxplot_list)]),
                           na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  sort_df <- sort_df[order(-sort_df$mean_balacc, -sort_df$mean_auc, na.last = TRUE), ]
  pdf(combined_pdf_path, width = 14, height = 10)
  for (g in sort_df$g) print(boxplot_list[[g]])
  dev.off()
  message("Saved ", length(boxplot_list), " gene plots to ", combined_pdf_path)
  
  # Calculate permutation-corrected p-values for the best RF model and elastic net for each gene
  perm_res_rf   <- permutation_pvalue(p_rf,   y, B_perm = B_perm)
  perm_res_enet <- permutation_pvalue(p_enet, y, B_perm = B_perm)
  perm_res_glm  <- permutation_pvalue(p_glm,  y, B_perm = B_perm)
  
  # ── Build gene summary ────────────────────────────────────────────────────
  gene_summary_rows <- lapply(gene_names, function(g) {
    p_rf   <- pred_store_rf[[g]]
    p_enet <- pred_store_enet[[g]]
    p_glm  <- pred_store_glm[[g]]
    if (is.null(p_rf)) return(NULL)
    rf_check <- find_opt_threshold(p_rf, y)
    if (is.na(rf_check$balacc)) return(NULL)
    
    rf_opt   <- find_opt_threshold(p_rf,   y)
    enet_opt <- find_opt_threshold(p_enet, y)
    glm_opt  <- find_opt_threshold(p_glm,  y)
    
    sym            <- hgnc_symbols[match(g, ensembl_ids)]; sym <- if (!is.na(sym)) sym else g
    balacc_gain    <- rf_opt$balacc - null_balacc
    
    chr_annot  <- if (!is.null(gene_chr_map)) gene_chr_map[g] else NA_character_
    on_sex_chr <- !is.na(chr_annot) & chr_annot %in% c("X", "chrX", "Y", "chrY")
    
    data.frame(
      ensembl_id        = g,
      gene_symbol       = sym,
      chromosome        = ifelse(is.na(chr_annot), "unknown", chr_annot),
      on_sex_chr        = on_sex_chr,
      # RF
      auc_rf            = round(auc_vec_rf[g],    3),
      ci_low_rf         = round(ci_low_rf[g],     3),
      ci_high_rf        = round(ci_high_rf[g],    3),
      balacc_opt_rf     = round(rf_opt$balacc,    3),
      sensitivity_rf    = round(rf_opt$sensitivity, 3),
      specificity_rf    = round(rf_opt$specificity, 3),
      best_node         = best_node_vec[g],
      # glmnet
      auc_enet          = round(auc_vec_enet[g],  3),
      ci_low_enet       = round(ci_low_enet[g],   3),
      ci_high_enet      = round(ci_high_enet[g],  3),
      balacc_opt_enet   = round(enet_opt$balacc,  3),
      sensitivity_enet  = round(enet_opt$sensitivity, 3),
      specificity_enet  = round(enet_opt$specificity, 3),
      best_alpha        = best_alpha_vec[g],
      # glm
      auc_glm           = round(auc_vec_glm[g],   3),
      ci_low_glm        = round(ci_low_glm[g],    3),
      ci_high_glm       = round(ci_high_glm[g],   3),
      balacc_opt_glm    = round(glm_opt$balacc,   3),
      sensitivity_glm   = round(glm_opt$sensitivity, 3),
      specificity_glm   = round(glm_opt$specificity, 3),
      # Permutation Results
      perm_pval_rf   = round(perm_res_rf$p_value,   4),
      perm_pval_enet = round(perm_res_enet$p_value, 4),
      perm_pval_glm  = round(perm_res_glm$p_value,  4),
      # Null comparison (descriptive)
      null_balacc       = round(null_balacc,      3),
      balacc_gain       = round(balacc_gain,      3),
      n_total           = length(y),
      n_AD              = sum(y == 1),
      n_CTRL            = sum(y == 0),
      stringsAsFactors  = FALSE
    )
  })
  
  gene_summary_df <- do.call(rbind, Filter(Negate(is.null), gene_summary_rows))
  
  gene_summary_df$pval_bh_rf   <- p.adjust(gene_summary_df$perm_pval_rf,   method = "BH")
  gene_summary_df$pval_bh_enet <- p.adjust(gene_summary_df$perm_pval_enet, method = "BH")
  gene_summary_df$pval_bh_glm  <- p.adjust(gene_summary_df$perm_pval_glm,  method = "BH")
  # Sig flag: significant in at least one model
  gene_summary_df$sig_any_bh   <- (gene_summary_df$pval_bh_rf   < 0.05 |
                                     gene_summary_df$pval_bh_enet < 0.05 |
                                     gene_summary_df$pval_bh_glm  < 0.05)
  
  # ── Rank by mean BalAcc across all three models ───────────────────────────
  # This is the primary sort key: genes that perform consistently well across
  # all methods rank higher than those that are strong in only one model.
  gene_summary_df$mean_balacc <- rowMeans(
    gene_summary_df[, c("balacc_opt_rf", "balacc_opt_enet", "balacc_opt_glm")],
    na.rm = TRUE
  )
  # Secondary sort: best single-model BalAcc (captures genes with one standout model)
  gene_summary_df$max_balacc <- pmax(
    gene_summary_df$balacc_opt_rf,
    gene_summary_df$balacc_opt_enet,
    gene_summary_df$balacc_opt_glm,
    na.rm = TRUE
  )
  # Tertiary: mean AUC (for ties)
  gene_summary_df$mean_auc <- rowMeans(
    gene_summary_df[, c("auc_rf", "auc_enet", "auc_glm")],
    na.rm = TRUE
  )
  gene_summary_df <- gene_summary_df[
    order(-gene_summary_df$mean_balacc,
          -gene_summary_df$max_balacc,
          -gene_summary_df$mean_auc), ]
  
  write.csv(gene_summary_df,
            file.path(output_dir, paste0("per_gene_summary_rf_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  # ── per_gene_df_sorted: used for PDF plot ordering and downstream ─────────
  per_gene_df        <- gene_summary_df
  per_gene_df_sorted <- per_gene_df[!is.na(per_gene_df$mean_balacc), ]
  # Already sorted above by mean_balacc; AUC is retained as tiebreaker
  # Strength bands based on mean BalAcc — primary metric
  per_gene_df_sorted$strength_rf <- cut(
    per_gene_df_sorted$mean_balacc,
    breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels = c("Poor (<0.6)", "Fair (0.6-0.7)", "Good (0.7-0.8)",
               "Strong (0.8-0.9)", "Excellent (>0.9)"),
    include.lowest = TRUE
  )
  
  # ── Top-20 by mean BalAcc ─────────────────────────────────────────────────
  top20            <- head(per_gene_df_sorted, 20)
  top20_sex_chr    <- top20[!is.na(top20$on_sex_chr) & top20$on_sex_chr, ]
  
  # Gene order on y-axis: ascending mean_balacc so best gene is at top
  gene_order <- top20$gene_symbol[order(top20$mean_balacc)]   # low → high
  n_genes    <- length(gene_order)
  
  # Dodge offsets: RF top, glmnet middle, glm bottom within each gene row
  # Row height = 1 unit; models spaced ±0.25 around centre
  model_offset <- c("RF" = 0.25, "glmnet" = 0, "glm" = -0.25)
  model_colors <- c("RF" = "#e41a1c", "glmnet" = "#377eb8", "glm" = "#888888")
  model_shapes <- c("RF" = 16,       "glmnet" = 17,          "glm" = 15)
  
  # Build long data frame with explicit numeric y positions
  make_long_row <- function(gs, model_name, auc, ci_lo, ci_hi, balacc) {
    data.frame(
      gene_symbol = gs,
      model       = model_name,
      auc         = auc,
      ci_low      = ci_lo,
      ci_high     = ci_hi,
      balacc_opt  = balacc,
      y_centre    = match(gs, gene_order),          # integer rank (1 = worst shown)
      y_pos       = match(gs, gene_order) + model_offset[model_name],
      stringsAsFactors = FALSE
    )
  }
  
  top20_long <- do.call(rbind, lapply(seq_len(nrow(top20)), function(i) {
    gs <- top20$gene_symbol[i]
    rbind(
      make_long_row(gs, "RF",
                    top20$auc_rf[i],   top20$ci_low_rf[i],   top20$ci_high_rf[i],
                    top20$balacc_opt_rf[i]),
      make_long_row(gs, "glmnet",
                    top20$auc_enet[i], top20$ci_low_enet[i], top20$ci_high_enet[i],
                    top20$balacc_opt_enet[i]),
      make_long_row(gs, "glm",
                    top20$auc_glm[i],  top20$ci_low_glm[i],  top20$ci_high_glm[i],
                    top20$balacc_opt_glm[i])
    )
  }))
  top20_long$model <- factor(top20_long$model, levels = c("RF", "glmnet", "glm"))
  
  # Mean BalAcc reference line per gene (plotted as a thin grey segment)
  top20_mean <- data.frame(
    gene_symbol = top20$gene_symbol,
    mean_balacc = top20$mean_balacc,
    y_centre    = match(top20$gene_symbol, gene_order),
    stringsAsFactors = FALSE
  )
  
  # Sex-chromosome gene y positions and labels
  sex_chr_y      <- if (nrow(top20_sex_chr) > 0)
    match(top20_sex_chr$gene_symbol, gene_order) else numeric(0)
  sex_chr_labels <- if (nrow(top20_sex_chr) > 0)
    top20_sex_chr$chromosome else character(0)
  
  p_biomarker <- ggplot(top20_long,
                        aes(y = y_pos, color = model, shape = model)) +
    
    # Horizontal reference lines at each gene's integer y (subtle grid)
    geom_hline(data = top20_mean, aes(yintercept = y_centre),
               color = "grey92", linewidth = 0.4, inherit.aes = FALSE) +
    
    # Mean BalAcc tick: thin dark vertical segment per gene
    geom_segment(data = top20_mean,
                 aes(x = mean_balacc, xend = mean_balacc,
                     y = y_centre - 0.42, yend = y_centre + 0.42),
                 color = "grey40", linewidth = 0.5, linetype = "solid",
                 inherit.aes = FALSE) +
    
    # AUC CI bars
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.10, linewidth = 0.55) +
    
    # AUC point
    geom_point(aes(x = auc), size = 3.2) +
    
    # BalAcc(opt) open diamond — same colour as model, no fill
    geom_point(aes(x = balacc_opt), size = 2.8, shape = 5) +
    
    # Sex-chromosome genes: small purple chromosome label (informational, not a penalty)
    {if (length(sex_chr_y) > 0)
      annotate("label", x = 0.405, y = sex_chr_y,
               label = sex_chr_labels,
               size = 2.4, color = "purple", fill = "white",
               label.padding = unit(0.10, "lines"), label.size = 0.25)
    } +
    
    scale_color_manual(values = model_colors, name = "Model") +
    scale_shape_manual(values = model_shapes, name = "Model") +
    
    # y-axis: one label per gene at integer positions
    scale_y_continuous(
      breaks = seq_len(n_genes),
      labels = gene_order,
      expand = expansion(add = 0.6)
    ) +
    
    geom_vline(xintercept = 0.5, linetype = "dashed",  color = "grey55", linewidth = 0.5) +
    geom_vline(xintercept = 0.7, linetype = "dotted",  color = "orange", linewidth = 0.6) +
    geom_vline(xintercept = 0.8, linetype = "dotted",  color = "red",    linewidth = 0.6) +
    
    annotate("text", x = 0.435, y = 0.05,
             label = paste0("● AUC + 95% CI   ◇ BalAcc(opt)   ",
                            "| = mean BalAcc across models   purple = sex chromosome gene"),
             hjust = 0, vjust = 0, size = 2.6, color = "grey35") +
    
    scale_x_continuous(limits = c(0.4, 1.02), breaks = seq(0.4, 1.0, 0.1)) +
    
    labs(
      title    = paste0("Top 20 AD Biomarkers — RF vs glmnet vs glm — ", label),
      subtitle = paste0(
        "Ranked by mean BalAcc(opt) across all three models",
        "  |  Null (age+sex) BalAcc: ", round(null_balacc * 100, 1), "%",
        "  |  balacc_gain = gene increment over demographic baseline"
      ),
      x     = "Metric value  (AUC ● with CI bars;  BalAcc ◇)",
      y     = "Gene  (ascending mean BalAcc — best at top)"
    ) +
    
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 10),
      plot.subtitle   = element_text(size = 8, color = "grey30"),
      panel.grid.major.y = element_blank(),   # suppress default grid; we drew our own
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom",
      legend.key.size    = unit(0.9, "lines")
    )
  
  ggsave(file.path(output_dir, paste0("top_biomarkers_rf_", label, "_", today, ".pdf")),
         plot   = p_biomarker,
         width  = 10,
         height = max(6, n_genes * 0.52 + 2),
         units  = "in",
         device = cairo_pdf)
  
  write.csv(per_gene_df_sorted,
            file.path(output_dir, paste0("per_gene_rf_summary_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  # ── Multi-gene RF ─────────────────────────────────────────────────────────
  cat("\nFitting multi-gene RF (LOOCV)...\n")
  sex_numeric  <- ifelse(toupper(meta_use$Sex) == "M", 1, 0)
  X_multi      <- cbind(X, Age = as.numeric(meta_use$Age), Sex = sex_numeric)
  
  # Tune min.node.size once on full data with inner CV
  cat("  Tuning multi-gene RF...\n")
  df_multi_full    <- as.data.frame(X_multi)
  tune_multi       <- tune_rf(df_multi_full, y,
                              node_grid    = node_grid,
                              nfolds_inner = min(nfolds_inner, n),
                              num.trees    = num.trees)
  best_node_multi  <- tune_multi$min.node.size
  cat("  Best min.node.size (multi-gene):", best_node_multi, "\n")
  
  pred_multi_rf <- numeric(n)
  importance_list <- list()
  
  for (i in 1:n) {
    train_idx   <- setdiff(1:n, i)
    df_tr_multi <- as.data.frame(X_multi[train_idx, , drop = FALSE])
    df_te_multi <- as.data.frame(X_multi[i,          , drop = FALSE])
    y_tr        <- y[train_idx]
    
    if (length(unique(y_tr)) < 2) { pred_multi_rf[i] <- 0.5; next }
    df_tr_multi$y <- factor(y_tr, levels = c(0, 1))
    
    fit_multi <- tryCatch(
      ranger::ranger(
        y             ~ .,
        data          = df_tr_multi,
        num.trees     = num.trees,
        mtry          = max(1, floor(sqrt(ncol(df_tr_multi) - 1))),
        min.node.size = best_node_multi,
        probability   = TRUE,
        importance    = "permutation",
        num.threads   = 1
      ),
      error = function(e) NULL
    )
    if (is.null(fit_multi)) { pred_multi_rf[i] <- 0.5; next }
    
    pred_multi_rf[i] <- tryCatch(
      predict(fit_multi, data = df_te_multi)$predictions[, "1"],
      error = function(e) 0.5
    )
    importance_list[[i]] <- fit_multi$variable.importance
  }
  
  # Aggregate permutation importance across folds (mean, SD)
  imp_mat <- do.call(rbind, Filter(Negate(is.null), importance_list))
  imp_df  <- data.frame(
    feature       = colnames(imp_mat),
    importance    = colMeans(imp_mat, na.rm = TRUE),
    importance_sd = apply(imp_mat, 2, sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  imp_df <- imp_df[!imp_df$feature %in% c("Age", "Sex"), ]
  imp_df$gene_symbol <- hgnc_symbols[match(imp_df$feature, ensembl_ids)]
  imp_df$gene_symbol[is.na(imp_df$gene_symbol)] <- imp_df$feature[is.na(imp_df$gene_symbol)]
  imp_df <- imp_df[order(-imp_df$importance), ]
  
  write.csv(imp_df,
            file.path(output_dir, paste0("rf_importance_multigene_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  # Multi-gene RF stats
  roc_multi_rf  <- tryCatch(roc(response = y, predictor = pred_multi_rf, quiet = TRUE),
                            error = function(e) NULL)
  auc_multi_rf  <- if (!is.null(roc_multi_rf)) as.numeric(auc(roc_multi_rf)) else NA
  boot_multi_rf <- bootstrap_auc_ci(y, pred_multi_rf, B = B_boot, conf = 0.95)
  rf_multi_opt  <- find_opt_threshold(pred_multi_rf, y)
  
  cat("  Multi-gene RF AUC:", round(auc_multi_rf, 3),
      "| BalAcc(opt):", round(rf_multi_opt$balacc * 100, 1), "%\n")
  
  # ── Importance plot ───────────────────────────────────────────────────────
  top_imp <- head(imp_df, 20)
  top_imp$ci_low  <- top_imp$importance - 1.96 * top_imp$importance_sd
  top_imp$ci_high <- top_imp$importance + 1.96 * top_imp$importance_sd
  
  p_importance <- ggplot(top_imp,
                         aes(x = importance,
                             y = reorder(gene_symbol, importance))) +
    geom_col(aes(fill = importance > 0), alpha = 0.7, width = 0.6) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.3, linewidth = 0.7, color = "grey30") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue"),
                      guide = "none") +
    labs(
      title    = paste0("RF Permutation Importance — Multi-gene — ", label),
      subtitle = paste0(
        "AUC: ", round(auc_multi_rf, 3),
        " [", round(boot_multi_rf$ci_low, 3), ", ", round(boot_multi_rf$ci_high, 3), "]",
        "  BalAcc(opt)=", round(rf_multi_opt$balacc * 100, 1), "%",
        "  Sens=", round(rf_multi_opt$sensitivity * 100, 1), "%",
        "  Spec=", round(rf_multi_opt$specificity * 100, 1), "%\n",
        "Error bars = ±1.96 SD across LOOCV folds   best_node=", best_node_multi
      ),
      x = "Mean permutation importance", y = "Gene"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8, color = "grey30"),
          legend.position = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))
  
  ggsave(file.path(output_dir, paste0("rf_importance_multigene_", label, "_", today, ".pdf")),
         plot   = p_importance,
         width  = 7,
         height = max(4, nrow(top_imp) * 0.45),
         units  = "in",
         device = cairo_pdf)
  
  # ── Bootstrap histogram (multi-gene RF) ───────────────────────────────────
  boot_df_rf <- data.frame(iteration = seq_along(boot_multi_rf$auc_boot),
                           auc = boot_multi_rf$auc_boot)
  rf_balacc_strength <- cut(rf_multi_opt$balacc, breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
                            labels = c("Poor", "Fair", "Good", "Strong", "Excellent"),
                            include.lowest = TRUE)
  
  p_boot_rf <- ggplot(boot_df_rf, aes(x = auc)) +
    geom_histogram(color = "black", fill = "#e41a1c", alpha = 0.7, bins = 30) +
    geom_vline(xintercept = auc_multi_rf,          color = "darkred",  linewidth = 1.1) +
    geom_vline(xintercept = boot_multi_rf$ci_low,  color = "darkred",  linetype = "dashed") +
    geom_vline(xintercept = boot_multi_rf$ci_high, color = "darkred",  linetype = "dashed") +
    geom_vline(xintercept = rf_multi_opt$balacc,   color = "tomato",   linetype = "dashed") +
    geom_vline(xintercept = 0.7,                   color = "orange",   linetype = "dotted") +
    geom_vline(xintercept = 0.8,                   color = "darkred",  linetype = "dotted") +
    annotate("text", x = rf_multi_opt$balacc, y = Inf,
             label = paste0("BalAcc(opt)=", round(rf_multi_opt$balacc, 3),
                            " (", rf_balacc_strength, ")"),
             vjust = 2, hjust = -0.1, color = "tomato", size = 3.5) +
    annotate("text", x = auc_multi_rf, y = Inf,
             label = paste0("AUC=", round(auc_multi_rf, 3)),
             vjust = 4, hjust = -0.1, color = "darkred", size = 3.2) +
    theme_bw(base_size = 10) +
    labs(title = paste0("Multi-gene RF LOOCV — ", label,
                        "\nBalAcc(opt): ", round(rf_multi_opt$balacc * 100, 1), "%",
                        " (", rf_balacc_strength, ")",
                        "  |  Sens: ", round(rf_multi_opt$sensitivity * 100, 1), "%",
                        "  |  Spec: ", round(rf_multi_opt$specificity * 100, 1), "%",
                        "\nAUC: ", round(boot_multi_rf$mean, 3),
                        " CI [", round(boot_multi_rf$ci_low, 3), ", ",
                        round(boot_multi_rf$ci_high, 3), "]",
                        "  best_node=", best_node_multi),
         x = "Bootstrap AUC", y = "Count")
  
  ggsave(file.path(output_dir, paste0("rf_multigene_bootstrap_hist_", label, "_", today, ".png")),
         p_boot_rf, width = 6, height = 4, dpi = 300)
  
  cat("RF analysis done:", label, "\n")
  cat("  Multi-gene RF AUC:", round(auc_multi_rf, 3),
      "CI [", round(boot_multi_rf$ci_low, 3), ",", round(boot_multi_rf$ci_high, 3), "]\n")
  
  list(
    per_gene_df         = per_gene_df_sorted,
    gene_summary_df     = gene_summary_df,
    auc_multi_rf        = auc_multi_rf,
    boot_mean_rf        = boot_multi_rf$mean,
    boot_ci_low_rf      = boot_multi_rf$ci_low,
    boot_ci_high_rf     = boot_multi_rf$ci_high,
    rf_importance       = imp_df,
    best_node_multi     = best_node_multi,
    top20               = top20          # includes mean_balacc, max_balacc, mean_auc columns
  )
}