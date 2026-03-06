#!/usr/bin/env Rscript
# RNApreprocessing_gene_kmer_analysisAndrobustness_pipeline.R
# Full Robustness Pipeline — HPC-Optimized + Parallel + Bayesian + Publication Ready

# Prevent broken rstan from loading — use cmdstanr throughout
options(brms.backend  = "cmdstanr")
options(mc.cores      = parallel::detectCores())
options(ranger.num.threads = 10)

# Intercept any attempt to load rstan
.rstan_blocked <- TRUE
if (!"cmdstanr" %in% rownames(installed.packages())) {
  stop("cmdstanr required — install with: install.packages('cmdstanr', repos='https://mc-stan.org/r-packages/')")
}

# install.packages("DESeq2") if needed

suppressPackageStartupMessages({
  library(openxlsx)
  library(data.table)
  library(tidyverse)
  library(brms)
  library(pROC)
  library(loo)
  library(future)
  library(future.apply)
  library(posterior)
  library(bayesplot)
  library(parallel)
  library(matrixStats)
  library(igraph)
  library(ComplexHeatmap)
  library(circlize)
  library(UpSetR)
  library(reshape2)
  library(stringi)
  library(stringr)
  library(progressr)
  library(ggseqlogo)
  library(Biostrings)
  library(DESeq2)          # varianceStabilizingTransformation, DESeqDataSetFromMatrix
  library(glmnet)          # cv.glmnet (elastic net)
  library(pROC)            # roc, auc
  library(ggplot2)         # all plotting
  library(dplyr)           # filter, mutate, arrange, etc.
  library(tidyr)           # unnest
  library(tibble)          # enframe
  library(stringr)         # str_squish
  library(igraph)          # graph_from_data_frame, cluster_louvain
  library(ggraph)          # ggraph, geom_edge_link, geom_node_point
  library(RColorBrewer)    # brewer.pal
  library(piano)           # runGSA, loadGSC, GSAsummaryTable
  library(DESeq2)
  source('run_analysis_rf.R')
  source('run_lodo_validation.R')
})

# Prevent dplyr from masking data.table functions
first <- data.table::first
last  <- data.table::last

# Setting memory usage per core for future.apply (adjust as needed for your HPC environment)
options(future.globals.maxSize = 10000 * 1024^2)  # 500MB per worker

# Setting output directory
output_dir <- "/home/thilsabeck/Documents/Gage2024/RNApreprocessing_gene_validation_files/"

# Ensure output_dir exists
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Loading multiple Gage2019 count versions for comparison of gene prediction accuracy across different count matrices and pipelines, I think these were done with GSNAP and with different genome versions
Gage2019_counts_AD_2019genome <- read.csv('/home/thilsabeck/Documents/Gage2019/mapping/paired_ensembl_trim/Fibroblast/merged_counts_AD_2019genome_SendMingchen.csv')
Gage2019_counts_noAD_2019genome <- read.csv('/home/thilsabeck/Documents/Gage2019/mapping/paired_ensembl_trim/Fibroblast/merged_counts_2019genome_noAD_SendMingchen.csv')
Gage2019_counts_2023genome_AD <- read.csv('/home/thilsabeck/Documents/Gage2019/mapping/paired_ensembl_trim/Fibroblast/merged_counts_AD_2023genome_SendMingchen.csv')
Gage2019_counts_2023genome_noAD <- read.csv('/home/thilsabeck/Documents/Gage2019/mapping/paired_ensembl_trim/Fibroblast/merged_counts_2023genome_noAD_SendMingchen.csv')
Gage2019_counts_2019genome_noAD <- read.csv('/home/thilsabeck/Documents/Gage2019/RNAseq Counts/SeqFibroblast_counts_reduced.csv')

# Loading additional Gage 2024 RNAseq sample metadata for potential use in GLMMs or diagnostics, I think this was done with STAR single-pass and star counter
Gage2024_metadata <- read.xlsx('/home/thilsabeck/Documents/Gage2024/RNAseq metadata_updatedAge_11-21-2024_reduced.xlsx', sheet = 1)
Gage_merged_metadata <- read.xlsx('/home/thilsabeck/Documents/Gage2019/mapping/Gage_merged_MetaData_updated_10-23-2024_AgeMatched_FastqCheck_10-6-2025_02-28-2026update.xlsx')
# Remove rows with NA in Sample.label of Gage_merged_metadata
Gage_merged_metadata <- Gage_merged_metadata[!is.na(Gage_merged_metadata$Sample.label), ]
# Create rownames based on Sample.label column
rownames(Gage2024_metadata) <- Gage2024_metadata$Sample.label
rownames(Gage_merged_metadata) <- Gage_merged_metadata$Sample.label
Gage2024_counts_b1 <- read.csv('/home/thilsabeck/Documents/Gage2024/Gage2024_paired_raw_counts_batch1.csv')
Gage2024_counts_b2 <- read.csv('/home/thilsabeck/Documents/Gage2024/Gage2024_paired_raw_counts_batch2.csv')
# Mingchen fibroblast merged and processed counts
Gage2024_Mingchen_counts <- read.csv('/home/thilsabeck/Documents/Gage2019/BatchCorrect_redo/Gage_Mingchen_fibroblast_merged_LimmaOnlyBatchCorrected_ForDiag_verify_11-06-24.csv')
Gage2024_fibroblast_counts <- read.csv('/home/thilsabeck/Documents/Gage2019/BatchCorrect_redo/Gage_fibroblast_merged_LimmaOnlyBatchCorrected_ForDiag_verify_11-06-24.csv')

# Merging two count tables and ensuring consistent formatting
Gage2024_counts <- full_join(Gage2024_counts_b1, Gage2024_counts_b2, by = "ensembl_id")
rownames(Gage2024_counts) <- Gage2024_counts$ensembl_id
Gage2024_counts <- Gage2024_counts[,-1]  # Remove redundant ensembl_id column after setting rownames

# Loading stable gene list from pipeline robustness analysis with 17 DEG sets
stable_genes <- fread("/home/thilsabeck/Documents/Gage2019/mapping/kmer_analysis_full/tables/gene_cross_deg_reproducibility_2026-02-26.csv")
### Load processed stable_genes for 17 DEG sets and bypass below processing until the save step
# Load stable_genes with hgnc_symbol and z_score to a new CSV for reference
stable_genes <- read.csv(file.path(output_dir, "gene_cross_deg_reproducibility_2026-02-26_stable_genes_with_hgnc_and_zscore_03-02-2026.csv"))
# Get hgnc_symbol for each ensembl_id in stable_genes using biomart or biomaRt package, and add as a new column
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = stable_genes$ensembl_id,
  mart = ensembl
)
# Mapping ensembl_id to hgnc_symbol in stable_genes
stable_genes <- merge(stable_genes, gene_info, by.x = "ensembl_id", by.y = "ensembl_gene_id", all.x = TRUE)

# Calculate z-score for cross-deg-fraction to identify significantly stable genes (e.g., z < -2)
stable_genes[, z_score := (cross_deg_fraction - mean(cross_deg_fraction, na.rm=TRUE)) / sd(cross_deg_fraction, na.rm=TRUE)]
stable_genes[, is_stable := abs(z_score) > 2]
cat("Loaded stable_genes with", nrow(stable_genes), "genes, of which",
    sum(stable_genes$is_stable, na.rm=TRUE), "are classified as stable (z < -2)\n")

# Save stable_genes with hgnc_symbol and z_score to a new CSV for reference
write.csv(stable_genes, file.path(output_dir, paste0("gene_cross_deg_reproducibility_2026-02-26_stable_genes_with_hgnc_and_zscore_", today, ".csv")), row.names = FALSE)
### Skip to hear if loading previously saved stable_genes with hgnc_symbol and z_score to bypass above processing steps

# Compare stable genes across different Diag groups using Gage2024_metadata Diag column and Sample.label
# Calculate mean values by gene in stable genes for genes with is_stable == TRUE by each Diag group
Gage2024_metadata$Diag <- factor(Gage2024_metadata$Diag, levels = c("CTRL", "AD", "Other"))
# for (gene in stable_genes[stable_genes$is_stable == TRUE,]$ensembl_id) {
#   cat(sprintf("Calculating mean counts for stable gene %s across Diag groups...\n", gene))
#   AD_mean <- Gage2024_counts %>%
#     filter(ensembl_id == gene) %>%
#     pivot_longer(-ensembl_id, names_to = "Sample", values_to = "Count") %>%
#     left_join(Gage2024_metadata, by = c("Sample" = "Sample.label")) %>%
#     filter(Diag == "AD") %>%
#     summarise(mean_count = mean(Count, na.rm=TRUE)) %>%
#     pull(mean_count)
#   CTRL_mean <- Gage2024_counts %>%
#     filter(ensembl_id == gene) %>%
#     pivot_longer(-ensembl_id, names_to = "Sample", values_to = "Count") %>%
#     left_join(Gage2024_metadata, by = c("Sample" = "Sample.label")) %>%
#     filter(Diag == "CTRL") %>%
#     summarise(mean_count = mean(Count, na.rm=TRUE)) %>%
#     pull(mean_count)
# }
# 
# stable_gene_means <- Gage2024_counts %>%
#   filter(ensembl_id %in% stable_genes[stable_genes$is_stable == TRUE,]$ensembl_id) %>%
#   pivot_longer(-ensembl_id, names_to = "Sample", values_to = "Count") %>%
#   left_join(Gage2024_metadata, by = c("Sample" = "Sample.label")) %>%
#   group_by(Diag) %>%
#   summarise(mean_count = mean(Count, na.rm=TRUE))

# Setup code to predict how well each stable_gene can predict Diag group using ROC AUC, comparing AD vs CTRL samples, and combination of multiple stable genes in this prediction accuracy
# Keep only genes present in both the counts and your ranked list
# 1. Sort by z-score (descending)
stable_genes <- stable_genes[order(stable_genes$z_score, decreasing = TRUE), ]

# 2. Keep only stable genes
stable_genes_stable <- subset(stable_genes, is_stable == TRUE)

# 3. Ranked Ensembl IDs (by z-score, among stable genes)
ranked_genes_ensemblID <- stable_genes_stable$ensembl_id
ranked_genes_hgnc_symbol <- stable_genes_stable$hgnc_symbol

# 4. Genes in common with the new count matrix for all counts matrix comparisons, we will use the same ranked list of genes but subset the count matrix to only those genes that are present in the ranked list and in the count matrix, and preserve the ranking from the original list
common_genes <- intersect(rownames(Gage2024_counts), ranked_genes)


# Preserve ranking from your original list
ranked_common_genes <- ranked_genes[ranked_genes %in% common_genes]

# Gage2024_counts_sub <- Gage2024_counts[ranked_common_genes, , drop = FALSE]
# 
# # Sample metadata with diagnosis labels
# # meta_new: data.frame with rownames = sample IDs, column Diag ∈ c("AD","CTRL")
# # Subset Gage2024_counts_sub by samples present in Gage2024_metadata and ensure order matches
# Gage2024_counts_sub <- Gage2024_counts_sub[, rownames(Gage2024_metadata), drop = TRUE]
# 
# dds <- DESeqDataSetFromMatrix(
#   countData = round(Gage2024_counts_sub),
#   colData   = Gage2024_metadata,
#   design    = ~ Diag
# )
# 
# # Filter very lowly expressed genes
# dds <- dds[rowSums(counts(dds)) > 1, ]
# 
# # Variance stabilizing transform
# vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
# expr <- assay(vsd)          # genes x samples
# expr <- t(expr)             # samples x genes (rows = samples, columns = genes)
# 
# # Binary outcome: 1 = AD, 0 = CTRL
# y <- ifelse(Gage2024_metadata$Diag == "AD", 1, 0)
# X <- as.matrix(expr)                 # samples x genes
# gene_names <- colnames(X)

### Performing above subsetting and normalization steps for each count matrix and then running the same per-gene and multi-gene prediction accuracy pipeline for each count matrix to compare how well stable genes can predict Diag group across different count matrices and pipelines, and also to compare how well a multi-gene elastic net model performs across different count matrices and pipelines, and to see if the same stable genes are selected in the multi-gene model across different count matrices and pipelines
normalize_counts <- function(counts_df, metadata, ranked_genes, stable_genes,
                             design_formula = "~ Diag + Age + Sex",
                             dataset_name = "dataset") {
  
  # ── 1. Set rownames from first column if needed ──────────────────────────
  if (!is.numeric(counts_df[, 1])) {
    rownames(counts_df) <- counts_df[, 1]
    counts_df <- counts_df[, -1]
  }
  
  # ── 2. Detect gene ID type and match to stable_genes ────────────────────
  sample_genes <- rownames(counts_df)
  
  is_ensembl <- grepl("^ENSG", sample_genes[1])
  
  if (is_ensembl) {
    common_stable <- stable_genes$ensembl_id[stable_genes$ensembl_id %in% sample_genes]
    ranked_common <- ranked_genes_ensemblID[ranked_genes_ensemblID %in% stable_genes$ensembl_id[stable_genes$ensembl_id %in% common_stable]]
    # Subset by ensembl
    counts_sub <- counts_df[ranked_common, , drop = FALSE]
  } else {
    common_stable <- stable_genes$hgnc_symbol[stable_genes$hgnc_symbol %in% sample_genes]
    ranked_common <- ranked_genes_hgnc_symbol[ranked_genes_hgnc_symbol %in% common_stable]
    # Subset by hgnc_symbol
    counts_sub <- counts_df[ranked_common, , drop = FALSE]
  }
  
  message(dataset_name, ": ", nrow(counts_sub), " genes matched, ",
          length(ranked_common), " in ranked list")
  
  # ── 3. Match samples to metadata ────────────────────────────────────────
  # Try matching by rownames, then by Sample.label or X1
  meta <- metadata
  
  if (!all(colnames(counts_sub) %in% rownames(meta))) {
    # Try stripping _counts suffix
    colnames(counts_sub) <- gsub("_counts$", "", colnames(counts_sub))
  }
  
  common_samples <- intersect(colnames(counts_sub), rownames(meta))
  
  if (length(common_samples) == 0) {
    # Try matching via Sample.label column
    rownames(meta) <- meta$Sample.label
    common_samples <- intersect(colnames(counts_sub), rownames(meta))
  }
  
  if (length(common_samples) == 0) {
    # Try matching via X1 column
    rownames(meta) <- gsub("_counts$", "", meta$X1)
    common_samples <- intersect(colnames(counts_sub), rownames(meta))
  }
  
  if (length(common_samples) == 0) {
    warning(dataset_name, ": No samples matched metadata, skipping")
    return(NULL)
  }
  
  message(dataset_name, ": ", length(common_samples), " samples matched metadata")
  
  counts_sub <- counts_sub[, common_samples, drop = FALSE]
  meta <- meta[common_samples, , drop = FALSE]
  
  # Remove constant columns from design (e.g. all-male)
  design_vars <- all.vars(as.formula(design_formula))
  design_vars <- design_vars[sapply(meta[, design_vars, drop = FALSE],
                                    function(x) length(unique(na.omit(x))) > 1)]
  
  if (length(design_vars) == 0) {
    warning(dataset_name, ": No valid design variables, skipping")
    return(NULL)
  }
  
  formula <- as.formula(paste("~", paste(design_vars, collapse = " + ")))
  message(dataset_name, ": Using formula: ", deparse(formula))
  
  # ── 4. DESeq2 normalization ──────────────────────────────────────────────
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData   = meta,
    design    = formula
  )
  dds <- dds[rowSums(counts(dds)) > 1, ]
  
  # Check if data looks pre-normalized (non-integer or negative values)
  raw_counts     <- counts(dds)
  is_prenormalized <- any(raw_counts < 0) ||
    mean(raw_counts == round(raw_counts)) < 0.95
  
  if (is_prenormalized) {
    # Pre-normalized data: use directly or log2 transform
    raw_range <- range(as.matrix(counts_sub), na.rm = TRUE)
    if (raw_range[1] < 0 || raw_range[2] < 100) {
      message(dataset_name, ": Data appears already log-scaled, using directly.")
      expr <- t(as.matrix(counts_sub))
    } else {
      message(dataset_name, ": Data appears pre-normalized, applying log2(x+1).")
      expr <- t(log2(as.matrix(counts_sub) + 1))
    }
  } else {
    # Raw counts: try VST, fall back to normTransform
    normalized <- tryCatch({
      message(dataset_name, ": Applying VST normalization.")
      varianceStabilizingTransformation(dds, blind = TRUE)
    }, error = function(e) {
      message(dataset_name, ": VST failed, falling back to normTransform.")
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersionsGeneEst(dds)
      dispersions(dds) <- mcols(dds)$dispGeneEst
      normTransform(dds)
    })
    
    expr <- t(assay(normalized))  # samples x genes
  }
  
  # ── 5. Preserve ranking ──────────────────────────────────────────────────
  ranked_common <- ranked_common[ranked_common %in% colnames(expr)]
  expr <- expr[, ranked_common, drop = FALSE]
  
  # Convert rownames to ensembl_id only if input was ensembl
  # If input was hgnc, keep hgnc and look up ensembl safely
  if (is_ensembl) {
    # Already ensembl rownames, keep as-is
    gene_names_out <- colnames(expr)
  } else {
    # Keep hgnc as primary, add ensembl as attribute
    hgnc_to_ens <- setNames(stable_genes$ensembl_id, stable_genes$hgnc_symbol)
    ensembl_ids <- hgnc_to_ens[colnames(expr)]
    
    # Only rename if all conversions succeeded
    if (all(!is.na(ensembl_ids))) {
      gene_names_out <- ensembl_ids
      colnames(expr) <- ensembl_ids
    } else {
      n_failed <- sum(is.na(ensembl_ids))
      message(dataset_name, ": ", n_failed, " genes could not be converted to ensembl_id, keeping hgnc symbols")
      gene_names_out <- colnames(expr)  # keep hgnc
    }
  }
  
  # ── 6. Outcome vector ────────────────────────────────────────────────────
  y <- ifelse(meta$Diag == "AD", 1, 0)
  
  return(list(
    X            = as.matrix(expr),
    y            = y,
    gene_names   = gene_names_out,
    metadata     = meta,
    dataset_name = dataset_name,
    is_ensembl   = is_ensembl
  ))
}

# ── Define all count sets ────────────────────────────────────────────────────
counts_list <- list(
  Gage2019_AD_2019genome    = Gage2019_counts_AD_2019genome,
  Gage2019_noAD_2019genome  = Gage2019_counts_noAD_2019genome,
  Gage2019_2023genome_AD    = Gage2019_counts_2023genome_AD,
  Gage2019_2023genome_noAD  = Gage2019_counts_2023genome_noAD,
  Gage2019_2019genome_noAD  = Gage2019_counts_2019genome_noAD,
  Gage2024_counts               = Gage2024_counts,
  Gage2024_Mingchen         = Gage2024_Mingchen_counts,
  Gage2024_fibroblast       = Gage2024_fibroblast_counts
)

# ── Run normalization on all count sets ──────────────────────────────────────
normalized_list <- lapply(names(counts_list), function(name) {
  message("\n=== Processing: ", name, " ===")
  tryCatch(
    normalize_counts(
      counts_df      = counts_list[[name]],
      metadata       = Gage_merged_metadata,
      ranked_genes   = ranked_genes,
      stable_genes   = stable_genes,
      design_formula = "~ Diag + Age + Sex",
      dataset_name   = name
    ),
    error = function(e) {
      warning(name, " failed: ", e$message)
      return(NULL)
    }
  )
})

names(normalized_list) <- names(counts_list)

# Remove failed datasets
normalized_list <- Filter(Negate(is.null), normalized_list)
message("\nSuccessfully normalized ", length(normalized_list), " datasets")

# library(pROC)
# library(caret)
# 
# set.seed(42)
# 
# K <- 5
# folds <- createFolds(y, k = K, list = TRUE, returnTrain = FALSE)
# 
# per_gene_results <- lapply(gene_names, function(g) {
#   aucs <- numeric(K)
#   
#   for (i in seq_len(K)) {
#     test_idx  <- folds[[i]]
#     train_idx <- setdiff(seq_along(y), test_idx)
#     
#     X_train <- X[train_idx, g, drop = FALSE]
#     X_test  <- X[test_idx,  g, drop = FALSE]
#     y_train <- y[train_idx]
#     y_test  <- y[test_idx]
#     
#     df_train <- data.frame(y = y_train, x = X_train[, 1])
#     fit <- glm(y ~ x, data = df_train, family = binomial())
#     
#     df_test <- data.frame(x = X_test[, 1])
#     p_hat <- predict(fit, newdata = df_test, type = "response")
#     
#     roc_obj <- roc(response = y_test, predictor = p_hat, quiet = TRUE)
#     aucs[i] <- as.numeric(auc(roc_obj))
#   }
#   
#   data.frame(
#     gene = g,
#     mean_auc = mean(aucs),
#     sd_auc   = sd(aucs),
#     n_splits = K
#   )
# })
# 
# per_gene_df <- do.call(rbind, per_gene_results)
# per_gene_df <- per_gene_df[order(-per_gene_df$mean_auc), ]
# head(per_gene_df, 10)
# 
# # pipeline from perplexity.ai
# ## =========================================
# ## 0. Setup
# ## =========================================
# 
# # Packages
# library(pROC)
# library(glmnet)
# library(ggplot2)
# library(reshape2)
# 
# # Data: expr (samples x genes), meta_new (Diag: AD/CTRL) ----------------------
# # Ensure rownames(meta_new) match rownames(expr)
# stopifnot(all(rownames(expr) %in% rownames(Gage2024_metadata)))
# Gage2024_metadata <- Gage2024_metadata[rownames(expr), , drop = FALSE]
# 
# y <- ifelse(Gage2024_metadata$Diag == "AD", 1, 0)
# X <- as.matrix(expr)
# gene_names <- colnames(X)
# n <- length(y)
# 
# set.seed(2026)
# n_iter <- 100
# 
# # Storage ---------------------------------------------------------------------
# per_gene_auc_mat <- matrix(
#   NA_real_,
#   nrow = length(gene_names),
#   ncol = n_iter,
#   dimnames = list(gene_names, paste0("iter_", 1:n_iter))
# )
# 
# multi_auc_vec <- numeric(n_iter)
# 
# # Elastic net alpha -----------------------------------------------------------
# alpha_val <- 0.5  # 0 = ridge, 1 = lasso, 0.5 = elastic net
# 
# ## =========================================
# ## 1. Main loop: 100 iterations
# ## =========================================
# 
# for (it in 1:n_iter) {
#   cat("Iteration", it, "...\n")
#   
#   ## -----------------------------
#   ## 1A. Per-gene models with LOOCV
#   ## -----------------------------
#   iter_gene_auc <- numeric(length(gene_names))
#   names(iter_gene_auc) <- gene_names
#   
#   for (g in gene_names) {
#     pred_all <- numeric(n)
#     
#     # LOOCV
#     for (i in 1:n) {
#       test_idx  <- i
#       train_idx <- setdiff(1:n, i)
#       
#       X_train <- X[train_idx, g, drop = FALSE]
#       X_test  <- X[test_idx,  g, drop = FALSE]
#       y_train <- y[train_idx]
#       
#       df_train <- data.frame(y = y_train, x = X_train[, 1])
#       fit <- glm(y ~ x, data = df_train, family = binomial())
#       
#       df_test <- data.frame(x = X_test[, 1])
#       p_hat <- predict(fit, newdata = df_test, type = "response")
#       
#       pred_all[test_idx] <- p_hat
#     }
#     
#     roc_g <- roc(response = y, predictor = pred_all, quiet = TRUE)
#     iter_gene_auc[g] <- as.numeric(auc(roc_g))
#   }
#   
#   per_gene_auc_mat[, it] <- iter_gene_auc
#   
#   
#   ## -----------------------------
#   ## 1B. Multi-gene elastic net with LOOCV + inner CV
#   ## -----------------------------
#   pred_all_multi <- numeric(n)
#   
#   for (i in 1:n) {
#     test_idx  <- i
#     train_idx <- setdiff(1:n, i)
#     
#     X_train <- X[train_idx, , drop = FALSE]
#     X_test  <- X[test_idx,  , drop = FALSE]
#     y_train <- y[train_idx]
#     
#     # 0 = unpenalized (always included), 1 = penalized (default)
#     penalty_factors <- c(
#       rep(1, ncol(X)),    # genes - penalized
#       rep(0, ncol(cov_mat))  # Age, Sex - unpenalized
#     )
#     
#     # Inner CV for lambda on training data only
#     cvfit <- cv.glmnet(
#       x              = X_train,
#       y              = y_train,
#       family         = "binomial",
#       alpha          = alpha_val,
#       nfolds         = min(5, length(y_train)),
#       type.measure   = "auc",
#       penalty.factor = penalty_factors
#     )
#     
#     lambda_best <- cvfit$lambda.1se
#     
#     p_hat <- predict(
#       cvfit,
#       newx = X_test,
#       s = lambda_best,
#       type = "response"
#     )[1, 1]
#     
#     pred_all_multi[test_idx] <- p_hat
#   }
#   
#   roc_multi <- roc(response = y, predictor = pred_all_multi, quiet = TRUE)
#   multi_auc_vec[it] <- as.numeric(auc(roc_multi))
#   
#   cat("  Multi-gene elastic net LOOCV AUC:",
#       round(multi_auc_vec[it], 3), "\n")
# }
# 
# 
# ## =========================================
# ## 2. Save numeric outputs
# ## =========================================
# 
# # Per-gene AUC across iterations ----------------------------------------------
# per_gene_auc_df <- as.data.frame(per_gene_auc_mat)
# per_gene_auc_df$gene <- rownames(per_gene_auc_df)
# 
# per_gene_auc_long <- melt(
#   per_gene_auc_df,
#   id.vars = "gene",
#   variable.name = "iteration",
#   value.name = "auc"
# )
# 
# write.csv(
#   per_gene_auc_df,
#   file = file.path(output_dir, "per_gene_auc_100iters_LOOCV_wide.csv"),
#   row.names = FALSE
# )
# 
# write.csv(
#   per_gene_auc_long,
#   file = file.path(output_dir, "per_gene_auc_100iters_LOOCV_long.csv"),
#   row.names = FALSE
# )
# 
# # Multigene AUC across iterations ---------------------------------------------
# multi_auc_df <- data.frame(
#   iteration = 1:n_iter,
#   auc       = multi_auc_vec
# )
# 
# write.csv(
#   multi_auc_df,
#   file = file.path(output_dir, "multigene_elasticnet_auc_100iters_LOOCV.csv"),
#   row.names = FALSE
# )
# 
# 
# ## =========================================
# ## 3. Plots: per-gene vs multigene AUC
# ## =========================================
# 
# # 3A. Per-gene AUC distribution (boxplot) -------------------------------------
# p_gene_box <- ggplot(per_gene_auc_long, aes(x = gene, y = auc)) +
#   geom_boxplot(outlier.size = 0.5) +
#   theme_bw(base_size = 10) +
#   theme(
#     axis.text.x = element_text(
#       angle = 90, hjust = 1, vjust = 0.5
#     )
#   ) +
#   labs(
#     title = "Per-gene ROC AUC across 100 LOOCV iterations",
#     x = "Gene",
#     y = "AUC"
#   )
# 
# ggsave(
#   filename = file.path(output_dir, "per_gene_auc_boxplot_LOOCV.png"),
#   plot = p_gene_box,
#   width = 10,
#   height = 5,
#   dpi = 300
# )
# 
# # 3B. Mean per-gene AUC (barplot, sorted) -------------------------------------
# gene_mean_auc <- aggregate(
#   auc ~ gene,
#   data = per_gene_auc_long,
#   FUN = mean
# )
# gene_mean_auc <- gene_mean_auc[order(-gene_mean_auc$auc), ]
# 
# p_gene_bar <- ggplot(gene_mean_auc, aes(x = reorder(gene, auc), y = auc)) +
#   geom_col() +
#   coord_flip() +
#   theme_bw(base_size = 10) +
#   labs(
#     title = "Mean per-gene ROC AUC (100 LOOCV iterations)",
#     x = "Gene",
#     y = "Mean AUC"
#   )
# 
# ggsave(
#   filename = file.path(output_dir, "per_gene_auc_mean_barplot_LOOCV.png"),
#   plot = p_gene_bar,
#   width = 6,
#   height = 8,
#   dpi = 300
# )
# 
# # 3C. Multi-gene AUC distribution (histogram) ---------------------------------
# p_multi <- ggplot(multi_auc_df, aes(x = auc)) +
#   geom_histogram(color = "black", fill = "steelblue", bins = 15) +
#   theme_bw(base_size = 10) +
#   labs(
#     title = "Multi-gene elastic net ROC AUC (100 LOOCV iterations)",
#     x = "AUC",
#     y = "Count"
#   )
# 
# ggsave(
#   filename = file.path(output_dir, "multigene_elasticnet_auc_hist_LOOCV.png"),
#   plot = p_multi,
#   width = 5,
#   height = 4,
#   dpi = 300
# )
# 
# # 3D. Compare per-gene vs multi-gene (median line) ----------------------------
# multi_auc_median <- median(multi_auc_vec)
# 
# p_compare <- ggplot(gene_mean_auc, aes(x = auc)) +
#   geom_histogram(color = "black", fill = "grey70", bins = 15) +
#   geom_vline(xintercept = multi_auc_median, color = "red", size = 1.1) +
#   theme_bw(base_size = 10) +
#   labs(
#     title = "Mean per-gene AUC vs multi-gene elastic net\nRed line = median multi-gene AUC",
#     x = "AUC",
#     y = "Number of genes"
#   )
# 
# ggsave(
#   filename = file.path(output_dir, "per_gene_vs_multigene_auc_compare_LOOCV.png"),
#   plot = p_compare,
#   width = 5,
#   height = 4,
#   dpi = 300
# )
# 
# 
# ## =========================================
# ## 4. Final elastic-net model on full data
# ## =========================================
# 
# cvfit_full <- cv.glmnet(
#   x = X,
#   y = y,
#   family = "binomial",
#   alpha = alpha_val,
#   nfolds = min(5, n),
#   type.measure = "auc"
# )
# 
# lambda_best_full <- cvfit_full$lambda.min # Trying with less aggressive penalizer compared to lambda.1se
# 
# # Coefficients at chosen lambda
# coef_full <- coef(cvfit_full, s = lambda_best_full)
# coef_vec  <- as.numeric(coef_full)
# names(coef_vec) <- rownames(coef_full)
# 
# nz <- which(coef_vec != 0)
# nz_genes <- setdiff(names(coef_vec)[nz], "(Intercept)")
# 
# final_model_df <- data.frame(
#   feature = names(coef_vec),
#   coefficient = coef_vec
# )
# 
# write.csv(
#   final_model_df,
#   file = file.path(output_dir, "final_elasticnet_coefficients_full_data.csv"),
#   row.names = FALSE
# )
# 
# cat("Final elastic net fit on full data:\n")
# cat("  alpha:", alpha_val, "\n")
# cat("  lambda.1se:", lambda_best_full, "\n")
# cat("  Non-zero genes:", paste(nz_genes, collapse = ", "), "\n")

############ Version of above gene stability and prediction accuracy pipeline that loops over multiple counts matrices and performs bootstrapping per gene and for multi-gene models to get confidence intervals on AUC estimates 2-28-2026

##
## Helper functions
##
# LOOCV elastic net, returns vector of predictions for all samples
run_loocv_elastic_net <- function(X, y, meta_use, alpha_val = 0.5) {
  n <- length(y)
  pred_all <- numeric(n)
  
  # Create covariate matrix (Sex as numeric, Age as numeric)
  sex_numeric <- ifelse(toupper(meta_use$Sex) == "M", 1, 0)
  cov_mat <- cbind(
    Age = as.numeric(meta_use$Age),
    Sex = sex_numeric
  )
  
  # Bind covariates to gene expression
  X_full <- cbind(X, cov_mat)
  
  for (i in 1:n) {
    test_idx  <- i
    train_idx <- setdiff(1:n, i)
    
    X_train <- X_full[train_idx, , drop = FALSE]
    X_test  <- X_full[test_idx,  , drop = FALSE]
    y_train <- y[train_idx]
    
    # 0 = unpenalized (always included), 1 = penalized (default)
    penalty_factors <- c(
      rep(1, ncol(X)),    # genes - penalized
      rep(0, ncol(cov_mat))  # Age, Sex - unpenalized
    )
    
    cvfit <- cv.glmnet(
      x              = X_train,
      y              = y_train,
      family         = "binomial",
      alpha          = alpha_val,
      nfolds         = min(5, length(y_train)),
      type.measure   = "auc",
      penalty.factor = penalty_factors
    )
    
    p_hat <- predict(
      cvfit,
      newx = X_test,
      s    = cvfit$lambda.1se,
      type = "response"
    )[1, 1]
    
    pred_all[test_idx] <- p_hat
  }
  
  pred_all
}

run_loocv_per_gene <- function(X, y, metadata, outdir, group_label, today) {
  gene_names <- colnames(X)
  n          <- length(y)
  auc_vec    <- numeric(length(gene_names))
  names(auc_vec) <- gene_names
  
  for (g in gene_names) {
    pred_all <- numeric(n)
    
    # Skip genes with no variance
    if (var(X[, g]) == 0) {
      auc_vec[g] <- NA
      next
    }
    
    for (i in 1:n) {
      train_idx <- setdiff(1:n, i)
      
      df_train <- data.frame(
        y   = y[train_idx],
        x   = X[train_idx, g],
        sex = as.factor(metadata$Sex[train_idx]),
        age = as.numeric(metadata$Age[train_idx])
      )
      df_test <- data.frame(
        x   = X[i, g],
        sex = as.factor(metadata$Sex[i]),
        age = as.numeric(metadata$Age[i])
      )
      
      # Drop sex from model if only one level in training fold
      has_sex_variance <- length(unique(df_train$sex)) > 1
      formula_i <- if (has_sex_variance) y ~ x + sex + age else y ~ x + age
      
      fit <- tryCatch(
        glm(formula_i, data = df_train, family = binomial()),
        error   = function(e) NULL,
        warning = function(w) suppressWarnings(
          glm(formula_i, data = df_train, family = binomial())
        )
      )
      
      if (is.null(fit)) { pred_all[i] <- 0.5; next }
      
      pred_all[i] <- tryCatch(
        predict(fit, newdata = df_test, type = "response"),
        error = function(e) 0.5
      )
    }
    
    roc_g <- tryCatch(
      roc(response = y, predictor = pred_all, quiet = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(roc_g)) { auc_vec[g] <- NA; next }
    
    auc_vec[g] <- as.numeric(auc(roc_g))
    
    # Box plot with sex as point shape and age as point color
    plot_df <- data.frame(
      expression = as.numeric(X[, g]),
      group      = as.factor(metadata[["Diag"]]),
      sex        = as.factor(metadata$Sex),
      age        = as.numeric(metadata$Age),
      sample     = rownames(metadata)
    )
    
    p_box <- ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(
        aes(shape = sex, color = age),
        width = 0.2, size = 2.5, alpha = 0.9
      ) +
      scale_shape_manual(values = c("M" = 16, "F" = 17),
                         na.value = 15) +
      scale_color_gradient(low = "lightblue", high = "darkred",
                           name = "Age") +
      labs(
        title = paste0(g, " | AUC: ", round(auc_vec[g], 3),
                       " (adj. age + sex)"),
        x     = "Diagnosis",
        y     = "Expression",
        fill  = "Diagnosis",
        shape = "Sex"
      ) +
      theme_bw() +
      theme(
        plot.title      = element_text(face = "bold", size = 12),
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
    
    ggsave(
      filename = file.path(outdir, paste0(group_label, "_", g,
                                          "_expression_by_Diag_", today, ".pdf")),
      plot   = p_box,
      width  = 7,
      height = 5,
      device = cairo_pdf
    )
  }
  
  auc_vec
}

# Bootstrap CI for AUC, given labels y and predictions p_hat
bootstrap_auc_ci <- function(y, p_hat, B = 1000, conf = 0.95) {
  n        <- length(y)
  auc_boot <- numeric(B)
  
  idx_ad   <- which(y == 1)
  idx_ctrl <- which(y == 0)
  
  for (b in 1:B) {
    # Stratified resample to ensure both classes present
    idx <- c(
      sample(idx_ad,   size = length(idx_ad),   replace = TRUE),
      sample(idx_ctrl, size = length(idx_ctrl), replace = TRUE)
    )
    
    roc_b <- tryCatch(
      roc(response = y[idx], predictor = p_hat[idx], quiet = TRUE),
      error = function(e) NULL
    )
    
    auc_boot[b] <- if (!is.null(roc_b)) as.numeric(auc(roc_b)) else NA
  }
  
  auc_boot <- auc_boot[!is.na(auc_boot)]
  alpha    <- (1 - conf) / 2
  
  list(
    mean     = mean(auc_boot, na.rm = TRUE),
    ci_low   = as.numeric(quantile(auc_boot, probs = alpha)),
    ci_high  = as.numeric(quantile(auc_boot, probs = 1 - alpha)),
    auc_boot = auc_boot
  )
}

find_opt_threshold <- function(pred_probs, y, metric = "balanced_acc") {
  thresholds    <- seq(0.01, 0.99, by = 0.01)
  sensitivities <- numeric(length(thresholds))
  specificities <- numeric(length(thresholds))
  
  scores <- sapply(seq_along(thresholds), function(j) {
    thresh     <- thresholds[j]
    pred_class <- ifelse(pred_probs >= thresh, 1, 0)
    
    tp <- sum(pred_class == 1 & y == 1)
    fp <- sum(pred_class == 1 & y == 0)
    tn <- sum(pred_class == 0 & y == 0)
    fn <- sum(pred_class == 0 & y == 1)
    
    sensitivities[j] <<- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    specificities[j] <<- ifelse((tn + fp) == 0, 0, tn / (tn + fp))
    
    if (metric == "f1") {
      precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
      recall    <- sensitivities[j]
      ifelse((precision + recall) == 0, 0,
             2 * precision * recall / (precision + recall))
    } else if (metric == "balanced_acc") {
      (sensitivities[j] + specificities[j]) / 2
    } else if (metric == "mcc") {
      denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
      ifelse(denom == 0, 0, (tp * tn - fp * fn) / denom)
    } else {
      stop("Unknown metric. Choose from: 'f1', 'balanced_acc', 'mcc'")
    }
  })
  
  best_idx <- which.max(scores)
  
  list(
    threshold      = thresholds[best_idx],
    score          = scores[best_idx],
    balacc         = if (metric == "balanced_acc") scores[best_idx] else NA,
    f1             = if (metric == "f1")           scores[best_idx] else NA,
    mcc            = if (metric == "mcc")          scores[best_idx] else NA,
    sensitivity    = sensitivities[best_idx],
    specificity    = specificities[best_idx],
    metric         = metric,
    all_thresholds = thresholds,
    all_scores     = scores,
    all_sens       = sensitivities,
    all_spec       = specificities
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER: permutation_pvalue
# Computes empirical p-value for a gene's observed BalAcc against a null
# distribution obtained by permuting y B_perm times.
# Uses the same find_opt_threshold() logic for consistency.
# ═══════════════════════════════════════════════════════════════════════════════
permutation_pvalue <- function(pred_obs, y, B_perm = 1000, metric = "balanced_acc") {
  obs_score <- find_opt_threshold(pred_obs, y, metric = metric)$score
  
  null_scores <- vapply(seq_len(B_perm), function(b) {
    y_perm <- sample(y)
    find_opt_threshold(pred_obs, y_perm, metric = metric)$score
  }, numeric(1))
  
  # One-sided: proportion of permutations >= observed
  p_val <- (sum(null_scores >= obs_score) + 1) / (B_perm + 1)
  
  list(
    p_value    = p_val,
    obs_score  = obs_score,
    null_mean  = mean(null_scores),
    null_sd    = sd(null_scores),
    null_scores = null_scores
  )
}

# # Settings
# output_dir <- "/home/thilsabeck/Documents/Gage2024/RNApreprocessing_gene_validation_files/AgeSex_confounds_residuals_accuracyFit_LODO/"
# if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# alpha_val <- 0.5
# B_boot    <- 1000
# today     <- format(Sys.Date(), "%m-%d-%Y")
# 
# results_summary <- list()
# 
# for (dname in names(normalized_list)) {
#   cat("=== Dataset:", dname, "===\n")
#   
#   nd       <- normalized_list[[dname]]
#   X        <- nd$X
#   y        <- nd$y
#   meta_use <- nd$metadata
#   
#   # Map gene names to both ensembl_id and hgnc_symbol
#   if (nd$gene_names[1] %in% stable_genes$ensembl_id) {
#     # Input is ensembl — map to hgnc
#     ensembl_ids    <- nd$gene_names
#     hgnc_symbols   <- stable_genes$hgnc_symbol[match(nd$gene_names, stable_genes$ensembl_id)]
#   } else if (nd$gene_names[1] %in% stable_genes$hgnc_symbol) {
#     # Input is hgnc — map to ensembl
#     hgnc_symbols   <- nd$gene_names
#     ensembl_ids    <- stable_genes$ensembl_id[match(nd$gene_names, stable_genes$hgnc_symbol)]
#   } else {
#     warning(dname, ": Gene names do not match either ensembl_id or hgnc_symbol in stable_genes, keeping original names")
#     hgnc_symbols <- nd$gene_names
#     ensembl_ids  <- rep(NA, length(nd$gene_names))
#   }
#   
#   stopifnot(nrow(X) == length(y))
#   stopifnot(all(rownames(X) == rownames(meta_use)))
#   
#   # 1. LOOCV per-gene AUCs + bootstrap CIs
#   gene_names   <- colnames(X)
#   gene_symbols <- hgnc_symbols[match(gene_names, ensembl_ids)]
#   n            <- length(y)
#   auc_vec      <- numeric(length(gene_names))
#   ci_low_vec   <- numeric(length(gene_names))
#   ci_high_vec  <- numeric(length(gene_names))
#   names(auc_vec) <- names(ci_low_vec) <- names(ci_high_vec) <- gene_names
#   gene_n = 0
#   # LOOCV per-gene AUCs - corrected to control for sex and age
#   for (g in gene_names) {
#     gene_n <- gene_n + 1
#     pred_all <- numeric(n)
#     
#     # Skip genes with no variance
#     if (var(X[, g]) == 0) {
#       auc_vec[g] <- NA; ci_low_vec[g] <- NA; ci_high_vec[g] <- NA
#       next
#     }
#     
#     for (i in 1:n) {
#       train_idx <- setdiff(1:n, i)
#       
#       df_train <- data.frame(
#         y   = y[train_idx],
#         x   = X[train_idx, g],
#         sex = as.factor(meta_use$Sex[train_idx]),
#         age = as.numeric(meta_use$Age[train_idx])
#       )
#       df_test <- data.frame(
#         x   = X[i, g],
#         sex = as.factor(meta_use$Sex[i]),
#         age = as.numeric(meta_use$Age[i])
#       )
#       
#       # Check variance in training covariates
#       if (length(unique(df_train$sex)) < 2) {
#         # Only one sex in training — drop sex from model
#         fit <- tryCatch(
#           glm(y ~ x + age, data = df_train, family = binomial()),
#           error   = function(e) NULL,
#           warning = function(w) suppressWarnings(
#             glm(y ~ x + age, data = df_train, family = binomial())
#           )
#         )
#       } else {
#         fit <- tryCatch(
#           glm(y ~ x + sex + age, data = df_train, family = binomial()),
#           error   = function(e) NULL,
#           warning = function(w) suppressWarnings(
#             glm(y ~ x + sex + age, data = df_train, family = binomial())
#           )
#         )
#       }
#       
#       if (is.null(fit)) { pred_all[i] <- 0.5; next }
#       
#       pred_all[i] <- tryCatch(
#         predict(fit, newdata = df_test, type = "response"),
#         error = function(e) 0.5
#       )
#     }
#     
#     roc_g <- tryCatch(
#       roc(response = y, predictor = pred_all, quiet = TRUE),
#       error = function(e) NULL
#     )
#     
#     if (is.null(roc_g)) {
#       auc_vec[g] <- NA; ci_low_vec[g] <- NA; ci_high_vec[g] <- NA; next
#     }
#     
#     auc_vec[g]     <- as.numeric(auc(roc_g))
#     ci_g           <- bootstrap_auc_ci(y = y, p_hat = pred_all, B = B_boot, conf = 0.95)
#     ci_low_vec[g]  <- ci_g$ci_low
#     ci_high_vec[g] <- ci_g$ci_high
#     
#     # Box plot with age and sex annotation
#     plot_df <- data.frame(
#       expression = as.numeric(X[, g]),
#       group      = as.factor(meta_use[["Diag"]]),
#       sex        = meta_use$Sex,
#       age        = meta_use$Age,
#       sample     = rownames(meta_use)
#     )
#     
#     p_box <- ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
#       geom_boxplot(outlier.shape = NA, alpha = 0.6) +
#       geom_jitter(aes(shape = sex), width = 0.15, size = 2, alpha = 0.8) +
#       labs(
#         title = paste0(hgnc_symbols[gene_n], ":", g,
#                        " | AUC: ", round(auc_vec[g], 3),
#                        " [", round(ci_low_vec[g], 3),
#                        ", ", round(ci_high_vec[g], 3), "]",
#                        " (adjusted for age + sex)"),
#         x    = "Diagnosis",
#         y    = "Expression",
#         fill = "Diagnosis",
#         shape = "Sex"
#       ) +
#       theme_bw() +
#       theme(
#         plot.title      = element_text(face = "bold", size = 11),
#         axis.text.x     = element_text(angle = 45, hjust = 1)
#       )
#     
#     ggsave(
#       filename = file.path(output_dir, paste0(dname, "_", hgnc_symbols[gene_n], "_",g,
#                                               "_expression_by_Diag_", today, ".pdf")),
#       plot = p_box, width = 6, height = 5, device = cairo_pdf
#     )
#   }
#   
#   per_gene_df <- data.frame(
#     ensembl_id  = ensembl_ids,
#     gene_symbol = hgnc_symbols,
#     auc         = as.numeric(auc_vec),
#     ci_low      = as.numeric(ci_low_vec),
#     ci_high     = as.numeric(ci_high_vec)
#   )
#   
#   # Filter to meaningful AUCs and sort
#   per_gene_df_plot <- per_gene_df[!is.na(per_gene_df$auc), ]
#   per_gene_df_plot <- per_gene_df_plot[order(-per_gene_df_plot$auc), ]
#   
#   # Classify biomarker strength
#   per_gene_df_plot$strength <- cut(
#     per_gene_df_plot$auc,
#     breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
#     labels = c("Poor (<0.6)", "Fair (0.6-0.7)", "Good (0.7-0.8)",
#                "Strong (0.8-0.9)", "Excellent (>0.9)"),
#     include.lowest = TRUE
#   )
#   
#   # Top 20 biomarker dot plot with CI
#   top20 <- head(per_gene_df_plot, 20)
#   
#   p_biomarker <- ggplot(top20,
#                         aes(x = auc,
#                             y = reorder(gene_symbol, auc),
#                             color = strength)) +
#     geom_point(size = 4) +
#     geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3, linewidth = 0.7) +
#     geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
#     geom_vline(xintercept = 0.7, linetype = "dotted", color = "orange") +
#     geom_vline(xintercept = 0.8, linetype = "dotted", color = "red") +
#     scale_color_manual(values = c(
#       "Poor (<0.6)"      = "grey60",
#       "Fair (0.6-0.7)"   = "steelblue",
#       "Good (0.7-0.8)"   = "darkgreen",
#       "Strong (0.8-0.9)" = "orange",
#       "Excellent (>0.9)" = "red"
#     )) +
#     scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
#     labs(
#       title    = paste0("Top AD Biomarkers - ", dname,
#                         "\n(LOOCV AUC adjusted for age + sex)"),
#       x        = "AUC (95% bootstrap CI)",
#       y        = "Gene",
#       color    = "Biomarker Strength"
#     ) +
#     theme_bw(base_size = 11) +
#     theme(
#       plot.title   = element_text(face = "bold", size = 12),
#       legend.position = "bottom"
#     )
#   
#   ggsave(
#     file.path(output_dir, paste0("top_biomarkers_", dname, "_", today, ".pdf")),
#     p_biomarker, width = 8, height = 7, device = cairo_pdf
#   )
#   
#   # Also save a summary table of top biomarkers
#   write.csv(
#     head(per_gene_df_plot, 20),
#     file      = file.path(output_dir, paste0("top20_biomarkers_", dname, "_", today, ".csv")),
#     row.names = FALSE
#   )
#   
#   write.csv(
#     per_gene_df,
#     file      = file.path(output_dir, paste0("per_gene_LOOCV_auc_", dname, "_", today, ".csv")),
#     row.names = FALSE
#   )
#   
#   # 2. LOOCV elastic net (sex + age unpenalized, already handled in run_loocv_elastic_net)
#   pred_multi <- run_loocv_elastic_net(X, y, meta_use = meta_use, alpha_val = alpha_val)
#   roc_multi  <- roc(response = y, predictor = pred_multi, quiet = TRUE)
#   auc_multi  <- as.numeric(auc(roc_multi))
#   
#   # 3. Bootstrap CI for multigene AUC
#   boot_res <- bootstrap_auc_ci(y, pred_multi, B = B_boot, conf = 0.95)
#   
#   boot_df <- data.frame(
#     iteration = 1:B_boot,
#     auc       = boot_res$auc_boot
#   )
#   write.csv(
#     boot_df,
#     file      = file.path(output_dir, paste0("multigene_bootstrap_auc_", dname, "_", today, ".csv")),
#     row.names = FALSE
#   )
#   
#   # Extract elastic net coefficients (averaged across LOOCV folds for stability)
#   # Refit on full data to get representative coefficients for interpretation
#   sex_numeric <- ifelse(toupper(meta_use$Sex) == "M", 1, 0)
#   cov_mat     <- cbind(Age = as.numeric(meta_use$Age), Sex = sex_numeric)
#   X_full      <- cbind(X, cov_mat)
#   penalty_factors <- c(rep(1, ncol(X)), rep(0, ncol(cov_mat)))
#   
#   # Fit elastic net and lasso at both lambdas
#   cvfit_enet  <- cv.glmnet(
#     x              = X_full,
#     y              = y,
#     family         = "binomial",
#     alpha          = alpha_val,   # elastic net
#     nfolds         = min(5, length(y)),
#     type.measure   = "auc",
#     penalty.factor = penalty_factors
#   )
#   
#   cvfit_lasso <- cv.glmnet(
#     x              = X_full,
#     y              = y,
#     family         = "binomial",
#     alpha          = 1,           # lasso
#     nfolds         = min(5, length(y)),
#     type.measure   = "auc",
#     penalty.factor = penalty_factors
#   )
#   
#   # Report selection at both lambdas for both models
#   cat("\nElastic net (alpha=", alpha_val, "):\n")
#   cat("  lambda.min selected:", sum(coef(cvfit_enet, s = "lambda.min")[-1] != 0), "genes\n")
#   cat("  lambda.1se selected:", sum(coef(cvfit_enet, s = "lambda.1se")[-1] != 0), "genes\n")
#   cat("Lasso (alpha=1):\n")
#   cat("  lambda.min selected:", sum(coef(cvfit_lasso, s = "lambda.min")[-1] != 0), "genes\n")
#   cat("  lambda.1se selected:", sum(coef(cvfit_lasso, s = "lambda.1se")[-1] != 0), "genes\n")
#   # Save CV curve plots
#   pdf(file.path(output_dir, paste0("glmnet_CV_curves_", dname, "_", today, ".pdf")),
#       width = 10, height = 5)
#   par(mfrow = c(1, 2))
#   plot(cvfit_enet,  main = paste("Elastic net -", dname))
#   plot(cvfit_lasso, main = paste("Lasso -", dname))
#   dev.off()
#   
#   # Choose best model: prefer lambda.min, prefer lasso if enet selects 0
#   best_fit    <- if (sum(coef(cvfit_enet,  s = "lambda.min")[-1] != 0) > 0) cvfit_enet  else cvfit_lasso
#   best_lambda <- "lambda.min"
#   coef_full   <- coef(best_fit, s = best_lambda)
#   
#   coef_df <- data.frame(
#     feature = rownames(coef_full)[-1],
#     coef    = as.numeric(coef_full)[-1]
#   )
#   coef_df <- coef_df[coef_df$coef != 0, ]
#   coef_df <- coef_df[!coef_df$feature %in% c("Age", "Sex"), ]
#   coef_df$gene_symbol <- stable_genes$hgnc_symbol[match(coef_df$feature, stable_genes$ensembl_id)]
#   coef_df$gene_symbol[is.na(coef_df$gene_symbol)] <- coef_df$feature[is.na(coef_df$gene_symbol)]
#   coef_df <- coef_df[order(-abs(coef_df$coef)), ]
#   
#   cat("Final selected genes:", nrow(coef_df), "\n")
#   cat("Selected:", paste(coef_df$gene_symbol, collapse = ", "), "\n")
#   
#   # Collinearity plot for top 20 genes by single-gene AUC
#   top20_genes <- head(per_gene_df_sorted$ensembl_id[
#     per_gene_df_sorted$ensembl_id %in% colnames(X)], 20)
#   top20_symbols <- stable_genes$hgnc_symbol[match(top20_genes, stable_genes$ensembl_id)]
#   top20_symbols[is.na(top20_symbols)] <- top20_genes[is.na(top20_symbols)]
#   
#   if (length(top20_genes) > 1) {
#     X_top20 <- X[, top20_genes, drop = FALSE]
#     colnames(X_top20) <- top20_symbols
#     cor_mat <- cor(X_top20, use = "pairwise.complete.obs")
#     
#     # Heatmap
#     p_cor_heat <- ggplot(
#       reshape2::melt(cor_mat),
#       aes(x = Var1, y = Var2, fill = value)
#     ) +
#       geom_tile() +
#       geom_text(aes(label = round(value, 2)), size = 2.5) +
#       scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
#                            midpoint = 0, limits = c(-1, 1)) +
#       labs(
#         title = paste0("Top gene collinearity - ", dname),
#         x = NULL, y = NULL, fill = "Pearson r"
#       ) +
#       theme_bw(base_size = 9) +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1))
#     
#     ggsave(
#       file.path(output_dir, paste0("collinearity_heatmap_", dname, "_", today, ".pdf")),
#       p_cor_heat,
#       width  = max(6, length(top20_genes) * 0.5),
#       height = max(5, length(top20_genes) * 0.5),
#       device = cairo_pdf
#     )
#     
#     # Clustermap version using pheatmap for cleaner clustering
#     pheatmap::pheatmap(
#       cor_mat,
#       color            = colorRampPalette(c("steelblue", "white", "tomato"))(100),
#       breaks           = seq(-1, 1, length.out = 101),
#       display_numbers  = TRUE,
#       number_format    = "%.2f",
#       fontsize_number  = 7,
#       main             = paste0("Top gene correlation - ", dname),
#       filename         = file.path(output_dir,
#                                    paste0("collinearity_clustered_", dname, "_", today, ".pdf")),
#       width  = max(6, length(top20_genes) * 0.4),
#       height = max(5, length(top20_genes) * 0.4)
#     )
#     
#     # Also plot collinearity for elastic net selected genes if any
#     if (nrow(coef_df) > 1) {
#       sel_genes   <- coef_df$feature[coef_df$feature %in% colnames(X)]
#       sel_symbols <- coef_df$gene_symbol[coef_df$feature %in% colnames(X)]
#       X_sel       <- X[, sel_genes, drop = FALSE]
#       colnames(X_sel) <- sel_symbols
#       cor_sel     <- cor(X_sel, use = "pairwise.complete.obs")
#       
#       pheatmap::pheatmap(
#         cor_sel,
#         color            = colorRampPalette(c("steelblue", "white", "tomato"))(100),
#         breaks           = seq(-1, 1, length.out = 101),
#         display_numbers  = TRUE,
#         number_format    = "%.2f",
#         fontsize_number  = 8,
#         main             = paste0("Selected gene correlation - ", dname),
#         filename         = file.path(output_dir,
#                                      paste0("selected_gene_correlation_", dname, "_", today, ".pdf")),
#         width  = max(4, nrow(coef_df) * 0.5),
#         height = max(4, nrow(coef_df) * 0.5)
#       )
#     }
#   }
#   
#   
#   # Map ensembl to symbol for selected genes
#   coef_df$gene_symbol <- stable_genes$hgnc_symbol[match(coef_df$feature, stable_genes$ensembl_id)]
#   coef_df$gene_symbol[is.na(coef_df$gene_symbol)] <- coef_df$feature[is.na(coef_df$gene_symbol)]
#   coef_df <- coef_df[order(-abs(coef_df$coef)), ]
#   
#   write.csv(
#     coef_df,
#     file      = file.path(output_dir, paste0("elastic_net_coefficients_", dname, "_", today, ".csv")),
#     row.names = FALSE
#   )
#   
#   # 4A. Per-gene AUC barplot with CI error bars
#   per_gene_df_sorted <- per_gene_df[!is.na(per_gene_df$auc), ]
#   per_gene_df_sorted <- per_gene_df_sorted[order(-per_gene_df_sorted$auc), ]
#   
#   p_gene_bar <- ggplot(per_gene_df_sorted,
#                        aes(x = reorder(gene_symbol, auc), y = auc)) +
#     geom_col(fill = "steelblue", alpha = 0.8) +
#     geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
#                   width = 0.3, color = "black") +
#     geom_hline(yintercept = 0.7, linetype = "dotted", color = "orange") +
#     geom_hline(yintercept = 0.8, linetype = "dotted", color = "red") +
#     coord_flip() +
#     theme_bw(base_size = 10) +
#     labs(
#       title    = paste("Per-gene LOOCV AUC (adj. age + sex) -", dname),
#       x        = "Gene",
#       y        = "AUC (95% CI)"
#     )
#   
#   ggsave(
#     file.path(output_dir, paste0("per_gene_LOOCV_auc_bar_", dname, "_", today, ".png")),
#     p_gene_bar, width = 6, height = 8, dpi = 300
#   )
#   
#   # 4B. Multigene bootstrap AUC histogram with strength annotations
#   auc_strength <- cut(
#     auc_multi,
#     breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
#     labels = c("Poor", "Fair", "Good", "Strong", "Excellent"),
#     include.lowest = TRUE
#   )
#   
#   p_boot <- ggplot(boot_df, aes(x = auc)) +
#     geom_histogram(color = "black", fill = "steelblue", bins = 30) +
#     geom_vline(xintercept = auc_multi,        color = "red",    linewidth = 1.1) +
#     geom_vline(xintercept = boot_res$ci_low,  color = "red",    linetype = "dashed") +
#     geom_vline(xintercept = boot_res$ci_high, color = "red",    linetype = "dashed") +
#     geom_vline(xintercept = 0.7,              color = "orange", linetype = "dotted") +
#     geom_vline(xintercept = 0.8,              color = "darkred",linetype = "dotted") +
#     annotate("text", x = auc_multi, y = Inf,
#              label = paste0("AUC=", round(auc_multi, 3), " (", auc_strength, ")"),
#              vjust = 2, hjust = -0.1, color = "red", size = 3.5) +
#     theme_bw(base_size = 10) +
#     labs(
#       title = paste0("Multigene LOOCV AUC (adj. age + sex) - ", dname,
#                      "\nMean: ", round(boot_res$mean, 3),
#                      " CI [", round(boot_res$ci_low, 3),
#                      ", ", round(boot_res$ci_high, 3), "]"),
#       x = "AUC",
#       y = "Count"
#     )
#   
#   ggsave(
#     file.path(output_dir, paste0("multigene_bootstrap_auc_hist_", dname, "_", today, ".png")),
#     p_boot, width = 6, height = 4, dpi = 300
#   )
#   
#   # 4C. Per-gene vs multigene comparison
#   p_compare <- ggplot(per_gene_df_sorted, aes(x = auc)) +
#     geom_histogram(color = "black", fill = "grey70", bins = 15) +
#     geom_vline(xintercept = auc_multi,        color = "red",    linewidth = 1.1) +
#     geom_vline(xintercept = boot_res$ci_low,  color = "red",    linetype = "dashed") +
#     geom_vline(xintercept = boot_res$ci_high, color = "red",    linetype = "dashed") +
#     geom_vline(xintercept = 0.5,              color = "grey40", linetype = "dashed") +
#     annotate("text", x = auc_multi, y = Inf,
#              label = paste0("Multigene\nAUC=", round(auc_multi, 3)),
#              vjust = 2, hjust = -0.1, color = "red", size = 3) +
#     theme_bw(base_size = 10) +
#     labs(
#       title = paste("Per-gene vs multigene AUC (adj. age + sex) -", dname),
#       x     = "AUC",
#       y     = "Number of genes"
#     )
#   
#   ggsave(
#     file.path(output_dir, paste0("per_gene_vs_multigene_auc_", dname, "_", today, ".png")),
#     p_compare, width = 5, height = 4, dpi = 300
#   )
#   
#   # 4D. Elastic net selected gene coefficients — most interpretable biomarker plot
#   if (nrow(coef_df) > 0) {
#     coef_df$direction <- ifelse(coef_df$coef > 0, "Higher in AD", "Lower in AD")
#     
#     p_coef <- ggplot(coef_df,
#                      aes(x = coef,
#                          y = reorder(gene_symbol, abs(coef)),
#                          fill = direction)) +
#       geom_col(alpha = 0.85) +
#       geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
#       scale_fill_manual(values = c("Higher in AD" = "tomato",
#                                    "Lower in AD"  = "steelblue")) +
#       labs(
#         title = paste0("Elastic Net Selected Genes - ", dname,
#                        "\n(LOOCV AUC: ", round(auc_multi, 3),
#                        " [", round(boot_res$ci_low, 3),
#                        ", ", round(boot_res$ci_high, 3), "]",
#                        " | adj. age + sex)"),
#         x    = "Coefficient (effect size)",
#         y    = "Gene",
#         fill = "Direction"
#       ) +
#       theme_bw(base_size = 11) +
#       theme(
#         plot.title      = element_text(face = "bold", size = 11),
#         legend.position = "bottom"
#       )
#     
#     ggsave(
#       file.path(output_dir, paste0("elastic_net_coef_", dname, "_", today, ".pdf")),
#       p_coef, width = 7, height = max(4, nrow(coef_df) * 0.35), device = cairo_pdf
#     )
#   }
#   
#   # 5. Store summary
#   results_summary[[dname]] <- list(
#     per_gene_df      = per_gene_df,
#     auc_multi        = auc_multi,
#     boot_mean        = boot_res$mean,
#     boot_ci_low      = boot_res$ci_low,
#     boot_ci_high     = boot_res$ci_high,
#     elastic_net_genes = coef_df
#   )
#   
#   cat("Dataset", dname, "done.\n")
#   cat("  Multigene LOOCV AUC (adj. age + sex):", round(auc_multi, 3), "\n")
#   cat("  Bootstrap AUC mean:", round(boot_res$mean, 3),
#       "CI [", round(boot_res$ci_low, 3), ",", round(boot_res$ci_high, 3), "]\n")
#   cat("  Elastic net selected", nrow(coef_df), "genes\n")
# }
# 
# # Save overall summary
# summary_df <- do.call(rbind, lapply(names(results_summary), function(dname) {
#   s <- results_summary[[dname]]
#   data.frame(
#     dataset   = dname,
#     auc_multi = s$auc_multi,
#     boot_mean = s$boot_mean,
#     ci_low    = s$boot_ci_low,
#     ci_high   = s$boot_ci_high
#   )
# }))
# 
# write.csv(
#   summary_df,
#   file      = file.path(output_dir, paste0("multidataset_auc_summary_", today, ".csv")),
#   row.names = FALSE
# )

### Modified version of above to also run the analysis after removing sex-related genes and to save the single gene plots in a single pdf 03-02-2026
### Commenting out below run_analysis and replacing with run_analysis_with_glm.R to make more modular and easier to track version changes on 3-5-2026
source("run_analysis_with_glm.R")
# run_analysis <- function(X, y, meta_use, dname, label, ensembl_ids, hgnc_symbols,
#                          output_dir, today, alpha_val, B_boot) {
#   
#   gene_names   <- colnames(X)
#   gene_symbols <- hgnc_symbols[match(gene_names, ensembl_ids)]
#   n            <- length(y)
#   auc_vec      <- numeric(length(gene_names))
#   ci_low_vec   <- numeric(length(gene_names))
#   ci_high_vec  <- numeric(length(gene_names))
#   names(auc_vec) <- names(ci_low_vec) <- names(ci_high_vec) <- gene_names
#   
#   pred_store   <- list()
#   boxplot_list <- list()
#   gene_n       <- 0
#   
#   for (g in gene_names) {
#     gene_n   <- gene_n + 1
#     pred_all <- numeric(n)
#     
#     if (var(X[, g]) == 0) {
#       auc_vec[g] <- NA; ci_low_vec[g] <- NA; ci_high_vec[g] <- NA
#       next
#     }
#     
#     for (i in 1:n) {
#       train_idx <- setdiff(1:n, i)
#       
#       # Build design matrix: gene + age + sex (sex unpenalized via penalty.factor)
#       sex_train <- ifelse(toupper(meta_use$Sex[train_idx]) == "M", 1, 0)
#       sex_test  <- ifelse(toupper(meta_use$Sex[i])         == "M", 1, 0)
#       
#       X_train <- cbind(
#         gene = X[train_idx, g],
#         age  = as.numeric(meta_use$Age[train_idx]),
#         sex  = sex_train
#       )
#       X_test <- matrix(c(
#         X[i, g],
#         as.numeric(meta_use$Age[i]),
#         sex_test
#       ), nrow = 1, dimnames = list(NULL, c("gene", "age", "sex")))
#       
#       y_train <- y[train_idx]
#       
#       # Skip if zero variance in gene or outcome is single class
#       if (var(X_train[, "gene"]) == 0 || length(unique(y_train)) < 2) {
#         pred_all[i] <- 0.5; next
#       }
#       
#       # Penalize gene only; age + sex are unpenalized covariates
#       pf <- c(gene = 1, age = 0, sex = 0)
#       
#       fit_cv <- tryCatch(
#         cv.glmnet(
#           x              = X_train,
#           y              = y_train,
#           family         = "binomial",
#           alpha          = alpha_val,
#           penalty.factor = pf,
#           nfolds         = min(5, length(y_train)),
#           type.measure   = "class",   # optimise accuracy
#           standardize    = TRUE
#         ),
#         error = function(e) NULL
#       )
#       
#       if (is.null(fit_cv)) { pred_all[i] <- 0.5; next }
#       
#       pred_all[i] <- tryCatch(
#         as.numeric(predict(fit_cv, newx = X_test,
#                            s = "lambda.min", type = "response")),
#         error = function(e) 0.5
#       )
#     }
#     
#     # Store completed pred_all after inner loop
#     pred_store[[g]] <- pred_all
#     
#     roc_g <- tryCatch(
#       roc(response = y, predictor = pred_all, quiet = TRUE),
#       error = function(e) NULL
#     )
#     
#     if (is.null(roc_g)) {
#       auc_vec[g] <- NA; ci_low_vec[g] <- NA; ci_high_vec[g] <- NA; next
#     }
#     
#     auc_vec[g]     <- as.numeric(auc(roc_g))
#     ci_g           <- bootstrap_auc_ci(y = y, p_hat = pred_all, B = B_boot, conf = 0.95)
#     ci_low_vec[g]  <- ci_g$ci_low
#     ci_high_vec[g] <- ci_g$ci_high
#     
#     sym <- if (!is.na(gene_symbols[gene_n])) gene_symbols[gene_n] else g
#     
#     # ── Optimal threshold from Youden's J ────────────────────────────────
#     best_coords <- tryCatch(
#       coords(roc_g, "best", best.method = "youden",
#              ret = c("threshold", "sensitivity", "specificity")),
#       error = function(e) data.frame(threshold = 0.5, sensitivity = NA, specificity = NA)
#     )
#     opt_thresh <- if (nrow(best_coords) > 1) best_coords$threshold[1] else best_coords$threshold
#     opt_thresh <- if (is.na(opt_thresh) || length(opt_thresh) == 0) 0.5 else opt_thresh
#     
#     # ── Stats at optimal threshold ────────────────────────────────────────
#     sensitivity_opt  <- sum((pred_all >= opt_thresh) & (y == 1)) / max(1, sum(y == 1))
#     specificity_opt  <- sum((pred_all <  opt_thresh) & (y == 0)) / max(1, sum(y == 0))
#     ppv_opt          <- sum((pred_all >= opt_thresh) & (y == 1)) / max(1, sum(pred_all >= opt_thresh))
#     npv_opt          <- sum((pred_all <  opt_thresh) & (y == 0)) / max(1, sum(pred_all <  opt_thresh))
#     balanced_acc_opt <- (sensitivity_opt + specificity_opt) / 2
#     correct_opt      <- ifelse(pred_all >= opt_thresh, "AD", "CTRL") == ifelse(y == 1, "AD", "CTRL")
#     
#     # ── Stats at fixed 0.5 threshold ──────────────────────────────────────
#     sensitivity_05   <- sum((pred_all >= 0.5) & (y == 1)) / max(1, sum(y == 1))
#     specificity_05   <- sum((pred_all <  0.5) & (y == 0)) / max(1, sum(y == 0))
#     balanced_acc_05  <- (sensitivity_05 + specificity_05) / 2
#     correct_05       <- ifelse(pred_all >= 0.5, "AD", "CTRL") == ifelse(y == 1, "AD", "CTRL")
#     
#     # ── Adjusted expression residuals ─────────────────────────────────────
#     df_full <- data.frame(
#       x   = as.numeric(X[, g]),
#       sex = as.factor(meta_use$Sex),
#       age = as.numeric(meta_use$Age)
#     )
#     has_sex_variance_full <- length(unique(df_full$sex)) > 1
#     lm_cov <- tryCatch(
#       if (has_sex_variance_full) lm(x ~ sex + age, data = df_full)
#       else                       lm(x ~ age,       data = df_full),
#       error = function(e) NULL
#     )
#     adj_expression <- if (!is.null(lm_cov)) residuals(lm_cov) else df_full$x - mean(df_full$x)
#     
#     # ── Build plot data frame ─────────────────────────────────────────────
#     plot_df <- data.frame(
#       raw_expression = as.numeric(X[, g]),
#       adj_expression = as.numeric(adj_expression),
#       predicted_prob = pred_all,
#       actual_diag    = ifelse(y == 1, "AD", "CTRL"),
#       group          = as.factor(meta_use[["Diag"]]),
#       sex            = as.factor(meta_use$Sex),
#       age            = as.numeric(meta_use$Age),
#       sample         = rownames(meta_use),
#       stringsAsFactors = FALSE
#     )
#     plot_df$predicted_opt <- ifelse(plot_df$predicted_prob >= opt_thresh, "AD", "CTRL")
#     plot_df$predicted_05  <- ifelse(plot_df$predicted_prob >= 0.5,        "AD", "CTRL")
#     plot_df$correct_opt   <- plot_df$predicted_opt == plot_df$actual_diag
#     plot_df$correct_05    <- plot_df$predicted_05  == plot_df$actual_diag
#     
#     # ── Panel A: Raw expression ───────────────────────────────────────────
#     p_raw <- ggplot(plot_df, aes(x = group, y = raw_expression, fill = group)) +
#       geom_boxplot(outlier.shape = NA, alpha = 0.6) +
#       geom_jitter(aes(shape = sex, color = age), width = 0.2, size = 2.5, alpha = 0.9) +
#       scale_shape_manual(values = c("M" = 16, "F" = 17), na.value = 15) +
#       scale_color_gradient(low = "lightblue", high = "darkred", name = "Age") +
#       labs(title = "A: Raw expression",
#            x = "Diagnosis", y = "Expression",
#            fill = "Diagnosis", shape = "Sex") +
#       theme_bw(base_size = 10) +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1))
#     
#     # ── Panel B: Adjusted expression ─────────────────────────────────────
#     p_adj <- ggplot(plot_df, aes(x = group, y = adj_expression, fill = group)) +
#       geom_boxplot(outlier.shape = NA, alpha = 0.6) +
#       geom_jitter(aes(shape = sex, color = age), width = 0.2, size = 2.5, alpha = 0.9) +
#       scale_shape_manual(values = c("M" = 16, "F" = 17), na.value = 15) +
#       scale_color_gradient(low = "lightblue", high = "darkred", name = "Age") +
#       geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#       labs(title = "B: Adjusted expression (sex + age residuals)",
#            x = "Diagnosis", y = "Adjusted Expression",
#            fill = "Diagnosis", shape = "Sex") +
#       theme_bw(base_size = 10) +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1))
#     
#     # ── Panel C: Predicted probability per sample ─────────────────────────
#     p_pred <- ggplot(plot_df,
#                      aes(x     = reorder(sample, predicted_prob),
#                          y     = predicted_prob,
#                          color = actual_diag,
#                          shape = correct_opt)) +
#       annotate("rect", xmin = -Inf, xmax = Inf, ymin = opt_thresh, ymax = 1.0,
#                fill = "lightgreen", alpha = 0.08) +
#       annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.0, ymax = opt_thresh,
#                fill = "lightyellow", alpha = 0.08) +
#       geom_hline(yintercept = opt_thresh, linetype = "dashed",
#                  color = "darkgreen", linewidth = 0.8) +
#       geom_hline(yintercept = 0.5, linetype = "dotted",
#                  color = "grey40", linewidth = 0.5) +
#       geom_point(size = 3, alpha = 0.9) +
#       annotate("text", x = 1, y = opt_thresh + 0.03,
#                label = paste0("Optimal threshold: ", round(opt_thresh, 2)),
#                color = "darkgreen", size = 3, hjust = 0) +
#       scale_color_manual(values = c("AD" = "tomato", "CTRL" = "steelblue"),
#                          name = "Actual Diagnosis") +
#       scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
#                          name   = "Prediction (optimal threshold)",
#                          labels = c("TRUE" = "Correct", "FALSE" = "Wrong")) +
#       scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
#       labs(title = "C: LOOCV predicted P(AD) per sample",
#            x     = "Sample (ordered by predicted probability)",
#            y     = "Predicted P(AD)") +
#       theme_bw(base_size = 10) +
#       theme(
#         axis.text.x     = element_text(angle = 90, hjust = 1, size = 7),
#         plot.title      = element_text(face = "bold", size = 10),
#         legend.position = "bottom"
#       )
#     
#     # ── Panel D: Dual confusion matrix (optimal vs 0.5) ───────────────────
#     make_conf_df <- function(predicted, thresh_label) {
#       df <- data.frame(
#         actual_diag    = ifelse(y == 1, "AD", "CTRL"),
#         predicted_diag = predicted,
#         stringsAsFactors = FALSE
#       ) %>%
#         group_by(actual_diag, predicted_diag) %>%
#         summarize(n = n(), .groups = "drop") %>%
#         mutate(correct   = actual_diag == predicted_diag,
#                threshold = thresh_label)
#       
#       all_combos <- expand.grid(
#         actual_diag    = c("AD", "CTRL"),
#         predicted_diag = c("AD", "CTRL"),
#         threshold      = thresh_label,
#         stringsAsFactors = FALSE
#       )
#       df <- merge(all_combos, df, all.x = TRUE)
#       df$n[is.na(df$n)]           <- 0
#       df$correct[is.na(df$correct)] <- FALSE
#       df
#     }
#     
#     conf_both <- bind_rows(
#       make_conf_df(plot_df$predicted_opt,
#                    paste0("Optimal (", round(opt_thresh, 2), ")\n",
#                           "BalAcc=", round(balanced_acc_opt * 100, 1), "%")),
#       make_conf_df(plot_df$predicted_05,
#                    paste0("Fixed (0.5)\n",
#                           "BalAcc=", round(balanced_acc_05 * 100, 1), "%"))
#     )
#     
#     p_conf <- ggplot(conf_both,
#                      aes(x = predicted_diag, y = actual_diag,
#                          fill = correct, label = n)) +
#       geom_tile(alpha = 0.7) +
#       geom_text(size = 7, fontface = "bold") +
#       facet_wrap(~ threshold) +
#       scale_fill_manual(values = c("TRUE" = "lightgreen",
#                                    "FALSE" = "lightsalmon"),
#                         guide = "none") +
#       labs(title = "D: Confusion matrix",
#            x = "Predicted", y = "Actual") +
#       theme_bw(base_size = 10) +
#       theme(plot.title  = element_text(face = "bold"),
#             strip.text  = element_text(face = "bold", size = 8))
#     
#     # ── Combine all 4 panels ──────────────────────────────────────────────
#     p_combined <- cowplot::plot_grid(
#       p_raw, p_adj, p_pred, p_conf,
#       ncol = 2, rel_heights = c(1, 1.3)
#     )
#     
#     title_grob <- cowplot::ggdraw() +
#       cowplot::draw_label(
#         paste0(sym, " (", g, ")",
#                "  |  AUC: ",       round(auc_vec[g], 3),
#                " [", round(ci_low_vec[g], 3), ", ", round(ci_high_vec[g], 3), "]",
#                "  |  BalAcc(opt): ", round(balanced_acc_opt * 100, 1), "%",
#                "  |  BalAcc(0.5): ", round(balanced_acc_05  * 100, 1), "%",
#                "  |  Sens: ",        round(sensitivity_opt  * 100, 1), "%",
#                "  |  Spec: ",        round(specificity_opt  * 100, 1), "%",
#                "  |  ", label),
#         fontface = "bold", size = 8, x = 0.01, hjust = 0
#       )
#     
#     boxplot_list[[g]] <- cowplot::plot_grid(
#       title_grob, p_combined,
#       ncol = 1, rel_heights = c(0.04, 1)
#     )
#     
#     cat("  Gene", gene_n, "/", length(gene_names), ":", sym,
#         "| AUC:", round(auc_vec[g], 3),
#         "| BalAcc(opt):", round(balanced_acc_opt * 100, 1), "%\n")
#   }
#   
#   # ── Save all gene plots to one combined PDF ───────────────────────────────
#   combined_pdf_path <- file.path(output_dir,
#                                  paste0(label, "_all_gene_plots_", today, ".pdf"))
#   message("Saving combined gene plot PDF: ", combined_pdf_path)
#   # Order plots by decreasing balanced_acc_opt, then decreasing AUC
#   plot_order <- gene_summary_df$ensembl_id  # gene_summary_df already sorted by balanced_acc_opt
#   plot_order_valid <- plot_order[plot_order %in% names(boxplot_list)]
#   # Any genes in boxplot_list but not in gene_summary_df (e.g. NA AUC) go last
#   remainder <- setdiff(names(boxplot_list), plot_order_valid)
#   
#   pdf(combined_pdf_path, width = 14, height = 10)
#   for (g in c(plot_order_valid, remainder)) print(boxplot_list[[g]])
#   dev.off()
#   message("Saved ", length(boxplot_list), " gene plots to ", combined_pdf_path)
#   
#   # ── Build gene summary with optimal threshold stats ───────────────────────
#   gene_summary_rows <- lapply(gene_names, function(g) {
#     pred_all <- pred_store[[g]]
#     if (is.null(pred_all) || is.na(auc_vec[g])) return(NULL)
#     
#     roc_g <- tryCatch(
#       roc(response = y, predictor = pred_all, quiet = TRUE),
#       error = function(e) NULL
#     )
#     if (is.null(roc_g)) return(NULL)
#     
#     best_coords <- tryCatch(
#       coords(roc_g, "best", best.method = "youden",
#              ret = c("threshold", "sensitivity", "specificity")),
#       error = function(e) data.frame(threshold = 0.5, sensitivity = NA, specificity = NA)
#     )
#     opt_thresh <- if (nrow(best_coords) > 1) best_coords$threshold[1] else best_coords$threshold
#     opt_thresh <- if (is.na(opt_thresh) || length(opt_thresh) == 0) 0.5 else opt_thresh
#     
#     sensitivity_opt  <- sum((pred_all >= opt_thresh) & (y == 1)) / max(1, sum(y == 1))
#     specificity_opt  <- sum((pred_all <  opt_thresh) & (y == 0)) / max(1, sum(y == 0))
#     ppv_opt          <- sum((pred_all >= opt_thresh) & (y == 1)) / max(1, sum(pred_all >= opt_thresh))
#     npv_opt          <- sum((pred_all <  opt_thresh) & (y == 0)) / max(1, sum(pred_all <  opt_thresh))
#     balanced_acc_opt <- (sensitivity_opt + specificity_opt) / 2
#     correct_opt      <- ifelse(pred_all >= opt_thresh, "AD", "CTRL") == ifelse(y == 1, "AD", "CTRL")
#     
#     sensitivity_05   <- sum((pred_all >= 0.5) & (y == 1)) / max(1, sum(y == 1))
#     specificity_05   <- sum((pred_all <  0.5) & (y == 0)) / max(1, sum(y == 0))
#     balanced_acc_05  <- (sensitivity_05 + specificity_05) / 2
#     correct_05       <- ifelse(pred_all >= 0.5, "AD", "CTRL") == ifelse(y == 1, "AD", "CTRL")
#     
#     sym <- hgnc_symbols[match(g, ensembl_ids)]
#     sym <- if (!is.na(sym)) sym else g
#     
#     data.frame(
#       ensembl_id       = g,
#       gene_symbol      = sym,
#       auc              = round(auc_vec[g],          3),
#       ci_low           = round(ci_low_vec[g],       3),
#       ci_high          = round(ci_high_vec[g],      3),
#       opt_threshold    = round(opt_thresh,          3),
#       sensitivity_opt  = round(sensitivity_opt,     3),
#       specificity_opt  = round(specificity_opt,     3),
#       ppv_opt          = round(ppv_opt,             3),
#       npv_opt          = round(npv_opt,             3),
#       balanced_acc_opt = round(balanced_acc_opt,    3),
#       accuracy_opt     = round(mean(correct_opt),   3),
#       n_correct_opt    = sum(correct_opt),
#       sensitivity_05   = round(sensitivity_05,      3),
#       specificity_05   = round(specificity_05,      3),
#       balanced_acc_05  = round(balanced_acc_05,     3),
#       accuracy_05      = round(mean(correct_05),    3),
#       n_correct_05     = sum(correct_05),
#       n_total          = length(y),
#       n_AD             = sum(y == 1),
#       n_CTRL           = sum(y == 0),
#       stringsAsFactors = FALSE
#     )
#   })
#   
#   gene_summary_df <- do.call(rbind, Filter(Negate(is.null), gene_summary_rows))
#   gene_summary_df <- gene_summary_df[order(-gene_summary_df$balanced_acc_opt), ]
#   
#   write.csv(
#     gene_summary_df,
#     file.path(output_dir, paste0("per_gene_summary_", label, "_", today, ".csv")),
#     row.names = FALSE
#   )
#   
#   # ── Build per_gene_df ─────────────────────────────────────────────────────
#   # Pull per-gene stats directly from gene_summary_df which is already correctly computed
#   # Do NOT rebuild from auc_vec/ci_low_vec/ci_high_vec — those vectors are indexed by
#   # gene name position and can misalign after merge with ensembl_ids/hgnc_symbols vectors.
#   per_gene_df <- gene_summary_df[, c("ensembl_id", "gene_symbol", "auc", "ci_low", "ci_high",
#                                      "balanced_acc_opt", "balanced_acc_05", "opt_threshold",
#                                      "sensitivity_opt", "specificity_opt")]
#   
#   # Sort by balanced accuracy
#   per_gene_df_sorted <- per_gene_df[!is.na(per_gene_df$auc), ]
#   per_gene_df_sorted <- per_gene_df_sorted[order(-per_gene_df_sorted$balanced_acc_opt), ]
#   stopifnot(all(per_gene_df_sorted$auc == per_gene_df_sorted$auc[
#     match(per_gene_df_sorted$ensembl_id, per_gene_df_sorted$ensembl_id)]))
#   per_gene_df_sorted$strength <- cut(
#     per_gene_df_sorted$auc,
#     breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
#     labels = c("Poor (<0.6)", "Fair (0.6-0.7)", "Good (0.7-0.8)",
#                "Strong (0.8-0.9)", "Excellent (>0.9)"),
#     include.lowest = TRUE
#   )
#   
#   # ── Top biomarker dot plot: AUC + balanced accuracy ───────────────────────
#   top20 <- head(per_gene_df_sorted, 20)
#   
#   p_biomarker <- ggplot(top20, aes(y = reorder(gene_symbol, balanced_acc_opt))) +
#     geom_point(aes(x = auc, color = "AUC"), size = 4) +
#     geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, color = "AUC"), height = 0.25, linewidth = 0.7) +
#     geom_point(aes(x = balanced_acc_opt, color = "Balanced Accuracy (optimal thresh)"),
#                size = 4, shape = 18) +
#     geom_point(aes(x = balanced_acc_05,  color = "Balanced Accuracy (0.5 thresh)"),
#                size = 3, shape = 15) +
#     geom_vline(xintercept = 0.5, linetype = "dashed",  color = "grey50") +
#     geom_vline(xintercept = 0.7, linetype = "dotted",  color = "orange") +
#     geom_vline(xintercept = 0.8, linetype = "dotted",  color = "red") +
#     scale_color_manual(values = c(
#       "AUC"                                  = "steelblue",
#       "Balanced Accuracy (optimal thresh)"   = "tomato",
#       "Balanced Accuracy (0.5 thresh)"       = "salmon"
#     )) +
#     scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
#     labs(
#       title = paste0("Top AD Biomarkers - ", label,
#                      "\n(circles=AUC±CI | diamonds=BalAcc optimal | squares=BalAcc@0.5)"),
#       x     = "Metric value",
#       y     = "Gene (sorted by balanced accuracy)",
#       color = "Metric"
#     ) +
#     theme_bw(base_size = 11) +
#     theme(plot.title      = element_text(face = "bold", size = 10),
#           legend.position = "bottom")
#   
#   ggsave(file.path(output_dir, paste0("top_biomarkers_", label, "_", today, ".pdf")),
#          p_biomarker, width = 9, height = 7, device = cairo_pdf)
#   
#   write.csv(head(per_gene_df_sorted, 20),
#             file.path(output_dir, paste0("top20_biomarkers_", label, "_", today, ".csv")),
#             row.names = FALSE)
#   write.csv(per_gene_df,
#             file.path(output_dir, paste0("per_gene_LOOCV_auc_", label, "_", today, ".csv")),
#             row.names = FALSE)
#   
#   # ── Elastic net ───────────────────────────────────────────────────────────
#   pred_multi <- run_loocv_elastic_net(X, y, meta_use = meta_use, alpha_val = alpha_val)
#   roc_multi  <- roc(response = y, predictor = pred_multi, quiet = TRUE)
#   auc_multi  <- as.numeric(auc(roc_multi))
#   boot_res   <- bootstrap_auc_ci(y, pred_multi, B = B_boot, conf = 0.95)
#   
#   boot_df <- data.frame(iteration = 1:B_boot, auc = boot_res$auc_boot)
#   write.csv(boot_df,
#             file.path(output_dir, paste0("multigene_bootstrap_auc_", label, "_", today, ".csv")),
#             row.names = FALSE)
#   
#   # ── Elastic net balanced accuracy stats (needed for all downstream plots) ──
#   roc_enet_obj <- tryCatch(
#     roc(response = y, predictor = pred_multi, quiet = TRUE),
#     error = function(e) NULL
#   )
#   
#   if (!is.null(roc_enet_obj)) {
#     bc_enet <- tryCatch(
#       coords(roc_enet_obj, "best", best.method = "youden",
#              ret = c("threshold", "sensitivity", "specificity")),
#       error = function(e) data.frame(threshold = 0.5, sensitivity = NA, specificity = NA)
#     )
#     enet_opt_thresh <- if (nrow(bc_enet) > 1) bc_enet$threshold[1] else bc_enet$threshold
#     enet_opt_thresh <- if (is.na(enet_opt_thresh) || length(enet_opt_thresh) == 0) 0.5 else enet_opt_thresh
#   } else {
#     enet_opt_thresh <- 0.5
#   }
#   
#   enet_sens_opt   <- sum((pred_multi >= enet_opt_thresh) & (y == 1)) / max(1, sum(y == 1))
#   enet_spec_opt   <- sum((pred_multi <  enet_opt_thresh) & (y == 0)) / max(1, sum(y == 0))
#   enet_balacc_opt <- (enet_sens_opt + enet_spec_opt) / 2
#   enet_sens_05    <- sum((pred_multi >= 0.5) & (y == 1)) / max(1, sum(y == 1))
#   enet_spec_05    <- sum((pred_multi <  0.5) & (y == 0)) / max(1, sum(y == 0))
#   enet_balacc_05  <- (enet_sens_05 + enet_spec_05) / 2
#   
#   sex_numeric     <- ifelse(toupper(meta_use$Sex) == "M", 1, 0)
#   cov_mat         <- cbind(Age = as.numeric(meta_use$Age), Sex = sex_numeric)
#   X_full          <- cbind(X, cov_mat)
#   penalty_factors <- c(rep(1, ncol(X)), rep(0, ncol(cov_mat)))
#   
#   cvfit_enet <- cv.glmnet(
#     x = X_full, y = y, family = "binomial", alpha = alpha_val,
#     nfolds = min(5, length(y)), type.measure = "class",   # optimise misclassification → accuracy
#     penalty.factor = penalty_factors
#   )
#   cvfit_lasso <- cv.glmnet(
#     x = X_full, y = y, family = "binomial", alpha = 1,
#     nfolds = min(5, length(y)), type.measure = "class",
#     penalty.factor = penalty_factors
#   )
#   
#   cat("\n", label, "- Elastic net (alpha=", alpha_val, "):\n")
#   cat("  lambda.min:", sum(coef(cvfit_enet,  s = "lambda.min")[-1] != 0), "genes\n")
#   cat("  lambda.1se:", sum(coef(cvfit_enet,  s = "lambda.1se")[-1] != 0), "genes\n")
#   cat("Lasso:\n")
#   cat("  lambda.min:", sum(coef(cvfit_lasso, s = "lambda.min")[-1] != 0), "genes\n")
#   cat("  lambda.1se:", sum(coef(cvfit_lasso, s = "lambda.1se")[-1] != 0), "genes\n")
#   
#   pdf(file.path(output_dir, paste0("glmnet_CV_curves_", label, "_", today, ".pdf")),
#       width = 10, height = 5)
#   par(mfrow = c(1, 2))
#   plot(cvfit_enet,  main = paste("Elastic net -", label))
#   plot(cvfit_lasso, main = paste("Lasso -", label))
#   dev.off()
#   
#   best_fit  <- if (sum(coef(cvfit_enet, s = "lambda.min")[-1] != 0) > 0) cvfit_enet else cvfit_lasso
#   coef_full <- coef(best_fit, s = "lambda.min")
#   coef_df   <- data.frame(
#     feature = rownames(coef_full)[-1],
#     coef    = as.numeric(coef_full)[-1]
#   )
#   coef_df <- coef_df[coef_df$coef != 0, ]
#   coef_df <- coef_df[!coef_df$feature %in% c("Age", "Sex"), ]
#   coef_df$gene_symbol <- stable_genes$hgnc_symbol[match(coef_df$feature, stable_genes$ensembl_id)]
#   coef_df$gene_symbol[is.na(coef_df$gene_symbol)] <- coef_df$feature[is.na(coef_df$gene_symbol)]
#   coef_df <- coef_df[order(-abs(coef_df$coef)), ]
#   
#   cat("Final selected genes:", nrow(coef_df), "\n")
#   cat("Selected:", paste(coef_df$gene_symbol, collapse = ", "), "\n")
#   write.csv(coef_df,
#             file.path(output_dir, paste0("elastic_net_coefficients_", label, "_", today, ".csv")),
#             row.names = FALSE)
#   
#   # ── Collinearity plots ────────────────────────────────────────────────────
#   top20_genes <- head(per_gene_df_sorted$ensembl_id[
#     per_gene_df_sorted$ensembl_id %in% colnames(X)], 20)
#   top20_sym <- stable_genes$hgnc_symbol[match(top20_genes, stable_genes$ensembl_id)]
#   top20_sym[is.na(top20_sym)] <- top20_genes[is.na(top20_sym)]
#   
#   if (length(top20_genes) > 1) {
#     X_top20 <- X[, top20_genes, drop = FALSE]
#     colnames(X_top20) <- top20_sym
#     cor_mat <- cor(X_top20, use = "pairwise.complete.obs")
#     
#     pheatmap::pheatmap(
#       cor_mat,
#       color           = colorRampPalette(c("steelblue", "white", "tomato"))(100),
#       breaks          = seq(-1, 1, length.out = 101),
#       display_numbers = TRUE, number_format = "%.2f", fontsize_number = 7,
#       main            = paste0("Top gene correlation - ", label),
#       filename        = file.path(output_dir,
#                                   paste0("collinearity_clustered_", label, "_", today, ".pdf")),
#       width  = max(6, length(top20_genes) * 0.4),
#       height = max(5, length(top20_genes) * 0.4)
#     )
#     
#     if (nrow(coef_df) > 1) {
#       sel_genes <- coef_df$feature[coef_df$feature %in% colnames(X)]
#       sel_sym   <- coef_df$gene_symbol[coef_df$feature %in% colnames(X)]
#       X_sel     <- X[, sel_genes, drop = FALSE]
#       colnames(X_sel) <- sel_sym
#       cor_sel   <- cor(X_sel, use = "pairwise.complete.obs")
#       pheatmap::pheatmap(
#         cor_sel,
#         color           = colorRampPalette(c("steelblue", "white", "tomato"))(100),
#         breaks          = seq(-1, 1, length.out = 101),
#         display_numbers = TRUE, number_format = "%.2f", fontsize_number = 8,
#         main            = paste0("Selected gene correlation - ", label),
#         filename        = file.path(output_dir,
#                                     paste0("selected_gene_correlation_", label, "_", today, ".pdf")),
#         width  = max(4, nrow(coef_df) * 0.5),
#         height = max(4, nrow(coef_df) * 0.5)
#       )
#     }
#   }
#   
#   # ── Summary plots ─────────────────────────────────────────────────────────
#   p_gene_bar <- ggplot(per_gene_df_sorted,
#                        aes(x = reorder(gene_symbol, balanced_acc_opt), y = balanced_acc_opt)) +
#     geom_col(fill = "tomato", alpha = 0.8) +
#     geom_point(aes(y = auc), color = "steelblue", size = 2) +
#     geom_hline(yintercept = 0.7, linetype = "dotted", color = "orange") +
#     geom_hline(yintercept = 0.8, linetype = "dotted", color = "red") +
#     coord_flip() +
#     theme_bw(base_size = 10) +
#     labs(title = paste("Per-gene balanced accuracy (bars) + AUC (dots) -", label),
#          subtitle = "Sorted by balanced accuracy at optimal threshold",
#          x = "Gene", y = "Metric value")
#   
#   ggsave(file.path(output_dir, paste0("per_gene_LOOCV_auc_bar_", label, "_", today, ".png")),
#          p_gene_bar, width = 6, height = 8, dpi = 300)
#   
#   auc_strength <- cut(auc_multi, breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
#                       labels = c("Poor", "Fair", "Good", "Strong", "Excellent"),
#                       include.lowest = TRUE)
#   
#   p_boot <- ggplot(boot_df, aes(x = auc)) +
#     geom_histogram(color = "black", fill = "steelblue", bins = 30) +
#     geom_vline(xintercept = auc_multi,        color = "red",     linewidth = 1.1) +
#     geom_vline(xintercept = boot_res$ci_low,  color = "red",     linetype = "dashed") +
#     geom_vline(xintercept = boot_res$ci_high, color = "red",     linetype = "dashed") +
#     geom_vline(xintercept = enet_balacc_opt,  color = "tomato",  linetype = "dashed", linewidth = 1.0) +
#     geom_vline(xintercept = enet_balacc_05,   color = "salmon",  linetype = "dashed", linewidth = 0.8) +
#     geom_vline(xintercept = 0.7,              color = "orange",  linetype = "dotted") +
#     geom_vline(xintercept = 0.8,              color = "darkred", linetype = "dotted") +
#     annotate("text", x = auc_multi, y = Inf,
#              label = paste0("AUC=", round(auc_multi, 3), " (", auc_strength, ")"),
#              vjust = 2, hjust = -0.1, color = "red", size = 3.5) +
#     annotate("text", x = enet_balacc_opt, y = Inf,
#              label = paste0("BalAcc(opt thresh=", round(enet_opt_thresh, 2), ")=",
#                             round(enet_balacc_opt, 3)),
#              vjust = 4, hjust = -0.05, color = "tomato", size = 3.2) +
#     annotate("text", x = enet_balacc_05, y = Inf,
#              label = paste0("BalAcc(@0.5)=", round(enet_balacc_05, 3)),
#              vjust = 6, hjust = -0.05, color = "salmon", size = 3.2) +
#     theme_bw(base_size = 10) +
#     labs(title = paste0("Multigene LOOCV — ", label,
#                         "\nAUC: ", round(boot_res$mean, 3),
#                         " CI [", round(boot_res$ci_low, 3), ", ", round(boot_res$ci_high, 3), "]",
#                         "  |  BalAcc(opt): ", round(enet_balacc_opt * 100, 1), "%",
#                         "  |  Sens: ",        round(enet_sens_opt   * 100, 1), "%",
#                         "  |  Spec: ",        round(enet_spec_opt   * 100, 1), "%",
#                         "  |  BalAcc(0.5): ", round(enet_balacc_05  * 100, 1), "%"),
#          x = "Bootstrap AUC", y = "Count")
#   ggsave(file.path(output_dir, paste0("multigene_bootstrap_auc_hist_", label, "_", today, ".png")),
#          p_boot, width = 6, height = 4, dpi = 300)
#   
#   # ── Structured comparison: elastic net vs top per-gene results ────────────
#   # Compute balanced accuracy for the multigene elastic net predictions
#   enet_opt_thresh  <- tryCatch({
#     roc_enet  <- roc(response = y, predictor = pred_multi, quiet = TRUE)
#     bc        <- coords(roc_enet, "best", best.method = "youden",
#                         ret = c("threshold", "sensitivity", "specificity"))
#     if (nrow(bc) > 1) bc$threshold[1] else bc$threshold
#   }, error = function(e) 0.5)
#   enet_opt_thresh  <- if (is.na(enet_opt_thresh) || length(enet_opt_thresh) == 0) 0.5 else enet_opt_thresh
#   
#   # Build comparison table: top 10 per-gene + elastic net row
#   top10_compare <- head(gene_summary_df, 10)[, c("gene_symbol", "auc", "ci_low", "ci_high",
#                                                  "balanced_acc_opt", "balanced_acc_05",
#                                                  "sensitivity_opt", "specificity_opt")]
#   top10_compare$model_type <- "Single gene (LOOCV logistic)"
#   
#   enet_row <- data.frame(
#     gene_symbol      = paste0("Elastic net (", nrow(coef_df), " genes)"),
#     auc              = round(auc_multi,            3),
#     ci_low           = round(enet_boot_res$ci_low, 3),   # was boot_res$ci_low
#     ci_high          = round(enet_boot_res$ci_high,3),   # was boot_res$ci_high
#     balanced_acc_opt = round(enet_balacc_opt,      3),
#     balanced_acc_05  = round(enet_balacc_05,       3),
#     sensitivity_opt  = round(enet_sens_opt,        3),
#     specificity_opt  = round(enet_spec_opt,        3),
#     model_type       = "Elastic net (LOOCV)",
#     stringsAsFactors = FALSE
#   )
#   
#   comparison_df <- rbind(top10_compare, enet_row)
#   write.csv(comparison_df,
#             file.path(output_dir, paste0("comparison_singlegene_vs_enet_", label, "_", today, ".csv")),
#             row.names = FALSE)
#   
#   # Compute enet bootstrap CI for the comparison plot
#   enet_boot_res <- bootstrap_auc_ci(y = y, p_hat = pred_multi, B = B_boot, conf = 0.95)
#   
#   comparison_df$is_enet <- comparison_df$model_type == "Elastic net (LOOCV)"
#   
#   p_compare <- ggplot(comparison_df,
#                       aes(y     = reorder(gene_symbol, balanced_acc_opt),
#                           color = model_type)) +
#     # AUC point + CI for all rows
#     geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
#                    height = 0.3, linewidth = 0.7) +
#     geom_point(aes(x = auc, shape = "AUC"), size = 4) +
#     # Balanced accuracy optimal threshold — filled diamond (shape 23)
#     geom_point(aes(x = balanced_acc_opt, shape = "BalAcc (optimal thresh)"),
#                size = 4, alpha = 0.85) +
#     # Balanced accuracy at 0.5 — filled square (shape 22)
#     geom_point(aes(x = balanced_acc_05, shape = "BalAcc (0.5 thresh)"),
#                size = 3, alpha = 0.7) +
#     scale_shape_manual(
#       name   = "Metric",
#       values = c(
#         "AUC"                    = 16,   # circle
#         "BalAcc (optimal thresh)"= 23,   # filled diamond
#         "BalAcc (0.5 thresh)"    = 22    # filled square
#       )
#     ) +
#     scale_color_manual(
#       name   = "Model type",
#       values = c(
#         "Single gene (LOOCV logistic)" = "steelblue",
#         "Elastic net (LOOCV)"          = "tomato"
#       )
#     ) +
#     geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
#     geom_vline(xintercept = 0.7, linetype = "dotted", color = "orange") +
#     geom_vline(xintercept = 0.8, linetype = "dotted", color = "red") +
#     scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
#     labs(
#       title    = paste0("Single-gene vs elastic net: AUC & balanced accuracy — ", label),
#       subtitle = paste0(
#         "Circles+CI=AUC | Diamonds=BalAcc(optimal) | Squares=BalAcc(@0.5)\n",
#         "Enet: AUC=",          round(auc_multi, 3),
#         " [",                  round(enet_boot_res$ci_low, 3), ", ",
#         round(enet_boot_res$ci_high, 3), "]",
#         "  BalAcc(opt)=",      round(enet_balacc_opt * 100, 1), "%",
#         "  Sens=",             round(enet_sens_opt   * 100, 1), "%",
#         "  Spec=",             round(enet_spec_opt   * 100, 1), "%"
#       ),
#       x = "Metric value", y = "Model"
#     ) +
#     theme_bw(base_size = 11) +
#     theme(
#       plot.title    = element_text(face = "bold", size = 10),
#       plot.subtitle = element_text(size = 8, color = "grey30"),
#       legend.position = "bottom"
#     )
#   
#   ggsave(file.path(output_dir, paste0("comparison_singlegene_vs_enet_", label, "_", today, ".pdf")),
#          p_compare,
#          width  = 9,
#          height = max(5, nrow(comparison_df) * 0.45),
#          device = cairo_pdf)
#   # ── Bootstrap CIs for elastic net coefficients ───────────────────────────
#   set.seed(42)
#   coef_boot_list <- lapply(1:B_boot, function(b) {
#     idx_b   <- sample(1:length(y), replace = TRUE)
#     y_b     <- y[idx_b]
#     X_full_b <- X_full[idx_b, , drop = FALSE]
#     
#     # Skip degenerate bootstrap samples
#     if (length(unique(y_b)) < 2) return(NULL)
#     
#     fit_b <- tryCatch(
#       cv.glmnet(x = X_full_b, y = y_b, family = "binomial",
#                 alpha = alpha_val, nfolds = min(5, length(y_b)),
#                 type.measure = "class", penalty.factor = penalty_factors,
#                 standardize = TRUE),
#       error = function(e) NULL
#     )
#     if (is.null(fit_b)) return(NULL)
#     
#     cf <- coef(fit_b, s = "lambda.min")
#     data.frame(
#       feature = rownames(cf)[-1],
#       coef    = as.numeric(cf)[-1],
#       stringsAsFactors = FALSE
#     )
#   })
#   
#   coef_boot_df <- do.call(rbind, Filter(Negate(is.null), coef_boot_list))
#   
#   # Summarise: mean, CI, selection frequency per feature
#   coef_boot_summary <- coef_boot_df %>%
#     group_by(feature) %>%
#     summarise(
#       coef_mean      = mean(coef,        na.rm = TRUE),
#       coef_ci_low    = quantile(coef, 0.025, na.rm = TRUE),
#       coef_ci_high   = quantile(coef, 0.975, na.rm = TRUE),
#       selection_freq = mean(coef != 0,   na.rm = TRUE),
#       .groups = "drop"
#     )
#   
#   # Merge CIs into coef_df (which has the point estimates from full-data fit)
#   coef_df <- merge(coef_df, coef_boot_summary, by = "feature", all.x = TRUE)
#   coef_df <- coef_df[order(-abs(coef_df$coef)), ]
#   
#   write.csv(coef_df,
#             file.path(output_dir, paste0("elastic_net_coefficients_", label, "_", today, ".csv")),
#             row.names = FALSE)
#   
#   if (nrow(coef_df) > 0) {
#     coef_df$direction <- ifelse(coef_df$coef > 0, "Higher in AD", "Lower in AD")
#     
#     p_coef <- ggplot(coef_df,
#                      aes(x = coef, y = reorder(gene_symbol, abs(coef)),
#                          fill = direction)) +
#       geom_col(alpha = 0.6, width = 0.6) +
#       # Bootstrap 95% CI error bars
#       geom_errorbarh(aes(xmin = coef_ci_low, xmax = coef_ci_high),
#                      height = 0.3, linewidth = 0.7, color = "grey30") +
#       geom_point(aes(x = coef, color = direction), size = 2.5) +
#       # Selection frequency as text label
#       geom_text(aes(x = ifelse(coef >= 0, coef_ci_high, coef_ci_low),
#                     label = paste0(round(selection_freq * 100), "%")),
#                 hjust = ifelse(coef_df$coef >= 0, -0.15, 1.15),
#                 size = 3, color = "grey30") +
#       geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
#       scale_fill_manual(values  = c("Higher in AD" = "tomato",  "Lower in AD" = "steelblue")) +
#       scale_color_manual(values = c("Higher in AD" = "tomato",  "Lower in AD" = "steelblue"),
#                          guide  = "none") +
#       labs(
#         title    = paste0("Elastic Net Selected Genes — ", label),
#         subtitle = paste0(
#           "AUC: ",        round(auc_multi, 3),
#           " [",           round(boot_res$ci_low, 3), ", ", round(boot_res$ci_high, 3), "]",
#           "  BalAcc(opt)=", round(enet_balacc_opt * 100, 1), "%",
#           "  Sens=",      round(enet_sens_opt * 100, 1), "%",
#           "  Spec=",      round(enet_spec_opt * 100, 1), "%\n",
#           "Bars=bootstrap 95% CI | % label=selection frequency across ", B_boot, " bootstraps"
#         ),
#         x    = "Coefficient (full-data fit)",
#         y    = "Gene",
#         fill = "Direction"
#       ) +
#       theme_bw(base_size = 11) +
#       theme(
#         plot.title    = element_text(face = "bold", size = 10),
#         plot.subtitle = element_text(size = 8, color = "grey30"),
#         legend.position = "bottom"
#       ) +
#       # Extra x margin to fit selection frequency labels
#       scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))
#     
#     ggsave(file.path(output_dir, paste0("elastic_net_coef_", label, "_", today, ".pdf")),
#            p_coef,
#            width  = 7,
#            height = max(4, nrow(coef_df) * 0.45),
#            device = cairo_pdf)
#   }
#   
#   cat("Dataset", label, "done.\n")
#   cat("  Multigene LOOCV AUC:", round(auc_multi, 3), "\n")
#   cat("  Bootstrap mean:", round(boot_res$mean, 3),
#       "CI [", round(boot_res$ci_low, 3), ",", round(boot_res$ci_high, 3), "]\n")
#   cat("  Elastic net selected", nrow(coef_df), "genes\n")
#   
#   list(
#     per_gene_df       = per_gene_df,
#     gene_summary_df   = gene_summary_df,
#     auc_multi         = auc_multi,
#     boot_mean         = boot_res$mean,
#     boot_ci_low       = boot_res$ci_low,
#     boot_ci_high      = boot_res$ci_high,
#     elastic_net_genes = coef_df
#   )
# }

# ── Run ───────────────────────────────────────────────────────────────────────
# ── Configuration ─────────────────────────────────────────────────────────────
output_dir <- "/home/thilsabeck/Documents/Gage2024/RNApreprocessing_gene_validation_files/AgeSex_confounds_residuals_accuracyFit_LODO_permed/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

alpha_val <- 0.5
B_boot    <- 1000
B_perm <- 1000   # permutations per gene — increase to 5000 for publication
today     <- format(Sys.Date(), "%m-%d-%Y")

# ── Chromosome annotation map ──────────────────────────────────────────────────
# Built from stable_genes if a chromosome column is available.
# Used to: (1) annotate genes in summaries, (2) identify sex-chr genes for
# the "nosexgenes" run without relying on a hardcoded symbol regex.
if ("chromosome" %in% names(stable_genes) || "chr" %in% names(stable_genes)) {
  chr_col     <- if ("chromosome" %in% names(stable_genes)) "chromosome" else "chr"
  gene_chr_map <- setNames(stable_genes[[chr_col]], stable_genes$ensembl_id)
  cat("Chromosome map built from stable_genes:", length(gene_chr_map), "entries\n")
} else {
  library(biomaRt)
  mart <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
  chr_df <- getBM(attributes=c('ensembl_gene_id','chromosome_name'),
                  filters='ensembl_gene_id', values=stable_genes$ensembl_id, mart=mart)
  stable_genes$chromosome <- chr_df$chromosome_name[ match(stable_genes$ensembl_id, chr_df$ensembl_gene_id)]
  gene_chr_map <- setNames(stable_genes$chromosome, stable_genes$ensembl_id)
  cat("No chromosome column found in stable_genes — chr annotation will be NA\n",
      "Consider adding it via biomaRt:\n",
      "  library(biomaRt)\n",
      "  mart <- useMart('ensembl', dataset='hsapiens_gene_ensembl')\n",
      "  chr_df <- getBM(attributes=c('ensembl_gene_id','chromosome_name'),\n",
      "                  filters='ensembl_gene_id', values=stable_genes$ensembl_id, mart=mart)\n",
      "  stable_genes$chromosome <- chr_df$chromosome_name[",
      "match(stable_genes$ensembl_id, chr_df$ensembl_gene_id)]\n")
}

# Sex-chromosome gene identification:
# Primary: use chromosome annotation from gene_chr_map (X/Y)
# Fallback: regex on symbol names (covers cases where chr annotation is missing)
sex_chr_symbol_regex <- "XIST|DDX3Y|USP9Y|RPS4Y1|ZFY|TSIX|EIF1AY|KDM5D|NLGN4Y|TMSB4Y|UTY|PRKY|PCDH11Y"

is_sex_chr_gene <- function(ensembl_ids_vec, hgnc_symbols_vec, gene_chr_map) {
  if (!is.null(gene_chr_map)) {
    chrs <- gene_chr_map[ensembl_ids_vec]
    chrs[is.na(chrs)] <- ""
    chr_flag <- chrs %in% c("X", "chrX", "Y", "chrY")
  } else {
    chr_flag <- rep(FALSE, length(ensembl_ids_vec))
  }
  # Union with regex fallback for any genes not in the map
  sym_flag <- grepl(sex_chr_symbol_regex, hgnc_symbols_vec, ignore.case = FALSE)
  chr_flag | sym_flag
}

# ── RF tuning grids ────────────────────────────────────────────────────────────
node_grid    <- c(1, 3, 5, 10)
alpha_grid   <- c(0, 0.25, 0.5, 0.75, 1.0)
nfolds_inner <- 5

# ── Per-dataset loop ───────────────────────────────────────────────────────────
results_summary  <- list()
dataset_list_lodo <- list()   # populated for LODO after loop

for (dname in names(normalized_list)) {
  cat("\n══════════════════════════════════════\n")
  cat("Dataset:", dname, "\n")
  cat("══════════════════════════════════════\n")
  
  nd       <- normalized_list[[dname]]
  X        <- nd$X
  y        <- nd$y
  meta_use <- nd$metadata
  
  # ── Resolve gene ID type ─────────────────────────────────────────────────
  if (nd$gene_names[1] %in% stable_genes$ensembl_id) {
    ensembl_ids  <- nd$gene_names
    hgnc_symbols <- stable_genes$hgnc_symbol[match(nd$gene_names, stable_genes$ensembl_id)]
  } else if (nd$gene_names[1] %in% stable_genes$hgnc_symbol) {
    hgnc_symbols <- nd$gene_names
    ensembl_ids  <- stable_genes$ensembl_id[match(nd$gene_names, stable_genes$hgnc_symbol)]
  } else {
    warning(dname, ": Gene names unmatched, keeping original")
    hgnc_symbols <- nd$gene_names
    ensembl_ids  <- rep(NA_character_, length(nd$gene_names))
  }
  
  stopifnot(nrow(X) == length(y))
  stopifnot(all(rownames(X) == rownames(meta_use)))
  
  # ── Sex-chromosome mask ───────────────────────────────────────────────────
  sex_mask <- is_sex_chr_gene(
    ensembl_ids_vec  = ensembl_ids,
    hgnc_symbols_vec = hgnc_symbols,
    gene_chr_map     = gene_chr_map
  )
  X_nosex   <- X[, !sex_mask, drop = FALSE]
  ens_nosex <- ensembl_ids[!sex_mask]
  hgn_nosex <- hgnc_symbols[!sex_mask]
  cat("  Sex-chr genes identified:", sum(sex_mask),
      " (", paste(hgnc_symbols[sex_mask], collapse = ", "), ")\n")
  
  # Store for LODO (use all-genes version)
  dataset_list_lodo[[dname]] <- list(X = X, y = y, meta_use = meta_use)
  
  # ════════════════════════════════════════════════════════════════════════════
  # Run 1: glmnet + glm — all genes
  # ════════════════════════════════════════════════════════════════════════════
  cat("\n── Run 1: glmnet + glm | all genes ──\n")
  results_summary[[paste0(dname, "_allgenes")]] <- run_analysis(
    X            = X,
    y            = y,
    meta_use     = meta_use,
    dname        = dname,
    label        = paste0(dname, "_allgenes"),
    ensembl_ids  = ensembl_ids,
    hgnc_symbols = hgnc_symbols,
    output_dir   = output_dir,
    today        = today,
    alpha_val    = alpha_val,
    B_boot       = B_boot,
    B_perm = B_perm,
    gene_chr_map = gene_chr_map
  )
  
  # ════════════════════════════════════════════════════════════════════════════
  # Run 2: glmnet + glm — sex-chromosome genes removed
  # ════════════════════════════════════════════════════════════════════════════
  if (ncol(X_nosex) == 0) {
    warning(dname, ": No genes remaining after sex-chr removal — skipping Run 2")
  } else {
    cat("\n── Run 2: glmnet + glm | sex-chr genes removed (", sum(sex_mask), "removed) ──\n")
    results_summary[[paste0(dname, "_nosexgenes")]] <- run_analysis(
      X            = X_nosex,
      y            = y,
      meta_use     = meta_use,
      dname        = dname,
      label        = paste0(dname, "_nosexgenes"),
      ensembl_ids  = ens_nosex,
      hgnc_symbols = hgn_nosex,
      output_dir   = output_dir,
      today        = today,
      alpha_val    = alpha_val,
      B_boot       = B_boot,
      B_perm = B_perm,
      gene_chr_map = gene_chr_map
    )
  }
  
  # ════════════════════════════════════════════════════════════════════════════
  # Run 3: RF + glmnet + glm — all genes
  # ════════════════════════════════════════════════════════════════════════════
  cat("\n── Run 3: RF + glmnet + glm | all genes ──\n")
  results_summary[[paste0(dname, "_rf_allgenes")]] <- run_analysis_rf(
    X            = X,
    y            = y,
    meta_use     = meta_use,
    dname        = dname,
    label        = paste0(dname, "_rf_allgenes"),
    ensembl_ids  = ensembl_ids,
    hgnc_symbols = hgnc_symbols,
    output_dir   = output_dir,
    today        = today,
    alpha_val    = alpha_val,
    B_boot       = B_boot,
    B_perm = B_perm,
    gene_chr_map = gene_chr_map,
    num.trees    = 500,
    node_grid    = node_grid,
    alpha_grid   = alpha_grid,
    nfolds_inner = nfolds_inner
  )
  
  # ════════════════════════════════════════════════════════════════════════════
  # Run 4: RF + glmnet + glm — sex-chromosome genes removed
  # ════════════════════════════════════════════════════════════════════════════
  if (ncol(X_nosex) == 0) {
    warning(dname, ": No genes remaining after sex-chr removal — skipping Run 4")
  } else {
    cat("\n── Run 4: RF + glmnet + glm | sex-chr genes removed (", sum(sex_mask), "removed) ──\n")
    results_summary[[paste0(dname, "_rf_nosexgenes")]] <- run_analysis_rf(
      X            = X_nosex,
      y            = y,
      meta_use     = meta_use,
      dname        = dname,
      label        = paste0(dname, "_rf_nosexgenes"),
      ensembl_ids  = ens_nosex,
      hgnc_symbols = hgn_nosex,
      output_dir   = output_dir,
      today        = today,
      alpha_val    = alpha_val,
      B_boot       = B_boot,
      B_perm = B_perm,
      gene_chr_map = gene_chr_map,
      num.trees    = 500,
      node_grid    = node_grid,
      alpha_grid   = alpha_grid,
      nfolds_inner = nfolds_inner
    )
  }
  
  cat("\nDataset", dname, "complete.\n")
}

# ── LODO cross-dataset validation ──────────────────────────────────────────────
# Only meaningful with ≥ 3 datasets. Skipped silently if fewer are present.
if (length(dataset_list_lodo) >= 3) {
  cat("\n══════════════════════════════════════\n")
  cat("Running LODO cross-dataset validation\n")
  cat("Datasets:", paste(names(dataset_list_lodo), collapse = ", "), "\n")
  cat("══════════════════════════════════════\n")
  
  lodo_results <- run_lodo_validation(
    dataset_list = dataset_list_lodo,
    ensembl_ids  = ensembl_ids,    # from last iteration — all datasets share same gene set
    hgnc_symbols = hgnc_symbols,
    output_dir   = output_dir,
    today        = today,
    alpha_val    = alpha_val,
    B_boot       = B_boot,
    gene_chr_map = gene_chr_map,
    num.trees    = 500,
    node_grid    = node_grid,
    alpha_grid   = alpha_grid,
    nfolds_inner = nfolds_inner
  )
  results_summary[["lodo"]] <- lodo_results
  
} else {
  cat("\nSkipping LODO: fewer than 3 datasets available (",
      length(dataset_list_lodo), ")\n")
  lodo_results <- NULL
}

# ── Cross-dataset summary table ────────────────────────────────────────────────
# Flatten results_summary into a single comparison row per run
cross_summary_rows <- lapply(names(results_summary), function(run_name) {
  res <- results_summary[[run_name]]
  if (is.null(res) || run_name == "lodo") return(NULL)
  
  # Handle both run_analysis and run_analysis_rf return structures
  is_rf <- !is.null(res$auc_multi_rf)
  
  data.frame(
    run             = run_name,
    model_type      = if (is_rf) "RF+glmnet+glm" else "glmnet+glm",
    auc_multi       = round(if (is_rf) res$auc_multi_rf    else res$auc_multi,    3),
    boot_ci_low     = round(if (is_rf) res$boot_ci_low_rf  else res$boot_ci_low,  3),
    boot_ci_high    = round(if (is_rf) res$boot_ci_high_rf else res$boot_ci_high, 3),
    n_genes_input   = if (!is.null(res$per_gene_df)) nrow(res$per_gene_df) else NA_integer_,
    n_genes_elastic = if (is_rf) nrow(res$rf_importance)
    else       nrow(res$elastic_net_genes),
    stringsAsFactors = FALSE
  )
})
cross_summary_df <- do.call(rbind, Filter(Negate(is.null), cross_summary_rows))

write.csv(cross_summary_df,
          file.path(output_dir, paste0("cross_run_summary_", today, ".csv")),
          row.names = FALSE)

cat("\n══════════════════════════════════════\n")
cat("Pipeline complete.\n")
cat("Runs completed:", nrow(cross_summary_df), "\n")
cat("Output directory:", output_dir, "\n")
print(cross_summary_df)

# ── Cross-dataset per-gene summary ───────────────────────────────────────────

# ══════════════════════════════════════════════════════════════════════════════
# downstream_analysis.R
#
# Cross-run comparisons, gene-level summaries, and visualizations.
# Run after the main pipeline loop (pipeline_main.R).
#
# Primary metric throughout: BalAcc(opt)
# Secondary metric (shown in outputs, not used for sorting): AUC
#
# Handles both run_analysis (glmnet+glm) and run_analysis_rf (RF+glmnet+glm)
# return structures transparently via helper extract_run_metrics().
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(tidyr)

# ── Helper: extract standardised scalar metrics from any result object ─────────
extract_run_metrics <- function(run_name, res) {
  if (is.null(res)) return(NULL)
  is_rf <- !is.null(res$auc_multi_rf)
  
  data.frame(
    run              = run_name,
    dataset          = sub("_rf_allgenes|_rf_nosexgenes|_allgenes|_nosexgenes", "", run_name),
    sex_removed      = grepl("nosexgenes", run_name),
    model_type       = if (is_rf) "RF+glmnet+glm" else "glmnet+glm",
    # Primary: BalAcc
    balacc_multi     = round(if (is_rf) res$balacc_multi_rf  else res$balacc_multi,    3),
    sens_multi       = round(if (is_rf) res$sens_multi_rf    else res$sens_multi,      3),
    spec_multi       = round(if (is_rf) res$spec_multi_rf    else res$spec_multi,      3),
    null_balacc      = round(res$null_balacc,                                          3),
    balacc_gain      = round((if (is_rf) res$balacc_multi_rf else res$balacc_multi) -
                               res$null_balacc,                                        3),
    # Secondary: AUC (bootstrap CI)
    auc_multi        = round(if (is_rf) res$auc_multi_rf     else res$auc_multi,       3),
    boot_mean        = round(if (is_rf) res$boot_mean_rf     else res$boot_mean,       3),
    boot_ci_low      = round(if (is_rf) res$boot_ci_low_rf   else res$boot_ci_low,     3),
    boot_ci_high     = round(if (is_rf) res$boot_ci_high_rf  else res$boot_ci_high,    3),
    # Gene counts
    n_genes_input    = if (!is.null(res$per_gene_df)) nrow(res$per_gene_df) else NA_integer_,
    n_genes_selected = if (is_rf) nrow(res$rf_importance)
    else       nrow(res$elastic_net_genes),
    stringsAsFactors = FALSE
  )
}

# ── Helper: extract per-gene rows from any result object ──────────────────────
# Returns a unified long-format data frame regardless of run type.
# Columns present in both run types are always included;
# RF-specific columns (balacc_opt_rf etc.) are NA for glmnet-only runs.
extract_gene_rows <- function(run_name, res) {
  if (is.null(res) || is.null(res$per_gene_df)) return(NULL)
  df     <- res$per_gene_df
  is_rf  <- !is.null(res$auc_multi_rf)
  
  # Unified primary BalAcc column:
  # glmnet runs:  balanced_acc_opt (glmnet single-gene)
  # RF runs:      mean_balacc      (mean of RF/glmnet/glm per gene)
  if (is_rf) {
    df$primary_balacc  <- df$mean_balacc
    df$primary_auc     <- if ("auc_rf" %in% names(df)) df$auc_rf else NA_real_
    df$primary_sens    <- if ("sensitivity_rf" %in% names(df)) df$sensitivity_rf else NA_real_
    df$primary_spec    <- if ("specificity_rf" %in% names(df)) df$specificity_rf else NA_real_
  } else {
    df$primary_balacc  <- df$balanced_acc_opt
    df$primary_auc     <- df$auc
    df$primary_sens    <- df$sensitivity_opt
    df$primary_spec    <- df$specificity_opt
  }
  
  df$run         <- run_name
  df$dataset     <- sub("_rf_allgenes|_rf_nosexgenes|_allgenes|_nosexgenes", "", run_name)
  df$sex_removed <- grepl("nosexgenes", run_name)
  df$model_type  <- if (is_rf) "RF+glmnet+glm" else "glmnet+glm"
  df$null_balacc <- res$null_balacc
  df
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Cross-run scalar summary
# ══════════════════════════════════════════════════════════════════════════════
run_metrics_rows <- lapply(names(results_summary), function(nm) {
  if (nm == "lodo") return(NULL)
  extract_run_metrics(nm, results_summary[[nm]])
})
cross_run_df <- do.call(rbind, Filter(Negate(is.null), run_metrics_rows))
rownames(cross_run_df) <- NULL

write.csv(cross_run_df,
          file.path(output_dir, paste0("cross_run_summary_", today, ".csv")),
          row.names = FALSE)

cat("\n══════════════════════════════════════\n")
cat("Cross-run summary\n")
cat("══════════════════════════════════════\n")
print(cross_run_df[, c("run", "model_type", "balacc_multi", "sens_multi",
                       "spec_multi", "balacc_gain", "auc_multi",
                       "n_genes_input", "n_genes_selected")])

# ══════════════════════════════════════════════════════════════════════════════
# 2. Cross-run scalar comparison plot
#    One row per run; primary = BalAcc(opt); secondary = AUC
# ══════════════════════════════════════════════════════════════════════════════
cross_run_long <- tidyr::pivot_longer(
  cross_run_df,
  cols      = c(balacc_multi, auc_multi, null_balacc),
  names_to  = "metric",
  values_to = "value"
)
cross_run_long$metric <- factor(cross_run_long$metric,
                                levels = c("balacc_multi", "auc_multi", "null_balacc"),
                                labels = c("BalAcc(opt)", "AUC", "Null BalAcc"))
cross_run_long$run_short <- gsub(paste(unique(cross_run_long$dataset),
                                       collapse = "|"), "", cross_run_long$run)
cross_run_long$run_short <- gsub("^_+|_+$", "", cross_run_long$run_short)

p_cross_run <- ggplot(cross_run_df,
                      aes(y = reorder(run, balacc_multi), color = model_type)) +
  # AUC bootstrap CI bar (secondary — lighter)
  geom_errorbarh(aes(xmin = boot_ci_low, xmax = boot_ci_high),
                 height = 0.25, linewidth = 0.5, alpha = 0.5) +
  # AUC point (open circle)
  geom_point(aes(x = auc_multi), size = 3, shape = 1, stroke = 1.0) +
  # BalAcc point (filled circle — primary)
  geom_point(aes(x = balacc_multi), size = 4, shape = 16) +
  # Sens (upward triangle) and Spec (downward triangle)
  geom_point(aes(x = sens_multi), size = 2.5, shape = 24, fill = "white") +
  geom_point(aes(x = spec_multi), size = 2.5, shape = 25, fill = "white") +
  # Null BalAcc tick per run
  geom_point(aes(x = null_balacc), size = 2, shape = 3, color = "grey50") +
  geom_vline(xintercept = 0.5, linetype = "solid",  color = "grey80", linewidth = 0.4) +
  geom_vline(xintercept = 0.7, linetype = "dotted", color = "orange", linewidth = 0.6) +
  geom_vline(xintercept = 0.8, linetype = "dotted", color = "red",    linewidth = 0.6) +
  scale_color_manual(values = c("RF+glmnet+glm" = "#e41a1c",
                                "glmnet+glm"    = "#377eb8"),
                     name = "Model type") +
  scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
  annotate("text", x = 0.41, y = 0.5,
           label = "● BalAcc(opt)   ○ AUC+CI   △ Sens   ▽ Spec   + Null BalAcc",
           hjust = 0, size = 2.8, color = "grey35") +
  labs(
    title    = "Multi-gene model performance — all runs",
    subtitle = "Primary: BalAcc(opt) ●   Secondary: AUC ○   Sorted by BalAcc",
    x        = "Metric value",
    y        = "Run"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", size = 10),
        plot.subtitle   = element_text(size = 8, color = "grey30"),
        legend.position = "bottom")

ggsave(file.path(output_dir, paste0("cross_run_comparison_", today, ".pdf")),
       p_cross_run,
       width  = 9,
       height = max(4, nrow(cross_run_df) * 0.45 + 1.5),
       units  = "in", device = cairo_pdf)

# ══════════════════════════════════════════════════════════════════════════════
# 3. Per-gene long table — pool across all runs
# ══════════════════════════════════════════════════════════════════════════════
all_gene_rows <- lapply(names(results_summary), function(nm) {
  if (nm == "lodo") return(NULL)
  extract_gene_rows(nm, results_summary[[nm]])
})
all_gene_stats <- do.call(rbind, Filter(Negate(is.null), all_gene_rows))
rownames(all_gene_stats) <- NULL

# Drop genes where primary BalAcc is NA (e.g. zero-variance genes)
all_gene_stats <- all_gene_stats[!is.na(all_gene_stats$primary_balacc), ]

# ══════════════════════════════════════════════════════════════════════════════
# 4. Cross-dataset gene summary
#    One row per gene; aggregated across all runs for that gene
# ══════════════════════════════════════════════════════════════════════════════
cross_gene_summary <- all_gene_stats %>%
  group_by(gene_symbol, ensembl_id, chromosome, on_sex_chr) %>%
  summarise(
    # Coverage
    n_runs              = n(),
    n_datasets          = n_distinct(dataset),
    n_allgenes_runs     = sum(!sex_removed),
    n_nosexgenes_runs   = sum(sex_removed),
    datasets_seen       = paste(sort(unique(dataset)), collapse = "; "),
    
    # Primary: BalAcc — mean, SD, consistency
    mean_balacc         = round(mean(primary_balacc,  na.rm = TRUE), 3),
    median_balacc       = round(median(primary_balacc, na.rm = TRUE), 3),
    sd_balacc           = round(sd(primary_balacc,    na.rm = TRUE), 3),
    min_balacc          = round(min(primary_balacc,   na.rm = TRUE), 3),
    max_balacc          = round(max(primary_balacc,   na.rm = TRUE), 3),
    mean_null_balacc    = round(mean(null_balacc,     na.rm = TRUE), 3),
    mean_balacc_gain    = round(mean(primary_balacc - null_balacc, na.rm = TRUE), 3),
    
    # BalAcc threshold counts
    n_balacc_gt_06      = sum(primary_balacc >= 0.6, na.rm = TRUE),
    n_balacc_gt_07      = sum(primary_balacc >= 0.7, na.rm = TRUE),
    n_balacc_gt_08      = sum(primary_balacc >= 0.8, na.rm = TRUE),
    pct_balacc_gt_06    = round(mean(primary_balacc >= 0.6, na.rm = TRUE) * 100, 1),
    pct_balacc_gt_07    = round(mean(primary_balacc >= 0.7, na.rm = TRUE) * 100, 1),
    pct_balacc_gt_08    = round(mean(primary_balacc >= 0.8, na.rm = TRUE) * 100, 1),
    runs_balacc_gt_07   = paste(run[primary_balacc >= 0.7], collapse = "; "),
    
    # Sens / Spec
    mean_sens           = round(mean(primary_sens, na.rm = TRUE), 3),
    mean_spec           = round(mean(primary_spec, na.rm = TRUE), 3),
    sd_sens             = round(sd(primary_sens,   na.rm = TRUE), 3),
    sd_spec             = round(sd(primary_spec,   na.rm = TRUE), 3),
    
    # Secondary: AUC
    mean_auc            = round(mean(primary_auc,    na.rm = TRUE), 3),
    median_auc          = round(median(primary_auc,  na.rm = TRUE), 3),
    sd_auc              = round(sd(primary_auc,      na.rm = TRUE), 3),
    n_auc_gt_07         = sum(primary_auc >= 0.7,    na.rm = TRUE),
    pct_auc_gt_07       = round(mean(primary_auc >= 0.7, na.rm = TRUE) * 100, 1),
    
    .groups = "drop"
  ) %>%
  # Robustness score: mean_balacc scaled by consistency, penalised for variance
  mutate(
    robustness_score = round(
      mean_balacc * (pct_balacc_gt_07 / 100) * (1 - sd_balacc),
      3
    ),
    strength = cut(
      mean_balacc,
      breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
      labels = c("Poor (<0.6)", "Fair (0.6–0.7)", "Good (0.7–0.8)",
                 "Strong (0.8–0.9)", "Excellent (>0.9)"),
      include.lowest = TRUE
    ),
    consistency = cut(
      pct_balacc_gt_07,
      breaks = c(-Inf, 0, 25, 50, 75, 100),
      labels = c("Never", "Rarely", "Sometimes", "Usually", "Always"),
      include.lowest = TRUE
    )
  ) %>%
  arrange(desc(robustness_score), desc(mean_balacc), sd_balacc)

write.csv(cross_gene_summary,
          file.path(output_dir, paste0("cross_dataset_gene_summary_", today, ".csv")),
          row.names = FALSE)

cat("\n══════════════════════════════════════\n")
cat("Cross-dataset gene summary\n")
cat("══════════════════════════════════════\n")
cat("Unique genes evaluated:   ", nrow(cross_gene_summary), "\n")
cat("Always good (BalAcc≥0.7): ",
    sum(cross_gene_summary$consistency == "Always" &
          cross_gene_summary$n_runs > 1), "\n")
cat("Always strong (BalAcc≥0.8):",
    sum(cross_gene_summary$n_balacc_gt_08 == cross_gene_summary$n_runs &
          cross_gene_summary$n_runs > 1), "\n")
cat("Top 10 most robust biomarkers:\n")
print(head(cross_gene_summary[, c("gene_symbol", "mean_balacc", "sd_balacc",
                                  "pct_balacc_gt_07", "robustness_score",
                                  "mean_sens", "mean_spec", "consistency")], 10),
      digits = 3)

# ══════════════════════════════════════════════════════════════════════════════
# 5. Plot A: Robustness bubble plot
#    x = % runs with BalAcc≥0.7 (consistency), y = mean BalAcc, size = n_runs
# ══════════════════════════════════════════════════════════════════════════════
top30_robust <- cross_gene_summary %>%
  filter(n_runs >= 2) %>%
  slice_head(n = 30)

p_robust <- ggplot(top30_robust,
                   aes(x     = pct_balacc_gt_07,
                       y     = mean_balacc,
                       size  = n_runs,
                       color = strength,
                       label = gene_symbol)) +
  geom_point(alpha = 0.8) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 25, seed = 42) +
  geom_hline(yintercept = 0.7,  linetype = "dotted", color = "orange", linewidth = 0.6) +
  geom_hline(yintercept = 0.8,  linetype = "dotted", color = "red",    linewidth = 0.6) +
  geom_vline(xintercept = 50,   linetype = "dashed",  color = "grey60", linewidth = 0.5) +
  geom_vline(xintercept = 100,  linetype = "dashed",  color = "darkgreen", linewidth = 0.5) +
  scale_color_manual(values = c(
    "Poor (<0.6)"      = "grey60",
    "Fair (0.6–0.7)"   = "steelblue",
    "Good (0.7–0.8)"   = "darkgreen",
    "Strong (0.8–0.9)" = "orange",
    "Excellent (>0.9)" = "red"
  ), name = "Mean BalAcc strength") +
  scale_size_continuous(range = c(3, 10), name = "N runs") +
  scale_x_continuous(limits = c(-2, 107), breaks = seq(0, 100, 25)) +
  scale_y_continuous(limits = c(0.45, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  labs(
    title    = "Gene biomarker robustness across all runs",
    subtitle = "Top-right = high mean BalAcc + consistent across all runs\nSize = number of runs evaluated",
    x        = "% runs with BalAcc(opt) ≥ 0.7  (consistency)",
    y        = "Mean BalAcc(opt) across runs"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", size = 10),
        plot.subtitle   = element_text(size = 8, color = "grey30"),
        legend.position = "bottom")

ggsave(file.path(output_dir, paste0("cross_dataset_robustness_bubble_", today, ".pdf")),
       p_robust, width = 11, height = 9, units = "in", device = cairo_pdf)

# ══════════════════════════════════════════════════════════════════════════════
# 6. Plot B: BalAcc heatmap — top 30 genes × all runs
#    Cells show BalAcc(opt); AUC shown as column annotation
# ══════════════════════════════════════════════════════════════════════════════
top30_genes <- cross_gene_summary %>%
  filter(n_runs >= 2) %>%
  slice_head(n = 30) %>%
  pull(gene_symbol)

heat_df <- all_gene_stats %>%
  filter(gene_symbol %in% top30_genes) %>%
  dplyr::select(gene_symbol, run, primary_balacc, primary_auc)

# BalAcc matrix
heat_balacc <- heat_df %>%
  dplyr::select(gene_symbol, run, primary_balacc) %>%
  tidyr::pivot_wider(names_from = run, values_from = primary_balacc, values_fill = NA)
heat_mat <- as.matrix(heat_balacc[, -1])
rownames(heat_mat) <- heat_balacc$gene_symbol
# Order rows by mean BalAcc descending
heat_mat <- heat_mat[order(-rowMeans(heat_mat, na.rm = TRUE)), ]

# Column annotation: model type and sex_removed flag
col_meta <- cross_run_df[match(colnames(heat_mat), cross_run_df$run),
                         c("model_type", "sex_removed", "balacc_multi")]
col_annot_df <- data.frame(
  Model       = col_meta$model_type,
  SexRemoved  = ifelse(col_meta$sex_removed, "Yes", "No"),
  row.names   = colnames(heat_mat)
)

# Row annotation: chromosome, consistency
row_meta <- cross_gene_summary[match(rownames(heat_mat), cross_gene_summary$gene_symbol), ]
row_annot_df <- data.frame(
  Consistency = as.character(row_meta$consistency),
  SexChr      = ifelse(row_meta$on_sex_chr, "Yes", "No"),
  row.names   = rownames(heat_mat)
)

pheatmap::pheatmap(
  heat_mat,
  color           = colorRampPalette(c("#d73027", "white", "#4575b4"))(100),
  breaks          = seq(0.4, 1.0, length.out = 101),
  display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6,
  na_col          = "grey90",
  cluster_rows    = FALSE,   # already sorted by mean BalAcc
  cluster_cols    = TRUE,
  annotation_row  = row_annot_df,
  annotation_col  = col_annot_df,
  annotation_colors = list(
    Model       = c("RF+glmnet+glm" = "#e41a1c", "glmnet+glm" = "#377eb8"),
    SexRemoved  = c("Yes" = "orchid", "No" = "grey90"),
    Consistency = c("Always"    = "#1a9641",
                    "Usually"   = "#a6d96a",
                    "Sometimes" = "#ffffbf",
                    "Rarely"    = "#fdae61",
                    "Never"     = "#d7191c"),
    SexChr      = c("Yes" = "purple", "No" = "grey90")
  ),
  main     = paste0("BalAcc(opt) per gene × run — top 30 by robustness score"),
  filename = file.path(output_dir, paste0("cross_dataset_balacc_heatmap_", today, ".pdf")),
  width    = max(8, ncol(heat_mat) * 0.9 + 3),
  height   = max(7, nrow(heat_mat) * 0.35 + 2)
)

# ══════════════════════════════════════════════════════════════════════════════
# 7. Plot C: Top 20 genes ranked by robustness score
#    Bar = robustness_score; error bars = mean_balacc ± sd_balacc
#    Text = mean BalAcc, Sens, Spec
# ══════════════════════════════════════════════════════════════════════════════
top20_robust <- cross_gene_summary %>%
  filter(n_runs >= 2) %>%
  slice_head(n = 20)

p_rank <- ggplot(top20_robust,
                 aes(x    = robustness_score,
                     y    = reorder(gene_symbol, robustness_score),
                     fill = strength)) +
  geom_col(alpha = 0.75, width = 0.65) +
  # Mean BalAcc ± SD error bars
  geom_errorbarh(aes(xmin = pmax(0, mean_balacc - sd_balacc),
                     xmax = pmin(1, mean_balacc + sd_balacc)),
                 height = 0.35, color = "grey20", linewidth = 0.7) +
  # Mean BalAcc point
  geom_point(aes(x = mean_balacc), size = 3, shape = 16, color = "grey10") +
  # Annotation: BalAcc | Sens | Spec | n_runs
  geom_text(aes(x = 0.01,
                label = paste0("BalAcc=", sprintf("%.2f", mean_balacc),
                               "  Sens=", sprintf("%.0f%%", mean_sens * 100),
                               "  Spec=", sprintf("%.0f%%", mean_spec * 100),
                               "  n=", n_runs)),
            hjust = 0, size = 2.7, color = "grey15") +
  # Consistency colour band at right edge
  geom_text(aes(x     = max(top20_robust$robustness_score) * 1.05,
                label = consistency),
            hjust = 0, size = 2.5, color = "grey30") +
  scale_fill_manual(values = c(
    "Poor (<0.6)"      = "grey70",
    "Fair (0.6–0.7)"   = "steelblue",
    "Good (0.7–0.8)"   = "darkgreen",
    "Strong (0.8–0.9)" = "orange",
    "Excellent (>0.9)" = "red"
  ), name = "Strength") +
  scale_x_continuous(limits = c(0, max(top20_robust$robustness_score) * 1.35),
                     breaks = seq(0, 1, 0.1)) +
  labs(
    title    = "Top 20 most robust AD biomarkers — all runs",
    subtitle = paste0("Robustness = mean BalAcc × consistency × (1 − SD)",
                      "  |  ● = mean BalAcc  |  bars = ±1 SD"),
    x        = "Robustness score",
    y        = "Gene  (ascending robustness)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", size = 10),
        plot.subtitle   = element_text(size = 8, color = "grey30"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")

ggsave(file.path(output_dir, paste0("top_robust_biomarkers_ranked_", today, ".pdf")),
       p_rank, width = 11, height = max(6, nrow(top20_robust) * 0.45 + 2),
       units = "in", device = cairo_pdf)

# ══════════════════════════════════════════════════════════════════════════════
# 8. Plot D: Sensitivity vs Specificity trade-off — top 20 robust genes
# ══════════════════════════════════════════════════════════════════════════════
top20_ss <- top20_robust %>% filter(!is.na(mean_sens) & !is.na(mean_spec))

if (nrow(top20_ss) > 0) {
  p_sensspec <- ggplot(top20_ss,
                       aes(x     = mean_spec,
                           y     = mean_sens,
                           color = strength,
                           size  = mean_balacc,
                           label = gene_symbol)) +
    geom_point(alpha = 0.85) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20, seed = 42) +
    # Error bars: ±SD across runs
    geom_errorbar(aes(ymin = pmax(0, mean_sens - sd_sens),
                      ymax = pmin(1, mean_sens + sd_sens)),
                  width = 0.01, linewidth = 0.4, alpha = 0.5) +
    geom_errorbarh(aes(xmin = pmax(0, mean_spec - sd_spec),
                       xmax = pmin(1, mean_spec + sd_spec)),
                   height = 0.01, linewidth = 0.4, alpha = 0.5) +
    geom_hline(yintercept = 0.7, linetype = "dotted", color = "grey50") +
    geom_vline(xintercept = 0.7, linetype = "dotted", color = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey75") +
    scale_color_manual(values = c(
      "Poor (<0.6)"      = "grey60",
      "Fair (0.6–0.7)"   = "steelblue",
      "Good (0.7–0.8)"   = "darkgreen",
      "Strong (0.8–0.9)" = "orange",
      "Excellent (>0.9)" = "red"
    ), name = "Strength") +
    scale_size_continuous(range = c(3, 9), name = "Mean BalAcc") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title    = "Mean Sensitivity vs Specificity — top robust biomarkers",
      subtitle = "Error bars = ±1 SD across runs   Diagonal = balanced trade-off",
      x        = "Mean Specificity",
      y        = "Mean Sensitivity"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title      = element_text(face = "bold", size = 10),
          plot.subtitle   = element_text(size = 8, color = "grey30"),
          legend.position = "bottom",
          aspect.ratio    = 1)
  
  ggsave(file.path(output_dir, paste0("sensitivity_specificity_tradeoff_", today, ".pdf")),
         p_sensspec, width = 8, height = 8, units = "in", device = cairo_pdf)
}

# ══════════════════════════════════════════════════════════════════════════════
# 9. Plot E: allgenes vs nosexgenes BalAcc delta per gene
#    Shows the effect of removing sex-chromosome genes on each gene's BalAcc
# ══════════════════════════════════════════════════════════════════════════════
delta_df <- all_gene_stats %>%
  group_by(gene_symbol, dataset, model_type) %>%
  summarise(
    balacc_all    = mean(primary_balacc[!sex_removed], na.rm = TRUE),
    balacc_nosex  = mean(primary_balacc[sex_removed],  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(balacc_all) & !is.na(balacc_nosex)) %>%
  mutate(delta_balacc = round(balacc_nosex - balacc_all, 3))

if (nrow(delta_df) > 0) {
  top_delta_genes <- delta_df %>%
    group_by(gene_symbol) %>%
    summarise(mean_delta = mean(abs(delta_balacc), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_delta)) %>%
    slice_head(n = 25) %>%
    pull(gene_symbol)
  
  p_delta <- ggplot(
    delta_df %>% filter(gene_symbol %in% top_delta_genes),
    aes(x     = delta_balacc,
        y     = reorder(gene_symbol, delta_balacc, FUN = mean),
        color = delta_balacc > 0,
        shape = model_type)
  ) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.6) +
    geom_jitter(size = 3, height = 0.15, alpha = 0.85) +
    scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "tomato"),
                       labels = c("TRUE" = "Improved", "FALSE" = "Worsened"),
                       name   = "Effect of sex-chr removal") +
    scale_shape_manual(values = c("RF+glmnet+glm" = 16, "glmnet+glm" = 17),
                       name   = "Model type") +
    geom_vline(xintercept =  0.05, linetype = "dotted", color = "darkgreen", linewidth = 0.4) +
    geom_vline(xintercept = -0.05, linetype = "dotted", color = "tomato",    linewidth = 0.4) +
    labs(
      title    = "BalAcc change after removing sex-chromosome genes",
      subtitle = "Positive = gene performs better after sex-chr removal  |  Each point = one run",
      x        = "ΔBalAcc(opt)  [nosex − allgenes]",
      y        = "Gene  (ordered by mean ΔBalAcc)"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title      = element_text(face = "bold", size = 10),
          plot.subtitle   = element_text(size = 8, color = "grey30"),
          legend.position = "bottom")
  
  ggsave(file.path(output_dir, paste0("sex_removal_delta_balacc_", today, ".pdf")),
         p_delta,
         width  = 9,
         height = max(6, length(top_delta_genes) * 0.38 + 2),
         units  = "in", device = cairo_pdf)
}

# ══════════════════════════════════════════════════════════════════════════════
# 10. If LODO results exist: append LODO BalAcc to cross_gene_summary
# ══════════════════════════════════════════════════════════════════════════════
if (!is.null(lodo_results) && !is.null(lodo_results$lodo_gene_scores)) {
  lgs <- lodo_results$lodo_gene_scores
  
  # Identify mean_balacc_lodo column (set in run_lodo_validation)
  if ("mean_balacc_lodo" %in% names(lgs)) {
    lodo_balacc <- lgs[, c("gene_symbol", "mean_balacc_lodo", "mean_auc_lodo")]
    cross_gene_summary <- merge(cross_gene_summary, lodo_balacc,
                                by = "gene_symbol", all.x = TRUE)
    
    # Re-rank incorporating LODO BalAcc: composite = 0.6 * within-run + 0.4 * LODO
    cross_gene_summary$composite_score <- round(
      0.6 * cross_gene_summary$robustness_score +
        0.4 * ifelse(is.na(cross_gene_summary$mean_balacc_lodo), 0,
                     cross_gene_summary$mean_balacc_lodo),
      3
    )
    cross_gene_summary <- cross_gene_summary %>%
      arrange(desc(composite_score), desc(robustness_score))
    
    write.csv(cross_gene_summary,
              file.path(output_dir, paste0("cross_dataset_gene_summary_with_lodo_", today, ".csv")),
              row.names = FALSE)
    
    cat("\nLODO BalAcc appended to cross-gene summary.\n")
    cat("Composite score = 0.6 × within-run robustness + 0.4 × mean LODO BalAcc\n")
    cat("Top 10 by composite score:\n")
    print(head(cross_gene_summary[, c("gene_symbol", "robustness_score",
                                      "mean_balacc_lodo", "composite_score",
                                      "consistency")], 10), digits = 3)
    
    # ── Plot F: within-run robustness vs LODO BalAcc scatter ──────────────
    lodo_scatter_df <- cross_gene_summary %>%
      filter(!is.na(mean_balacc_lodo) & n_runs >= 2) %>%
      slice_head(n = 40)
    
    p_lodo_scatter <- ggplot(lodo_scatter_df,
                             aes(x     = robustness_score,
                                 y     = mean_balacc_lodo,
                                 color = strength,
                                 size  = n_datasets,
                                 label = gene_symbol)) +
      geom_point(alpha = 0.85) +
      ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20, seed = 42) +
      geom_hline(yintercept = 0.6, linetype = "dotted", color = "orange",    linewidth = 0.5) +
      geom_hline(yintercept = 0.7, linetype = "dotted", color = "red",       linewidth = 0.5) +
      geom_vline(xintercept = 0.4, linetype = "dashed",  color = "grey60",   linewidth = 0.4) +
      scale_color_manual(values = c(
        "Poor (<0.6)"      = "grey60",
        "Fair (0.6–0.7)"   = "steelblue",
        "Good (0.7–0.8)"   = "darkgreen",
        "Strong (0.8–0.9)" = "orange",
        "Excellent (>0.9)" = "red"
      ), name = "Strength") +
      scale_size_continuous(range = c(3, 8), name = "N datasets") +
      labs(
        title    = "Within-run robustness vs LODO generalisation",
        subtitle = "Top-right = performs well both within-dataset and on held-out datasets",
        x        = "Within-run robustness score",
        y        = "Mean LODO BalAcc (held-out datasets)"
      ) +
      theme_bw(base_size = 11) +
      theme(plot.title      = element_text(face = "bold", size = 10),
            plot.subtitle   = element_text(size = 8, color = "grey30"),
            legend.position = "bottom")
    
    ggsave(file.path(output_dir, paste0("robustness_vs_lodo_scatter_", today, ".pdf")),
           p_lodo_scatter, width = 9, height = 8, units = "in", device = cairo_pdf)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Final console summary
# ══════════════════════════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════\n")
cat("Downstream analysis complete.\n")
cat("Output directory:", output_dir, "\n")
cat("Files written:\n")
cat("  cross_run_summary_", today, ".csv\n", sep = "")
cat("  cross_run_comparison_", today, ".pdf\n", sep = "")
cat("  cross_dataset_gene_summary_", today, ".csv\n", sep = "")
cat("  cross_dataset_robustness_bubble_", today, ".pdf\n", sep = "")
cat("  cross_dataset_balacc_heatmap_", today, ".pdf\n", sep = "")
cat("  top_robust_biomarkers_ranked_", today, ".pdf\n", sep = "")
cat("  sensitivity_specificity_tradeoff_", today, ".pdf\n", sep = "")
cat("  sex_removal_delta_balacc_", today, ".pdf\n", sep = "")
if (!is.null(lodo_results))
  cat("  cross_dataset_gene_summary_with_lodo_", today, ".csv\n",
      "  robustness_vs_lodo_scatter_", today, ".pdf\n", sep = "")
cat("══════════════════════════════════════\n")