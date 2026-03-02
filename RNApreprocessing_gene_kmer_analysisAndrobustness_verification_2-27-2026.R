#!/usr/bin/env Rscript
# RNApreprocessing_gene_kmer_analysisAndrobustness_pipeline.R
# Full Robustness Pipeline — HPC-Optimized + Parallel + Bayesian + Publication Ready

# Prevent broken rstan from loading — use cmdstanr throughout
options(brms.backend  = "cmdstanr")
options(mc.cores      = parallel::detectCores())

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
rownames(Gage2024_metadata) <- Gage2024_metadata$Sample.label
rownames(Gage_merged_metadata) <- Gage_merged_metadata$Sample.label
# Remove rows with NA in Sample.label of Gage_merged_metadata
Gage_merged_metadata <- Gage_merged_metadata[!is.na(Gage_merged_metadata$Sample.label), ]
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

# Compare stable genes across different Diag groups using Gage2024_metadata Diag column and Sample.label
# Calculate mean values by gene in stable genes for genes with is_stable == TRUE by each Diag group
Gage2024_metadata$Diag <- factor(Gage2024_metadata$Diag, levels = c("CTRL", "AD", "Other"))
for (gene in stable_genes[stable_genes$is_stable == TRUE,]$ensembl_id) {
  cat(sprintf("Calculating mean counts for stable gene %s across Diag groups...\n", gene))
  AD_mean <- Gage2024_counts %>%
    filter(ensembl_id == gene) %>%
    pivot_longer(-ensembl_id, names_to = "Sample", values_to = "Count") %>%
    left_join(Gage2024_metadata, by = c("Sample" = "Sample.label")) %>%
    filter(Diag == "AD") %>%
    summarise(mean_count = mean(Count, na.rm=TRUE)) %>%
    pull(mean_count)
  CTRL_mean <- Gage2024_counts %>%
    filter(ensembl_id == gene) %>%
    pivot_longer(-ensembl_id, names_to = "Sample", values_to = "Count") %>%
    left_join(Gage2024_metadata, by = c("Sample" = "Sample.label")) %>%
    filter(Diag == "CTRL") %>%
    summarise(mean_count = mean(Count, na.rm=TRUE)) %>%
    pull(mean_count)
}
stable_gene_means <- Gage2024_counts %>%
  filter(ensembl_id %in% stable_genes[stable_genes$is_stable == TRUE,]$ensembl_id) %>%
  pivot_longer(-ensembl_id, names_to = "Sample", values_to = "Count") %>%
  left_join(Gage2024_metadata, by = c("Sample" = "Sample.label")) %>%
  group_by(Diag) %>%
  summarise(mean_count = mean(Count, na.rm=TRUE))

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

Gage2024_counts_sub <- Gage2024_counts[ranked_common_genes, , drop = FALSE]

# Sample metadata with diagnosis labels
# meta_new: data.frame with rownames = sample IDs, column Diag ∈ c("AD","CTRL")
# Subset Gage2024_counts_sub by samples present in Gage2024_metadata and ensure order matches
Gage2024_counts_sub <- Gage2024_counts_sub[, rownames(Gage2024_metadata), drop = TRUE]

dds <- DESeqDataSetFromMatrix(
  countData = round(Gage2024_counts_sub),
  colData   = Gage2024_metadata,
  design    = ~ Diag
)

# Filter very lowly expressed genes
dds <- dds[rowSums(counts(dds)) > 1, ]

# Variance stabilizing transform
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
expr <- assay(vsd)          # genes x samples
expr <- t(expr)             # samples x genes (rows = samples, columns = genes)

# Binary outcome: 1 = AD, 0 = CTRL
y <- ifelse(Gage2024_metadata$Diag == "AD", 1, 0)
X <- as.matrix(expr)                 # samples x genes
gene_names <- colnames(X)

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

library(pROC)
library(caret)

set.seed(42)

K <- 5
folds <- createFolds(y, k = K, list = TRUE, returnTrain = FALSE)

per_gene_results <- lapply(gene_names, function(g) {
  aucs <- numeric(K)
  
  for (i in seq_len(K)) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_along(y), test_idx)
    
    X_train <- X[train_idx, g, drop = FALSE]
    X_test  <- X[test_idx,  g, drop = FALSE]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    df_train <- data.frame(y = y_train, x = X_train[, 1])
    fit <- glm(y ~ x, data = df_train, family = binomial())
    
    df_test <- data.frame(x = X_test[, 1])
    p_hat <- predict(fit, newdata = df_test, type = "response")
    
    roc_obj <- roc(response = y_test, predictor = p_hat, quiet = TRUE)
    aucs[i] <- as.numeric(auc(roc_obj))
  }
  
  data.frame(
    gene = g,
    mean_auc = mean(aucs),
    sd_auc   = sd(aucs),
    n_splits = K
  )
})

per_gene_df <- do.call(rbind, per_gene_results)
per_gene_df <- per_gene_df[order(-per_gene_df$mean_auc), ]
head(per_gene_df, 10)

# pipeline from perplexity.ai
## =========================================
## 0. Setup
## =========================================

# Packages
library(pROC)
library(glmnet)
library(ggplot2)
library(reshape2)

# Data: expr (samples x genes), meta_new (Diag: AD/CTRL) ----------------------
# Ensure rownames(meta_new) match rownames(expr)
stopifnot(all(rownames(expr) %in% rownames(Gage2024_metadata)))
Gage2024_metadata <- Gage2024_metadata[rownames(expr), , drop = FALSE]

y <- ifelse(Gage2024_metadata$Diag == "AD", 1, 0)
X <- as.matrix(expr)
gene_names <- colnames(X)
n <- length(y)

set.seed(2026)
n_iter <- 100

# Storage ---------------------------------------------------------------------
per_gene_auc_mat <- matrix(
  NA_real_,
  nrow = length(gene_names),
  ncol = n_iter,
  dimnames = list(gene_names, paste0("iter_", 1:n_iter))
)

multi_auc_vec <- numeric(n_iter)

# Elastic net alpha -----------------------------------------------------------
alpha_val <- 0.5  # 0 = ridge, 1 = lasso, 0.5 = elastic net

## =========================================
## 1. Main loop: 100 iterations
## =========================================

for (it in 1:n_iter) {
  cat("Iteration", it, "...\n")
  
  ## -----------------------------
  ## 1A. Per-gene models with LOOCV
  ## -----------------------------
  iter_gene_auc <- numeric(length(gene_names))
  names(iter_gene_auc) <- gene_names
  
  for (g in gene_names) {
    pred_all <- numeric(n)
    
    # LOOCV
    for (i in 1:n) {
      test_idx  <- i
      train_idx <- setdiff(1:n, i)
      
      X_train <- X[train_idx, g, drop = FALSE]
      X_test  <- X[test_idx,  g, drop = FALSE]
      y_train <- y[train_idx]
      
      df_train <- data.frame(y = y_train, x = X_train[, 1])
      fit <- glm(y ~ x, data = df_train, family = binomial())
      
      df_test <- data.frame(x = X_test[, 1])
      p_hat <- predict(fit, newdata = df_test, type = "response")
      
      pred_all[test_idx] <- p_hat
    }
    
    roc_g <- roc(response = y, predictor = pred_all, quiet = TRUE)
    iter_gene_auc[g] <- as.numeric(auc(roc_g))
  }
  
  per_gene_auc_mat[, it] <- iter_gene_auc
  
  
  ## -----------------------------
  ## 1B. Multi-gene elastic net with LOOCV + inner CV
  ## -----------------------------
  pred_all_multi <- numeric(n)
  
  for (i in 1:n) {
    test_idx  <- i
    train_idx <- setdiff(1:n, i)
    
    X_train <- X[train_idx, , drop = FALSE]
    X_test  <- X[test_idx,  , drop = FALSE]
    y_train <- y[train_idx]
    
    # 0 = unpenalized (always included), 1 = penalized (default)
    penalty_factors <- c(
      rep(1, ncol(X)),    # genes - penalized
      rep(0, ncol(cov_mat))  # Age, Sex - unpenalized
    )
    
    # Inner CV for lambda on training data only
    cvfit <- cv.glmnet(
      x              = X_train,
      y              = y_train,
      family         = "binomial",
      alpha          = alpha_val,
      nfolds         = min(5, length(y_train)),
      type.measure   = "auc",
      penalty.factor = penalty_factors
    )
    
    lambda_best <- cvfit$lambda.1se
    
    p_hat <- predict(
      cvfit,
      newx = X_test,
      s = lambda_best,
      type = "response"
    )[1, 1]
    
    pred_all_multi[test_idx] <- p_hat
  }
  
  roc_multi <- roc(response = y, predictor = pred_all_multi, quiet = TRUE)
  multi_auc_vec[it] <- as.numeric(auc(roc_multi))
  
  cat("  Multi-gene elastic net LOOCV AUC:",
      round(multi_auc_vec[it], 3), "\n")
}


## =========================================
## 2. Save numeric outputs
## =========================================

# Per-gene AUC across iterations ----------------------------------------------
per_gene_auc_df <- as.data.frame(per_gene_auc_mat)
per_gene_auc_df$gene <- rownames(per_gene_auc_df)

per_gene_auc_long <- melt(
  per_gene_auc_df,
  id.vars = "gene",
  variable.name = "iteration",
  value.name = "auc"
)

write.csv(
  per_gene_auc_df,
  file = file.path(output_dir, "per_gene_auc_100iters_LOOCV_wide.csv"),
  row.names = FALSE
)

write.csv(
  per_gene_auc_long,
  file = file.path(output_dir, "per_gene_auc_100iters_LOOCV_long.csv"),
  row.names = FALSE
)

# Multigene AUC across iterations ---------------------------------------------
multi_auc_df <- data.frame(
  iteration = 1:n_iter,
  auc       = multi_auc_vec
)

write.csv(
  multi_auc_df,
  file = file.path(output_dir, "multigene_elasticnet_auc_100iters_LOOCV.csv"),
  row.names = FALSE
)


## =========================================
## 3. Plots: per-gene vs multigene AUC
## =========================================

# 3A. Per-gene AUC distribution (boxplot) -------------------------------------
p_gene_box <- ggplot(per_gene_auc_long, aes(x = gene, y = auc)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5
    )
  ) +
  labs(
    title = "Per-gene ROC AUC across 100 LOOCV iterations",
    x = "Gene",
    y = "AUC"
  )

ggsave(
  filename = file.path(output_dir, "per_gene_auc_boxplot_LOOCV.png"),
  plot = p_gene_box,
  width = 10,
  height = 5,
  dpi = 300
)

# 3B. Mean per-gene AUC (barplot, sorted) -------------------------------------
gene_mean_auc <- aggregate(
  auc ~ gene,
  data = per_gene_auc_long,
  FUN = mean
)
gene_mean_auc <- gene_mean_auc[order(-gene_mean_auc$auc), ]

p_gene_bar <- ggplot(gene_mean_auc, aes(x = reorder(gene, auc), y = auc)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(
    title = "Mean per-gene ROC AUC (100 LOOCV iterations)",
    x = "Gene",
    y = "Mean AUC"
  )

ggsave(
  filename = file.path(output_dir, "per_gene_auc_mean_barplot_LOOCV.png"),
  plot = p_gene_bar,
  width = 6,
  height = 8,
  dpi = 300
)

# 3C. Multi-gene AUC distribution (histogram) ---------------------------------
p_multi <- ggplot(multi_auc_df, aes(x = auc)) +
  geom_histogram(color = "black", fill = "steelblue", bins = 15) +
  theme_bw(base_size = 10) +
  labs(
    title = "Multi-gene elastic net ROC AUC (100 LOOCV iterations)",
    x = "AUC",
    y = "Count"
  )

ggsave(
  filename = file.path(output_dir, "multigene_elasticnet_auc_hist_LOOCV.png"),
  plot = p_multi,
  width = 5,
  height = 4,
  dpi = 300
)

# 3D. Compare per-gene vs multi-gene (median line) ----------------------------
multi_auc_median <- median(multi_auc_vec)

p_compare <- ggplot(gene_mean_auc, aes(x = auc)) +
  geom_histogram(color = "black", fill = "grey70", bins = 15) +
  geom_vline(xintercept = multi_auc_median, color = "red", size = 1.1) +
  theme_bw(base_size = 10) +
  labs(
    title = "Mean per-gene AUC vs multi-gene elastic net\nRed line = median multi-gene AUC",
    x = "AUC",
    y = "Number of genes"
  )

ggsave(
  filename = file.path(output_dir, "per_gene_vs_multigene_auc_compare_LOOCV.png"),
  plot = p_compare,
  width = 5,
  height = 4,
  dpi = 300
)


## =========================================
## 4. Final elastic-net model on full data
## =========================================

cvfit_full <- cv.glmnet(
  x = X,
  y = y,
  family = "binomial",
  alpha = alpha_val,
  nfolds = min(5, n),
  type.measure = "auc"
)

lambda_best_full <- cvfit_full$lambda.min # Trying with less aggressive penalizer compared to lambda.1se

# Coefficients at chosen lambda
coef_full <- coef(cvfit_full, s = lambda_best_full)
coef_vec  <- as.numeric(coef_full)
names(coef_vec) <- rownames(coef_full)

nz <- which(coef_vec != 0)
nz_genes <- setdiff(names(coef_vec)[nz], "(Intercept)")

final_model_df <- data.frame(
  feature = names(coef_vec),
  coefficient = coef_vec
)

write.csv(
  final_model_df,
  file = file.path(output_dir, "final_elasticnet_coefficients_full_data.csv"),
  row.names = FALSE
)

cat("Final elastic net fit on full data:\n")
cat("  alpha:", alpha_val, "\n")
cat("  lambda.1se:", lambda_best_full, "\n")
cat("  Non-zero genes:", paste(nz_genes, collapse = ", "), "\n")

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

# Settings
output_dir <- "/home/thilsabeck/Documents/Gage2024/RNApreprocessing_gene_validation_files/AgeSex_confounds/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

alpha_val <- 0.5
B_boot    <- 1000
today     <- format(Sys.Date(), "%m-%d-%Y")

results_summary <- list()

for (dname in names(normalized_list)) {
  cat("=== Dataset:", dname, "===\n")
  
  nd       <- normalized_list[[dname]]
  X        <- nd$X
  y        <- nd$y
  meta_use <- nd$metadata
  
  # Map gene names to both ensembl_id and hgnc_symbol
  if (nd$gene_names[1] %in% stable_genes$ensembl_id) {
    # Input is ensembl — map to hgnc
    ensembl_ids    <- nd$gene_names
    hgnc_symbols   <- stable_genes$hgnc_symbol[match(nd$gene_names, stable_genes$ensembl_id)]
  } else if (nd$gene_names[1] %in% stable_genes$hgnc_symbol) {
    # Input is hgnc — map to ensembl
    hgnc_symbols   <- nd$gene_names
    ensembl_ids    <- stable_genes$ensembl_id[match(nd$gene_names, stable_genes$hgnc_symbol)]
  } else {
    warning(dname, ": Gene names do not match either ensembl_id or hgnc_symbol in stable_genes, keeping original names")
    hgnc_symbols <- nd$gene_names
    ensembl_ids  <- rep(NA, length(nd$gene_names))
  }
  
  stopifnot(nrow(X) == length(y))
  stopifnot(all(rownames(X) == rownames(meta_use)))
  
  # 1. LOOCV per-gene AUCs + bootstrap CIs
  gene_names   <- colnames(X)
  gene_symbols <- hgnc_symbols[match(gene_names, ensembl_ids)]
  n            <- length(y)
  auc_vec      <- numeric(length(gene_names))
  ci_low_vec   <- numeric(length(gene_names))
  ci_high_vec  <- numeric(length(gene_names))
  names(auc_vec) <- names(ci_low_vec) <- names(ci_high_vec) <- gene_names
  gene_n = 0
  # LOOCV per-gene AUCs - corrected to control for sex and age
  for (g in gene_names) {
    gene_n <- gene_n + 1
    pred_all <- numeric(n)
    
    # Skip genes with no variance
    if (var(X[, g]) == 0) {
      auc_vec[g] <- NA; ci_low_vec[g] <- NA; ci_high_vec[g] <- NA
      next
    }
    
    for (i in 1:n) {
      train_idx <- setdiff(1:n, i)
      
      df_train <- data.frame(
        y   = y[train_idx],
        x   = X[train_idx, g],
        sex = as.factor(meta_use$Sex[train_idx]),
        age = as.numeric(meta_use$Age[train_idx])
      )
      df_test <- data.frame(
        x   = X[i, g],
        sex = as.factor(meta_use$Sex[i]),
        age = as.numeric(meta_use$Age[i])
      )
      
      # Check variance in training covariates
      if (length(unique(df_train$sex)) < 2) {
        # Only one sex in training — drop sex from model
        fit <- tryCatch(
          glm(y ~ x + age, data = df_train, family = binomial()),
          error   = function(e) NULL,
          warning = function(w) suppressWarnings(
            glm(y ~ x + age, data = df_train, family = binomial())
          )
        )
      } else {
        fit <- tryCatch(
          glm(y ~ x + sex + age, data = df_train, family = binomial()),
          error   = function(e) NULL,
          warning = function(w) suppressWarnings(
            glm(y ~ x + sex + age, data = df_train, family = binomial())
          )
        )
      }
      
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
    
    if (is.null(roc_g)) {
      auc_vec[g] <- NA; ci_low_vec[g] <- NA; ci_high_vec[g] <- NA; next
    }
    
    auc_vec[g]     <- as.numeric(auc(roc_g))
    ci_g           <- bootstrap_auc_ci(y = y, p_hat = pred_all, B = B_boot, conf = 0.95)
    ci_low_vec[g]  <- ci_g$ci_low
    ci_high_vec[g] <- ci_g$ci_high
    
    # Box plot with age and sex annotation
    plot_df <- data.frame(
      expression = as.numeric(X[, g]),
      group      = as.factor(meta_use[["Diag"]]),
      sex        = meta_use$Sex,
      age        = meta_use$Age,
      sample     = rownames(meta_use)
    )
    
    p_box <- ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(aes(shape = sex), width = 0.15, size = 2, alpha = 0.8) +
      labs(
        title = paste0(hgnc_symbols[gene_n], ":", g,
                       " | AUC: ", round(auc_vec[g], 3),
                       " [", round(ci_low_vec[g], 3),
                       ", ", round(ci_high_vec[g], 3), "]",
                       " (adjusted for age + sex)"),
        x    = "Diagnosis",
        y    = "Expression",
        fill = "Diagnosis",
        shape = "Sex"
      ) +
      theme_bw() +
      theme(
        plot.title      = element_text(face = "bold", size = 11),
        axis.text.x     = element_text(angle = 45, hjust = 1)
      )
    
    ggsave(
      filename = file.path(output_dir, paste0(dname, "_", hgnc_symbols[gene_n], "_",g,
                                              "_expression_by_Diag_", today, ".pdf")),
      plot = p_box, width = 6, height = 5, device = cairo_pdf
    )
  }
  
  per_gene_df <- data.frame(
    ensembl_id  = ensembl_ids,
    gene_symbol = hgnc_symbols,
    auc         = as.numeric(auc_vec),
    ci_low      = as.numeric(ci_low_vec),
    ci_high     = as.numeric(ci_high_vec)
  )
  
  # Filter to meaningful AUCs and sort
  per_gene_df_plot <- per_gene_df[!is.na(per_gene_df$auc), ]
  per_gene_df_plot <- per_gene_df_plot[order(-per_gene_df_plot$auc), ]
  
  # Classify biomarker strength
  per_gene_df_plot$strength <- cut(
    per_gene_df_plot$auc,
    breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels = c("Poor (<0.6)", "Fair (0.6-0.7)", "Good (0.7-0.8)",
               "Strong (0.8-0.9)", "Excellent (>0.9)"),
    include.lowest = TRUE
  )
  
  # Top 20 biomarker dot plot with CI
  top20 <- head(per_gene_df_plot, 20)
  
  p_biomarker <- ggplot(top20,
                        aes(x = auc,
                            y = reorder(gene_symbol, auc),
                            color = strength)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0.7, linetype = "dotted", color = "orange") +
    geom_vline(xintercept = 0.8, linetype = "dotted", color = "red") +
    scale_color_manual(values = c(
      "Poor (<0.6)"      = "grey60",
      "Fair (0.6-0.7)"   = "steelblue",
      "Good (0.7-0.8)"   = "darkgreen",
      "Strong (0.8-0.9)" = "orange",
      "Excellent (>0.9)" = "red"
    )) +
    scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
    labs(
      title    = paste0("Top AD Biomarkers - ", dname,
                        "\n(LOOCV AUC adjusted for age + sex)"),
      x        = "AUC (95% bootstrap CI)",
      y        = "Gene",
      color    = "Biomarker Strength"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
  
  ggsave(
    file.path(output_dir, paste0("top_biomarkers_", dname, "_", today, ".pdf")),
    p_biomarker, width = 8, height = 7, device = cairo_pdf
  )
  
  # Also save a summary table of top biomarkers
  write.csv(
    head(per_gene_df_plot, 20),
    file      = file.path(output_dir, paste0("top20_biomarkers_", dname, "_", today, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    per_gene_df,
    file      = file.path(output_dir, paste0("per_gene_LOOCV_auc_", dname, "_", today, ".csv")),
    row.names = FALSE
  )
  
  # 2. LOOCV elastic net (sex + age unpenalized, already handled in run_loocv_elastic_net)
  pred_multi <- run_loocv_elastic_net(X, y, meta_use = meta_use, alpha_val = alpha_val)
  roc_multi  <- roc(response = y, predictor = pred_multi, quiet = TRUE)
  auc_multi  <- as.numeric(auc(roc_multi))
  
  # 3. Bootstrap CI for multigene AUC
  boot_res <- bootstrap_auc_ci(y, pred_multi, B = B_boot, conf = 0.95)
  
  boot_df <- data.frame(
    iteration = 1:B_boot,
    auc       = boot_res$auc_boot
  )
  write.csv(
    boot_df,
    file      = file.path(output_dir, paste0("multigene_bootstrap_auc_", dname, "_", today, ".csv")),
    row.names = FALSE
  )
  
  # Extract elastic net coefficients (averaged across LOOCV folds for stability)
  # Refit on full data to get representative coefficients for interpretation
  sex_numeric <- ifelse(toupper(meta_use$Sex) == "M", 1, 0)
  cov_mat     <- cbind(Age = as.numeric(meta_use$Age), Sex = sex_numeric)
  X_full      <- cbind(X, cov_mat)
  penalty_factors <- c(rep(1, ncol(X)), rep(0, ncol(cov_mat)))
  
  # Fit elastic net and lasso at both lambdas
  cvfit_enet  <- cv.glmnet(
    x              = X_full,
    y              = y,
    family         = "binomial",
    alpha          = alpha_val,   # elastic net
    nfolds         = min(5, length(y)),
    type.measure   = "auc",
    penalty.factor = penalty_factors
  )
  
  cvfit_lasso <- cv.glmnet(
    x              = X_full,
    y              = y,
    family         = "binomial",
    alpha          = 1,           # lasso
    nfolds         = min(5, length(y)),
    type.measure   = "auc",
    penalty.factor = penalty_factors
  )
  
  # Report selection at both lambdas for both models
  cat("\nElastic net (alpha=", alpha_val, "):\n")
  cat("  lambda.min selected:", sum(coef(cvfit_enet, s = "lambda.min")[-1] != 0), "genes\n")
  cat("  lambda.1se selected:", sum(coef(cvfit_enet, s = "lambda.1se")[-1] != 0), "genes\n")
  cat("Lasso (alpha=1):\n")
  cat("  lambda.min selected:", sum(coef(cvfit_lasso, s = "lambda.min")[-1] != 0), "genes\n")
  cat("  lambda.1se selected:", sum(coef(cvfit_lasso, s = "lambda.1se")[-1] != 0), "genes\n")
  # Save CV curve plots
  pdf(file.path(output_dir, paste0("glmnet_CV_curves_", dname, "_", today, ".pdf")),
      width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot(cvfit_enet,  main = paste("Elastic net -", dname))
  plot(cvfit_lasso, main = paste("Lasso -", dname))
  dev.off()
  
  # Choose best model: prefer lambda.min, prefer lasso if enet selects 0
  best_fit    <- if (sum(coef(cvfit_enet,  s = "lambda.min")[-1] != 0) > 0) cvfit_enet  else cvfit_lasso
  best_lambda <- "lambda.min"
  coef_full   <- coef(best_fit, s = best_lambda)
  
  coef_df <- data.frame(
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
  
  # Collinearity plot for top 20 genes by single-gene AUC
  top20_genes <- head(per_gene_df_sorted$ensembl_id[
    per_gene_df_sorted$ensembl_id %in% colnames(X)], 20)
  top20_symbols <- stable_genes$hgnc_symbol[match(top20_genes, stable_genes$ensembl_id)]
  top20_symbols[is.na(top20_symbols)] <- top20_genes[is.na(top20_symbols)]
  
  if (length(top20_genes) > 1) {
    X_top20 <- X[, top20_genes, drop = FALSE]
    colnames(X_top20) <- top20_symbols
    cor_mat <- cor(X_top20, use = "pairwise.complete.obs")
    
    # Heatmap
    p_cor_heat <- ggplot(
      reshape2::melt(cor_mat),
      aes(x = Var1, y = Var2, fill = value)
    ) +
      geom_tile() +
      geom_text(aes(label = round(value, 2)), size = 2.5) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                           midpoint = 0, limits = c(-1, 1)) +
      labs(
        title = paste0("Top gene collinearity - ", dname),
        x = NULL, y = NULL, fill = "Pearson r"
      ) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(
      file.path(output_dir, paste0("collinearity_heatmap_", dname, "_", today, ".pdf")),
      p_cor_heat,
      width  = max(6, length(top20_genes) * 0.5),
      height = max(5, length(top20_genes) * 0.5),
      device = cairo_pdf
    )
    
    # Clustermap version using pheatmap for cleaner clustering
    pheatmap::pheatmap(
      cor_mat,
      color            = colorRampPalette(c("steelblue", "white", "tomato"))(100),
      breaks           = seq(-1, 1, length.out = 101),
      display_numbers  = TRUE,
      number_format    = "%.2f",
      fontsize_number  = 7,
      main             = paste0("Top gene correlation - ", dname),
      filename         = file.path(output_dir,
                                   paste0("collinearity_clustered_", dname, "_", today, ".pdf")),
      width  = max(6, length(top20_genes) * 0.4),
      height = max(5, length(top20_genes) * 0.4)
    )
    
    # Also plot collinearity for elastic net selected genes if any
    if (nrow(coef_df) > 1) {
      sel_genes   <- coef_df$feature[coef_df$feature %in% colnames(X)]
      sel_symbols <- coef_df$gene_symbol[coef_df$feature %in% colnames(X)]
      X_sel       <- X[, sel_genes, drop = FALSE]
      colnames(X_sel) <- sel_symbols
      cor_sel     <- cor(X_sel, use = "pairwise.complete.obs")
      
      pheatmap::pheatmap(
        cor_sel,
        color            = colorRampPalette(c("steelblue", "white", "tomato"))(100),
        breaks           = seq(-1, 1, length.out = 101),
        display_numbers  = TRUE,
        number_format    = "%.2f",
        fontsize_number  = 8,
        main             = paste0("Selected gene correlation - ", dname),
        filename         = file.path(output_dir,
                                     paste0("selected_gene_correlation_", dname, "_", today, ".pdf")),
        width  = max(4, nrow(coef_df) * 0.5),
        height = max(4, nrow(coef_df) * 0.5)
      )
    }
  }
  
  
  # Map ensembl to symbol for selected genes
  coef_df$gene_symbol <- stable_genes$hgnc_symbol[match(coef_df$feature, stable_genes$ensembl_id)]
  coef_df$gene_symbol[is.na(coef_df$gene_symbol)] <- coef_df$feature[is.na(coef_df$gene_symbol)]
  coef_df <- coef_df[order(-abs(coef_df$coef)), ]
  
  write.csv(
    coef_df,
    file      = file.path(output_dir, paste0("elastic_net_coefficients_", dname, "_", today, ".csv")),
    row.names = FALSE
  )
  
  # 4A. Per-gene AUC barplot with CI error bars
  per_gene_df_sorted <- per_gene_df[!is.na(per_gene_df$auc), ]
  per_gene_df_sorted <- per_gene_df_sorted[order(-per_gene_df_sorted$auc), ]
  
  p_gene_bar <- ggplot(per_gene_df_sorted,
                       aes(x = reorder(gene_symbol, auc), y = auc)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0.3, color = "black") +
    geom_hline(yintercept = 0.7, linetype = "dotted", color = "orange") +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "red") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(
      title    = paste("Per-gene LOOCV AUC (adj. age + sex) -", dname),
      x        = "Gene",
      y        = "AUC (95% CI)"
    )
  
  ggsave(
    file.path(output_dir, paste0("per_gene_LOOCV_auc_bar_", dname, "_", today, ".png")),
    p_gene_bar, width = 6, height = 8, dpi = 300
  )
  
  # 4B. Multigene bootstrap AUC histogram with strength annotations
  auc_strength <- cut(
    auc_multi,
    breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels = c("Poor", "Fair", "Good", "Strong", "Excellent"),
    include.lowest = TRUE
  )
  
  p_boot <- ggplot(boot_df, aes(x = auc)) +
    geom_histogram(color = "black", fill = "steelblue", bins = 30) +
    geom_vline(xintercept = auc_multi,        color = "red",    linewidth = 1.1) +
    geom_vline(xintercept = boot_res$ci_low,  color = "red",    linetype = "dashed") +
    geom_vline(xintercept = boot_res$ci_high, color = "red",    linetype = "dashed") +
    geom_vline(xintercept = 0.7,              color = "orange", linetype = "dotted") +
    geom_vline(xintercept = 0.8,              color = "darkred",linetype = "dotted") +
    annotate("text", x = auc_multi, y = Inf,
             label = paste0("AUC=", round(auc_multi, 3), " (", auc_strength, ")"),
             vjust = 2, hjust = -0.1, color = "red", size = 3.5) +
    theme_bw(base_size = 10) +
    labs(
      title = paste0("Multigene LOOCV AUC (adj. age + sex) - ", dname,
                     "\nMean: ", round(boot_res$mean, 3),
                     " CI [", round(boot_res$ci_low, 3),
                     ", ", round(boot_res$ci_high, 3), "]"),
      x = "AUC",
      y = "Count"
    )
  
  ggsave(
    file.path(output_dir, paste0("multigene_bootstrap_auc_hist_", dname, "_", today, ".png")),
    p_boot, width = 6, height = 4, dpi = 300
  )
  
  # 4C. Per-gene vs multigene comparison
  p_compare <- ggplot(per_gene_df_sorted, aes(x = auc)) +
    geom_histogram(color = "black", fill = "grey70", bins = 15) +
    geom_vline(xintercept = auc_multi,        color = "red",    linewidth = 1.1) +
    geom_vline(xintercept = boot_res$ci_low,  color = "red",    linetype = "dashed") +
    geom_vline(xintercept = boot_res$ci_high, color = "red",    linetype = "dashed") +
    geom_vline(xintercept = 0.5,              color = "grey40", linetype = "dashed") +
    annotate("text", x = auc_multi, y = Inf,
             label = paste0("Multigene\nAUC=", round(auc_multi, 3)),
             vjust = 2, hjust = -0.1, color = "red", size = 3) +
    theme_bw(base_size = 10) +
    labs(
      title = paste("Per-gene vs multigene AUC (adj. age + sex) -", dname),
      x     = "AUC",
      y     = "Number of genes"
    )
  
  ggsave(
    file.path(output_dir, paste0("per_gene_vs_multigene_auc_", dname, "_", today, ".png")),
    p_compare, width = 5, height = 4, dpi = 300
  )
  
  # 4D. Elastic net selected gene coefficients — most interpretable biomarker plot
  if (nrow(coef_df) > 0) {
    coef_df$direction <- ifelse(coef_df$coef > 0, "Higher in AD", "Lower in AD")
    
    p_coef <- ggplot(coef_df,
                     aes(x = coef,
                         y = reorder(gene_symbol, abs(coef)),
                         fill = direction)) +
      geom_col(alpha = 0.85) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
      scale_fill_manual(values = c("Higher in AD" = "tomato",
                                   "Lower in AD"  = "steelblue")) +
      labs(
        title = paste0("Elastic Net Selected Genes - ", dname,
                       "\n(LOOCV AUC: ", round(auc_multi, 3),
                       " [", round(boot_res$ci_low, 3),
                       ", ", round(boot_res$ci_high, 3), "]",
                       " | adj. age + sex)"),
        x    = "Coefficient (effect size)",
        y    = "Gene",
        fill = "Direction"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title      = element_text(face = "bold", size = 11),
        legend.position = "bottom"
      )
    
    ggsave(
      file.path(output_dir, paste0("elastic_net_coef_", dname, "_", today, ".pdf")),
      p_coef, width = 7, height = max(4, nrow(coef_df) * 0.35), device = cairo_pdf
    )
  }
  
  # 5. Store summary
  results_summary[[dname]] <- list(
    per_gene_df      = per_gene_df,
    auc_multi        = auc_multi,
    boot_mean        = boot_res$mean,
    boot_ci_low      = boot_res$ci_low,
    boot_ci_high     = boot_res$ci_high,
    elastic_net_genes = coef_df
  )
  
  cat("Dataset", dname, "done.\n")
  cat("  Multigene LOOCV AUC (adj. age + sex):", round(auc_multi, 3), "\n")
  cat("  Bootstrap AUC mean:", round(boot_res$mean, 3),
      "CI [", round(boot_res$ci_low, 3), ",", round(boot_res$ci_high, 3), "]\n")
  cat("  Elastic net selected", nrow(coef_df), "genes\n")
}

# Save overall summary
summary_df <- do.call(rbind, lapply(names(results_summary), function(dname) {
  s <- results_summary[[dname]]
  data.frame(
    dataset   = dname,
    auc_multi = s$auc_multi,
    boot_mean = s$boot_mean,
    ci_low    = s$boot_ci_low,
    ci_high   = s$boot_ci_high
  )
}))

write.csv(
  summary_df,
  file      = file.path(output_dir, paste0("multidataset_auc_summary_", today, ".csv")),
  row.names = FALSE
)

### Modified version of above to also run the analysis after removing sex-related genes and to save the single gene plots in a single pdf 03-02-2026
output_dir <- "/home/thilsabeck/Documents/Gage2024/RNApreprocessing_gene_validation_files/AgeSex_confounds/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

alpha_val <- 0.5
B_boot    <- 1000
today     <- format(Sys.Date(), "%m-%d-%Y")

# Sex chromosome gene patterns to remove
sex_chr_patterns <- "XIST|DDX3Y|USP9Y|RPS4Y1|ZFY|TSIX|EIF1AY|KDM5D|NLGN4Y|TMSB4Y|UTY|PRKY|PCDH11Y"

results_summary <- list()

# Helper to run the full per-gene + elastic net analysis for one X matrix
run_analysis <- function(X, y, meta_use, dname, label, ensembl_ids, hgnc_symbols,
                         output_dir, today, alpha_val, B_boot) {
  
  gene_names   <- colnames(X)
  gene_symbols <- hgnc_symbols[match(gene_names, ensembl_ids)]
  n            <- length(y)
  auc_vec      <- numeric(length(gene_names))
  ci_low_vec   <- numeric(length(gene_names))
  ci_high_vec  <- numeric(length(gene_names))
  names(auc_vec) <- names(ci_low_vec) <- names(ci_high_vec) <- gene_names
  
  # Collect per-gene plots for combined PDF
  boxplot_list <- list()
  gene_n <- 0
  
  for (g in gene_names) {
    gene_n <- gene_n + 1
    pred_all <- numeric(n)
    
    if (var(X[, g]) == 0) {
      auc_vec[g] <- NA; ci_low_vec[g] <- NA; ci_high_vec[g] <- NA
      next
    }
    
    for (i in 1:n) {
      train_idx <- setdiff(1:n, i)
      
      df_train <- data.frame(
        y   = y[train_idx],
        x   = X[train_idx, g],
        sex = as.factor(meta_use$Sex[train_idx]),
        age = as.numeric(meta_use$Age[train_idx])
      )
      df_test <- data.frame(
        x   = X[i, g],
        sex = as.factor(meta_use$Sex[i]),
        age = as.numeric(meta_use$Age[i])
      )
      
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
    
    if (is.null(roc_g)) {
      auc_vec[g] <- NA; ci_low_vec[g] <- NA; ci_high_vec[g] <- NA; next
    }
    
    auc_vec[g]    <- as.numeric(auc(roc_g))
    ci_g          <- bootstrap_auc_ci(y = y, p_hat = pred_all, B = B_boot, conf = 0.95)
    ci_low_vec[g] <- ci_g$ci_low
    ci_high_vec[g]<- ci_g$ci_high
    
    # Box plot
    sym <- if (!is.na(gene_symbols[gene_n])) gene_symbols[gene_n] else g
    plot_df <- data.frame(
      expression = as.numeric(X[, g]),
      group      = as.factor(meta_use[["Diag"]]),
      sex        = as.factor(meta_use$Sex),
      age        = as.numeric(meta_use$Age),
      sample     = rownames(meta_use)
    )
    
    p_box <- ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(aes(shape = sex, color = age), width = 0.2, size = 2.5, alpha = 0.9) +
      scale_shape_manual(values = c("M" = 16, "F" = 17), na.value = 15) +
      scale_color_gradient(low = "lightblue", high = "darkred", name = "Age") +
      labs(
        title = paste0(sym, ":", g,
                       " | AUC: ", round(auc_vec[g], 3),
                       " [", round(ci_low_vec[g], 3),
                       ", ", round(ci_high_vec[g], 3), "]",
                       " (adj. age+sex)"),
        x = "Diagnosis", y = "Expression",
        fill = "Diagnosis", shape = "Sex"
      ) +
      theme_bw() +
      theme(
        plot.title      = element_text(face = "bold", size = 10),
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
    
    boxplot_list[[g]] <- p_box
  }
  
  # Save all gene boxplots in one combined PDF
  combined_pdf_path <- file.path(output_dir,
                                 paste0(label, "_all_gene_boxplots_", today, ".pdf"))
  pdf(combined_pdf_path, width = 7, height = 5)
  for (p in boxplot_list) print(p)
  dev.off()
  message(label, ": Saved ", length(boxplot_list), " gene plots to ", combined_pdf_path)
  
  # Build per_gene_df
  per_gene_df <- data.frame(
    ensembl_id  = ensembl_ids,
    gene_symbol = hgnc_symbols,
    auc         = as.numeric(auc_vec),
    ci_low      = as.numeric(ci_low_vec),
    ci_high     = as.numeric(ci_high_vec)
  )
  
  per_gene_df_sorted <- per_gene_df[!is.na(per_gene_df$auc), ]
  per_gene_df_sorted <- per_gene_df_sorted[order(-per_gene_df_sorted$auc), ]
  
  per_gene_df_sorted$strength <- cut(
    per_gene_df_sorted$auc,
    breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels = c("Poor (<0.6)", "Fair (0.6-0.7)", "Good (0.7-0.8)",
               "Strong (0.8-0.9)", "Excellent (>0.9)"),
    include.lowest = TRUE
  )
  
  top20 <- head(per_gene_df_sorted, 20)
  
  p_biomarker <- ggplot(top20, aes(x = auc, y = reorder(gene_symbol, auc), color = strength)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3) +
    geom_vline(xintercept = 0.5, linetype = "dashed",  color = "grey50") +
    geom_vline(xintercept = 0.7, linetype = "dotted",  color = "orange") +
    geom_vline(xintercept = 0.8, linetype = "dotted",  color = "red") +
    scale_color_manual(values = c(
      "Poor (<0.6)"      = "grey60",
      "Fair (0.6-0.7)"   = "steelblue",
      "Good (0.7-0.8)"   = "darkgreen",
      "Strong (0.8-0.9)" = "orange",
      "Excellent (>0.9)" = "red"
    )) +
    scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
    labs(
      title  = paste0("Top AD Biomarkers - ", label, "\n(LOOCV AUC adj. age + sex)"),
      x      = "AUC (95% bootstrap CI)",
      y      = "Gene",
      color  = "Biomarker Strength"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12), legend.position = "bottom")
  
  ggsave(file.path(output_dir, paste0("top_biomarkers_", label, "_", today, ".pdf")),
         p_biomarker, width = 8, height = 7, device = cairo_pdf)
  
  write.csv(head(per_gene_df_sorted, 20),
            file.path(output_dir, paste0("top20_biomarkers_", label, "_", today, ".csv")),
            row.names = FALSE)
  write.csv(per_gene_df,
            file.path(output_dir, paste0("per_gene_LOOCV_auc_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  # ── Elastic net ──────────────────────────────────────────────────────────
  pred_multi <- run_loocv_elastic_net(X, y, meta_use = meta_use, alpha_val = alpha_val)
  roc_multi  <- roc(response = y, predictor = pred_multi, quiet = TRUE)
  auc_multi  <- as.numeric(auc(roc_multi))
  boot_res   <- bootstrap_auc_ci(y, pred_multi, B = B_boot, conf = 0.95)
  
  boot_df <- data.frame(iteration = 1:B_boot, auc = boot_res$auc_boot)
  write.csv(boot_df,
            file.path(output_dir, paste0("multigene_bootstrap_auc_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  sex_numeric     <- ifelse(toupper(meta_use$Sex) == "M", 1, 0)
  cov_mat         <- cbind(Age = as.numeric(meta_use$Age), Sex = sex_numeric)
  X_full          <- cbind(X, cov_mat)
  penalty_factors <- c(rep(1, ncol(X)), rep(0, ncol(cov_mat)))
  
  cvfit_enet <- cv.glmnet(
    x = X_full, y = y, family = "binomial", alpha = alpha_val,
    nfolds = min(5, length(y)), type.measure = "auc",
    penalty.factor = penalty_factors
  )
  cvfit_lasso <- cv.glmnet(
    x = X_full, y = y, family = "binomial", alpha = 1,
    nfolds = min(5, length(y)), type.measure = "auc",
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
  
  coef_df <- data.frame(feature = rownames(coef_full)[-1], coef = as.numeric(coef_full)[-1])
  coef_df <- coef_df[coef_df$coef != 0 & !coef_df$feature %in% c("Age", "Sex"), ]
  coef_df$gene_symbol <- stable_genes$hgnc_symbol[match(coef_df$feature, stable_genes$ensembl_id)]
  coef_df$gene_symbol[is.na(coef_df$gene_symbol)] <- coef_df$feature[is.na(coef_df$gene_symbol)]
  coef_df <- coef_df[order(-abs(coef_df$coef)), ]
  
  cat("Final selected genes:", nrow(coef_df), "\n")
  cat("Selected:", paste(coef_df$gene_symbol, collapse = ", "), "\n")
  
  write.csv(coef_df,
            file.path(output_dir, paste0("elastic_net_coefficients_", label, "_", today, ".csv")),
            row.names = FALSE)
  
  # Collinearity plots
  top20_genes <- head(per_gene_df_sorted$ensembl_id[
    per_gene_df_sorted$ensembl_id %in% colnames(X)], 20)
  top20_sym <- stable_genes$hgnc_symbol[match(top20_genes, stable_genes$ensembl_id)]
  top20_sym[is.na(top20_sym)] <- top20_genes[is.na(top20_sym)]
  
  if (length(top20_genes) > 1) {
    X_top20 <- X[, top20_genes, drop = FALSE]
    colnames(X_top20) <- top20_sym
    cor_mat <- cor(X_top20, use = "pairwise.complete.obs")
    
    p_cor_heat <- ggplot(reshape2::melt(cor_mat), aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      geom_text(aes(label = round(value, 2)), size = 2.5) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                           midpoint = 0, limits = c(-1, 1)) +
      labs(title = paste0("Top gene collinearity - ", label),
           x = NULL, y = NULL, fill = "Pearson r") +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(output_dir, paste0("collinearity_heatmap_", label, "_", today, ".pdf")),
           p_cor_heat,
           width = max(6, length(top20_genes) * 0.5),
           height = max(5, length(top20_genes) * 0.5),
           device = cairo_pdf)
    
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
  
  # Summary plots
  p_gene_bar <- ggplot(per_gene_df_sorted, aes(x = reorder(gene_symbol, auc), y = auc)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.3, color = "black") +
    geom_hline(yintercept = 0.7, linetype = "dotted", color = "orange") +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "red") +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(title = paste("Per-gene LOOCV AUC (adj. age+sex) -", label),
         x = "Gene", y = "AUC (95% CI)")
  
  ggsave(file.path(output_dir, paste0("per_gene_LOOCV_auc_bar_", label, "_", today, ".png")),
         p_gene_bar, width = 6, height = 8, dpi = 300)
  
  auc_strength <- cut(auc_multi, breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
                      labels = c("Poor", "Fair", "Good", "Strong", "Excellent"),
                      include.lowest = TRUE)
  
  p_boot <- ggplot(boot_df, aes(x = auc)) +
    geom_histogram(color = "black", fill = "steelblue", bins = 30) +
    geom_vline(xintercept = auc_multi,        color = "red",     linewidth = 1.1) +
    geom_vline(xintercept = boot_res$ci_low,  color = "red",     linetype = "dashed") +
    geom_vline(xintercept = boot_res$ci_high, color = "red",     linetype = "dashed") +
    geom_vline(xintercept = 0.7,              color = "orange",  linetype = "dotted") +
    geom_vline(xintercept = 0.8,              color = "darkred", linetype = "dotted") +
    annotate("text", x = auc_multi, y = Inf,
             label = paste0("AUC=", round(auc_multi, 3), " (", auc_strength, ")"),
             vjust = 2, hjust = -0.1, color = "red", size = 3.5) +
    theme_bw(base_size = 10) +
    labs(title = paste0("Multigene LOOCV AUC (adj. age+sex) - ", label,
                        "\nMean: ", round(boot_res$mean, 3),
                        " CI [", round(boot_res$ci_low, 3),
                        ", ", round(boot_res$ci_high, 3), "]"),
         x = "AUC", y = "Count")
  
  ggsave(file.path(output_dir, paste0("multigene_bootstrap_auc_hist_", label, "_", today, ".png")),
         p_boot, width = 6, height = 4, dpi = 300)
  
  p_compare <- ggplot(per_gene_df_sorted, aes(x = auc)) +
    geom_histogram(color = "black", fill = "grey70", bins = 15) +
    geom_vline(xintercept = auc_multi,        color = "red",    linewidth = 1.1) +
    geom_vline(xintercept = boot_res$ci_low,  color = "red",    linetype = "dashed") +
    geom_vline(xintercept = boot_res$ci_high, color = "red",    linetype = "dashed") +
    geom_vline(xintercept = 0.5,              color = "grey40", linetype = "dashed") +
    annotate("text", x = auc_multi, y = Inf,
             label = paste0("Multigene\nAUC=", round(auc_multi, 3)),
             vjust = 2, hjust = -0.1, color = "red", size = 3) +
    theme_bw(base_size = 10) +
    labs(title = paste("Per-gene vs multigene AUC -", label),
         x = "AUC", y = "Number of genes")
  
  ggsave(file.path(output_dir, paste0("per_gene_vs_multigene_auc_", label, "_", today, ".png")),
         p_compare, width = 5, height = 4, dpi = 300)
  
  if (nrow(coef_df) > 0) {
    coef_df$direction <- ifelse(coef_df$coef > 0, "Higher in AD", "Lower in AD")
    p_coef <- ggplot(coef_df, aes(x = coef, y = reorder(gene_symbol, abs(coef)), fill = direction)) +
      geom_col(alpha = 0.85) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
      scale_fill_manual(values = c("Higher in AD" = "tomato", "Lower in AD" = "steelblue")) +
      labs(title = paste0("Elastic Net Selected Genes - ", label,
                          "\n(LOOCV AUC: ", round(auc_multi, 3),
                          " [", round(boot_res$ci_low, 3), ", ", round(boot_res$ci_high, 3), "]",
                          " | adj. age+sex)"),
           x = "Coefficient", y = "Gene", fill = "Direction") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom")
    
    ggsave(file.path(output_dir, paste0("elastic_net_coef_", label, "_", today, ".pdf")),
           p_coef, width = 7, height = max(4, nrow(coef_df) * 0.35), device = cairo_pdf)
  }
  
  cat("Dataset", label, "done.\n")
  cat("  Multigene LOOCV AUC:", round(auc_multi, 3), "\n")
  cat("  Bootstrap mean:", round(boot_res$mean, 3),
      "CI [", round(boot_res$ci_low, 3), ",", round(boot_res$ci_high, 3), "]\n")
  cat("  Elastic net selected", nrow(coef_df), "genes\n")
  
  list(
    per_gene_df       = per_gene_df,
    auc_multi         = auc_multi,
    boot_mean         = boot_res$mean,
    boot_ci_low       = boot_res$ci_low,
    boot_ci_high      = boot_res$ci_high,
    elastic_net_genes = coef_df
  )
}

# ── Main loop ────────────────────────────────────────────────────────────────
for (dname in names(normalized_list)) {
  cat("\n=== Dataset:", dname, "===\n")
  
  nd       <- normalized_list[[dname]]
  X        <- nd$X
  y        <- nd$y
  meta_use <- nd$metadata
  
  # Map gene IDs
  if (nd$gene_names[1] %in% stable_genes$ensembl_id) {
    ensembl_ids  <- nd$gene_names
    hgnc_symbols <- stable_genes$hgnc_symbol[match(nd$gene_names, stable_genes$ensembl_id)]
  } else if (nd$gene_names[1] %in% stable_genes$hgnc_symbol) {
    hgnc_symbols <- nd$gene_names
    ensembl_ids  <- stable_genes$ensembl_id[match(nd$gene_names, stable_genes$hgnc_symbol)]
  } else {
    warning(dname, ": Gene names unmatched, keeping original")
    hgnc_symbols <- nd$gene_names
    ensembl_ids  <- rep(NA, length(nd$gene_names))
  }
  
  stopifnot(nrow(X) == length(y))
  stopifnot(all(rownames(X) == rownames(meta_use)))
  
  # ── Run 1: All genes ─────────────────────────────────────────────────────
  cat("\n-- Run: all genes --\n")
  results_summary[[paste0(dname, "_allgenes")]] <- run_analysis(
    X           = X,
    y           = y,
    meta_use    = meta_use,
    dname       = dname,
    label       = paste0(dname, "_allgenes"),
    ensembl_ids = ensembl_ids,
    hgnc_symbols= hgnc_symbols,
    output_dir  = output_dir,
    today       = today,
    alpha_val   = alpha_val,
    B_boot      = B_boot
  )
  
  # ── Run 2: Sex chromosome genes removed ──────────────────────────────────
  sex_mask     <- grepl(sex_chr_patterns, colnames(X)) |
    grepl(sex_chr_patterns, hgnc_symbols[match(colnames(X), ensembl_ids)])
  X_nosex      <- X[, !sex_mask, drop = FALSE]
  ens_nosex    <- ensembl_ids[!sex_mask]
  hgnc_nosex   <- hgnc_symbols[!sex_mask]
  
  cat("\n-- Run: sex genes removed (", sum(sex_mask), "genes removed) --\n")
  
  if (ncol(X_nosex) == 0) {
    warning(dname, ": No genes remaining after sex chromosome removal, skipping")
  } else {
    results_summary[[paste0(dname, "_nosexgenes")]] <- run_analysis(
      X           = X_nosex,
      y           = y,
      meta_use    = meta_use,
      dname       = dname,
      label       = paste0(dname, "_nosexgenes"),
      ensembl_ids = ens_nosex,
      hgnc_symbols= hgnc_nosex,
      output_dir  = output_dir,
      today       = today,
      alpha_val   = alpha_val,
      B_boot      = B_boot
    )
  }
}

# ── Save overall summary ─────────────────────────────────────────────────────
summary_df <- do.call(rbind, lapply(names(results_summary), function(lname) {
  s <- results_summary[[lname]]
  data.frame(
    label     = lname,
    auc_multi = s$auc_multi,
    boot_mean = s$boot_mean,
    ci_low    = s$boot_ci_low,
    ci_high   = s$boot_ci_high,
    n_elastic = nrow(s$elastic_net_genes)
  )
}))

write.csv(summary_df,
          file.path(output_dir, paste0("multidataset_auc_summary_", today, ".csv")),
          row.names = FALSE)
