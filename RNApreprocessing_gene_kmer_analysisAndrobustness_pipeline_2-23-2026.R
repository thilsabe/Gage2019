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

suppressPackageStartupMessages({
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
  library(universalmotif)
  library(ggseqlogo)
  library(Biostrings)
  library(cmdstanr)
})

# Prevent dplyr from masking data.table functions
first <- data.table::first
last  <- data.table::last

# Setting memory usage per core for future.apply (adjust as needed for your HPC environment)
options(future.globals.maxSize = 10000 * 1024^2)  # 500MB per worker

# ============================================================
# 0. CONFIGURATION
# ============================================================

n_cores <- 16      # increase to 32 for full HPC runs
testrun <- FALSE   # set FALSE for full production run
handlers("txtprogressbar")

# Input paths
deg_sets_dir       <- "/cjc/data/Gage_align_bamFiles/count_DEG_analysis/WithConfounds/HOMER/overlap_DEGs"
kmer_deseq_dir     <- "/netapp/snl/scratch25/AHAllen_projects/Gage2019/kmerator_results/deseq2_results"
kmer_map_file      <- "kmer_rMATS_annotations/kmer_gene_map_from_gtf_20260222.csv"
kmer_features_file <- "kmer_features/kmer_splice_motif_position_features_20260222.csv"
gtf_file           <- "/home/thilsabeck/Documents/genomes/Reference/Homo_sapiens.GRCh38.112.gtf"

# Output paths
output_dir <- "/netapp/snl/scratch25/AHAllen_projects/Gage2019/FINAL_ALIGNER_KMER_ANALYSIS"
today      <- format(Sys.Date(), "%Y-%m-%d")

for (d in c(output_dir,
            file.path(output_dir, "figures"),
            file.path(output_dir, "tables"),
            file.path(output_dir, "models"),
            file.path(output_dir, "diagnostics"))) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}

# Known pipeline dimension labels
trimmers <- c("fastp","trimmomatic","trimgalore","untrimmed")
aligners  <- c("hisat2","star_single","star_twopass","gsnap","gsnap_2024genome")
counters  <- c("htseq","star","homer_exon","homer_exon_unique","homer_gene",
               "homer_gene_condensed","homer_gene_condensed_unique","homer_gene_unique")

# Sub-kmer sizes for motif enrichment
subkmer_sizes <- c(4L, 6L, 8L)

# Stratified background sampling cap — per trimmer, per LFC bin
# Total non-sig rows per trimmer = n_lfc_bins * bg_per_bin = 6 * 50000 = 300000
# Total across 4 trimmers = ~1.2M + all sig kmers
bg_per_bin  <- 50000L   # max non-sig kmers per LFC bin per trimmer
lfc_breaks  <- c(-Inf, -2, -1, 0, 1, 2, Inf)
n_lfc_bins  <- length(lfc_breaks) - 1L   # = 6

# ============================================================
# HPC PARALLEL CONFIG
# ============================================================

options(mc.cores = n_cores)
plan(multicore, workers = n_cores)

# ============================================================
# 1. LOAD REFERENCE FILES
# ============================================================
cat("Loading reference files...\n")
kmer_gene_map <- fread(kmer_map_file)
cat("  kmer_gene_map raw columns:", paste(names(kmer_gene_map), collapse=", "), "\n")

# Select only columns that exist
base_cols     <- c("kmer", "kmer_chrom", "kmer_start", "kmer_end", "ensembl_id")
opt_cols      <- c("gene_symbol")
keep_cols     <- c(base_cols, intersect(opt_cols, names(kmer_gene_map)))
kmer_gene_map <- kmer_gene_map[, keep_cols, with=FALSE]

# Add gene_symbol as NA if missing
if (!"gene_symbol" %in% names(kmer_gene_map)) {
  kmer_gene_map[, gene_symbol := NA_character_]
  cat("  [INFO] gene_symbol not in kmer_gene_map — filled with NA\n")
}
cat("  kmer_gene_map:", nrow(kmer_gene_map), "rows |",
    "gene_symbol present:", sum(!is.na(kmer_gene_map$gene_symbol)), "\n")

kmer_features <- fread(kmer_features_file)
cat("  kmer_features:", nrow(kmer_features), "rows\n")

# Collapse kmer_features to one row per kmer — done once globally
kmer_features_collapsed <- kmer_features[, .(
  gc_content        = mean(gc_content,    na.rm=TRUE),
  length            = .SD$length[1],
  chrom             = .SD$chrom[1],
  start             = .SD$start[1],
  end               = .SD$end[1],
  dist_to_tss       = min(dist_to_tss,    na.rm=TRUE),
  dist_to_tes       = min(dist_to_tes,    na.rm=TRUE),
  has_exact_overlap = any(has_exact_overlap),
  has_rMATS_overlap = any(has_rMATS_overlap),
  has_gene          = any(has_gene),
  n_genes_mapped    = uniqueN(ensembl_id[ensembl_id != "unassigned"]),
  ensembl_id        = ensembl_id[ensembl_id != "unassigned"][1],
  gene_symbol       = if ("gene_symbol" %in% names(.SD))
    gene_symbol[gene_symbol != "unassigned"][1]
  else NA_character_
), by=kmer]
cat("  kmer_features_collapsed:", nrow(kmer_features_collapsed), "unique kmers\n")

# ============================================================
# 2. PARSE FILENAME METADATA
# ============================================================

parse_metadata <- function(fname) {
  parts <- str_split(basename(fname), "_")[[1]]
  base  <- gsub("\\.csv$", "", basename(fname))
  
  # What dimension is VARYING across col 1 of this file
  grouping_type <- ifelse(grepl("aligners", fname), "aligner",
                          ifelse(grepl("trimmers", fname), "trimmer",
                                 ifelse(grepl("counters", fname), "counter", "unknown")))
  
  # Extract raw fixed value string from "within_<value>"
  within_idx  <- which(parts == "within")
  fixed_value <- if (length(within_idx) > 0 && within_idx + 1 <= length(parts)) {
    remaining <- paste(parts[(within_idx + 1):length(parts)], collapse="_")
    remaining <- gsub("\\.csv$", "", remaining)
    remaining <- sub("_\\d{2}-\\d{2}-\\d{4}$", "", remaining)
    remaining
  } else {
    NA_character_
  }
  
  matched_trimmer <- if (grouping_type != "trimmer") {
    trimmers[sapply(trimmers, function(t)
      grepl(paste0("(^|_)", t, "(_|$)"), fixed_value))]
  } else character(0)
  
  matched_aligner <- if (grouping_type != "aligner") {
    aligners[sapply(aligners, function(a)
      grepl(paste0("(^|_)", a, "(_|$)"), fixed_value))]
  } else character(0)
  
  matched_counter <- if (grouping_type != "counter") {
    counters[sapply(counters, function(c)
      grepl(paste0("(^|_)", c, "(_|$)"), fixed_value))]
  } else character(0)
  
  matched_trimmer <- sort(matched_trimmer, decreasing=TRUE)
  matched_aligner <- sort(matched_aligner, decreasing=TRUE)
  matched_counter <- sort(matched_counter, decreasing=TRUE)
  
  deg_list_trimmer <- if (length(matched_trimmer) > 0) matched_trimmer[1] else NA_character_
  deg_list_aligner <- if (length(matched_aligner) > 0) matched_aligner[1] else NA_character_
  deg_list_counter <- if (length(matched_counter) > 0) matched_counter[1] else NA_character_
  
  n_matched_dims <- (!is.na(deg_list_trimmer)) +
    (!is.na(deg_list_aligner)) +
    (!is.na(deg_list_counter))
  if (n_matched_dims == 0) {
    cat(sprintf("  [WARNING] parse_metadata: fixed_value '%s' did not match any known",
                fixed_value),
        "trimmer/aligner/counter in", basename(fname), "\n")
  }
  if (n_matched_dims > 1) {
    cat(sprintf("  [WARNING] parse_metadata: fixed_value '%s' matched multiple dimensions",
                fixed_value),
        sprintf("(trimmer=%s, aligner=%s, counter=%s) in %s\n",
                deg_list_trimmer, deg_list_aligner, deg_list_counter, basename(fname)))
  }
  
  fastq_trimmer <- deg_list_trimmer %||% NA_character_
  
  counter_match <- sort(counters[sapply(counters, grepl, x=base, fixed=TRUE)],
                        decreasing=TRUE)
  counter <- if (length(counter_match) > 0) counter_match[1] else "all"
  
  n_varying <- switch(grouping_type,
                      "aligner" = length(aligners),
                      "trimmer" = length(trimmers),
                      "counter" = length(counters),
                      1L
  )
  
  list(
    grouping_type    = grouping_type,
    fixed_value      = fixed_value,
    deg_list_trimmer = deg_list_trimmer,
    deg_list_aligner = deg_list_aligner,
    deg_list_counter = deg_list_counter,
    fastq_trimmer    = fastq_trimmer,
    counter          = counter,
    n_varying        = n_varying
  )
}

# ============================================================
# 3. LOAD DEG FILE (aligner = col 1, gene = last col)
# ============================================================

load_deg_aligners <- function(deg_file) {
  dt          <- fread(deg_file)
  aligner_col <- names(dt)[1]
  gene_col    <- names(dt)[ncol(dt)]
  unique(dt[, .(aligner = get(aligner_col), ensembl_id = get(gene_col))])
}

# ============================================================
# 4. MATCH KMER FILE
# ============================================================

match_kmer_file <- function(deg_file, fastq_trimmer=NULL) {
  new_name      <- sub("\\.csv$", "_AD_vs_CTRL.csv", basename(deg_file))
  new_name      <- gsub("-", "_", new_name)
  base_name     <- sub("_AD_vs_CTRL\\.csv$", ".csv",     new_name)
  filtered_name <- sub("_AD_vs_CTRL\\.csv$", ".csv.zst", new_name)
  
  trimmers_to_use <- if (is.null(fastq_trimmer) || all(is.na(fastq_trimmer))) {
    trimmers
  } else {
    fastq_trimmer
  }
  
  paths <- lapply(trimmers_to_use, function(t) {
    list(
      deseq2     = file.path(kmer_deseq_dir,
                             paste0("kmc_deseq2_",     t, "_", new_name)),
      postfilter = file.path(kmer_deseq_dir,
                             paste0("kmc_postfilter_", t, "_", base_name)),
      filtered   = file.path(kmer_deseq_dir,
                             paste0("kmc_filtered_",   t, "_", filtered_name))
    )
  })
  names(paths) <- trimmers_to_use
  return(paths)
}

# ============================================================
# 5. HELPER FUNCTIONS
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Kmer-level hypergeometric motif enrichment ----
motif_enrichment <- function(kmer_df, background_kmers, progressor=NULL) {
  if (nrow(kmer_df) == 0) return(data.table())
  
  obs       <- table(kmer_df$kmer)
  bg_counts <- table(background_kmers$kmer)
  total_bg  <- nrow(background_kmers)
  total_obs <- nrow(kmer_df)
  
  kmer_names <- names(obs)
  obs_counts <- as.numeric(obs)
  bg_matched <- as.numeric(bg_counts[kmer_names])
  bg_matched[is.na(bg_matched)] <- 1
  
  pvals <- phyper(
    q          = obs_counts - 1,
    m          = bg_matched,
    n          = total_bg - bg_matched,
    k          = total_obs,
    lower.tail = FALSE
  )
  
  if (!is.null(progressor)) progressor(amount = length(kmer_names))
  
  res_df <- data.table(
    kmer            = kmer_names,
    count           = obs_counts,
    bg_count        = bg_matched,
    expected        = total_obs * (bg_matched / total_bg),
    fold_enrichment = obs_counts / (total_obs * (bg_matched / total_bg)),
    pval            = pvals
  )
  res_df[, `-log10p` := -log10(pval + 1e-12)]
  res_df[, padj      := p.adjust(pval, method="BH")]
  setorder(res_df, pval)
  return(res_df)
}

# ---- Sub-kmer extraction ----
extract_subkmers <- function(kmers, k=6L) {
  chunks <- split(unique(kmers),
                  ceiling(seq_along(unique(kmers)) / 50000))
  rbindlist(future_lapply(chunks, function(chunk) {
    rbindlist(lapply(chunk, function(seq) {
      n <- nchar(seq) - k + 1L
      if (n < 1) return(NULL)
      data.table(parent_kmer = seq,
                 subkmer     = substring(seq, 1:n, k:nchar(seq)))
    }))
  }, future.seed=TRUE))
}

# ---- Sub-kmer hypergeometric enrichment ----
subkmer_enrichment <- function(foreground_kmers, background_kmers,
                               k=6L, progressor=NULL) {
  if (nrow(foreground_kmers) == 0) return(data.table())
  
  fg_sub <- extract_subkmers(foreground_kmers$kmer, k)
  bg_sub <- extract_subkmers(background_kmers$kmer, k)
  
  if (nrow(fg_sub) == 0 || nrow(bg_sub) == 0) return(data.table())
  
  obs       <- table(fg_sub$subkmer)
  bg_counts <- table(bg_sub$subkmer)
  total_bg  <- nrow(bg_sub)
  total_obs <- nrow(fg_sub)
  
  subkmer_names <- names(obs)
  obs_counts    <- as.numeric(obs)
  bg_matched    <- as.numeric(bg_counts[subkmer_names])
  bg_matched[is.na(bg_matched)] <- 1
  
  pvals <- phyper(
    q          = obs_counts - 1,
    m          = bg_matched,
    n          = total_bg - bg_matched,
    k          = total_obs,
    lower.tail = FALSE
  )
  
  if (!is.null(progressor)) progressor(amount = length(subkmer_names))
  
  res_df <- data.table(
    subkmer         = subkmer_names,
    k               = k,
    count           = obs_counts,
    bg_count        = bg_matched,
    expected        = total_obs * (bg_matched / total_bg),
    fold_enrichment = obs_counts / (total_obs * (bg_matched / total_bg)),
    pval            = pvals
  )
  res_df[, `-log10p` := -log10(pval + 1e-12)]
  res_df[, padj      := p.adjust(pval, method="BH")]
  setorder(res_df, pval)
  return(res_df)
}

# ---- universalmotif PWM-based enrichment ----
run_universalmotif_enrichment <- function(foreground_kmers, background_kmers,
                                          deg_out, deg_name, direction, today) {
  tryCatch({
    fg_seqs <- DNAStringSet(unique(foreground_kmers$kmer))
    bg_seqs <- DNAStringSet(unique(background_kmers$kmer))
    
    if (length(fg_seqs) < 5) {
      cat(sprintf("  [SKIP] universalmotif %s: too few unique sequences (%d)\n",
                  direction, length(fg_seqs)))
      return(invisible(NULL))
    }
    
    motifs <- tryCatch(
      universalmotif::create_motif(fg_seqs, type="PPM"),
      error = function(e) NULL
    )
    if (is.null(motifs)) return(invisible(NULL))
    
    enrich_res <- universalmotif::enrich_motifs(
      motifs,
      sequences     = fg_seqs,
      bkg.sequences = bg_seqs,
      RC            = TRUE,
      verbose       = 0
    )
    
    if (!is.null(enrich_res) && nrow(as.data.frame(enrich_res)) > 0) {
      enrich_dt <- as.data.table(as.data.frame(enrich_res))
      enrich_dt[, direction := direction]
      fwrite(enrich_dt,
             file.path(deg_out,
                       paste0("universalmotif_enrichment_", direction, "_", today, ".csv")))
      
      top_n <- min(10, nrow(enrich_dt))
      pdf(file.path(deg_out,
                    paste0("universalmotif_logo_", direction, "_", today, ".pdf")),
          width=12, height=3 * top_n)
      universalmotif::view_motifs(motifs[seq_len(top_n)])
      dev.off()
      
      cat(sprintf("  [OK] universalmotif %s: %d motifs saved\n", direction, nrow(enrich_dt)))
    }
  }, error = function(e) {
    cat(sprintf("  [WARNING] universalmotif failed (%s): %s\n", direction, e$message))
  })
}

# ---- Sequence logo plotting ----
plot_seqlogo_robust <- function(top_motifs, k_size, deg_out, today, deg_name) {
  up_motifs   <- top_motifs[direction == "up"   & k == k_size & padj < 0.05, subkmer]
  down_motifs <- top_motifs[direction == "down" & k == k_size & padj < 0.05, subkmer]
  
  has_up   <- length(up_motifs)   > 0
  has_down <- length(down_motifs) > 0
  
  if (!has_up && !has_down) {
    cat(sprintf("  [SKIP] No significant %d-mer motifs in either direction for %s\n",
                k_size, deg_name))
    return(invisible(NULL))
  }
  
  pdf(file.path(deg_out, paste0("subkmer_seqlogo_k", k_size, "_", today, ".pdf")),
      width=10, height=5 * (has_up + has_down))
  
  if (has_up) {
    tryCatch({
      print(ggseqlogo(up_motifs) +
              ggtitle(sprintf("Upregulated enriched %d-mers — %s (n=%d)",
                              k_size, deg_name, length(up_motifs))) +
              theme_classic())
    }, error = function(e) {
      cat(sprintf("  [WARNING] ggseqlogo failed (up, k=%d): %s\n", k_size, e$message))
      plot.new(); title(sprintf("ggseqlogo failed: %s", e$message))
    })
  } else {
    plot.new()
    title(sprintf("No significant upregulated %d-mers\n%s", k_size, deg_name))
  }
  
  if (has_down) {
    tryCatch({
      print(ggseqlogo(down_motifs) +
              ggtitle(sprintf("Downregulated enriched %d-mers — %s (n=%d)",
                              k_size, deg_name, length(down_motifs))) +
              theme_classic())
    }, error = function(e) {
      cat(sprintf("  [WARNING] ggseqlogo failed (down, k=%d): %s\n", k_size, e$message))
      plot.new(); title(sprintf("ggseqlogo failed: %s", e$message))
    })
  } else {
    plot.new()
    title(sprintf("No significant downregulated %d-mers\n%s", k_size, deg_name))
  }
  
  dev.off()
  cat(sprintf("  Saved seqlogo: k=%d, up=%d, down=%d significant motifs\n",
              k_size, length(up_motifs), length(down_motifs)))
}

# ---- Jaccard similarity matrix ----
compute_jaccard <- function(method_list) {
  mnames <- names(method_list)
  mat    <- matrix(0, nrow=length(mnames), ncol=length(mnames),
                   dimnames=list(mnames, mnames))
  for (i in seq_along(mnames)) {
    for (j in seq_along(mnames)) {
      a          <- method_list[[i]]
      b          <- method_list[[j]]
      inter      <- length(intersect(a, b))
      union_size <- length(unique(c(a, b)))
      mat[i, j]  <- if (union_size > 0) inter / union_size else 0
    }
  }
  return(mat)
}

# ============================================================
# GLMM ANALYSIS FUNCTION
# ============================================================
# FIX: All data.table column references inside the function use
# `dt[, col_vec, with=FALSE]` instead of `..col_vec` prefix.
# The `..` variable-lookup mechanism resolves to the *calling* frame,
# not the function frame, causing "object not found" errors inside
# functions. `with=FALSE` is the correct idiomatic replacement.
#
# FIX: Random effect is on `kmer` (unit of observation) not
# `ensembl_id`. In this context many kmers map to each gene, so
# using ensembl_id as the grouping factor collapses to effectively
# one level per gene, causing "grouping factors must have > 1
# sampled level" errors. kmer is the correct grouping factor —
# it accounts for kmer-level baseline detectability differences
# across pipeline methods.
# ============================================================
run_glmm_analysis <- function(dt, grouping_type, fixed_value,
                              deg_name, out_dir, today,
                              min_rows=100L) {
  
  if (nrow(dt) < min_rows) {
    cat(sprintf("  [SKIP] GLMM: too few rows (%d < %d)\n", nrow(dt), min_rows))
    return(invisible(NULL))
  }
  
  required_pkgs <- c("lme4", "broom.mixed")
  missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly=TRUE)]
  if (length(missing_pkgs) > 0) {
    cat(sprintf("  [SKIP] GLMM: missing packages: %s\n",
                paste(missing_pkgs, collapse=", ")))
    return(invisible(NULL))
  }
  
  suppressPackageStartupMessages({
    library(lme4)
    library(broom.mixed)
  })
  
  results <- list()
  
  vary_col <- switch(grouping_type,
                     "aligner" = "aligner",
                     "trimmer" = "deg_list_trimmer",
                     "counter" = "counter",
                     NULL
  )
  
  has_vary <- !is.null(vary_col) &&
    vary_col %in% names(dt) &&
    uniqueN(dt[[vary_col]][!is.na(dt[[vary_col]])]) > 1
  
  has_fastq_trimmer <- "fastq_trimmer" %in% names(dt) &&
    uniqueN(dt$fastq_trimmer[!is.na(dt$fastq_trimmer)]) > 1
  
  # FIX: report n unique levels, not n_rows which was always 119M
  cat(sprintf("  GLMM: vary_col=%s (n_levels=%d) | fastq_trimmer=%s (n_levels=%d) | grouping=%s\n",
              vary_col %||% "none",
              if (has_vary) uniqueN(dt[[vary_col]], na.rm=TRUE) else 0,
              ifelse(has_fastq_trimmer, "included", "excluded"),
              if (has_fastq_trimmer) uniqueN(dt$fastq_trimmer, na.rm=TRUE) else 0,
              grouping_type))
  
  # --- Model 1: Varying dimension only ---
  if (has_vary) {
    tryCatch({
      dt_model <- dt[!is.na(get(vary_col)) & !is.na(kmer) & !is.na(DEG_binary)]
      dt_model[, vary := as.factor(get(vary_col))]
      # FIX: use kmer as random effect — it is the unit of observation.
      # ensembl_id would collapse all rows of the same gene to one random intercept,
      # leaving only 1 level when a DEG set has kmers from a single gene region.
      dt_model[, kmer_f := as.factor(kmer)]
      
      fit <- glmer(
        DEG_binary ~ vary + (1 | kmer_f),
        data    = dt_model,
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))
      )
      coefs <- as.data.table(tidy(fit, effects="fixed"))
      coefs[, model       := paste0(vary_col, "_only")]
      coefs[, vary_col    := vary_col]
      coefs[, fixed_value := fixed_value]
      coefs[, deg_set     := deg_name]
      results[["vary_dim"]] <- coefs
      cat(sprintf("  GLMM %s_only: AIC=%.1f\n", vary_col, AIC(fit)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] GLMM vary_dim failed: %s\n", e$message))
    })
  }
  
  # --- Model 2: fastq_trimmer only ---
  if (has_fastq_trimmer) {
    tryCatch({
      dt_model <- dt[!is.na(fastq_trimmer) & !is.na(kmer) & !is.na(DEG_binary)]
      dt_model[, trimmer_f := as.factor(fastq_trimmer)]
      dt_model[, kmer_f    := as.factor(kmer)]
      
      fit <- glmer(
        DEG_binary ~ trimmer_f + (1 | kmer_f),
        data    = dt_model,
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))
      )
      coefs <- as.data.table(tidy(fit, effects="fixed"))
      coefs[, model       := "fastq_trimmer_only"]
      coefs[, vary_col    := "fastq_trimmer"]
      coefs[, fixed_value := fixed_value]
      coefs[, deg_set     := deg_name]
      results[["fastq_trimmer"]] <- coefs
      cat(sprintf("  GLMM fastq_trimmer_only: AIC=%.1f\n", AIC(fit)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] GLMM fastq_trimmer failed: %s\n", e$message))
    })
  }
  
  # --- Model 3: vary_dim + fastq_trimmer additive ---
  if (has_vary && has_fastq_trimmer) {
    tryCatch({
      dt_model <- dt[!is.na(get(vary_col)) & !is.na(fastq_trimmer) &
                       !is.na(kmer) & !is.na(DEG_binary)]
      dt_model[, vary      := as.factor(get(vary_col))]
      dt_model[, trimmer_f := as.factor(fastq_trimmer)]
      dt_model[, kmer_f    := as.factor(kmer)]
      
      fit <- glmer(
        DEG_binary ~ vary + trimmer_f + (1 | kmer_f),
        data    = dt_model,
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))
      )
      coefs <- as.data.table(tidy(fit, effects="fixed"))
      coefs[, model       := paste0(vary_col, "_plus_fastq_trimmer")]
      coefs[, vary_col    := vary_col]
      coefs[, fixed_value := fixed_value]
      coefs[, deg_set     := deg_name]
      results[["vary_plus_trimmer"]] <- coefs
      cat(sprintf("  GLMM %s+fastq_trimmer: AIC=%.1f\n", vary_col, AIC(fit)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] GLMM vary+trimmer additive failed: %s\n", e$message))
    })
  }
  
  # --- Model 4: vary_dim x fastq_trimmer interaction ---
  if (has_vary && has_fastq_trimmer) {
    tryCatch({
      dt_model <- dt[!is.na(get(vary_col)) & !is.na(fastq_trimmer) &
                       !is.na(kmer) & !is.na(DEG_binary)]
      dt_model[, vary      := as.factor(get(vary_col))]
      dt_model[, trimmer_f := as.factor(fastq_trimmer)]
      dt_model[, kmer_f    := as.factor(kmer)]
      
      fit <- glmer(
        DEG_binary ~ vary * trimmer_f + (1 | kmer_f),
        data    = dt_model,
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      )
      coefs <- as.data.table(tidy(fit, effects="fixed"))
      coefs[, model       := paste0(vary_col, "_x_fastq_trimmer")]
      coefs[, vary_col    := vary_col]
      coefs[, fixed_value := fixed_value]
      coefs[, deg_set     := deg_name]
      results[["vary_x_trimmer"]] <- coefs
      cat(sprintf("  GLMM %s x fastq_trimmer: AIC=%.1f\n", vary_col, AIC(fit)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] GLMM vary x fastq_trimmer interaction failed: %s\n",
                  e$message))
    })
  }
  
  # --- Model 5: kmer entropy ---
  if ("kmer_entropy" %in% names(dt) &&
      sum(!is.na(dt$kmer_entropy)) > min_rows) {
    tryCatch({
      dt_model <- dt[!is.na(kmer_entropy) & !is.na(kmer) & !is.na(DEG_binary)]
      dt_model[, kmer_f := as.factor(kmer)]
      
      fit <- glmer(
        DEG_binary ~ scale(kmer_entropy) + (1 | kmer_f),
        data    = dt_model,
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))
      )
      coefs <- as.data.table(tidy(fit, effects="fixed"))
      coefs[, model       := "entropy_vs_detection"]
      coefs[, fixed_value := fixed_value]
      coefs[, deg_set     := deg_name]
      results[["entropy"]] <- coefs
      cat(sprintf("  GLMM entropy: AIC=%.1f\n", AIC(fit)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] GLMM entropy failed: %s\n", e$message))
    })
  }
  
  # --- Model 6: GC content ---
  if ("gc_content" %in% names(dt) &&
      sum(!is.na(dt$gc_content)) > min_rows) {
    tryCatch({
      dt_model <- dt[!is.na(gc_content) & !is.na(kmer) & !is.na(DEG_binary)]
      dt_model[, kmer_f := as.factor(kmer)]
      
      fit <- glmer(
        DEG_binary ~ scale(gc_content) + (1 | kmer_f),
        data    = dt_model,
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))
      )
      coefs <- as.data.table(tidy(fit, effects="fixed"))
      coefs[, model       := "gc_content_vs_detection"]
      coefs[, fixed_value := fixed_value]
      coefs[, deg_set     := deg_name]
      results[["gc"]] <- coefs
      cat(sprintf("  GLMM gc_content: AIC=%.1f\n", AIC(fit)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] GLMM gc_content failed: %s\n", e$message))
    })
  }
  
  # --- Model 7: Full model ---
  # FIX: replace ..base_cols with explicit with=FALSE to avoid
  # data.table scoping failure inside function frames
  if (has_vary && "kmer_entropy" %in% names(dt) && "gc_content" %in% names(dt)) {
    tryCatch({
      base_col_vec <- c(vary_col, "kmer", "DEG_binary",
                        "kmer_entropy", "gc_content")
      if (has_fastq_trimmer) base_col_vec <- c(base_col_vec, "fastq_trimmer")
      base_col_vec <- intersect(base_col_vec, names(dt))
      
      # FIX: with=FALSE instead of ..base_col_vec
      dt_model <- dt[rowSums(is.na(dt[, base_col_vec, with=FALSE])) == 0]
      dt_model[, vary   := as.factor(get(vary_col))]
      dt_model[, kmer_f := as.factor(kmer)]
      
      formula_str <- if (has_fastq_trimmer) {
        dt_model[, trimmer_f := as.factor(fastq_trimmer)]
        paste0("DEG_binary ~ vary + trimmer_f + scale(kmer_entropy) +",
               " scale(gc_content) + vary:scale(kmer_entropy) +",
               " trimmer_f:scale(kmer_entropy) + (1 | kmer_f)")
      } else {
        paste0("DEG_binary ~ vary + scale(kmer_entropy) +",
               " scale(gc_content) + vary:scale(kmer_entropy) + (1 | kmer_f)")
      }
      
      fit <- glmer(
        as.formula(formula_str),
        data    = dt_model,
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      )
      coefs <- as.data.table(tidy(fit, effects="fixed"))
      coefs[, model       := "full_model"]
      coefs[, vary_col    := vary_col]
      coefs[, fixed_value := fixed_value]
      coefs[, deg_set     := deg_name]
      results[["full"]] <- coefs
      cat(sprintf("  GLMM full: AIC=%.1f | %d terms\n", AIC(fit), nrow(coefs)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] GLMM full model failed: %s\n", e$message))
    })
  }
  
  if (length(results) == 0) return(invisible(NULL))
  
  combined <- rbindlist(results, fill=TRUE)
  combined[, grouping_type := grouping_type]
  combined[!is.na(p.value) & term != "(Intercept)",
           padj_within_model := p.adjust(p.value, method="BH"),
           by=model]
  combined[!is.na(p.value) & term != "(Intercept)",
           padj_across_models := p.adjust(p.value, method="BH")]
  fwrite(combined,
         file.path(out_dir, paste0("glmm_results_", today, ".csv")))
  
  sig_results <- combined[!is.na(p.value) & term != "(Intercept)"]
  if (nrow(sig_results) > 0) {
    sig_results[, neg_log10p  := -log10(p.value + 1e-300)]
    sig_results[, significant := p.value < 0.05]
    
    p_glmm <- ggplot(sig_results,
                     aes(x=estimate, y=neg_log10p,
                         colour=significant, label=term)) +
      geom_point(size=2, alpha=0.7) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="red") +
      geom_vline(xintercept=0, linetype="dotted") +
      facet_wrap(~ model, scales="free") +
      theme_classic() +
      labs(title=paste0("GLMM results — ", deg_name),
           x="Effect size (log-odds)", y="-log10(p-value)",
           colour="p<0.05")
    ggsave(file.path(out_dir, paste0("glmm_volcano_", today, ".pdf")),
           p_glmm, width=14, height=10)
  }
  
  return(combined)
}

# ============================================================
# 6. GTF LOOKUP — load once before the loop
# ============================================================

symbol_lookup_file     <- "kmer_rMATS_annotations/gene_symbol_to_ensembl_lookup.csv"
transcript_lookup_file <- "kmer_rMATS_annotations/transcript_to_ensembl_lookup.csv"

if (file.exists(symbol_lookup_file) && file.exists(transcript_lookup_file)) {
  cat("Loading pre-saved GTF lookups...\n")
  gtf_lookup         <- fread(symbol_lookup_file)
  transcript_to_gene <- fread(transcript_lookup_file)
} else {
  cat("Parsing GTF for ID lookups:", gtf_file, "\n")
  gtf_raw <- fread(gtf_file, sep="\t", header=FALSE, skip=5,
                   col.names=c("chrom","source","feature","start","end",
                               "score","strand","frame","attributes"),
                   colClasses="character", fill=TRUE)
  gtf_raw <- gtf_raw[feature == "exon"]
  gtf_raw[, `:=`(
    ensembl_id    = sub("\\..*$", "", str_match(attributes, 'gene_id "([^"]+)"')[,2]),
    gene_symbol   = str_match(attributes, 'gene_name "([^"]+)"')[,2],
    transcript_id = sub("\\..*$", "", str_match(attributes, 'transcript_id "([^"]+)"')[,2])
  )]
  gtf_raw <- gtf_raw[grepl("^ENSG", ensembl_id)]
  
  gtf_lookup <- unique(gtf_raw[!is.na(gene_symbol), .(gene_symbol, ensembl_id)])
  gtf_lookup <- gtf_lookup[, .SD[1], by=gene_symbol]
  
  transcript_to_gene <- unique(gtf_raw[!is.na(transcript_id), .(transcript_id, ensembl_id)])
  transcript_to_gene <- transcript_to_gene[, .SD[1], by=transcript_id]
  
  fwrite(gtf_lookup,         symbol_lookup_file)
  fwrite(transcript_to_gene, transcript_lookup_file)
  cat("  Saved symbol lookup    :", nrow(gtf_lookup), "entries\n")
  cat("  Saved transcript lookup:", nrow(transcript_to_gene), "entries\n")
  rm(gtf_raw)
}

# ============================================================
# 7. MAIN LOADING + PER-DEG-SET ANALYSIS LOOP
# ============================================================

priority_deg_lists <- c(
  "all_DEGs_aligners_within_fastp_11-11-2025",
  "all_DEGs_aligners_within_trimmomatic_11-11-2025",
  "all_DEGs_aligners_within_trimgalore_11-11-2025",
  "all_DEGs_aligners_within_fastp_homer_gene_11-11-2025",
  "all_DEGs_aligners_within_fastp_htseq_11-11-2025",
  "all_DEGs_aligners_within_fastp_homer_exon_11-11-2025",
  "all_DEGs_aligners_within_trimmomatic_homer_gene_11-11-2025",
  "all_DEGs_aligners_within_trimmomatic_htseq_11-11-2025",
  "all_DEGs_aligners_within_trimmomatic_homer_exon_11-11-2025",
  "all_DEGs_aligners_within_trimgalore_homer_gene_11-11-2025",
  "all_DEGs_aligners_within_trimgalore_htseq_11-11-2025",
  "all_DEGs_aligners_within_trimgalore_homer_exon_11-11-2025",
  "all_DEGs_trimmers_within_gsnap_11-11-2025",
  "all_DEGs_trimmers_within_gsnap_2024genome_11-11-2025",
  "all_DEGs_trimmers_within_hisat2_11-11-2025",
  "all_DEGs_trimmers_within_star_single_11-11-2025",
  "all_DEGs_trimmers_within_star_twopass_11-11-2025"
)

deg_files <- list.files(deg_sets_dir, pattern="^all_DEGs_.*\\.csv$", full.names=TRUE)
cat("Found", length(deg_files), "DEG files\n")

priority_deg_files <- paste0(priority_deg_lists, ".csv")
deg_files <- deg_files[basename(deg_files) %in% priority_deg_files]
cat("After priority filter:", length(deg_files), "DEG files\n")

if (testrun) {
  deg_files <- deg_files[1:3]
  cat("TEST RUN with", length(deg_files), "DEG files\n")
}

# Accumulator lists
full_dataset           <- list()
gene_tracker           <- list()
all_summaries          <- list()
all_deg_dt             <- list()
kmers_by_method        <- list()
all_kmer_gene_lists    <- list()
all_kmer_aligner_lists <- list()
fastq_sensitivity_all  <- list()
motif_results_all      <- list()
subkmer_results_all    <- list()
glmm_results_all       <- list()

for (deg_file in deg_files) {
  
  meta     <- parse_metadata(deg_file)
  meta     <- lapply(meta, function(x) if (is.null(x)) NA_character_ else x)
  deg_name <- gsub("\\.csv$", "", basename(deg_file))
  
  cat("\n============================================================\n")
  cat("Processing:", deg_name, "\n")
  cat("  fastq_trimmer   :", meta$fastq_trimmer    %||% "NA", "\n")
  cat("  deg_list_trimmer:", meta$deg_list_trimmer  %||% "NA", "\n")
  cat("  deg_list_aligner:", meta$deg_list_aligner  %||% "NA", "\n")
  cat("  fixed_value     :", meta$fixed_value       %||% "NA", "\n")
  cat("  counter         :", meta$counter,           "\n")
  cat("  grouping_type   :", meta$grouping_type,     "\n")
  
  # ----------------------------------------------------------
  # Load DEG file + map IDs to ENSG
  # ----------------------------------------------------------
  deg_dt <- load_deg_aligners(deg_file)
  
  deg_dt[, id_type := fcase(
    grepl("^ENSG", ensembl_id), "ensembl_gene",
    grepl("^ENST", ensembl_id), "ensembl_transcript",
    default =                    "gene_symbol"
  )]
  cat("\nID type breakdown:\n")
  print(deg_dt[, .N, by=id_type])
  
  deg_dt[id_type == "ensembl_transcript",
         ensembl_id := transcript_to_gene$ensembl_id[
           match(ensembl_id, transcript_to_gene$transcript_id)]]
  deg_dt[id_type == "gene_symbol",
         ensembl_id := gtf_lookup$ensembl_id[
           match(ensembl_id, gtf_lookup$gene_symbol)]]
  deg_dt[, id_type := NULL]
  
  cat("  ENSG IDs after mapping:", sum(grepl("^ENSG", deg_dt$ensembl_id), na.rm=TRUE), "\n")
  cat("  Unresolved/NA         :", sum(!grepl("^ENSG", deg_dt$ensembl_id) |
                                         is.na(deg_dt$ensembl_id), na.rm=TRUE), "\n")
  
  # ----------------------------------------------------------
  # Store raw DEG gene list
  # ----------------------------------------------------------
  all_deg_dt[[deg_name]] <- deg_dt[, .(
    aligner          = aligner,
    ensembl_id       = ensembl_id,
    gene_symbol      = if ("gene_symbol" %in% names(.SD)) gene_symbol else NA_character_,
    deg_name         = deg_name,
    deg_list_trimmer = meta$deg_list_trimmer,
    counter          = meta$counter,
    grouping_type    = meta$grouping_type
  )]
  
  # ----------------------------------------------------------
  # KMER LOADING — STRATIFIED SAMPLING
  # FIX: replaces the previous flat fread of all 5M rows per
  # trimmer. The cartesian join of 4 × 5M rows against
  # kmer_gene_map exploded to ~120M rows and killed the process
  # (OOM, ~184GB RSS). Strategy:
  #   1. Keep ALL significant kmers (padj < 0.05) — these are
  #      the foreground for motif enrichment and GLMM DEG_binary=1.
  #   2. Stratified background: sample non-sig kmers across 6
  #      log2FoldChange bins so the full effect-size distribution
  #      is represented. This ensures GLMM has variance in
  #      DEG_binary (not just 1s) without loading all 5M rows.
  #   3. gc() after each trimmer to return memory before the next.
  # Total per DEG set: ~300K bg + n_sig << 20M original.
  # ----------------------------------------------------------
  kmer_dt_list <- list()
  trimmers_to_try <- trimmers
  
  cat("  Expected deseq2 path pattern:\n")
  for (t in trimmers) {
    test_path <- match_kmer_file(deg_file, t)
    cat(sprintf("    %s:\n      deseq2: %s\n", t, test_path[[t]]$deseq2))
  }
  cat("  Actual files in deseq2_results matching this deg:\n")
  deg_name_safe <- gsub("-", "_", gsub("\\.csv$", "", basename(deg_file)))
  actual <- Sys.glob(file.path(kmer_deseq_dir,
                               paste0("*", deg_name_safe, "*")))
  for (f in actual) cat(sprintf("    %s\n", f))
  
  kmer_paths <- match_kmer_file(deg_file, trimmers_to_try)
  
  for (fastq_trimmer in names(kmer_paths)) {
    paths <- kmer_paths[[fastq_trimmer]]
    
    kmer_file_match <- Sys.glob(sub("\\.csv$", "*.csv", paths$deseq2))
    
    if (length(kmer_file_match) == 0) {
      cat(sprintf("  [INFO] deseq2 not found for trimmer=%s, trying postfilter\n",
                  fastq_trimmer))
      kmer_file_match <- Sys.glob(sub("\\.csv$", "*.csv", paths$postfilter))
    }
    
    if (length(kmer_file_match) == 0) {
      cat(sprintf("  [SKIP] No kmer files found for trimmer=%s\n", fastq_trimmer))
      next
    }
    
    tryCatch({
      dt <- fread(kmer_file_match[1])
      dt[, fastq_trimmer    := fastq_trimmer]
      dt[, kmer_source_file := basename(kmer_file_match[1])]
      dt[, kmer_source_type := ifelse(grepl("deseq2",    basename(kmer_file_match[1])),
                                      "deseq2",
                                      ifelse(grepl("postfilter", basename(kmer_file_match[1])),
                                             "postfilter", "other"))]
      
      n_total <- nrow(dt)
      n_sig   <- sum(!is.na(dt$padj) & dt$padj < 0.05)
      
      # --- Keep ALL significant kmers ---
      sig_dt <- dt[!is.na(padj) & padj < 0.05]
      
      # --- Stratified background sample ---
      # Sample non-sig kmers proportionally across LFC bins so the
      # background represents the full effect-size distribution.
      # This preserves power for method comparison while capping memory.
      bg_dt <- dt[is.na(padj) | padj >= 0.05]
      
      if (nrow(bg_dt) > (bg_per_bin * n_lfc_bins)) {
        bg_dt[, lfc_bin := cut(log2FoldChange,
                               breaks = lfc_breaks,
                               labels = FALSE)]
        bg_dt[is.na(lfc_bin), lfc_bin := 3L]   # NA LFC → middle bin
        bg_sample <- bg_dt[, .SD[sample(.N, min(.N, bg_per_bin))], by=lfc_bin]
        bg_sample[, lfc_bin := NULL]
      } else {
        bg_sample <- bg_dt
      }
      
      dt_final <- rbindlist(list(sig_dt, bg_sample), fill=TRUE)
      cat(sprintf("  [OK] %s (%s): %d sig + %d bg sampled (from %d total, %.1f%% retained)\n",
                  fastq_trimmer,
                  dt$kmer_source_type[1],
                  nrow(sig_dt),
                  nrow(bg_sample),
                  n_total,
                  100 * nrow(dt_final) / n_total))
      
      kmer_dt_list[[fastq_trimmer]] <- dt_final
      rm(dt, sig_dt, bg_dt, bg_sample)
      gc()
      
    }, error = function(e) {
      cat(sprintf("  [WARNING] Failed to load kmer file for trimmer=%s: %s\n",
                  fastq_trimmer, e$message))
    })
  }
  
  if (length(kmer_dt_list) == 0) {
    cat(sprintf("  [SKIP] No kmer files found for %s — continuing to next DEG set\n",
                deg_name))
    next
  }
  
  cat(sprintf("  Loaded kmer data for %d/%d trimmer(s): %s\n",
              length(kmer_dt_list), length(trimmers),
              paste(names(kmer_dt_list), collapse=", ")))
  
  kmer_dt <- rbindlist(kmer_dt_list, fill=TRUE)
  rm(kmer_dt_list)
  gc()
  
  # ----------------------------------------------------------
  # Compute sequence features for all kmers
  # ----------------------------------------------------------
  all_kmers_dt <- data.table(kmer = unique(kmer_dt$kmer))
  all_kmers_dt[, `:=`(
    freq_A = stri_count_fixed(kmer, "A") / nchar(kmer),
    freq_C = stri_count_fixed(kmer, "C") / nchar(kmer),
    freq_G = stri_count_fixed(kmer, "G") / nchar(kmer),
    freq_T = stri_count_fixed(kmer, "T") / nchar(kmer)
  )]
  all_kmers_dt[, kmer_entropy := {
    h <- function(p) ifelse(p > 0, -p * log2(p), 0)
    h(freq_A) + h(freq_C) + h(freq_G) + h(freq_T)
  }]
  all_kmers_dt[, c("freq_A","freq_C","freq_G","freq_T") := NULL]
  all_kmers_dt[, low_complexity  := kmer_entropy < 1.5 |
                 grepl("A{10,}|C{10,}|G{10,}|T{10,}", kmer)]
  all_kmers_dt[, has_exon_overlap := kmer %in% kmer_gene_map$kmer]
  
  gene_counts  <- kmer_gene_map[, .(n_genes_mapped = uniqueN(ensembl_id)), by=kmer]
  all_kmers_dt <- merge(all_kmers_dt, gene_counts, by="kmer", all.x=TRUE)
  all_kmers_dt[is.na(n_genes_mapped), n_genes_mapped := 0L]
  
  kmer_coords_dt <- unique(kmer_gene_map[, .(kmer, kmer_chrom, kmer_start, kmer_end)],
                           by="kmer")
  kmer_pos_dt    <- kmer_features_collapsed[, .(kmer, dist_to_tss, dist_to_tes,
                                                gc_content, length)]
  all_kmers_dt   <- merge(all_kmers_dt, kmer_coords_dt, by="kmer", all.x=TRUE)
  all_kmers_dt   <- merge(all_kmers_dt, kmer_pos_dt,    by="kmer", all.x=TRUE)
  
  # ----------------------------------------------------------
  # Merge sequence features into kmer_dt
  # ----------------------------------------------------------
  kmer_dt <- merge(kmer_dt, all_kmers_dt, by="kmer", all.x=TRUE)
  rm(all_kmers_dt)
  gc()
  
  # ----------------------------------------------------------
  # CARTESIAN GENE JOIN — deferred and bounded
  # FIX: the original code joined ALL ~20M kmer rows against
  # kmer_gene_map (many-to-many), producing ~120M rows and
  # causing OOM kill at ~184GB RSS. The fix:
  #   1. Compute pre-join size estimate and warn if large.
  #   2. If expected rows > 10M, fall back to one-gene-per-kmer
  #      (primary gene) to cap memory.
  #   3. Full cartesian only when the join is small enough.
  # ----------------------------------------------------------
  n_kmers_to_join  <- uniqueN(kmer_dt$kmer)
  genes_per_kmer_n <- kmer_gene_map[kmer %in% kmer_dt$kmer,
                                    .N, by=kmer][, mean(N)]
  expected_join_rows <- n_kmers_to_join * genes_per_kmer_n
  cat(sprintf("  Pre-join size estimate: %d kmers × %.1f genes/kmer = ~%.0fM rows\n",
              n_kmers_to_join, genes_per_kmer_n, expected_join_rows / 1e6))
  
  if (expected_join_rows > 1e7) {
    cat("  [INFO] Large join expected — using primary gene per kmer to cap memory\n")
    kmer_gene_map_join <- kmer_gene_map[
      kmer %in% kmer_dt$kmer,
      .SD[which.max(!is.na(ensembl_id))],
      by=kmer]
  } else {
    kmer_gene_map_join <- kmer_gene_map[kmer %in% kmer_dt$kmer]
  }
  
  kmer_dt <- merge(kmer_dt,
                   kmer_gene_map_join[, .(kmer, ensembl_id, gene_symbol)],
                   by="kmer", all.x=TRUE,
                   allow.cartesian=(expected_join_rows <= 1e7))
  rm(kmer_gene_map_join)
  gc()
  
  # Resolve .x/.y conflicts — kmer_gene_map is authoritative
  if ("ensembl_id.y" %in% names(kmer_dt)) {
    kmer_dt[, ensembl_id  := fifelse(!is.na(ensembl_id.y),  ensembl_id.y,  ensembl_id.x)]
    kmer_dt[, gene_symbol := fifelse(!is.na(gene_symbol.y), gene_symbol.y, gene_symbol.x)]
    kmer_dt[, c("ensembl_id.x","ensembl_id.y","gene_symbol.x","gene_symbol.y") := NULL]
  }
  
  kmer_dt[, annotation_status := fcase(
    !is.na(ensembl_id) & ensembl_id != "unassigned", "gene_annotated",
    low_complexity == TRUE,                            "low_complexity",
    default =                                          "intergenic_unknown"
  )]
  
  cat("\nAnnotation status:\n")
  print(kmer_dt[, .N, by=annotation_status])
  cat("  gc_content non-NA  :", sum(!is.na(kmer_dt$gc_content)), "\n")
  cat("  kmer_start non-NA  :", sum(!is.na(kmer_dt$kmer_start)), "\n")
  cat("  dist_to_tss non-NA :", sum(!is.na(kmer_dt$dist_to_tss)), "\n")
  cat("  ensembl_id non-NA  :", sum(!is.na(kmer_dt$ensembl_id)), "\n")
  
  # ----------------------------------------------------------
  # Cross-reference with DEG genes per aligner → merged
  # ----------------------------------------------------------
  cat("  kmer_dt rows before DEG merge:", nrow(kmer_dt), "\n")
  cat("  deg_dt rows:", nrow(deg_dt), "\n")
  cat("  unique ensembl_id in kmer_dt:", uniqueN(kmer_dt$ensembl_id), "\n")
  cat("  unique ensembl_id in deg_dt :", uniqueN(deg_dt$ensembl_id), "\n")
  cat("  overlapping ensembl_ids     :",
      length(intersect(kmer_dt$ensembl_id, deg_dt$ensembl_id)), "\n")
  
  kmer_dt[, is_deg_gene  := ensembl_id %in% deg_dt$ensembl_id]
  kmer_dt[, n_aligners_called := {
    deg_dt[match(ensembl_id, deg_dt$ensembl_id), .N]
  }]
  
  cat("  kmers in DEG genes:", sum(kmer_dt$is_deg_gene, na.rm=TRUE), "\n")
  
  merged <- merge(
    kmer_dt[!is.na(ensembl_id) & is_deg_gene == TRUE],
    unique(deg_dt[, .(ensembl_id, aligner)]),
    by = "ensembl_id", allow.cartesian = TRUE
  )
  
  cat("  merged rows after aligner join:", nrow(merged), "\n")
  cat("  Full kmer_dt rows     :", nrow(kmer_dt), "\n")
  cat("  DEG-overlapping rows  :", nrow(merged), "\n")
  cat("  Unannotated kmers     :", sum(kmer_dt$annotation_status == "intergenic_unknown"), "\n")
  cat("  Low complexity kmers  :", sum(kmer_dt$annotation_status == "low_complexity"), "\n")
  
  if (nrow(merged) == 0) {
    cat("  [SKIP] No rows after merging kmer_dt with deg_dt for", deg_name, "\n")
    next
  }
  
  cat("  Assigning metadata to merged...\n")
  cat("  nrow(merged) before := :", nrow(merged), "\n")
  cat("  meta$fixed_value:", meta$fixed_value %||% "NULL", "\n")
  
  merged[, `:=`(
    grouping_type    = meta$grouping_type,
    fixed_value      = meta$fixed_value,
    deg_list_trimmer = meta$deg_list_trimmer,
    deg_list_aligner = meta$deg_list_aligner,
    deg_list_counter = meta$deg_list_counter,
    counter          = meta$counter,
    method           = paste(meta$grouping_type,
                             meta$fixed_value %||% "unknown",
                             meta$counter     %||% "all",
                             sep="_"),
    DEG_binary       = 1L,
    gene_specific    = ensembl_id %in% deg_dt$ensembl_id
  )]
  
  cat("  fixed_value assigned:", unique(merged$fixed_value), "\n")
  cat("  method assigned     :", unique(merged$method), "\n")
  cat("  NA methods:", sum(is.na(merged$method)), "\n")
  cat("  Unique methods:", uniqueN(merged$method), "\n")
  cat("  nrow merged:", nrow(merged), "\n")
  cat("  length unique kmer:", length(unique(merged$kmer)), "\n")
  
  merged[is.na(method), method := paste(meta$grouping_type,
                                        meta$fixed_value %||% "unknown",
                                        meta$counter     %||% "all",
                                        sep="_")]
  
  full_dataset[[deg_name]]           <- merged
  gene_tracker[[deg_name]]           <- unique(merged$ensembl_id)
  
  merged[, upset_group := paste(
    ifelse(!is.na(aligner)       & aligner       != "", aligner,       ""),
    ifelse(!is.na(fastq_trimmer) & fastq_trimmer != "", fastq_trimmer, ""),
    ifelse(!is.na(counter)       & counter       != "", counter,       ""),
    sep="_"
  )]
  merged[, upset_group := gsub("_+", "_", gsub("^_|_$", "", upset_group))]
  
  kmers_by_method[[deg_name]] <- lapply(
    split(merged$kmer, merged$upset_group), unique)
  kmers_by_method[[deg_name]] <- kmers_by_method[[deg_name]][
    sapply(kmers_by_method[[deg_name]], length) > 0]
  
  all_kmer_gene_lists[[deg_name]] <- lapply(
    split(merged$ensembl_id, merged$upset_group), unique)
  
  all_kmer_aligner_lists[[deg_name]] <- lapply(
    split(merged$ensembl_id, merged$aligner[!is.na(merged$aligner)]), unique)
  
  cat(sprintf("  Method sets for UpSet/discordant: %d groups\n",
              length(kmers_by_method[[deg_name]])))
  cat(sprintf("  Groups: %s\n",
              paste(names(kmers_by_method[[deg_name]]), collapse=", ")))
  
  all_kmer_aligner_lists[[deg_name]] <- lapply(
    split(merged$ensembl_id, merged$aligner), unique)
  
  deg_out <- file.path(output_dir, deg_name)
  dir.create(deg_out, recursive=TRUE, showWarnings=FALSE)
  
  kmers_dt <- merged
  
  n_varying <- fcase(
    meta$grouping_type == "aligner", length(aligners),
    meta$grouping_type == "trimmer", length(trimmers),
    meta$grouping_type == "counter", length(counters),
    default = 1L
  )
  
  # -------------------------------------------------------
  # 6a. FASTQ-TRIMMER SENSITIVITY vs DEG-LIST TRIMMER
  # -------------------------------------------------------
  if (meta$grouping_type == "aligner") {
    fastq_sens <- kmers_dt[, .(
      n_genes = uniqueN(ensembl_id),
      n_kmers = uniqueN(kmer),
      n_rows  = .N
    ), by=.(fixed_value, aligner)]
    fastq_sens[, genes_norm    := n_genes / max(n_genes, na.rm=TRUE)]
    fastq_sens[, kmers_norm    := n_kmers / max(n_kmers, na.rm=TRUE)]
    fastq_sens[, varying_value := aligner]
    fastq_sens[, grouping_type := "aligner"]
    fwrite(fastq_sens,
           file.path(deg_out, paste0("fastq_trimmer_sensitivity_vs_deglist_trimmer_",
                                     today, ".csv")))
    p_sens <- ggplot(fastq_sens, aes(x=aligner, y=genes_norm, fill=aligner)) +
      geom_col() + facet_wrap(~ fixed_value) + theme_classic() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      labs(title=paste0("Aligner sensitivity (fixed trimmer: ",
                        meta$fixed_value, ")\n", deg_name),
           y="Normalised gene recovery", x="Aligner")
    ggsave(file.path(deg_out, paste0("fastq_trimmer_sensitivity_", today, ".pdf")),
           p_sens, width=12, height=8)
    fastq_sensitivity_all[[deg_name]] <- fastq_sens
    
  } else if (meta$grouping_type == "trimmer") {
    trimmer_sens <- kmers_dt[, .(
      n_genes = uniqueN(ensembl_id),
      n_kmers = uniqueN(kmer),
      n_rows  = .N
    ), by=.(fixed_value, fastq_trimmer)]
    trimmer_sens[, genes_norm    := n_genes / max(n_genes, na.rm=TRUE)]
    trimmer_sens[, kmers_norm    := n_kmers / max(n_kmers, na.rm=TRUE)]
    trimmer_sens[, varying_value := fastq_trimmer]
    trimmer_sens[, grouping_type := "trimmer"]
    fwrite(trimmer_sens,
           file.path(deg_out, paste0("trimmer_sensitivity_vs_fixed_aligner_",
                                     today, ".csv")))
    p_sens <- ggplot(trimmer_sens, aes(x=fastq_trimmer, y=genes_norm, fill=fastq_trimmer)) +
      geom_col() + facet_wrap(~ fixed_value) + theme_classic() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      labs(title=paste0("Trimmer sensitivity (fixed aligner: ",
                        meta$deg_list_aligner, ")\n", deg_name),
           y="Normalised gene recovery", x="Trimmer")
    ggsave(file.path(deg_out, paste0("trimmer_sensitivity_", today, ".pdf")),
           p_sens, width=12, height=8)
    fastq_sensitivity_all[[deg_name]] <- trimmer_sens
    
  } else if (meta$grouping_type == "counter") {
    counter_sens <- kmers_dt[, .(
      n_genes = uniqueN(ensembl_id),
      n_kmers = uniqueN(kmer),
      n_rows  = .N
    ), by=.(fixed_value, counter)]
    counter_sens[, genes_norm    := n_genes / max(n_genes, na.rm=TRUE)]
    counter_sens[, kmers_norm    := n_kmers / max(n_kmers, na.rm=TRUE)]
    counter_sens[, varying_value := counter]
    counter_sens[, grouping_type := "counter"]
    fwrite(counter_sens,
           file.path(deg_out, paste0("counter_sensitivity_vs_fixed_dim_",
                                     today, ".csv")))
    fastq_sensitivity_all[[deg_name]] <- counter_sens
  }
  
  # -------------------------------------------------------
  # 6b. KMER STABILITY ACROSS FASTQ TRIMMERS
  # -------------------------------------------------------
  gene_symbol_lookup <- if ("gene_symbol" %in% names(kmers_dt)) {
    unique(kmers_dt[!is.na(gene_symbol) & gene_symbol != "unassigned",
                    .(ensembl_id, gene_symbol)])
  } else {
    data.table(ensembl_id=character(0), gene_symbol=character(0))
  }
  gene_symbol_lookup <- gene_symbol_lookup[, .(gene_symbol = gene_symbol[1]), by=ensembl_id]
  
  varying_col <- meta$grouping_type
  
  gene_stability <- kmers_dt[, .(
    n_varying_detected = uniqueN(get(varying_col))
  ), by=.(ensembl_id, fixed_value, counter)]
  gene_stability[, stability_index := n_varying_detected / meta$n_varying]
  
  gene_stability <- merge(gene_stability, gene_symbol_lookup, by="ensembl_id", all.x=TRUE)
  setcolorder(gene_stability, c("ensembl_id", "gene_symbol",
                                "fixed_value", "counter",
                                "n_varying_detected", "stability_index"))
  fwrite(gene_stability,
         file.path(deg_out, paste0("gene_stability_across_fastq_", today, ".csv")))
  
  kmer_stability <- kmers_dt[, .(
    n_varying_detected = uniqueN(get(varying_col))
  ), by=.(kmer, ensembl_id, fixed_value, counter)]
  kmer_stability[, stability_index := n_varying_detected / meta$n_varying]
  kmer_stability <- merge(kmer_stability, gene_symbol_lookup, by="ensembl_id", all.x=TRUE)
  setcolorder(kmer_stability, c("kmer", "ensembl_id", "gene_symbol",
                                "fixed_value", "counter",
                                "n_varying_detected", "stability_index"))
  fwrite(kmer_stability,
         file.path(deg_out, paste0("kmer_stability_across_fastq_", today, ".csv")))
  
  kmer_stability_per_gene <- kmer_stability[, .(
    n_stable_kmers      = sum(stability_index >= 1.0),
    n_partial_kmers     = sum(stability_index > 0 & stability_index < 1.0),
    n_total_kmers       = .N,
    mean_kmer_stability = mean(stability_index),
    max_kmer_stability  = max(stability_index),
    most_stable_kmer    = kmer[which.max(stability_index)],
    stable_kmers_list   = paste(kmer[stability_index >= 1.0], collapse="|")
  ), by=.(ensembl_id, fixed_value, counter)]
  
  gene_kmer_stability_merged <- merge(
    gene_stability,
    kmer_stability_per_gene,
    by    = c("ensembl_id", "fixed_value", "counter"),
    all.x = TRUE
  )
  
  gene_kmer_stability_merged[is.na(n_total_kmers),       n_total_kmers      := 0L]
  gene_kmer_stability_merged[is.na(n_stable_kmers),      n_stable_kmers     := 0L]
  gene_kmer_stability_merged[is.na(n_partial_kmers),     n_partial_kmers    := 0L]
  gene_kmer_stability_merged[is.na(mean_kmer_stability), mean_kmer_stability := 0]
  gene_kmer_stability_merged[is.na(max_kmer_stability),  max_kmer_stability  := 0]
  gene_kmer_stability_merged[is.na(stable_kmers_list),   stable_kmers_list  := ""]
  
  gene_kmer_stability_merged[, joint_stable         := stability_index >= 1.0 & n_stable_kmers > 0]
  gene_kmer_stability_merged[, joint_stability_score := stability_index * mean_kmer_stability]
  
  desired_cols <- c("ensembl_id", "gene_symbol",
                    "fixed_value", "counter",
                    "stability_index", "n_varying_detected",
                    "mean_kmer_stability", "max_kmer_stability",
                    "n_stable_kmers", "n_partial_kmers", "n_total_kmers",
                    "joint_stable", "joint_stability_score",
                    "most_stable_kmer", "stable_kmers_list")
  desired_cols <- desired_cols[desired_cols %in% names(gene_kmer_stability_merged)]
  setcolorder(gene_kmer_stability_merged, desired_cols)
  
  setorder(gene_kmer_stability_merged, -joint_stability_score)
  fwrite(gene_kmer_stability_merged,
         file.path(deg_out, paste0("gene_kmer_joint_stability_", today, ".csv")))
  
  cat(sprintf("  Joint stable genes (gene + kmer fully stable): %d\n",
              sum(gene_kmer_stability_merged$joint_stable, na.rm=TRUE)))
  
  p_stab <- ggplot(gene_stability,
                   aes(x=stability_index, fill=fixed_value)) +
    geom_histogram(binwidth=1/meta$n_varying, boundary=0, colour="white") +
    facet_wrap(~ counter, scales="free_y") +
    theme_classic() +
    labs(title=paste0("Gene stability across ", meta$grouping_type, "s\n", deg_name),
         x=paste0("Stability index (fraction of ", meta$grouping_type, "s detected)"),
         y="Number of genes",
         fill=paste0("Fixed ", switch(meta$grouping_type,
                                      "aligner" = "trimmer",
                                      "trimmer" = "aligner",
                                      "counter" = "fixed dim")))
  ggsave(file.path(deg_out, paste0("gene_stability_histogram_", today, ".pdf")),
         p_stab, width=10, height=6)
  
  # -------------------------------------------------------
  # 6c. MOTIF ENRICHMENT
  # -------------------------------------------------------
  n_sig <- sum(kmer_dt$padj < 0.05, na.rm=TRUE)
  cat(sprintf("\n  Significant kmers (padj<0.05): %d\n", n_sig))
  
  if (n_sig > 0) {
    foreground_up   <- kmer_dt[padj < 0.05 & log2FoldChange > 0]
    foreground_down <- kmer_dt[padj < 0.05 & log2FoldChange < 0]
    
    # Background: prefer postfilter file; fall back to non-sig rows in kmer_dt
    paths <- kmer_paths[[names(kmer_dt_list)[1]]]
    if (!exists("paths") || is.null(paths)) paths <- list()
    
    if (!is.null(paths$postfilter) && file.exists(paths$postfilter) &&
        !is.null(paths$deseq2)     && file.exists(paths$deseq2)) {
      cat("  Loading DESeq2 script outputs as background controls...\n")
      deseq2_res_bg    <- fread(paths$deseq2)
      background_kmers <- deseq2_res_bg[is.na(padj) | padj >= 0.05, kmer]
      postfilter_wide  <- fread(paths$postfilter)
      kmer_bg          <- postfilter_wide[kmer %in% background_kmers, .(kmer)]
      cat(sprintf("  Background: %d non-sig k-mers (postfilter, padj>=0.05 or NA)\n",
                  nrow(kmer_bg)))
    } else {
      cat("  [WARN] DESeq2 background outputs not found — using non-sig rows from kmer_dt\n")
      kmer_bg <- kmer_dt[is.na(padj) | padj >= 0.05, .(kmer)]
    }
    
    cat(sprintf("  Motif enrichment: %s (up=%d, down=%d kmers, bg=%d kmers)\n",
                deg_name, nrow(foreground_up), nrow(foreground_down), nrow(kmer_bg)))
    
    with_progress({
      p <- progressor(steps = 2)
      motif_res_up   <- motif_enrichment(foreground_up,   kmer_bg, progressor=p)
      motif_res_down <- motif_enrichment(foreground_down, kmer_bg, progressor=p)
    })
    motif_res_up[,   direction := "up"]
    motif_res_down[, direction := "down"]
    motif_res <- rbindlist(list(motif_res_up, motif_res_down), fill=TRUE)
    motif_res[, deg_set := deg_name]
    fwrite(motif_res,
           file.path(deg_out, paste0("motif_enrichment_", today, ".csv")))
    motif_results_all[[deg_name]] <- motif_res
    
    cat(sprintf("  Sub-kmer enrichment (k=%s) in parallel...\n",
                paste(subkmer_sizes, collapse=",")))
    
    with_progress({
      p <- progressor(steps = length(subkmer_sizes) * 2)
      
      subkmer_res_list <- future_lapply(subkmer_sizes, function(k_size) {
        p(message = sprintf("Sub-kmer k=%d up", k_size))
        sub_up   <- subkmer_enrichment(foreground_up,   kmer_bg, k=k_size)
        p(message = sprintf("Sub-kmer k=%d down", k_size))
        sub_down <- subkmer_enrichment(foreground_down, kmer_bg, k=k_size)
        if (nrow(sub_up)   > 0) sub_up[,   direction := "up"]
        if (nrow(sub_down) > 0) sub_down[, direction := "down"]
        rbindlist(list(sub_up, sub_down), fill=TRUE)
      }, future.seed=TRUE)
    })
    
    names(subkmer_res_list) <- as.character(subkmer_sizes)
    subkmer_results <- rbindlist(subkmer_res_list, fill=TRUE)
    subkmer_results[, deg_set := deg_name]
    fwrite(subkmer_results,
           file.path(deg_out, paste0("subkmer_motif_enrichment_", today, ".csv")))
    subkmer_results_all[[deg_name]] <- subkmer_results
    
    top_motifs <- subkmer_results[padj < 0.05][
      order(pval), .SD[1:min(.N, 20)], by=.(k, direction)]
    for (k_size in subkmer_sizes) {
      plot_seqlogo_robust(top_motifs, k_size, deg_out, today, deg_name)
    }
    
    cat("  Running universalmotif enrichment...\n")
    run_universalmotif_enrichment(foreground_up,   kmer_bg, deg_out, deg_name, "up",   today)
    run_universalmotif_enrichment(foreground_down, kmer_bg, deg_out, deg_name, "down", today)
    
  } else {
    cat("  [SKIP] No significant kmers for motif enrichment\n")
  }
  
  # -------------------------------------------------------
  # 6c-ii. METHOD-DISCORDANT KMER MOTIF ENRICHMENT
  # -------------------------------------------------------
  method_sets <- kmers_by_method[[deg_name]]
  method_sets <- method_sets[sapply(method_sets, length) > 0]
  
  if (length(method_sets) >= 2) {
    
    kmers_all_methods <- Reduce(intersect, method_sets)
    kmers_any_method  <- Reduce(union,     method_sets)
    kmers_discordant  <- setdiff(kmers_any_method, kmers_all_methods)
    
    cat(sprintf("  Method-discordant kmers: %d consensus | %d discordant | %d total\n",
                length(kmers_all_methods),
                length(kmers_discordant),
                length(kmers_any_method)))
    
    if (length(kmers_discordant) >= 10 && length(kmers_all_methods) >= 10) {
      
      discordant_dt <- kmer_dt[kmer %in% kmers_discordant]
      consensus_dt  <- kmer_dt[kmer %in% kmers_all_methods]
      
      motif_discordant <- motif_enrichment(discordant_dt, consensus_dt)
      motif_discordant[, context := "method_discordant_vs_consensus"]
      motif_discordant[, deg_set := deg_name]
      
      fwrite(motif_discordant,
             file.path(deg_out,
                       paste0("motif_method_discordant_vs_consensus_", today, ".csv")))
      
      with_progress({
        p <- progressor(steps = length(subkmer_sizes) * 2)
        
        disc_subkmer_list <- future_lapply(subkmer_sizes, function(k_size) {
          p(message = sprintf("Discordant sub-kmer k=%d", k_size))
          sub_disc <- subkmer_enrichment(discordant_dt, consensus_dt, k=k_size)
          sub_cons <- subkmer_enrichment(consensus_dt,  discordant_dt, k=k_size)
          if (nrow(sub_disc) > 0) sub_disc[, direction := "discordant"]
          if (nrow(sub_cons) > 0) sub_cons[, direction := "consensus_enriched"]
          p(message = sprintf("Consensus sub-kmer k=%d", k_size))
          rbindlist(list(sub_disc, sub_cons), fill=TRUE)
        }, future.seed=TRUE)
      })
      
      names(disc_subkmer_list) <- as.character(subkmer_sizes)
      disc_subkmer_results <- rbindlist(disc_subkmer_list, fill=TRUE)
      disc_subkmer_results[, deg_set := deg_name]
      
      fwrite(disc_subkmer_results,
             file.path(deg_out,
                       paste0("subkmer_method_discordant_vs_consensus_", today, ".csv")))
      
      top_disc_motifs <- disc_subkmer_results[padj < 0.05][
        order(pval), .SD[1:min(.N, 20)], by=.(k, direction)]
      for (k_size in subkmer_sizes) {
        if (nrow(top_disc_motifs[k == k_size]) > 0) {
          plot_seqlogo_robust(top_disc_motifs, k_size, deg_out, today,
                              paste0(deg_name, "_discordant"))
        }
      }
      
      cat(sprintf("  Discordant motif enrichment: %d significant kmers\n",
                  sum(motif_discordant$padj < 0.05, na.rm=TRUE)))
      
    } else {
      cat(sprintf("  [SKIP] Too few kmers for discordant enrichment (discordant=%d, consensus=%d)\n",
                  length(kmers_discordant), length(kmers_all_methods)))
    }
    
    cat("  Per-method specific motif enrichment...\n")
    
    method_specific_motifs <- rbindlist(future_lapply(
      names(method_sets), function(meth) {
        
        other_kmers  <- Reduce(union, method_sets[names(method_sets) != meth])
        unique_kmers <- setdiff(method_sets[[meth]], other_kmers)
        
        if (length(unique_kmers) < 10 || length(other_kmers) < 10) return(NULL)
        
        fg <- kmer_dt[kmer %in% unique_kmers]
        bg <- kmer_dt[kmer %in% other_kmers]
        
        res <- motif_enrichment(fg, bg)
        if (nrow(res) == 0) return(NULL)
        
        res[, method  := meth]
        res[, deg_set := deg_name]
        res[, n_unique_kmers := length(unique_kmers)]
        res
      }, future.seed=TRUE), fill=TRUE)
    
    if (!is.null(method_specific_motifs) && nrow(method_specific_motifs) > 0) {
      fwrite(method_specific_motifs,
             file.path(deg_out, paste0("motif_method_specific_", today, ".csv")))
      
      method_sig_summary <- method_specific_motifs[padj < 0.05, .(
        n_sig_motifs     = .N,
        top_motif        = kmer[which.min(pval)],
        mean_enrichment  = mean(fold_enrichment, na.rm=TRUE),
        n_unique_kmers   = first(n_unique_kmers)
      ), by=method]
      setorder(method_sig_summary, -n_sig_motifs)
      
      fwrite(method_sig_summary,
             file.path(deg_out, paste0("motif_method_specific_summary_", today, ".csv")))
      
      cat("  Per-method specific motif summary:\n")
      print(method_sig_summary)
      
    } else {
      cat("  [SKIP] No significant method-specific motifs found\n")
    }
    
    for (dim_col in c("aligner", "fastq_trimmer", "counter")) {
      if (!dim_col %in% names(kmer_dt)) next
      dim_vals <- unique(kmer_dt[[dim_col]])
      dim_vals <- dim_vals[!is.na(dim_vals)]
      if (length(dim_vals) < 2) next
      
      dim_sets <- lapply(dim_vals, function(v) {
        unique(kmer_dt[get(dim_col) == v & padj < 0.05, kmer])
      })
      names(dim_sets) <- dim_vals
      dim_sets <- dim_sets[sapply(dim_sets, length) > 0]
      if (length(dim_sets) < 2) next
      
      dim_consensus   <- Reduce(intersect, dim_sets)
      dim_discordant  <- setdiff(Reduce(union, dim_sets), dim_consensus)
      
      if (length(dim_discordant) < 10 || length(dim_consensus) < 10) next
      
      dim_fg <- kmer_dt[kmer %in% dim_discordant]
      dim_bg <- kmer_dt[kmer %in% dim_consensus]
      
      dim_motifs <- motif_enrichment(dim_fg, dim_bg)
      if (nrow(dim_motifs) == 0) next
      
      dim_motifs[, dimension := dim_col]
      dim_motifs[, deg_set   := deg_name]
      
      fwrite(dim_motifs,
             file.path(deg_out,
                       paste0("motif_discordant_by_", dim_col, "_", today, ".csv")))
      
      cat(sprintf("  %s-discordant motifs: %d sig | %d discordant kmers\n",
                  dim_col,
                  sum(dim_motifs$padj < 0.05, na.rm=TRUE),
                  length(dim_discordant)))
    }
    
  } else {
    cat("  [SKIP] Fewer than 2 non-empty method sets — skipping discordant enrichment\n")
  }
  
  # -------------------------------------------------------
  # 6d. WITHIN-DEG-SET REPRODUCIBILITY
  # -------------------------------------------------------
  by_cols <- intersect(c("ensembl_id", "gene_symbol"), names(kmers_dt))
  reproducibility <- kmers_dt[, .(N = .N), by=by_cols]
  reproducibility[, reproducibility_score := N / length(unique(kmers_dt$method))]
  fwrite(reproducibility,
         file.path(deg_out, paste0("cross_deg_reproducibility_", today, ".csv")))
  
  # -------------------------------------------------------
  # 6e. GENE-STRATIFIED KMER NETWORK + PERTURBATION ROBUSTNESS
  # -------------------------------------------------------
  net_edges <- kmers_dt[!is.na(ensembl_id), .(kmer, ensembl_id)]
  
  if (nrow(net_edges) > 0) {
    g <- graph_from_data_frame(net_edges, directed=FALSE)
    saveRDS(g, file.path(deg_out, paste0("gene_kmer_network_", today, ".rds")))
    
    network_robustness <- replicate(100, {
      remove <- sample(V(g)$name, floor(vcount(g) * 0.1))
      g2     <- delete_vertices(g, remove)
      length(components(g2)$no)
    })
    fwrite(data.table(n_components=network_robustness),
           file.path(deg_out, paste0("network_robustness_", today, ".csv")))
  } else {
    cat("  [SKIP] Network: no annotated kmer-gene edges\n")
  }
  
  # -------------------------------------------------------
  # 6f. METHOD INTERSECTION (UpSet)
  # -------------------------------------------------------
  upset_input <- kmers_by_method[[deg_name]]
  upset_input <- upset_input[sapply(upset_input, length) > 0]
  
  if (length(upset_input) >= 2) {
    tryCatch({
      pdf(file.path(deg_out, paste0("upset_methods_", today, ".pdf")),
          width=10, height=6)
      print(upset(fromList(upset_input),
                  nsets      = min(10, length(upset_input)),
                  order.by   = "freq",
                  decreasing = TRUE))
      dev.off()
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat(sprintf("  [WARNING] UpSet plot failed: %s\n", e$message))
    })
  } else {
    cat("  [SKIP] UpSet: fewer than 2 non-empty method sets\n")
  }
  
  # -------------------------------------------------------
  # 6g. GC-CONTENT DISTRIBUTION
  # -------------------------------------------------------
  if ("gc_content" %in% names(kmers_dt) &&
      sum(!is.na(kmers_dt$gc_content)) > 0) {
    p_gc <- ggplot(kmers_dt, aes(gc_content, fill=fixed_value)) +
      geom_density(alpha=0.4) +
      facet_wrap(~ fixed_value) +
      theme_classic() +
      labs(title=paste0("GC-content by FASTQ trimmer\n", deg_name))
    ggsave(file.path(deg_out, paste0("gc_density_", today, ".pdf")), p_gc)
  }
  
  # -------------------------------------------------------
  # 6h. SUMMARY TABLE
  # -------------------------------------------------------
  summary_dt <- kmers_dt[, .(
    n_sig         = .N,
    unique_genes  = uniqueN(ensembl_id),
    mean_gc       = mean(gc_content,  na.rm=TRUE),
    mean_baseMean = mean(baseMean,    na.rm=TRUE)
  ), by=.(method, grouping_type, fixed_value, counter)]
  
  fwrite(summary_dt,
         file.path(deg_out, paste0("method_summary_", today, ".csv")))
  all_summaries[[deg_name]] <- summary_dt
  
  # -------------------------------------------------------
  # 6i. PER-DEG-SET GLMM ANALYSIS
  # FIX: prepare a deduplicated input of unique kmer x method
  # rows before passing to GLMM. The full `kmers_dt` at this
  # point is the DEG-filtered merged table — still large due
  # to the DEG-gene x aligner cartesian. Collapsing to unique
  # kmer x pipeline dimension combinations gives the model the
  # correct unit of analysis (one binary detection observation
  # per kmer per method) at a tractable size.
  # -------------------------------------------------------
  cat("  Running per-DEG-set GLMM analysis...\n")
  
  # Build GLMM input: one row per unique kmer x pipeline combination
  # Include all pipeline dimensions + sequence features
  glmm_cols <- intersect(
    c("kmer", "ensembl_id", "DEG_binary",
      "aligner", "fastq_trimmer", "counter", "deg_list_trimmer",
      "kmer_entropy", "gc_content"),
    names(kmers_dt)
  )
  glmm_input <- unique(kmers_dt[, glmm_cols, with=FALSE])
  
  cat(sprintf("  GLMM input: %d unique kmer x method rows (from %d merged rows)\n",
              nrow(glmm_input), nrow(kmers_dt)))
  
  glmm_results_per_deg <- run_glmm_analysis(
    dt            = glmm_input,
    grouping_type = meta$grouping_type,
    fixed_value   = meta$fixed_value %||% "unknown",
    deg_name      = deg_name,
    out_dir       = deg_out,
    today         = today
  )
  
  if (!is.null(glmm_results_per_deg)) {
    glmm_results_all[[deg_name]] <- glmm_results_per_deg
    cat(sprintf("  GLMM: %d models fitted | %d significant terms (p<0.05)\n",
                uniqueN(glmm_results_per_deg$model),
                sum(glmm_results_per_deg$p.value < 0.05, na.rm=TRUE)))
  }
  
  cat("  [DONE]", deg_name, "\n")
}

# ============================================================
# POST-LOOP: COMBINE ALL RESULTS
# ============================================================

cat("\nCombining all DEG sets...\n")

reconstruct_fixed_value <- function(dt) {
  result <- rep("unknown", nrow(dt))
  if ("deg_list_trimmer" %in% names(dt)) {
    mask <- !is.na(dt$deg_list_trimmer) & dt$deg_list_trimmer != ""
    result[mask] <- dt$deg_list_trimmer[mask]
  }
  if ("deg_list_aligner" %in% names(dt)) {
    mask <- result == "unknown" & !is.na(dt$deg_list_aligner) & dt$deg_list_aligner != ""
    result[mask] <- dt$deg_list_aligner[mask]
  }
  if ("deg_list_counter" %in% names(dt)) {
    mask <- result == "unknown" & !is.na(dt$deg_list_counter) & dt$deg_list_counter != ""
    result[mask] <- dt$deg_list_counter[mask]
  }
  return(result)
}

analysis_dt <- rbindlist(full_dataset, fill=TRUE)

if (!"fixed_value" %in% names(analysis_dt)) {
  cat("[WARNING] fixed_value missing from analysis_dt — reconstructing\n")
  analysis_dt[, fixed_value := reconstruct_fixed_value(analysis_dt)]
} else if (any(is.na(analysis_dt$fixed_value))) {
  cat("[WARNING] fixed_value has NAs — filling\n")
  na_rows <- is.na(analysis_dt$fixed_value)
  analysis_dt[na_rows, fixed_value := reconstruct_fixed_value(analysis_dt[na_rows])]
}

cat("fixed_value distribution:\n")
print(analysis_dt[, .N, by=fixed_value])

# ============================================================
# CROSS-DEG ANALYSES
# ============================================================

all_deg_dt_combined <- rbindlist(all_deg_dt, fill=TRUE)

if (!"fixed_value" %in% names(all_deg_dt_combined)) {
  cat("[WARNING] fixed_value missing from all_deg_dt_combined — reconstructing\n")
  all_deg_dt_combined[, fixed_value := reconstruct_fixed_value(all_deg_dt_combined)]
} else if (any(is.na(all_deg_dt_combined$fixed_value))) {
  na_rows <- is.na(all_deg_dt_combined$fixed_value)
  all_deg_dt_combined[na_rows,
                      fixed_value := reconstruct_fixed_value(all_deg_dt_combined[na_rows])]
}

fwrite(all_deg_dt_combined,
       file.path(output_dir, "tables",
                 paste0("all_deg_tables_combined_", today, ".csv")))

gene_cross_deg <- all_deg_dt_combined[, .(
  n_deg_sets    = uniqueN(deg_name),
  n_aligners    = uniqueN(aligner),
  deg_sets_list = paste(unique(deg_name), collapse="|"),
  aligners_list = paste(unique(aligner),  collapse="|")
), by=ensembl_id]
gene_cross_deg[, cross_deg_fraction := n_deg_sets / length(deg_files)]

fwrite(gene_cross_deg,
       file.path(output_dir, "tables",
                 paste0("gene_cross_deg_reproducibility_", today, ".csv")))

gene_cross_deg_aligner <- all_deg_dt_combined[, .(
  n_deg_sets         = uniqueN(deg_name),
  cross_deg_fraction = uniqueN(deg_name) / length(deg_files)
), by=.(ensembl_id, aligner)]

fwrite(gene_cross_deg_aligner,
       file.path(output_dir, "tables",
                 paste0("gene_cross_deg_per_aligner_", today, ".csv")))

gene_cross_deg_counter <- all_deg_dt_combined[, .(
  n_deg_sets         = uniqueN(deg_name),
  cross_deg_fraction = uniqueN(deg_name) / length(deg_files)
), by=.(ensembl_id, counter)]

fwrite(gene_cross_deg_counter,
       file.path(output_dir, "tables",
                 paste0("gene_cross_deg_per_counter_", today, ".csv")))

gene_cross_deg_fixed <- all_deg_dt_combined[, .(
  n_deg_sets         = uniqueN(deg_name),
  cross_deg_fraction = uniqueN(deg_name) / length(deg_files)
), by=.(ensembl_id, fixed_value)]

fwrite(gene_cross_deg_fixed,
       file.path(output_dir, "tables",
                 paste0("gene_cross_deg_per_fixed_value_", today, ".csv")))

p_cross_deg <- ggplot(gene_cross_deg, aes(x=cross_deg_fraction)) +
  geom_histogram(bins=50, fill="steelblue", colour="white") +
  theme_classic() +
  labs(title="Gene cross-DEG reproducibility",
       x="Fraction of DEG sets gene appeared in",
       y="Number of genes")
ggsave(file.path(output_dir, "figures",
                 paste0("cross_deg_reproducibility_distribution_", today, ".pdf")),
       p_cross_deg, width=8, height=6)

p_cross_deg_aligner <- ggplot(gene_cross_deg_aligner,
                              aes(x=cross_deg_fraction, fill=aligner)) +
  geom_density(alpha=0.4) +
  theme_classic() +
  labs(title="Cross-DEG gene reproducibility by aligner",
       x="Fraction of DEG sets gene appeared in", y="Density")
ggsave(file.path(output_dir, "figures",
                 paste0("cross_deg_reproducibility_by_aligner_", today, ".pdf")),
       p_cross_deg_aligner, width=10, height=6)

p_cross_deg_fixed <- ggplot(gene_cross_deg_fixed,
                            aes(x=cross_deg_fraction, fill=fixed_value)) +
  geom_density(alpha=0.4) +
  theme_classic() +
  labs(title="Cross-DEG gene reproducibility by fixed dimension",
       x="Fraction of DEG sets gene appeared in", y="Density")
ggsave(file.path(output_dir, "figures",
                 paste0("cross_deg_reproducibility_by_fixed_value_", today, ".pdf")),
       p_cross_deg_fixed, width=10, height=6)

cat("\nCross-DEG summary:\n")
cat("  Total unique genes across all DEG sets:",
    uniqueN(all_deg_dt_combined$ensembl_id), "\n")
cat("  Genes in all DEG sets (fraction=1)    :",
    sum(gene_cross_deg$cross_deg_fraction == 1), "\n")
cat("  Genes in >50% of DEG sets             :",
    sum(gene_cross_deg$cross_deg_fraction > 0.5), "\n")

# ============================================================
# COMBINED GLMM ANALYSIS ACROSS ALL DEG SETS
# ============================================================

cat("\nRunning combined GLMM analysis across all DEG sets...\n")

if (nrow(analysis_dt) > 100 && requireNamespace("lme4", quietly=TRUE) &&
    requireNamespace("broom.mixed", quietly=TRUE)) {
  
  suppressPackageStartupMessages({
    library(lme4)
    library(broom.mixed)
  })
  
  glmm_combined_results <- list()
  
  glmm_dt <- analysis_dt[
    !is.na(DEG_binary) & !is.na(ensembl_id) & !is.na(deg_name),
    .(
      DEG_binary       = DEG_binary,
      ensembl_id       = ensembl_id,
      deg_name         = deg_name,
      grouping_type    = grouping_type,
      fixed_value      = fixed_value,
      aligner          = aligner,
      fastq_trimmer    = fastq_trimmer,
      deg_list_trimmer = if ("deg_list_trimmer" %in% names(analysis_dt))
        deg_list_trimmer else NA_character_,
      counter          = counter,
      kmer_entropy     = if ("kmer_entropy" %in% names(analysis_dt))
        kmer_entropy else NA_real_,
      gc_content       = if ("gc_content" %in% names(analysis_dt))
        gc_content else NA_real_
    )
  ]
  
  glmm_dt[, gene_f         := as.factor(ensembl_id)]
  glmm_dt[, deg_set_f      := as.factor(deg_name)]
  glmm_dt[, grouping_f     := as.factor(grouping_type)]
  glmm_dt[, fixed_value_f  := as.factor(fixed_value)]
  glmm_dt[, aligner_f      := as.factor(aligner)]
  glmm_dt[, fastq_trim_f   := as.factor(fastq_trimmer)]
  glmm_dt[, counter_f      := as.factor(counter)]
  
  has_aligner       <- uniqueN(glmm_dt$aligner_f,    na.rm=TRUE) > 1
  has_fastq_trimmer <- uniqueN(glmm_dt$fastq_trim_f, na.rm=TRUE) > 1
  has_counter       <- uniqueN(glmm_dt$counter_f,    na.rm=TRUE) > 1
  has_entropy       <- sum(!is.na(glmm_dt$kmer_entropy)) > 100
  has_gc            <- sum(!is.na(glmm_dt$gc_content))   > 100
  
  cat(sprintf("  Combined GLMM data: %d rows | %d genes | %d DEG sets\n",
              nrow(glmm_dt),
              uniqueN(glmm_dt$gene_f),
              uniqueN(glmm_dt$deg_set_f)))
  cat(sprintf("  Dimensions: aligner=%s | fastq_trimmer=%s | counter=%s\n",
              has_aligner, has_fastq_trimmer, has_counter))
  
  if (has_aligner) {
    tryCatch({
      fit_c1 <- glmer(
        DEG_binary ~ aligner_f + (1 | gene_f) + (1 | deg_set_f),
        data    = glmm_dt[!is.na(aligner_f)],
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      )
      coefs <- as.data.table(tidy(fit_c1, effects="fixed"))
      coefs[, model := "aligner_crossed_random"]
      glmm_combined_results[["c1_aligner"]] <- coefs
      cat(sprintf("  Combined GLMM aligner: AIC=%.1f\n", AIC(fit_c1)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] Combined GLMM aligner failed: %s\n", e$message))
    })
  }
  
  if (has_fastq_trimmer) {
    tryCatch({
      fit_c2 <- glmer(
        DEG_binary ~ fastq_trim_f + (1 | gene_f) + (1 | deg_set_f),
        data    = glmm_dt[!is.na(fastq_trim_f)],
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      )
      coefs <- as.data.table(tidy(fit_c2, effects="fixed"))
      coefs[, model := "fastq_trimmer_crossed_random"]
      glmm_combined_results[["c2_fastq_trimmer"]] <- coefs
      cat(sprintf("  Combined GLMM fastq_trimmer: AIC=%.1f\n", AIC(fit_c2)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] Combined GLMM fastq_trimmer failed: %s\n", e$message))
    })
  }
  
  if (has_counter) {
    tryCatch({
      fit_c3 <- glmer(
        DEG_binary ~ counter_f + (1 | gene_f) + (1 | deg_set_f),
        data    = glmm_dt[!is.na(counter_f)],
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      )
      coefs <- as.data.table(tidy(fit_c3, effects="fixed"))
      coefs[, model := "counter_crossed_random"]
      glmm_combined_results[["c3_counter"]] <- coefs
      cat(sprintf("  Combined GLMM counter: AIC=%.1f\n", AIC(fit_c3)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] Combined GLMM counter failed: %s\n", e$message))
    })
  }
  
  if (has_aligner && has_fastq_trimmer) {
    tryCatch({
      fit_c4 <- glmer(
        DEG_binary ~ aligner_f + fastq_trim_f + (1 | gene_f) + (1 | deg_set_f),
        data    = glmm_dt[!is.na(aligner_f) & !is.na(fastq_trim_f)],
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      )
      coefs <- as.data.table(tidy(fit_c4, effects="fixed"))
      coefs[, model := "aligner_plus_fastq_trimmer"]
      glmm_combined_results[["c4_aligner_plus_trimmer"]] <- coefs
      cat(sprintf("  Combined GLMM aligner+fastq_trimmer: AIC=%.1f\n", AIC(fit_c4)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] Combined GLMM aligner+fastq_trimmer failed: %s\n",
                  e$message))
    })
  }
  
  if (has_aligner && has_fastq_trimmer) {
    tryCatch({
      fit_c5 <- glmer(
        DEG_binary ~ aligner_f * fastq_trim_f + (1 | gene_f) + (1 | deg_set_f),
        data    = glmm_dt[!is.na(aligner_f) & !is.na(fastq_trim_f)],
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=3e5))
      )
      coefs <- as.data.table(tidy(fit_c5, effects="fixed"))
      coefs[, model := "aligner_x_fastq_trimmer"]
      glmm_combined_results[["c5_aligner_x_trimmer"]] <- coefs
      cat(sprintf("  Combined GLMM aligner x fastq_trimmer: AIC=%.1f\n", AIC(fit_c5)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] Combined GLMM aligner x fastq_trimmer failed: %s\n",
                  e$message))
    })
  }
  
  if (uniqueN(glmm_dt$grouping_f) > 1) {
    tryCatch({
      fit_c6 <- glmer(
        DEG_binary ~ grouping_f + (1 | gene_f) + (1 | deg_set_f),
        data    = glmm_dt,
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      )
      coefs <- as.data.table(tidy(fit_c6, effects="fixed"))
      coefs[, model := "grouping_type_effect"]
      glmm_combined_results[["c6_grouping_type"]] <- coefs
      cat(sprintf("  Combined GLMM grouping_type: AIC=%.1f\n", AIC(fit_c6)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] Combined GLMM grouping_type failed: %s\n", e$message))
    })
  }
  
  if (has_entropy && has_gc) {
    tryCatch({
      dt_seq <- glmm_dt[!is.na(kmer_entropy) & !is.na(gc_content)]
      fit_c7 <- glmer(
        DEG_binary ~ scale(kmer_entropy) + scale(gc_content) +
          (1 | gene_f) + (1 | deg_set_f),
        data    = dt_seq,
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
      )
      coefs <- as.data.table(tidy(fit_c7, effects="fixed"))
      coefs[, model := "entropy_gc_crossed_random"]
      glmm_combined_results[["c7_seq_features"]] <- coefs
      cat(sprintf("  Combined GLMM entropy+gc: AIC=%.1f\n", AIC(fit_c7)))
    }, error = function(e) {
      cat(sprintf("  [WARNING] Combined GLMM entropy+gc failed: %s\n", e$message))
    })
  }
  
  # Model C8: Full combined model
  # FIX: replace ..req_cols with with=FALSE for in-function data.table scoping
  tryCatch({
    fixed_terms <- character(0)
    if (has_aligner)       fixed_terms <- c(fixed_terms, "aligner_f")
    if (has_fastq_trimmer) fixed_terms <- c(fixed_terms, "fastq_trim_f")
    if (has_counter)       fixed_terms <- c(fixed_terms, "counter_f")
    if (has_entropy)       fixed_terms <- c(fixed_terms, "scale(kmer_entropy)")
    if (has_gc)            fixed_terms <- c(fixed_terms, "scale(gc_content)")
    if (has_aligner && has_entropy)
      fixed_terms <- c(fixed_terms, "aligner_f:scale(kmer_entropy)")
    if (has_fastq_trimmer && has_entropy)
      fixed_terms <- c(fixed_terms, "fastq_trim_f:scale(kmer_entropy)")
    
    if (length(fixed_terms) >= 2) {
      formula_str <- paste0(
        "DEG_binary ~ ",
        paste(fixed_terms, collapse=" + "),
        " + (1 | gene_f) + (1 | deg_set_f)"
      )
      cat(sprintf("  Combined GLMM full formula: %s\n", formula_str))
      
      req_cols <- c("DEG_binary", "gene_f", "deg_set_f",
                    if (has_aligner)       "aligner_f",
                    if (has_fastq_trimmer) "fastq_trim_f",
                    if (has_counter)       "counter_f",
                    if (has_entropy)       "kmer_entropy",
                    if (has_gc)            "gc_content")
      # FIX: with=FALSE instead of ..req_cols
      dt_full <- glmm_dt[rowSums(is.na(glmm_dt[, req_cols, with=FALSE])) == 0]
      
      fit_c8 <- glmer(
        as.formula(formula_str),
        data    = dt_full,
        family  = binomial(),
        control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=3e5))
      )
      coefs <- as.data.table(tidy(fit_c8, effects="fixed"))
      coefs[, model := "full_combined_model"]
      glmm_combined_results[["c8_full"]] <- coefs
      cat(sprintf("  Combined GLMM full: AIC=%.1f | %d terms\n",
                  AIC(fit_c8), nrow(coefs)))
      
      saveRDS(fit_c8,
              file.path(output_dir, "models",
                        paste0("glmm_full_combined_", today, ".rds")))
    }
  }, error = function(e) {
    cat(sprintf("  [WARNING] Combined GLMM full model failed: %s\n", e$message))
  })
  
  if (length(glmm_combined_results) > 0) {
    
    glmm_combined_dt <- rbindlist(glmm_combined_results, fill=TRUE)
    glmm_combined_dt[!is.na(p.value) & term != "(Intercept)",
                     padj_within_model := p.adjust(p.value, method="BH"),
                     by=model]
    glmm_combined_dt[!is.na(p.value) & term != "(Intercept)",
                     padj_global := p.adjust(p.value, method="BH")]
    glmm_combined_dt[!is.na(p.value) & term != "(Intercept)",
                     padj_bonferroni := p.adjust(p.value, method="bonferroni")]
    fwrite(glmm_combined_dt,
           file.path(output_dir, "tables",
                     paste0("glmm_combined_across_degsets_", today, ".csv")))
    
    glmm_combined_summary <- glmm_combined_dt[
      !is.na(p.value) & term != "(Intercept)", .(
        n_models_tested  = uniqueN(model),
        n_significant    = sum(p.value < 0.05),
        frac_significant = mean(p.value < 0.05),
        mean_estimate    = mean(estimate,  na.rm=TRUE),
        sd_estimate      = sd(estimate,    na.rm=TRUE),
        min_p            = min(p.value,    na.rm=TRUE)
      ), by=.(model, term)]
    setorder(glmm_combined_summary, model, min_p)
    
    fwrite(glmm_combined_summary,
           file.path(output_dir, "tables",
                     paste0("glmm_combined_summary_", today, ".csv")))
    
    cat("\nCombined GLMM summary — significant terms per model:\n")
    print(glmm_combined_summary[frac_significant > 0])
    
    sig_terms <- glmm_combined_dt[
      !is.na(padj_within_model) &
        term != "(Intercept)" &
        padj_within_model < 0.05]
    
    if (nrow(sig_terms) > 0) {
      sig_terms[, neg_log10p  := -log10(p.value + 1e-300)]
      sig_terms[, significant := p.value < 0.05]
      
      p_combined_forest <- ggplot(
        sig_terms,
        aes(x      = estimate,
            xmin   = estimate - 1.96 * std.error,
            xmax   = estimate + 1.96 * std.error,
            y      = term,
            colour = significant)) +
        geom_pointrange(size=0.4) +
        geom_vline(xintercept=0, linetype="dashed") +
        facet_wrap(~ model, scales="free_x") +
        scale_colour_manual(
          values = c("FALSE"="grey60", "TRUE"="steelblue"),
          labels = c("n.s.", "p<0.05")) +
        theme_classic() +
        theme(axis.text.y  = element_text(size=7),
              strip.text   = element_text(size=7)) +
        labs(title  = "Combined GLMM — pipeline dimension effects across all DEG sets",
             x      = "Effect size (log-odds) ± 95% CI",
             y      = "Term",
             colour = "Significant")
      ggsave(
        file.path(output_dir, "figures",
                  paste0("glmm_combined_forest_", today, ".pdf")),
        p_combined_forest,
        width  = 16,
        height = max(8, uniqueN(sig_terms$term) * 0.35)
      )
      
      dim_summary <- sig_terms[, .(
        mean_abs_effect = mean(abs(estimate), na.rm=TRUE),
        n_terms         = .N,
        n_sig           = sum(p.value < 0.05)
      ), by=model]
      
      p_variance <- ggplot(dim_summary,
                           aes(x=model, y=mean_abs_effect, fill=model)) +
        geom_col() +
        theme_classic() +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              legend.position="none") +
        labs(title="Mean absolute effect size by pipeline dimension model",
             x="Model", y="Mean |effect size| (log-odds)")
      ggsave(
        file.path(output_dir, "figures",
                  paste0("glmm_combined_variance_partition_", today, ".pdf")),
        p_variance, width=10, height=6
      )
    }
    
    # Combine per-DEG GLMM results for consistency comparison
    glmm_per_deg_dt <- if (length(glmm_results_all) > 0) {
      rbindlist(glmm_results_all, fill=TRUE)
    } else {
      NULL
    }
    
    if (!is.null(glmm_per_deg_dt) && nrow(glmm_per_deg_dt) > 0) {
      fwrite(glmm_per_deg_dt,
             file.path(output_dir, "tables",
                       paste0("glmm_per_degset_all_", today, ".csv")))
      
      glmm_consistency <- merge(
        glmm_per_deg_dt[term != "(Intercept)", .(
          per_deg_mean_estimate = mean(estimate,       na.rm=TRUE),
          per_deg_sd_estimate   = sd(estimate,         na.rm=TRUE),
          per_deg_frac_sig      = mean(p.value < 0.05, na.rm=TRUE),
          per_deg_frac_sig_adj  = mean(!is.na(padj_within_model) &
                                         padj_within_model < 0.05, na.rm=TRUE),
          n_deg_sets            = uniqueN(deg_set)
        ), by=.(model, term)],
        glmm_combined_dt[term != "(Intercept)", .(
          combined_estimate    = mean(estimate,         na.rm=TRUE),
          combined_p           = min(p.value,           na.rm=TRUE),
          combined_padj        = min(padj_within_model, na.rm=TRUE),
          combined_padj_global = min(padj_global,       na.rm=TRUE),
          combined_bonferroni  = min(padj_bonferroni,   na.rm=TRUE)
        ), by=.(model, term)],
        by  = c("model", "term"),
        all = TRUE
      )
      
      glmm_consistency[, consistent_sig := per_deg_frac_sig > 0.5 &
                         combined_padj < 0.05]
      glmm_consistency[, consistent_sig_strict := per_deg_frac_sig_adj > 0.5 &
                         combined_padj_global < 0.05]
      glmm_consistency[, direction_consistent := sign(per_deg_mean_estimate) ==
                         sign(combined_estimate)]
      glmm_consistency[, consistency_score :=
                         per_deg_frac_sig *
                         as.integer(direction_consistent) *
                         (-log10(combined_padj + 1e-300) /
                            max(-log10(combined_padj + 1e-300), na.rm=TRUE))]
      
      setorder(glmm_consistency, -consistency_score)
      
      fwrite(glmm_consistency,
             file.path(output_dir, "tables",
                       paste0("glmm_consistency_perDEG_vs_combined_", today, ".csv")))
      
      cat("\nGLMM consistency — top terms (consistent_sig == TRUE):\n")
      print(glmm_consistency[consistent_sig == TRUE,
                             .(model, term, per_deg_frac_sig, combined_padj,
                               combined_padj_global, direction_consistent,
                               consistency_score)])
      
      cat(sprintf("\n  Consistent terms (lenient) : %d\n",
                  sum(glmm_consistency$consistent_sig,        na.rm=TRUE)))
      cat(sprintf("  Consistent terms (strict)  : %d\n",
                  sum(glmm_consistency$consistent_sig_strict, na.rm=TRUE)))
    } else {
      cat("[INFO] No per-DEG GLMM results to compare against combined model\n")
    }
    
  } else {
    cat("[INFO] No combined GLMM models converged\n")
  }
  
} else {
  cat("[INFO] Skipping combined GLMM — insufficient data or missing packages\n")
}

# ============================================================
# COMBINE MOTIF RESULTS
# ============================================================

if (length(motif_results_all) > 0) {
  motif_results_dt <- rbindlist(motif_results_all, fill=TRUE)
  fwrite(motif_results_dt,
         file.path(output_dir, "tables",
                   paste0("motif_enrichment_all_degsets_", today, ".csv")))
}

if (length(subkmer_results_all) > 0) {
  subkmer_results_dt <- rbindlist(subkmer_results_all, fill=TRUE)
  fwrite(subkmer_results_dt,
         file.path(output_dir, "tables",
                   paste0("subkmer_enrichment_all_degsets_", today, ".csv")))
  
  top_global <- subkmer_results_dt[padj < 0.05][
    order(pval), .SD[1:min(.N, 50)], by=.(k, direction)]
  fwrite(top_global,
         file.path(output_dir, "tables",
                   paste0("top_subkmer_motifs_global_", today, ".csv")))
}

disc_motif_files <- list.files(output_dir, pattern="motif_method_discordant_vs_consensus_",
                               recursive=TRUE, full.names=TRUE)
if (length(disc_motif_files) > 0) {
  disc_motifs_all <- rbindlist(lapply(disc_motif_files, fread), fill=TRUE)
  fwrite(disc_motifs_all,
         file.path(output_dir, "tables",
                   paste0("motif_method_discordant_all_degsets_", today, ".csv")))
  
  disc_motifs_global <- disc_motifs_all[padj < 0.05, .(
    n_deg_sets       = uniqueN(deg_set),
    mean_enrichment  = mean(fold_enrichment, na.rm=TRUE),
    min_padj         = min(padj, na.rm=TRUE),
    contexts         = paste(unique(context), collapse="|")
  ), by=kmer]
  setorder(disc_motifs_global, -n_deg_sets, min_padj)
  fwrite(disc_motifs_global,
         file.path(output_dir, "tables",
                   paste0("motif_discordant_global_summary_", today, ".csv")))
  cat("\nTop globally discordant motifs:\n")
  print(head(disc_motifs_global, 10))
}

for (dim_col in c("aligner", "fastq_trimmer", "counter")) {
  dim_files <- list.files(output_dir,
                          pattern=paste0("motif_discordant_by_", dim_col, "_"),
                          recursive=TRUE, full.names=TRUE)
  if (length(dim_files) == 0) next
  
  dim_motifs_all <- rbindlist(lapply(dim_files, fread), fill=TRUE)
  fwrite(dim_motifs_all,
         file.path(output_dir, "tables",
                   paste0("motif_discordant_by_", dim_col, "_all_degsets_", today, ".csv")))
  
  cat(sprintf("\n%s-discordant motifs across all DEG sets: %d significant\n",
              dim_col, sum(dim_motifs_all$padj < 0.05, na.rm=TRUE)))
}

# ============================================================
# 8. GLOBAL FASTQ-TRIMMER SENSITIVITY SUMMARY
# ============================================================

fastq_sensitivity_dt <- rbindlist(fastq_sensitivity_all, fill=TRUE, idcol="deg_set")

global_fastq_sens <- fastq_sensitivity_dt[, .(
  mean_genes_norm = mean(genes_norm, na.rm=TRUE),
  mean_kmers_norm = mean(kmers_norm, na.rm=TRUE),
  sd_genes_norm   = sd(genes_norm,   na.rm=TRUE),
  n_deg_sets      = uniqueN(deg_set)
), by=.(grouping_type, fixed_value, varying_value)]

fwrite(global_fastq_sens,
       file.path(output_dir, "tables",
                 paste0("global_fastq_sensitivity_summary_", today, ".csv")))

p_global_sens <- ggplot(global_fastq_sens,
                        aes(x=varying_value, y=mean_genes_norm,
                            ymin=mean_genes_norm - sd_genes_norm,
                            ymax=mean_genes_norm + sd_genes_norm,
                            colour=varying_value)) +
  geom_pointrange() +
  facet_grid(grouping_type ~ fixed_value, scales="free_x") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        strip.text.y=element_text(angle=0)) +
  labs(title="Global sensitivity by varying dimension and fixed context",
       y="Mean normalised gene recovery",
       x="Varying dimension value",
       colour="Varying value")
ggsave(file.path(output_dir, "figures",
                 paste0("global_fastq_sensitivity_", today, ".pdf")),
       p_global_sens, width=14, height=8)

# ============================================================
# 9. GENE ROBUSTNESS STABILITY CURVES
# ============================================================

gene_freq <- table(unlist(gene_tracker))
gene_robustness_global <- data.table(
  ensembl_id = names(gene_freq),
  freq        = as.numeric(gene_freq)
)
gene_robustness_global[, robustness := freq / length(gene_tracker)]

fwrite(gene_robustness_global,
       file.path(output_dir, "tables",
                 paste0("gene_robustness_global_", today, ".csv")))

png(file.path(output_dir, "figures",
              paste0("GeneStabilityCurve_", today, ".png")),
    width=2400, height=2000, res=600)
print(
  ggplot(gene_robustness_global, aes(x=robustness)) +
    geom_density(fill="steelblue", alpha=0.5) +
    theme_classic() +
    labs(title="Gene robustness across all DEG sets",
         x="Robustness (fraction of DEG sets gene appeared in)")
)
dev.off()

# ============================================================
# 10. CROSS-DEG CONSENSUS & ROBUSTNESS RANKING
# ============================================================

all_summaries_dt <- rbindlist(all_summaries, idcol="deg_set")

consensus <- all_summaries_dt[, .(
  total_sig   = sum(n_sig),
  total_genes = sum(unique_genes),
  avg_gc      = mean(mean_gc,       na.rm=TRUE),
  avg_expr    = mean(mean_baseMean, na.rm=TRUE)
), by=.(method, counter, fixed_value)]

consensus[, robustness_score :=
            as.numeric(scale(total_genes)) -
            abs(as.numeric(scale(avg_gc))) -
            abs(as.numeric(scale(avg_expr)))]

fwrite(consensus,
       file.path(output_dir, "diagnostics",
                 paste0("ROBUSTNESS_RANKING_", today, ".csv")))

# ============================================================
# 11. BUILD FULL GENE x PIPELINE COMBINATION MATRIX
# ============================================================

required_cols <- c("ensembl_id", "aligner", "counter", "grouping_type",
                   "fixed_value", "DEG_binary")
missing_cols  <- setdiff(required_cols, names(analysis_dt))
if (length(missing_cols) > 0) {
  cat(sprintf("[WARNING] Missing columns in analysis_dt: %s\n",
              paste(missing_cols, collapse=", ")))
}

all_ensembl   <- unique(analysis_dt$ensembl_id[!is.na(analysis_dt$ensembl_id)])
all_aligners  <- unique(analysis_dt$aligner[!is.na(analysis_dt$aligner)])
all_trimmers  <- unique(analysis_dt$fastq_trimmer[!is.na(analysis_dt$fastq_trimmer)])
all_counters  <- unique(analysis_dt$counter[!is.na(analysis_dt$counter)])
all_fixed     <- unique(analysis_dt$fixed_value[!is.na(analysis_dt$fixed_value)])
all_groupings <- unique(analysis_dt$grouping_type[!is.na(analysis_dt$grouping_type)])

cat(sprintf("  Grid dimensions:\n"))
cat(sprintf("    ensembl_id    : %d unique genes\n",    length(all_ensembl)))
cat(sprintf("    aligner       : %d (%s)\n",
            length(all_aligners), paste(all_aligners, collapse=", ")))
cat(sprintf("    fastq_trimmer : %d (%s)\n",
            length(all_trimmers), paste(all_trimmers, collapse=", ")))
cat(sprintf("    counter       : %d (%s)\n",
            length(all_counters), paste(all_counters, collapse=", ")))
cat(sprintf("    fixed_value   : %d (%s)\n",
            length(all_fixed), paste(all_fixed, collapse=", ")))

full_grid <- CJ(
  ensembl_id    = all_ensembl,
  aligner       = all_aligners,
  fastq_trimmer = all_trimmers,
  counter       = all_counters
)

join_cols <- intersect(
  c("ensembl_id", "aligner", "fastq_trimmer", "counter"),
  names(analysis_dt)
)
detected <- unique(analysis_dt[!is.na(ensembl_id), join_cols, with=FALSE])
detected[, DEG_binary := 1L]

full_dt <- merge(full_grid, detected,
                 by = join_cols,
                 all.x = TRUE)
full_dt[is.na(DEG_binary), DEG_binary := 0L]

grouping_map <- analysis_dt[, .(
  grouping_type = grouping_type[1],
  fixed_value   = fixed_value[1]
), by=.(aligner, fastq_trimmer, counter)]
grouping_map <- unique(grouping_map)

full_dt <- merge(full_dt, grouping_map,
                 by  = intersect(c("aligner", "fastq_trimmer", "counter"),
                                 names(full_dt)),
                 all.x = TRUE)

cat(sprintf("  full_dt: %d rows | %d genes\n",
            nrow(full_dt), uniqueN(full_dt$ensembl_id)))
cat(sprintf("  DEG_binary=1: %d (%.1f%%) | DEG_binary=0: %d (%.1f%%)\n",
            sum(full_dt$DEG_binary == 1L),
            100 * mean(full_dt$DEG_binary == 1L),
            sum(full_dt$DEG_binary == 0L),
            100 * mean(full_dt$DEG_binary == 0L)))

n_detected   <- nrow(detected)
n_in_full    <- sum(full_dt$DEG_binary == 1L)
if (n_detected != n_in_full) {
  cat(sprintf("  [WARNING] Detected=%d but DEG_binary=1 rows=%d — check join cols\n",
              n_detected, n_in_full))
} else {
  cat(sprintf("  [OK] All %d detected DEG combinations flagged correctly\n",
              n_detected))
}

gene_feature_cols <- intersect(
  c("ensembl_id", "gc_content", "kmer_entropy", "length",
    "has_exon_overlap", "gene_symbol"),
  names(analysis_dt)
)
if (length(gene_feature_cols) > 1) {
  gene_features <- analysis_dt[, lapply(.SD, function(x) {
    if (is.numeric(x)) mean(x, na.rm=TRUE) else x[!is.na(x)][1]
  }), by=ensembl_id, .SDcols=setdiff(gene_feature_cols, "ensembl_id")]
  
  full_dt <- merge(full_dt, gene_features, by="ensembl_id", all.x=TRUE)
  cat(sprintf("  Merged %d feature columns into full_dt\n",
              length(gene_feature_cols) - 1))
}

cat(sprintf("Full gene x pipeline matrix: %d rows | DEG_binary=1: %d (%.1f%%)\n",
            nrow(full_dt),
            sum(full_dt$DEG_binary),
            100 * mean(full_dt$DEG_binary)))

# ============================================================
# 12. KMER GLOBAL STABILITY
# ============================================================

kmer_global_stability <- analysis_dt[, .(
  n_deg_sets       = uniqueN(method),
  n_fastq_trimmers = uniqueN(fastq_trimmer[!is.na(fastq_trimmer)]),
  n_counters       = uniqueN(counter[!is.na(counter)]),
  n_aligners       = uniqueN(aligner[!is.na(aligner)]),
  n_fixed_values   = uniqueN(fixed_value[!is.na(fixed_value)]),
  mean_baseMean    = mean(baseMean, na.rm=TRUE)
), by=.(kmer, ensembl_id)]

kmer_global_stability[, kmer_stability_score :=
                        (n_fastq_trimmers / length(trimmers)) *
                        (n_counters       / length(counters)) *
                        (n_aligners       / length(aligners)) *
                        (n_deg_sets       / length(deg_files))]

setorder(kmer_global_stability, -kmer_stability_score)
fwrite(kmer_global_stability,
       file.path(output_dir, "tables",
                 paste0("kmer_global_stability_", today, ".tsv")))

gene_kmer_stability <- kmer_global_stability[, .(
  n_stable_kmers      = sum(kmer_stability_score >= 0.5),
  mean_kmer_stability = mean(kmer_stability_score),
  max_kmer_stability  = max(kmer_stability_score)
), by=ensembl_id]
fwrite(gene_kmer_stability,
       file.path(output_dir, "tables",
                 paste0("gene_kmer_stability_", today, ".tsv")))

# ============================================================
# 13. FINALIZED BIOMARKER OUTPUT
# ============================================================

biomarker_dt <- merge(gene_robustness_global, gene_kmer_stability,
                      by="ensembl_id", all=TRUE)

gene_detection <- full_dt[, .(
  detection_rate       = mean(DEG_binary),
  n_pipelines_detected = sum(DEG_binary),
  n_pipelines_total    = .N
), by=ensembl_id]
biomarker_dt <- merge(biomarker_dt, gene_detection, by="ensembl_id", all=TRUE)

gene_expr <- kmer_global_stability[, .(mean_baseMean=mean(mean_baseMean, na.rm=TRUE)),
                                   by=ensembl_id]
biomarker_dt <- merge(biomarker_dt, gene_expr, by="ensembl_id", all.x=TRUE)

biomarker_dt <- merge(biomarker_dt,
                      gene_cross_deg[, .(ensembl_id, cross_deg_fraction, n_aligners)],
                      by="ensembl_id", all.x=TRUE)
biomarker_dt[is.na(cross_deg_fraction), cross_deg_fraction := 0]

biomarker_dt[, biomarker_score :=
               as.numeric(scale(detection_rate))      +
               as.numeric(scale(mean_kmer_stability)) +
               as.numeric(scale(robustness))          +
               as.numeric(scale(cross_deg_fraction))]

setorder(biomarker_dt, -biomarker_score)

biomarker_dt[, tier := fcase(
  biomarker_score >  1, "Tier1_HighConfidence",
  biomarker_score >= 0, "Tier2_Moderate",
  biomarker_score <  0, "Tier3_LowConfidence"
)]

fwrite(biomarker_dt,
       file.path(output_dir, "tables",
                 paste0("FINALIZED_BIOMARKERS_", today, ".tsv")))
cat("Biomarker tier summary:\n")
print(biomarker_dt[, .N, by=tier])

png(file.path(output_dir, "figures",
              paste0("biomarker_score_distribution_", today, ".png")),
    width=2400, height=1800, res=300)
print(ggplot(biomarker_dt, aes(x=biomarker_score, fill=tier)) +
        geom_histogram(bins=60, colour="white") +
        theme_classic() +
        labs(title="Finalized biomarker score distribution",
             x="Biomarker score", y="Number of genes"))
dev.off()

png(file.path(output_dir, "figures",
              paste0("detection_vs_kmer_stability_", today, ".png")),
    width=2400, height=2000, res=300)
print(ggplot(biomarker_dt,
             aes(x=mean_kmer_stability, y=detection_rate, colour=tier)) +
        geom_point(alpha=0.5, size=1) +
        theme_classic() +
        labs(title="Detection rate vs kmer stability",
             x="Mean kmer stability score", y="Pipeline detection rate"))
dev.off()

png(file.path(output_dir, "figures",
              paste0("cross_deg_vs_detection_", today, ".png")),
    width=2400, height=2000, res=300)
print(ggplot(biomarker_dt,
             aes(x=cross_deg_fraction, y=detection_rate, colour=tier)) +
        geom_point(alpha=0.5, size=1) +
        theme_classic() +
        labs(title="Cross-DEG reproducibility vs pipeline detection rate",
             x="Cross-DEG fraction (from raw DEG lists)",
             y="Pipeline detection rate"))
dev.off()

# ============================================================
# 14. BAYESIAN HIERARCHICAL MODEL
# ============================================================

formula_full     <- NULL
formula_additive <- NULL
brms_data        <- NULL
fit_full         <- NULL
loo_full         <- NULL

if (!requireNamespace("brms", quietly=TRUE)) {
  cat("[SKIP] brms not available — skipping Bayesian hierarchical model\n")
} else {
  
  suppressPackageStartupMessages(library(brms))
  options(brms.backend = "cmdstanr")
  cat(sprintf("  brms backend: %s\n", getOption("brms.backend")))
  
  has_aligner  <- "aligner"       %in% names(full_dt) &
    uniqueN(full_dt$aligner[!is.na(full_dt$aligner)])             > 1
  has_trimmer  <- "fastq_trimmer" %in% names(full_dt) &
    uniqueN(full_dt$fastq_trimmer[!is.na(full_dt$fastq_trimmer)]) > 1
  has_counter  <- "counter"       %in% names(full_dt) &
    uniqueN(full_dt$counter[!is.na(full_dt$counter)])             > 1
  has_entropy  <- "kmer_entropy"  %in% names(full_dt) &
    sum(!is.na(full_dt$kmer_entropy) & is.finite(full_dt$kmer_entropy)) > 100
  has_gc       <- "gc_content"    %in% names(full_dt) &
    sum(!is.na(full_dt$gc_content)   & is.finite(full_dt$gc_content))   > 100
  has_length   <- "length"        %in% names(full_dt) &
    sum(!is.na(full_dt$length)       & is.finite(full_dt$length))       > 100
  
  cat("  Bayesian model dimensions:\n")
  cat(sprintf("    aligner      : %s (n=%d)\n", has_aligner,
              if (has_aligner) uniqueN(full_dt$aligner,       na.rm=TRUE) else 0))
  cat(sprintf("    fastq_trimmer: %s (n=%d)\n", has_trimmer,
              if (has_trimmer) uniqueN(full_dt$fastq_trimmer, na.rm=TRUE) else 0))
  cat(sprintf("    counter      : %s (n=%d)\n", has_counter,
              if (has_counter) uniqueN(full_dt$counter,       na.rm=TRUE) else 0))
  cat(sprintf("    kmer_entropy : %s\n", has_entropy))
  cat(sprintf("    gc_content   : %s\n", has_gc))
  cat(sprintf("    length       : %s\n", has_length))
  
  deg_frac <- mean(full_dt$DEG_binary, na.rm=TRUE)
  cat(sprintf("  DEG_binary: %d zeros | %d ones (%.1f%% positive)\n",
              sum(full_dt$DEG_binary == 0L, na.rm=TRUE),
              sum(full_dt$DEG_binary == 1L, na.rm=TRUE),
              100 * deg_frac))
  
  if (deg_frac == 0 || deg_frac == 1) {
    cat(sprintf("[SKIP] Bayesian model — no variance in DEG_binary (fraction=%.3f)\n",
                deg_frac))
  } else {
    
    model_cols <- unique(c(
      "DEG_binary", "ensembl_id",
      if (has_entropy) "kmer_entropy",
      if (has_gc)      "gc_content",
      if (has_length)  "length",
      if (has_aligner) "aligner",
      if (has_trimmer) "fastq_trimmer",
      if (has_counter) "counter"
    ))
    model_cols <- intersect(model_cols, names(full_dt))
    
    brms_data <- as.data.frame(full_dt[, model_cols, with=FALSE])
    
    if (has_aligner) brms_data$aligner       <- as.factor(brms_data$aligner)
    if (has_trimmer) brms_data$fastq_trimmer  <- as.factor(brms_data$fastq_trimmer)
    if (has_counter) brms_data$counter        <- as.factor(brms_data$counter)
    brms_data$ensembl_id <- as.factor(brms_data$ensembl_id)
    
    num_cols <- model_cols[sapply(brms_data[model_cols], is.numeric)]
    for (col in num_cols) {
      n_inf <- sum(!is.finite(brms_data[[col]]), na.rm=TRUE)
      if (n_inf > 0) {
        cat(sprintf("    %s: %d non-finite values replaced with NA\n", col, n_inf))
        brms_data[[col]][!is.finite(brms_data[[col]])] <- NA
      }
    }
    
    cat("  NAs per column before filtering:\n")
    for (col in model_cols) {
      n_na <- sum(is.na(brms_data[[col]]))
      if (n_na > 0) cat(sprintf("    %s: %d NAs\n", col, n_na))
    }
    
    brms_data <- brms_data[complete.cases(brms_data), ]
    
    zero_var_cols <- num_cols[sapply(num_cols, function(col) {
      v <- var(brms_data[[col]], na.rm=TRUE)
      is.na(v) || v == 0
    })]
    if (length(zero_var_cols) > 0) {
      cat(sprintf("  [WARNING] Zero-variance columns: %s — excluded from model\n",
                  paste(zero_var_cols, collapse=", ")))
      if ("kmer_entropy" %in% zero_var_cols) has_entropy <- FALSE
      if ("gc_content"   %in% zero_var_cols) has_gc      <- FALSE
      if ("length"       %in% zero_var_cols) has_length  <- FALSE
    }
    
    cat(sprintf("  brms data after cleaning: %d rows | DEG_binary=1: %d (%.1f%%)\n",
                nrow(brms_data),
                sum(brms_data$DEG_binary),
                100 * mean(brms_data$DEG_binary)))
    
    fixed_terms <- character(0)
    if (has_entropy) fixed_terms <- c(fixed_terms, "scale(kmer_entropy)")
    if (has_gc)      fixed_terms <- c(fixed_terms, "scale(gc_content)")
    if (has_length)  fixed_terms <- c(fixed_terms, "scale(length)")
    if (has_aligner) fixed_terms <- c(fixed_terms, "aligner")
    if (has_trimmer) fixed_terms <- c(fixed_terms, "fastq_trimmer")
    if (has_counter) fixed_terms <- c(fixed_terms, "counter")
    if (has_aligner && has_trimmer)
      fixed_terms <- c(fixed_terms, "aligner:fastq_trimmer")
    if (has_aligner && has_entropy)
      fixed_terms <- c(fixed_terms, "aligner:scale(kmer_entropy)")
    if (has_trimmer && has_entropy)
      fixed_terms <- c(fixed_terms, "fastq_trimmer:scale(kmer_entropy)")
    
    additive_terms <- fixed_terms[!grepl(":", fixed_terms)]
    
    formula_full     <- paste0("DEG_binary ~ ",
                               paste(fixed_terms,     collapse=" + "),
                               " + (1 | ensembl_id)")
    formula_additive <- paste0("DEG_binary ~ ",
                               paste(additive_terms, collapse=" + "),
                               " + (1 | ensembl_id)")
    
    cat(sprintf("  Full formula     : %s\n", formula_full))
    cat(sprintf("  Additive formula : %s\n", formula_additive))
    
    priors <- c(
      prior(normal(0, 2),   class="Intercept"),
      prior(normal(0, 1),   class="b"),
      prior(exponential(1), class="sd")
    )
    
    dir.create(file.path(output_dir, "models"), showWarnings=FALSE)
    brms_backend <- if (requireNamespace("cmdstanr", quietly=TRUE)) "cmdstanr" else "rstan"
    fit_full     <- NULL
    
    cat("  Fitting full Bayesian model (with interactions)...\n")
    tryCatch({
      fit_full <- brm(
        formula = as.formula(formula_full),
        data    = brms_data,
        family  = bernoulli(),
        prior   = priors,
        chains  = 4,
        cores   = min(4L, n_cores),
        iter    = 4000,
        warmup  = 1000,
        seed    = 42,
        backend = brms_backend,
        file    = file.path(output_dir, "models",
                            paste0("brms_full_", today))
      )
      cat("  Full model fitted successfully\n")
    }, error = function(e) {
      cat(sprintf("  [WARNING] Full model failed: %s\n", e$message))
    })
    
    if (is.null(fit_full)) {
      cat("  Fitting additive fallback model (no interactions)...\n")
      tryCatch({
        fit_full <<- brm(
          formula = as.formula(formula_additive),
          data    = brms_data,
          family  = bernoulli(),
          prior   = priors,
          chains  = 4,
          cores   = min(4L, n_cores),
          iter    = 4000,
          warmup  = 1000,
          seed    = 42,
          backend = brms_backend,
          file    = file.path(output_dir, "models",
                              paste0("brms_additive_", today))
        )
        cat("  Additive model fitted successfully\n")
      }, error = function(e) {
        cat(sprintf("  [ERROR] Additive model also failed: %s\n", e$message))
      })
    }
    
    if (is.null(fit_full) && (has_entropy || has_gc)) {
      cat("  Fitting minimal fallback model (sequence features only)...\n")
      minimal_terms   <- additive_terms[!grepl("aligner|fastq_trimmer|counter",
                                               additive_terms)]
      formula_minimal <- paste0("DEG_binary ~ ",
                                paste(minimal_terms, collapse=" + "),
                                " + (1 | ensembl_id)")
      cat(sprintf("  Minimal formula: %s\n", formula_minimal))
      tryCatch({
        fit_full <<- brm(
          formula = as.formula(formula_minimal),
          data    = brms_data,
          family  = bernoulli(),
          prior   = priors,
          chains  = 4,
          cores   = min(4L, n_cores),
          iter    = 4000,
          warmup  = 1000,
          seed    = 42,
          backend = brms_backend,
          file    = file.path(output_dir, "models",
                              paste0("brms_minimal_", today))
        )
        cat("  Minimal model fitted successfully\n")
      }, error = function(e) {
        cat(sprintf("  [ERROR] Minimal model also failed: %s\n", e$message))
      })
    }
    
    if (!is.null(fit_full)) {
      
      saveRDS(fit_full,
              file.path(output_dir, "models",
                        paste0("brms_full_", today, ".rds")))
      
      tryCatch({
        cat("\n  Bayesian model summary:\n")
        suppressPackageStartupMessages(library(posterior))
        draws <- as_draws_df(fit_full)
        summ  <- summarise_draws(draws,
                                 mean, sd,
                                 ~quantile(.x, probs=c(0.025, 0.975)),
                                 ~rhat(.x),
                                 ~ess_bulk(.x))
        print(summ[grepl("^b_", summ$variable), ])
      }, error = function(e) {
        cat(sprintf("  [WARNING] Summary failed: %s\n", e$message))
      })
      
      brms_coefs <- as.data.table(
        fixef(fit_full, summary=TRUE), keep.rownames="term")
      
      expected_names <- c("Estimate", "Est.Error", "Q2.5", "Q97.5")
      actual_names   <- names(brms_coefs)[names(brms_coefs) != "term"]
      cat(sprintf("  fixef columns: %s\n", paste(actual_names, collapse=", ")))
      
      if (all(expected_names %in% names(brms_coefs))) {
        setnames(brms_coefs,
                 c("Estimate", "Est.Error", "Q2.5", "Q97.5"),
                 c("estimate", "std_error", "ci_lower", "ci_upper"))
      } else {
        setnames(brms_coefs,
                 actual_names,
                 c("estimate", "std_error", "ci_lower", "ci_upper")[seq_along(actual_names)])
      }
      
      brms_coefs[, significant  := sign(ci_lower) == sign(ci_upper)]
      brms_coefs[, model        := "brms"]
      brms_coefs[, abs_estimate := abs(estimate)]
      setorder(brms_coefs, -abs_estimate)
      brms_coefs[, abs_estimate := NULL]
      
      fwrite(brms_coefs,
             file.path(output_dir, "tables",
                       paste0("brms_coefficients_", today, ".csv")))
      
      cat(sprintf("\n  Significant terms (95%% CI excludes 0): %d\n",
                  sum(brms_coefs$significant, na.rm=TRUE)))
      print(brms_coefs[significant == TRUE,
                       .(term, estimate, ci_lower, ci_upper)])
      
      p_brms <- ggplot(
        brms_coefs[term != "Intercept"],
        aes(x      = estimate,
            xmin   = ci_lower,
            xmax   = ci_upper,
            y      = reorder(term, estimate),
            colour = significant)) +
        geom_pointrange() +
        geom_vline(xintercept=0, linetype="dashed") +
        scale_colour_manual(
          values = c("FALSE"="grey60", "TRUE"="steelblue"),
          labels = c("CI includes 0", "CI excludes 0")) +
        theme_classic() +
        theme(axis.text.y=element_text(size=8)) +
        labs(title  = "Bayesian hierarchical model — pipeline dimension effects",
             x      = "Posterior estimate (log-odds) ± 95% CI",
             y      = "Term",
             colour = "Significant")
      ggsave(file.path(output_dir, "figures",
                       paste0("brms_forest_plot_", today, ".pdf")),
             p_brms, width=12,
             height=max(6, nrow(brms_coefs[term != "Intercept"]) * 0.4))
      
      tryCatch({
        ce_plots <- conditional_effects(fit_full)
        pdf(file.path(output_dir, "figures",
                      paste0("brms_conditional_effects_", today, ".pdf")),
            width=10, height=6)
        for (ce_name in names(ce_plots)) {
          print(plot(ce_plots, plot=FALSE)[[ce_name]] +
                  theme_classic() +
                  labs(title=paste("Conditional effect:", ce_name)))
        }
        dev.off()
      }, error = function(e) {
        cat(sprintf("  [WARNING] Conditional effects failed: %s\n", e$message))
        if (dev.cur() > 1) dev.off()
      })
      
      tryCatch({
        loo_full <- loo(fit_full)
        cat("\n  LOO cross-validation:\n")
        print(loo_full)
        saveRDS(loo_full,
                file.path(output_dir, "models",
                          paste0("brms_loo_", today, ".rds")))
      }, error = function(e) {
        cat(sprintf("  [WARNING] LOO failed: %s\n", e$message))
      })
      
    } else {
      cat("[SKIP] All Bayesian model attempts failed — check data diagnostics above\n")
    }
  }
}

# ============================================================
# 15. MODEL COMPARISON (WAIC)
# ============================================================

if (exists("fit_full") && !is.null(fit_full)) {
  tryCatch({
    waic_full <- waic(fit_full)
    saveRDS(waic_full, file.path(output_dir, "models",
                                 paste0("waic_full_", today, ".rds")))
  }, error = function(e) cat(sprintf("[WARNING] WAIC failed: %s\n", e$message)))
} else {
  cat("[SKIP] Section 15 — fit_full not available\n")
}

# ============================================================
# 16. K-FOLD CROSS-VALIDATION
# ============================================================

if (exists("fit_full") && !is.null(fit_full)) {
  tryCatch({
    kfold_res <- kfold(fit_full, K=5)
    saveRDS(kfold_res, file.path(output_dir, "models",
                                 paste0("kfold_", today, ".rds")))
  }, error = function(e) cat(sprintf("[WARNING] K-fold failed: %s\n", e$message)))
} else {
  cat("[SKIP] Section 16 — fit_full not available\n")
}

# ============================================================
# 17. LEAVE-ONE-ALIGNER-OUT VALIDATION
# ============================================================

if (exists("fit_full") && !is.null(fit_full) && exists("formula_full")) {
  unique_aligners <- unique(full_dt$aligner[!is.na(full_dt$aligner)])
  brms_backend    <- if (requireNamespace("cmdstanr", quietly=TRUE)) "cmdstanr" else "rstan"
  loo_align <- lapply(unique_aligners, function(a) {
    tryCatch({
      brm(as.formula(formula_full),
          data    = brms_data[brms_data$aligner != a, ],
          family  = bernoulli(),
          chains  = 2,
          cores   = min(4L, n_cores),
          iter    = 2000,
          seed    = 42,
          backend = brms_backend)
    }, error = function(e) {
      cat(sprintf("  [WARNING] LOO aligner %s failed: %s\n", a, e$message))
      NULL
    })
  })
  names(loo_align) <- unique_aligners
  loo_align <- loo_align[!sapply(loo_align, is.null)]
  saveRDS(loo_align, file.path(output_dir, "models",
                               paste0("leave_one_aligner_out_", today, ".rds")))
} else {
  cat("[SKIP] Section 17 — fit_full or formula_full not available\n")
}

# ============================================================
# 18. ROC / AUC PER ALIGNER
# ============================================================

if (exists("fit_full") && !is.null(fit_full)) {
  tryCatch({
    suppressPackageStartupMessages(library(pROC))
    fitted_probs <- fitted(fit_full)[, "Estimate"]
    full_dt[seq_len(nrow(brms_data)), fitted := fitted_probs]
    
    unique_aligners <- unique(full_dt$aligner[!is.na(full_dt$aligner)])
    roc_list <- lapply(unique_aligners, function(a) {
      sub <- full_dt[aligner == a & !is.na(fitted)]
      if (nrow(sub) < 10) return(NULL)
      roc(sub$DEG_binary, sub$fitted, quiet=TRUE)
    })
    names(roc_list) <- unique_aligners
    roc_list <- roc_list[!sapply(roc_list, is.null)]
    
    pdf(file.path(output_dir, "figures",
                  paste0("ROC_comparison_", today, ".pdf")),
        width=8, height=6)
    cols <- RColorBrewer::brewer.pal(min(8, length(roc_list)), "Set1")
    for (i in seq_along(roc_list)) {
      if (i == 1) plot(roc_list[[i]],  col=cols[i], main="ROC — per aligner")
      else        lines(roc_list[[i]], col=cols[i])
    }
    legend("bottomright", names(roc_list), col=cols, lwd=2, cex=0.8)
    dev.off()
    
    auc_dt <- data.table(
      aligner = names(roc_list),
      auc     = sapply(roc_list, function(r) as.numeric(auc(r)))
    )
    fwrite(auc_dt, file.path(output_dir, "tables",
                             paste0("roc_auc_per_aligner_", today, ".csv")))
    cat("  AUC per aligner:\n")
    print(auc_dt)
  }, error = function(e) {
    cat(sprintf("[WARNING] ROC/AUC failed: %s\n", e$message))
  })
} else {
  cat("[SKIP] Section 18 — fit_full not available\n")
}

# ============================================================
# 19. VARIANCE DECOMPOSITION
# ============================================================

if (exists("fit_full") && !is.null(fit_full)) {
  tryCatch({
    vc    <- VarCorr(fit_full)
    vc_dt <- as.data.table(as.data.frame(vc))
    fwrite(vc_dt, file.path(output_dir, "tables",
                            paste0("variance_decomposition_", today, ".csv")))
    saveRDS(vc, file.path(output_dir, "models",
                          paste0("variance_components_", today, ".rds")))
    cat("\n  Variance components:\n")
    print(vc)
  }, error = function(e) {
    cat(sprintf("[WARNING] Variance decomposition failed: %s\n", e$message))
  })
} else {
  cat("[SKIP] Section 19 — fit_full not available\n")
}

# ============================================================
# 20. SHAP-STYLE VARIANCE ATTRIBUTION
# ============================================================

if (exists("fit_full") && !is.null(fit_full)) {
  tryCatch({
    suppressPackageStartupMessages({
      library(posterior)
      library(matrixStats)
    })
    posterior_draws <- as_draws_df(fit_full)
    draw_matrix     <- as.matrix(posterior_draws[, grep("^b_", names(posterior_draws))])
    var_components  <- colVars(draw_matrix)
    var_dt <- data.table(
      parameter = names(var_components),
      variance  = var_components,
      sd        = sqrt(var_components)
    )
    setorder(var_dt, -variance)
    fwrite(var_dt, file.path(output_dir, "tables",
                             paste0("variance_attribution_", today, ".tsv")),
           sep="\t")
    cat("\n  Top variance components:\n")
    print(head(var_dt, 10))
  }, error = function(e) {
    cat(sprintf("[WARNING] Variance attribution failed: %s\n", e$message))
  })
} else {
  cat("[SKIP] Section 20 — fit_full not available\n")
}

# ============================================================
# 21. HIERARCHICAL CLUSTERING
# ============================================================

if (nrow(full_dt) > 0) {
  tryCatch({
    clust_cols <- intersect(c("aligner", "fastq_trimmer", "counter"),
                            names(full_dt))
    full_dt[, combo := do.call(paste, c(.SD, sep="_")), .SDcols=clust_cols]
    
    clust_mat <- dcast(full_dt, ensembl_id ~ combo,
                       value.var="DEG_binary", fill=0)
    clust_matrix <- as.matrix(clust_mat[, -"ensembl_id", with=FALSE])
    rownames(clust_matrix) <- clust_mat$ensembl_id
    
    dist_mat <- dist(clust_matrix, method="binary")
    hc       <- hclust(dist_mat, method="ward.D2")
    
    pdf(file.path(output_dir, "figures",
                  paste0("hierarchical_dendrogram_", today, ".pdf")),
        width=12, height=8)
    plot(hc, labels=clust_mat$ensembl_id, cex=0.5,
         main="Gene clustering by pipeline detection pattern")
    dev.off()
    
    k         <- min(8L, nrow(clust_matrix) - 1L)
    clusters  <- cutree(hc, k=k)
    cluster_dt <- data.table(ensembl_id=names(clusters), cluster=clusters)
    fwrite(cluster_dt, file.path(output_dir, "tables",
                                 paste0("gene_clusters_", today, ".csv")))
  }, error = function(e) {
    cat(sprintf("[WARNING] Hierarchical clustering failed: %s\n", e$message))
  })
} else {
  cat("[SKIP] Section 21 — full_dt empty\n")
}

# ============================================================
# 22. POSTERIOR PREDICTIVE CHECKS
# ============================================================

if (exists("fit_full") && !is.null(fit_full)) {
  tryCatch({
    dir.create(file.path(output_dir, "diagnostics"), showWarnings=FALSE)
    pdf(file.path(output_dir, "diagnostics",
                  paste0("posterior_predictive_", today, ".pdf")),
        width=10, height=6)
    print(pp_check(fit_full, ndraws=100))
    dev.off()
  }, error = function(e) {
    cat(sprintf("[WARNING] Posterior predictive check failed: %s\n", e$message))
    if (dev.cur() > 1) dev.off()
  })
} else {
  cat("[SKIP] Section 22 — fit_full not available\n")
}

# ============================================================
# 23. INTERACTION VISUALIZATION
# ============================================================

if (exists("fit_full") && !is.null(fit_full)) {
  tryCatch({
    model_terms     <- rownames(fixef(fit_full))
    has_interaction <- any(grepl("aligner.*fastq_trimmer|fastq_trimmer.*aligner",
                                 model_terms))
    if (has_interaction) {
      pdf(file.path(output_dir, "figures",
                    paste0("interaction_plots_", today, ".pdf")),
          width=12, height=8)
      print(plot(conditional_effects(
        fit_full, effects="aligner:fastq_trimmer"), plot=FALSE)[[1]] +
          theme_classic() +
          labs(title="Aligner x FASTQ trimmer interaction"))
      dev.off()
    } else {
      cat("[SKIP] Section 23 — no aligner:fastq_trimmer interaction in model\n")
    }
  }, error = function(e) {
    cat(sprintf("[WARNING] Interaction plot failed: %s\n", e$message))
    if (dev.cur() > 1) dev.off()
  })
} else {
  cat("[SKIP] Section 23 — fit_full not available\n")
}

# ============================================================
# 24. PDF ROBUSTNESS REPORT WITH METHOD SIMILARITY HEATMAPS
# ============================================================

tryCatch({
  dir.create(file.path(output_dir, "figures"), showWarnings=FALSE)
  pdf(file.path(output_dir,
                paste0("Robustness_Report_with_Heatmaps_", today, ".pdf")),
      width=14, height=10)
  
  if (exists("all_summaries_dt") && !is.null(all_summaries_dt) &&
      nrow(all_summaries_dt) > 0 &&
      all(c("method", "deg_set", "unique_genes") %in% names(all_summaries_dt))) {
    heatmap_mat <- reshape2::acast(
      all_summaries_dt, method ~ deg_set, value.var="unique_genes")
    heatmap_mat[is.na(heatmap_mat)] <- 0
    draw(Heatmap(heatmap_mat,
                 name         = "Unique Genes",
                 column_title = "Unique genes per method x DEG set",
                 col          = colorRamp2(
                   c(0, max(heatmap_mat, na.rm=TRUE)),
                   c("white", "steelblue"))))
  }
  
  if (exists("consensus") && !is.null(consensus) &&
      all(c("total_genes", "method") %in% names(consensus))) {
    barplot(consensus$total_genes,
            names.arg = consensus$method,
            las=2, main="Total genes per method", col="steelblue")
  }
  
  if (exists("all_kmer_gene_lists") && length(all_kmer_gene_lists) > 0) {
    for (dn in names(all_kmer_gene_lists)) {
      jac_mat <- compute_jaccard(all_kmer_gene_lists[[dn]])
      if (!is.null(jac_mat) && nrow(jac_mat) > 1) {
        draw(Heatmap(jac_mat,
                     name         = "Jaccard",
                     column_title = paste0("Method similarity: ", dn),
                     col          = colorRamp2(c(0,1),
                                               c("white","darkgreen"))))
      }
    }
  }
  
  dev.off()
}, error = function(e) {
  cat(sprintf("[WARNING] Robustness report failed: %s\n", e$message))
  if (dev.cur() > 1) dev.off()
})

# ============================================================
# 25. ZSTD COMPRESSION
# ============================================================

compress_outputs <- function(path) {
  files <- list.files(path, full.names=TRUE, recursive=FALSE,
                      pattern="\\.(csv|tsv|rds)$")
  files <- files[!grepl("\\.zst$", files)]
  cat(sprintf("  Compressing %d files in %s\n", length(files), path))
  for (f in files) {
    tryCatch({
      ret <- system(paste("zstd -19 -f --rm", shQuote(f)), intern=FALSE)
      if (ret != 0) cat(sprintf("  [WARNING] zstd failed for: %s\n", f))
    }, error = function(e) {
      cat(sprintf("  [WARNING] Compression error for %s: %s\n", f, e$message))
    })
  }
}

compress_outputs(file.path(output_dir, "tables"))
compress_outputs(file.path(output_dir, "models"))

# ============================================================
# 26. CLEAN SHUTDOWN
# ============================================================
cat("Shutting down parallel workers...\n")
tryCatch({
  future::plan(sequential)
  gc()
}, error = function(e) {
  cat(sprintf("[INFO] Parallel shutdown (non-critical): %s\n", e$message))
})

cat("\n============================================================\n")
cat("Pipeline completed successfully.\n")
cat(sprintf("Completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("============================================================\n")
cat("Key outputs:\n")
cat(sprintf("  Finalized biomarkers            : %s\n",
            file.path(output_dir, "tables", paste0("FINALIZED_BIOMARKERS_",          today, ".tsv"))))
cat(sprintf("  All DEG tables combined         : %s\n",
            file.path(output_dir, "tables", paste0("all_deg_tables_combined_",       today, ".csv"))))
cat(sprintf("  Gene cross-DEG reproducibility  : %s\n",
            file.path(output_dir, "tables", paste0("gene_cross_deg_reproducibility_",today, ".csv"))))
cat(sprintf("  Cross-DEG per aligner           : %s\n",
            file.path(output_dir, "tables", paste0("gene_cross_deg_per_aligner_",    today, ".csv"))))
cat(sprintf("  Cross-DEG per counter           : %s\n",
            file.path(output_dir, "tables", paste0("gene_cross_deg_per_counter_",    today, ".csv"))))
cat(sprintf("  Cross-DEG per fixed value       : %s\n",
            file.path(output_dir, "tables", paste0("gene_cross_deg_per_fixed_value_",today, ".csv"))))
cat(sprintf("  Kmer global stability           : %s\n",
            file.path(output_dir, "tables", paste0("kmer_global_stability_",         today, ".tsv"))))
cat(sprintf("  Gene kmer stability             : %s\n",
            file.path(output_dir, "tables", paste0("gene_kmer_stability_",           today, ".tsv"))))
cat(sprintf("  Motif enrichment (all DEG sets) : %s\n",
            file.path(output_dir, "tables", paste0("motif_enrichment_all_degsets_",  today, ".csv"))))
cat(sprintf("  Sub-kmer enrichment (all DEG)   : %s\n",
            file.path(output_dir, "tables", paste0("subkmer_enrichment_all_degsets_",today, ".csv"))))
cat(sprintf("  Global fastq sensitivity        : %s\n",
            file.path(output_dir, "tables", paste0("global_fastq_sensitivity_summary_",today, ".csv"))))
cat(sprintf("  GLMM combined results           : %s\n",
            file.path(output_dir, "tables", paste0("glmm_combined_across_degsets_",  today, ".csv"))))
cat(sprintf("  GLMM consistency                : %s\n",
            file.path(output_dir, "tables", paste0("glmm_consistency_perDEG_vs_combined_",today, ".csv"))))
if (!is.null(fit_full)) {
  cat(sprintf("  Bayesian model                  : %s\n",
              file.path(output_dir, "models", paste0("brms_full_",                   today, ".rds"))))
  cat(sprintf("  Bayesian coefficients           : %s\n",
              file.path(output_dir, "tables", paste0("brms_coefficients_",           today, ".csv"))))
}
cat(sprintf("  Robustness ranking              : %s\n",
            file.path(output_dir, "diagnostics", paste0("ROBUSTNESS_RANKING_",       today, ".csv"))))
