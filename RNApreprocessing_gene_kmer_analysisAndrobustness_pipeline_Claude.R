### "Final" Pipeline 2-20-2026
# Saved as RNApreprocessing_gene_kmer_analysisAndrobustness_pipeline.R
# Run with:
# Rscript RNApreprocessing_gene_kmer_analysisAndrobustness_pipeline.R

#!/usr/bin/env Rscript

# ============================================================
# FULL ROBUSTNESS PIPELINE
# HPC-Optimized + Parallel + Bayesian + Publication Ready
# ============================================================

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
})

# ============================================================
# 0. CONFIGURATION
# ============================================================

n_cores <- 8   # increase to 32 for full HPC runs

# Input directories and files
deg_sets_dir       <- "/cjc/data/Gage_align_bamFiles/count_DEG_analysis/WithConfounds/HOMER/overlap_DEGs/DEG_sets"
kmer_deseq_dir     <- "/netapp/snl/scratch25/AHAllen_projects/Gage2019/kmerator_results/deseq2_results"
kmer_map_file      <- "kmer_rMATS_annotations/kmer_gene_map_from_gtf_20260212.csv"
kmer_features_file <- "kmer_features/kmer_splice_motif_position_features_20260212.csv"

# Output directories
base_output_dir <- "/netapp/snl/scratch25/AHAllen_projects/Gage2019/kmer_analysis_outputs"
output_dir      <- "/netapp/snl/scratch25/AHAllen_projects/Gage2019/FINAL_ALIGNER_KMER_ANALYSIS"

today <- format(Sys.Date(), "%Y-%m-%d")

dir.create(base_output_dir,                        recursive=TRUE, showWarnings=FALSE)
dir.create(output_dir,                             recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output_dir, "figures"),       showWarnings=FALSE)
dir.create(file.path(output_dir, "tables"),        showWarnings=FALSE)
dir.create(file.path(output_dir, "models"),        showWarnings=FALSE)
dir.create(file.path(output_dir, "diagnostics"),   showWarnings=FALSE)

# Known pipeline dimension labels — used for filename parsing and for
# building the complete gene × pipeline combination matrix downstream
trimmers <- c("fastp","trimmomatic","trimgalore","untrimmed")
aligners  <- c("hisat2","star_single","star_twopass","gsnap","gsnap_2024genome")
counters  <- c("htseq","star","homer_exon","homer_exon_unique","homer_gene",
               "homer_gene_condensed","homer_gene_condensed_unique","homer_gene_unique")

# ============================================================
# HPC PARALLEL CONFIG
# ============================================================

options(mc.cores = n_cores)
plan(multicore, workers = n_cores)

# ============================================================
# 1. LOAD REFERENCE FILES
# ============================================================

kmer_gene_map <- fread(kmer_map_file)
kmer_features <- fread(kmer_features_file)

# ============================================================
# 2. PARSE FILENAME METADATA
# ============================================================
# Expected filename pattern (underscore-delimited):
#   <prefix>_<something>_<fastq_trimmer>_..._within_<deg_list_trimmer>_..._<counter>...csv
#
# fastq_trimmer    : trimmer applied to the FASTQ before alignment/counting
# deg_list_trimmer : trimmer whose DEG list was used as the overlap reference
#                    (token immediately after "within" in filename)
# grouping_type    : which pipeline axis is being varied in this file
# counter          : first token matching the known counters vector

parse_metadata <- function(fname) {
  parts <- str_split(basename(fname), "_")[[1]]

  fastq_trimmer <- parts[3]

  grouping_type <- ifelse(grepl("aligners", fname), "aligner",
                   ifelse(grepl("trimmers", fname), "trimmer",
                   ifelse(grepl("counters", fname), "counter", "unknown")))

  within_idx       <- which(parts == "within")
  deg_list_trimmer <- if (length(within_idx) > 0 && within_idx + 1 <= length(parts)) {
    parts[within_idx + 1]
  } else {
    NA_character_
  }

  # Match against the known counters vector
  counter_match <- parts[parts %in% counters]
  counter       <- if (length(counter_match) > 0) counter_match[1] else "none"

  list(
    fastq_trimmer    = fastq_trimmer,
    grouping_type    = grouping_type,
    deg_list_trimmer = deg_list_trimmer,
    counter          = counter
  )
}

# ============================================================
# 3. LOAD DEG FILE (aligner = col 1, gene = last col)
# ============================================================

load_deg_aligners <- function(deg_file) {
  dt          <- fread(deg_file)
  aligner_col <- names(dt)[1]
  gene_col    <- names(dt)[ncol(dt)]
  deg_long    <- dt[, .(
    aligner    = get(aligner_col),
    ensembl_id = get(gene_col)
  )]
  unique(deg_long)
}

# ============================================================
# 4. MATCH KMER FILE
# ============================================================

match_kmer_file <- function(deg_file) {
  file.path(kmer_deseq_dir, basename(deg_file))
}

# ============================================================
# 5. HELPER FUNCTIONS
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

motif_enrichment <- function(kmer_df, background_kmers) {
  obs       <- table(kmer_df$kmer)
  bg_counts <- table(background_kmers$kmer)
  total_bg  <- length(background_kmers$kmer)
  total_obs <- nrow(kmer_df)
  res <- lapply(names(obs), function(k) {
    c(kmer  = k,
      count = obs[k],
      pval  = phyper(obs[k] - 1,
                     bg_counts[k] %||% 1,
                     total_bg,
                     total_obs,
                     lower.tail = FALSE))
  })
  res_df <- as.data.table(do.call(rbind, res))
  res_df[, `-log10p` := -log10(as.numeric(pval) + 1e-12)]
  return(res_df[order(as.numeric(pval))])
}

compute_jaccard <- function(method_list) {
  mnames <- names(method_list)
  mat    <- matrix(0,
                   nrow     = length(mnames),
                   ncol     = length(mnames),
                   dimnames = list(mnames, mnames))
  for (i in seq_along(mnames)) {
    for (j in seq_along(mnames)) {
      a          <- method_list[[i]]
      b          <- method_list[[j]]
      inter      <- length(intersect(a, b))
      union_size <- length(unique(c(a, b)))
      mat[i, j]  <- ifelse(union_size > 0, inter / union_size, 0)
    }
  }
  return(mat)
}

# ============================================================
# 6. MAIN LOADING + PER-DEG-SET ANALYSIS LOOP
# ============================================================

deg_files           <- list.files(deg_sets_dir, pattern="\\.csv$", full.names=TRUE)
full_dataset        <- list()
gene_tracker        <- list()
all_summaries       <- list()
kmers_by_method     <- list()
all_kmer_gene_lists <- list()
fastq_sensitivity_all <- list()

for (deg_file in deg_files) {

  meta     <- parse_metadata(deg_file)
  deg_name <- gsub("\\.csv$", "", basename(deg_file))
  cat("Processing:", deg_name, "\n")
  cat("  fastq_trimmer   :", meta$fastq_trimmer,    "\n")
  cat("  deg_list_trimmer:", meta$deg_list_trimmer,  "\n")
  cat("  counter         :", meta$counter,           "\n")
  cat("  grouping_type   :", meta$grouping_type,     "\n")

  deg_dt    <- load_deg_aligners(deg_file)
  kmer_file <- match_kmer_file(deg_file)
  if (!file.exists(kmer_file)) {
    cat("  [SKIP] kmer file not found:", kmer_file, "\n")
    next
  }

  kmer_dt <- fread(kmer_file)
  kmer_dt <- merge(kmer_dt, kmer_gene_map, by="kmer", all.x=TRUE)
  kmer_dt <- merge(kmer_dt, kmer_features, by="kmer", all.x=TRUE)

  # Cross-reference with DEG genes per aligner
  merged <- merge(kmer_dt, deg_dt, by="ensembl_id", allow.cartesian=TRUE)

  merged[, `:=`(
    fastq_trimmer    = meta$fastq_trimmer,
    deg_list_trimmer = meta$deg_list_trimmer,
    counter          = meta$counter,
    grouping_type    = meta$grouping_type,
    method           = deg_name,
    DEG_binary       = 1L
  )]

  full_dataset[[deg_name]]        <- merged
  gene_tracker[[deg_name]]        <- unique(merged$ensembl_id)
  kmers_by_method[[deg_name]]     <- unique(merged$kmer)
  all_kmer_gene_lists[[deg_name]] <- split(merged$ensembl_id, merged$aligner)

  deg_out <- file.path(base_output_dir, deg_name)
  dir.create(deg_out, recursive=TRUE, showWarnings=FALSE)

  kmers_dt <- merged

  # -------------------------------------------------------
  # 6a. FASTQ-TRIMMER SENSITIVITY vs DEG-LIST TRIMMER
  #
  # For each deg_list_trimmer context, measures how many genes
  # and kmers are recovered when the FASTQ was trimmed by each
  # fastq_trimmer. Reveals whether FASTQ trimming interacts with
  # the trimmer context of the DEG list being queried.
  # -------------------------------------------------------
  fastq_sens <- kmers_dt[, .(
    n_genes = uniqueN(ensembl_id),
    n_kmers = uniqueN(kmer),
    n_rows  = .N
  ), by=.(deg_list_trimmer, fastq_trimmer, aligner)]

  fastq_sens[, genes_norm := n_genes / max(n_genes, na.rm=TRUE),
             by=deg_list_trimmer]
  fastq_sens[, kmers_norm := n_kmers / max(n_kmers, na.rm=TRUE),
             by=deg_list_trimmer]

  fwrite(fastq_sens,
         file.path(deg_out,
                   paste0("fastq_trimmer_sensitivity_vs_deglist_trimmer_", today, ".csv")))

  p_sens <- ggplot(fastq_sens,
                   aes(x=fastq_trimmer, y=genes_norm, fill=fastq_trimmer)) +
    geom_col() +
    facet_grid(aligner ~ deg_list_trimmer, scales="free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(title = paste0("FASTQ-trimmer sensitivity per DEG-list trimmer\n", deg_name),
         y = "Normalised gene recovery", x = "FASTQ trimmer")
  ggsave(file.path(deg_out, paste0("fastq_trimmer_sensitivity_", today, ".pdf")),
         p_sens, width=12, height=8)

  fastq_sensitivity_all[[deg_name]] <- fastq_sens

  # -------------------------------------------------------
  # 6b. KMER STABILITY ACROSS FASTQ TRIMMERS
  #
  # Within a fixed DEG list (deg_list_trimmer + counter),
  # stable kmers/genes appear regardless of which fastq_trimmer
  # was used. Stability index = fraction of trimmers detected in.
  # -------------------------------------------------------
  gene_stability <- kmers_dt[, .(
    n_trimmers_detected = uniqueN(fastq_trimmer)
  ), by=.(ensembl_id, deg_list_trimmer, counter)]
  gene_stability[, stability_index := n_trimmers_detected / length(trimmers)]
  fwrite(gene_stability,
         file.path(deg_out, paste0("gene_stability_across_fastq_", today, ".csv")))

  kmer_stability <- kmers_dt[, .(
    n_trimmers_detected = uniqueN(fastq_trimmer)
  ), by=.(kmer, ensembl_id, deg_list_trimmer, counter)]
  kmer_stability[, stability_index := n_trimmers_detected / length(trimmers)]
  fwrite(kmer_stability,
         file.path(deg_out, paste0("kmer_stability_across_fastq_", today, ".csv")))

  p_stab <- ggplot(gene_stability,
                   aes(x=stability_index, fill=deg_list_trimmer)) +
    geom_histogram(binwidth=1/length(trimmers), boundary=0, colour="white") +
    facet_wrap(~ counter, scales="free_y") +
    theme_classic() +
    labs(title = paste0("Gene stability across FASTQ trimmers\n", deg_name),
         x = "Stability index (fraction of trimmers detected)",
         y = "Number of genes")
  ggsave(file.path(deg_out, paste0("gene_stability_histogram_", today, ".pdf")),
         p_stab, width=10, height=6)

  # -------------------------------------------------------
  # 6c. MOTIF ENRICHMENT
  # -------------------------------------------------------
  if ("gene_specific" %in% names(kmers_dt)) {
    motif_res <- motif_enrichment(kmers_dt[gene_specific == TRUE], kmers_dt)
    fwrite(motif_res,
           file.path(deg_out, paste0("motif_enrichment_", today, ".csv")))
  }

  # -------------------------------------------------------
  # 6d. CROSS-DEG REPRODUCIBILITY
  # -------------------------------------------------------
  reproducibility <- kmers_dt[, .(N = .N), by=.(ensembl_id)]
  reproducibility[, reproducibility_score := N / length(unique(kmers_dt$method))]
  fwrite(reproducibility,
         file.path(deg_out, paste0("cross_deg_reproducibility_", today, ".csv")))

  # -------------------------------------------------------
  # 6e. GENE-STRATIFIED K-MER NETWORK + PERTURBATION ROBUSTNESS
  # -------------------------------------------------------
  g <- graph_from_data_frame(
    kmers_dt[!is.na(ensembl_id), .(kmer, ensembl_id)],
    directed = FALSE
  )
  saveRDS(g, file.path(deg_out, paste0("gene_kmer_network_", today, ".rds")))

  network_robustness <- replicate(100, {
    remove <- sample(V(g)$name, floor(vcount(g) * 0.1))
    g2     <- delete_vertices(g, remove)
    length(components(g2)$no)
  })
  fwrite(data.table(n_components = network_robustness),
         file.path(deg_out, paste0("network_robustness_", today, ".csv")))

  # -------------------------------------------------------
  # 6f. METHOD INTERSECTION (UpSet)
  # -------------------------------------------------------
  pdf(file.path(deg_out, paste0("upset_methods_", today, ".pdf")),
      width=10, height=6)
  upset(fromList(kmers_by_method))
  dev.off()

  # -------------------------------------------------------
  # 6g. GC-CONTENT DISTRIBUTION
  # -------------------------------------------------------
  if ("gc_content" %in% names(kmers_dt)) {
    p_gc <- ggplot(kmers_dt, aes(gc_content, fill=fastq_trimmer)) +
      geom_density(alpha=0.4) +
      facet_wrap(~ deg_list_trimmer) +
      theme_classic() +
      labs(title = paste0("GC-content by FASTQ trimmer\n", deg_name))
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
  ), by=.(method, fastq_trimmer, deg_list_trimmer, counter)]

  fwrite(summary_dt,
         file.path(deg_out, paste0("method_summary_", today, ".csv")))
  all_summaries[[deg_name]] <- summary_dt
}

# Combine all DEG sets into single analysis table
analysis_dt <- rbindlist(full_dataset, fill=TRUE)

# ============================================================
# 7. GLOBAL FASTQ-TRIMMER SENSITIVITY SUMMARY
#
# Across all DEG sets: for each deg_list_trimmer context, which
# fastq_trimmer recovers the most genes/kmers consistently?
# ============================================================

fastq_sensitivity_dt <- rbindlist(fastq_sensitivity_all, idcol="deg_set")

global_fastq_sens <- fastq_sensitivity_dt[, .(
  mean_genes_norm = mean(genes_norm, na.rm=TRUE),
  mean_kmers_norm = mean(kmers_norm, na.rm=TRUE),
  sd_genes_norm   = sd(genes_norm,   na.rm=TRUE)
), by=.(deg_list_trimmer, fastq_trimmer)]

fwrite(global_fastq_sens,
       file.path(base_output_dir,
                 paste0("global_fastq_sensitivity_summary_", today, ".csv")))

p_global_sens <- ggplot(global_fastq_sens,
                         aes(x=fastq_trimmer, y=mean_genes_norm,
                             ymin=mean_genes_norm - sd_genes_norm,
                             ymax=mean_genes_norm + sd_genes_norm,
                             colour=fastq_trimmer)) +
  geom_pointrange() +
  facet_wrap(~ deg_list_trimmer) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title = "Global FASTQ-trimmer sensitivity by DEG-list trimmer",
       y = "Mean normalised gene recovery", x = "FASTQ trimmer")
ggsave(file.path(base_output_dir,
                 paste0("global_fastq_sensitivity_", today, ".pdf")),
       p_global_sens, width=12, height=6)

# ============================================================
# 8. GENE ROBUSTNESS STABILITY CURVES
# ============================================================

gene_freq <- table(unlist(gene_tracker))

gene_robustness_global <- data.table(
  ensembl_id = names(gene_freq),
  freq        = as.numeric(gene_freq)
)
gene_robustness_global[, robustness := freq / length(gene_tracker)]

png(file.path(output_dir, paste0("GeneStabilityCurve_", today, ".png")),
    width=2400, height=2000, res=600)
print(
  ggplot(gene_robustness_global, aes(x=robustness)) +
    geom_density(fill="steelblue", alpha=0.5) +
    theme_classic() +
    labs(title = "Gene robustness across all DEG sets",
         x = "Robustness (fraction of DEG sets gene appeared in)")
)
dev.off()

# ============================================================
# 9. CROSS-DEG CONSENSUS & ROBUSTNESS RANKING
# ============================================================

all_summaries_dt <- rbindlist(all_summaries, idcol="deg_set")

consensus <- all_summaries_dt[, .(
  total_sig   = sum(n_sig),
  total_genes = sum(unique_genes),
  avg_gc      = mean(mean_gc,       na.rm=TRUE),
  avg_expr    = mean(mean_baseMean, na.rm=TRUE)
), by=.(method, counter, deg_list_trimmer)]

consensus[, robustness_score :=
            as.numeric(scale(total_genes)) -
            abs(as.numeric(scale(avg_gc))) -
            abs(as.numeric(scale(avg_expr)))]

fwrite(consensus,
       file.path(base_output_dir,
                 paste0("ROBUSTNESS_RANKING_", today, ".csv")))

# ============================================================
# 10. BUILD FULL GENE × PIPELINE COMBINATION MATRIX
#     All dimensions inferred from parsed data — no external
#     metadata file required.
# ============================================================

all_genes    <- unique(analysis_dt$ensembl_id)
all_aligners <- unique(analysis_dt$aligner)
all_trimrs   <- unique(analysis_dt$fastq_trimmer)
all_counters <- unique(analysis_dt$counter)

full_dt <- CJ(
  ensembl_id = all_genes,
  aligner    = all_aligners,
  trimmer    = all_trimrs,
  counter    = all_counters,
  unique     = TRUE
)
full_dt[, DEG_binary := 0L]

full_dt[
  analysis_dt,
  on = .(ensembl_id, aligner, trimmer = fastq_trimmer, counter),
  DEG_binary := 1L
]

# Attach any gene-level features from analysis_dt
reserved_cols <- c("ensembl_id","aligner","fastq_trimmer","counter",
                   "DEG_binary","method","grouping_type","deg_list_trimmer","kmer")
gene_meta_cols <- setdiff(names(analysis_dt), reserved_cols)
gene_meta      <- unique(analysis_dt[, c("ensembl_id", gene_meta_cols), with=FALSE])
full_dt        <- merge(full_dt, gene_meta, by="ensembl_id", all.x=TRUE)

# ============================================================
# 11. KMER STABILITY — GLOBAL (across all DEG sets)
#
# A kmer is globally stable if it is detected across many
# fastq_trimmer × counter × aligner × deg_set combinations.
# ============================================================

kmer_global_stability <- analysis_dt[, .(
  n_deg_sets          = uniqueN(method),
  n_fastq_trimmers    = uniqueN(fastq_trimmer),
  n_counters          = uniqueN(counter),
  n_aligners          = uniqueN(aligner),
  n_deg_list_trimmers = uniqueN(deg_list_trimmer),
  mean_baseMean       = mean(baseMean, na.rm=TRUE)
), by=.(kmer, ensembl_id)]

# Composite stability score: product of fractional coverage per dimension
kmer_global_stability[, kmer_stability_score :=
  (n_fastq_trimmers    / length(trimmers))  *
  (n_counters          / length(counters))  *
  (n_aligners          / length(aligners))  *
  (n_deg_sets          / length(deg_files))]

setorder(kmer_global_stability, -kmer_stability_score)

fwrite(kmer_global_stability,
       file.path(output_dir, "tables",
                 paste0("kmer_global_stability_", today, ".tsv")))

# Gene-level aggregation of kmer stability
gene_kmer_stability <- kmer_global_stability[, .(
  n_stable_kmers      = sum(kmer_stability_score >= 0.5),
  mean_kmer_stability = mean(kmer_stability_score),
  max_kmer_stability  = max(kmer_stability_score)
), by=ensembl_id]

fwrite(gene_kmer_stability,
       file.path(output_dir, "tables",
                 paste0("gene_kmer_stability_", today, ".tsv")))

# ============================================================
# 12. FINALIZED BIOMARKER OUTPUT
#
# Integrates:
#   - global gene robustness (fraction of DEG sets)
#   - kmer stability (fraction of pipeline combinations)
#   - pipeline detection rate across full binary matrix
#   - mean expression
# ============================================================

biomarker_dt <- merge(
  gene_robustness_global,
  gene_kmer_stability,
  by  = "ensembl_id",
  all = TRUE
)

gene_detection <- full_dt[, .(
  detection_rate       = mean(DEG_binary),
  n_pipelines_detected = sum(DEG_binary),
  n_pipelines_total    = .N
), by=ensembl_id]

biomarker_dt <- merge(biomarker_dt, gene_detection, by="ensembl_id", all=TRUE)

gene_expr <- kmer_global_stability[,
  .(mean_baseMean = mean(mean_baseMean, na.rm=TRUE)),
  by=ensembl_id]
biomarker_dt <- merge(biomarker_dt, gene_expr, by="ensembl_id", all.x=TRUE)

# Composite biomarker score
biomarker_dt[, biomarker_score :=
  as.numeric(scale(detection_rate))      +
  as.numeric(scale(mean_kmer_stability)) +
  as.numeric(scale(robustness))]

setorder(biomarker_dt, -biomarker_score)

biomarker_dt[, tier := fcase(
  biomarker_score >  1, "Tier1_HighConfidence",
  biomarker_score >= 0, "Tier2_Moderate",
  biomarker_score <  0, "Tier3_LowConfidence"
)]

fwrite(biomarker_dt,
       file.path(output_dir, "tables",
                 paste0("FINALIZED_BIOMARKERS_", today, ".tsv")))

cat("Biomarker summary:\n")
print(biomarker_dt[, .N, by=tier])

png(file.path(output_dir, "figures",
              paste0("biomarker_score_distribution_", today, ".png")),
    width=2400, height=1800, res=300)
print(
  ggplot(biomarker_dt, aes(x=biomarker_score, fill=tier)) +
    geom_histogram(bins=60, colour="white") +
    theme_classic() +
    labs(title = "Finalized biomarker score distribution",
         x = "Biomarker score", y = "Number of genes")
)
dev.off()

png(file.path(output_dir, "figures",
              paste0("detection_vs_kmer_stability_", today, ".png")),
    width=2400, height=2000, res=300)
print(
  ggplot(biomarker_dt,
         aes(x=mean_kmer_stability, y=detection_rate, colour=tier)) +
    geom_point(alpha=0.5, size=1) +
    theme_classic() +
    labs(title = "Detection rate vs kmer stability",
         x = "Mean kmer stability score",
         y = "Pipeline detection rate")
)
dev.off()

# ============================================================
# 13. BAYESIAN HIERARCHICAL MODEL
# ============================================================

formula_full <- bf(
  DEG_binary ~
    kmer_entropy + GC + length +
    aligner + trimmer + counter +
    aligner:trimmer +
    aligner:counter +
    trimmer:counter +
    aligner:trimmer:counter +
    (1 | ensembl_id),
  family = bernoulli()
)

fit_full <- brm(
  formula_full,
  data    = full_dt,
  chains  = 4,
  cores   = n_cores,
  iter    = 4000,
  backend = "cmdstanr"
)

saveRDS(fit_full, file.path(output_dir, "models/full_model.rds"))

# ============================================================
# 14. MODEL COMPARISON (LOO / WAIC)
# ============================================================

loo_full  <- loo(fit_full)
waic_full <- waic(fit_full)

saveRDS(loo_full,  file.path(output_dir, "models/loo_full.rds"))
saveRDS(waic_full, file.path(output_dir, "models/waic_full.rds"))

# ============================================================
# 15. K-FOLD CROSS-VALIDATION
# ============================================================

kfold_res <- kfold(fit_full, K=10)
saveRDS(kfold_res, file.path(output_dir, "models/kfold.rds"))

# ============================================================
# 16. LEAVE-ONE-ALIGNER-OUT VALIDATION
# ============================================================

unique_aligners <- unique(full_dt$aligner)

loo_align <- future_lapply(unique_aligners, function(a) {
  sub_dt <- full_dt[aligner != a]
  brm(formula_full, data=sub_dt, chains=2, cores=2, iter=2000, backend="cmdstanr")
})
names(loo_align) <- unique_aligners

saveRDS(loo_align, file.path(output_dir, "models/leave_one_aligner_out.rds"))

# ============================================================
# 17. ROC / AUC PER ALIGNER
# ============================================================

fitted_probs      <- fitted(fit_full)[, 1]
full_dt[, fitted := fitted_probs]

roc_list <- lapply(unique_aligners, function(a) {
  sub <- full_dt[aligner == a]
  roc(sub$DEG_binary, sub$fitted)
})
names(roc_list) <- unique_aligners

pdf(file.path(output_dir, "figures/ROC_comparison.pdf"))
for (a in unique_aligners) plot(roc_list[[a]], main=paste("ROC —", a))
dev.off()

# ============================================================
# 18. VARIANCE DECOMPOSITION
# ============================================================

vc <- VarCorr(fit_full)
saveRDS(vc, file.path(output_dir, "tables/variance_decomposition.rds"))

# ============================================================
# 19. SHAP-STYLE VARIANCE ATTRIBUTION
# ============================================================

posterior_draws <- as_draws_df(fit_full)
var_components  <- colVars(as.matrix(posterior_draws))

fwrite(
  data.table(parameter = names(var_components), variance = var_components),
  file.path(output_dir, "tables/variance_attribution.tsv")
)

# ============================================================
# 20. HIERARCHICAL CLUSTERING
# ============================================================

clust_mat <- dcast(full_dt,
                   ensembl_id ~ aligner + trimmer + counter,
                   value.var = "DEG_binary",
                   fill      = 0)
dist_mat <- dist(as.matrix(clust_mat[, -1]))
hc       <- hclust(dist_mat)

pdf(file.path(output_dir, "figures/hierarchical_dendrogram.pdf"))
plot(hc, labels=FALSE, main="Gene clustering by pipeline detection pattern")
dev.off()

# ============================================================
# 21. POSTERIOR PREDICTIVE CHECKS
# ============================================================

pdf(file.path(output_dir, "diagnostics/posterior_predictive.pdf"))
pp_check(fit_full)
dev.off()

# ============================================================
# 22. INTERACTION VISUALIZATION
# ============================================================

pdf(file.path(output_dir, "figures/interaction_plots.pdf"))
plot(conditional_effects(fit_full, effects="aligner:trimmer"))
dev.off()

# ============================================================
# 23. PDF ROBUSTNESS REPORT WITH METHOD SIMILARITY HEATMAPS
# ============================================================

pdf(file.path(base_output_dir,
              paste0("Robustness_Report_with_Heatmaps_", today, ".pdf")),
    width=14, height=10)

heatmap_mat <- reshape2::acast(all_summaries_dt,
                               method ~ deg_set,
                               value.var = "unique_genes")
draw(Heatmap(heatmap_mat,
             name         = "Unique Genes",
             column_title = "Unique genes per method x DEG set",
             col          = colorRamp2(
               c(min(heatmap_mat, na.rm=TRUE), max(heatmap_mat, na.rm=TRUE)),
               c("white","steelblue"))))

barplot(consensus$total_genes,
        names.arg = consensus$method,
        las  = 2,
        main = "Total genes per method",
        col  = "steelblue")

for (deg_name in names(all_kmer_gene_lists)) {
  jac_mat <- compute_jaccard(all_kmer_gene_lists[[deg_name]])
  draw(Heatmap(jac_mat,
               name         = "Jaccard",
               column_title = paste0("Method similarity: ", deg_name),
               col          = colorRamp2(c(0, 1), c("white","darkgreen"))))
}

dev.off()

# ============================================================
# 24. ZSTD COMPRESSION
# ============================================================

compress_files <- function(path) {
  files <- list.files(path, full.names=TRUE, recursive=FALSE)
  future_lapply(files, function(f) system(paste("zstd -19 -f", shQuote(f))))
}

compress_files(file.path(output_dir, "tables"))
compress_files(file.path(output_dir, "models"))

# ============================================================
# DONE
# ============================================================

cat("\nPipeline completed successfully.\n")
cat("Key outputs:\n")
cat("  Finalized biomarkers :",
    file.path(output_dir, "tables", paste0("FINALIZED_BIOMARKERS_",       today, ".tsv")), "\n")
cat("  Kmer global stability:",
    file.path(output_dir, "tables", paste0("kmer_global_stability_",      today, ".tsv")), "\n")
cat("  Gene kmer stability  :",
    file.path(output_dir, "tables", paste0("gene_kmer_stability_",        today, ".tsv")), "\n")
cat("  Global fastq sens    :",
    file.path(base_output_dir,      paste0("global_fastq_sensitivity_summary_", today, ".csv")), "\n")
cat("  Robustness ranking   :",
    file.path(base_output_dir,      paste0("ROBUSTNESS_RANKING_",         today, ".csv")), "\n")
