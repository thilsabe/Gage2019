### TEST VERSION — RNApreprocessing_gene_kmer_analysisAndrobustness_pipeline.R
# Generates synthetic data to verify parsing, merging, analyses, and output paths
# without requiring HPC paths, real files, or Bayesian model runs.
#
# Run with:
#   Rscript TEST_pipeline.R
#
# All outputs go to: ./TEST_OUTPUT/

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(igraph)
  library(ComplexHeatmap)
  library(circlize)
  library(UpSetR)
  library(reshape2)
  library(matrixStats)
})

cat("============================================================\n")
cat("  PIPELINE TEST — synthetic data\n")
cat("============================================================\n\n")

# ============================================================
# 0. TEST CONFIGURATION
# ============================================================

n_cores <- 2

# Override all HPC paths with local test directories
deg_sets_dir       <- "./TEST_INPUT/deg_sets"
kmer_deseq_dir     <- "./TEST_INPUT/kmer_deseq"
kmer_map_file      <- "./TEST_INPUT/kmer_gene_map.csv"
kmer_features_file <- "./TEST_INPUT/kmer_features.csv"
base_output_dir    <- "./TEST_OUTPUT/kmer_analysis_outputs"
output_dir         <- "./TEST_OUTPUT/FINAL_ALIGNER_KMER_ANALYSIS"

today <- format(Sys.Date(), "%Y-%m-%d")

# Create all output directories
for (d in c(
  deg_sets_dir, kmer_deseq_dir,
  base_output_dir, output_dir,
  file.path(output_dir, "figures"),
  file.path(output_dir, "tables"),
  file.path(output_dir, "models"),
  file.path(output_dir, "diagnostics")
)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

# Pipeline dimension labels (small subsets for speed)
trimmers <- c("fastp","trimmomatic","trimgalore","untrimmed")
aligners  <- c("hisat2","star_single","star_twopass")
counters  <- c("htseq","homer_gene","homer_gene_condensed_unique")

cat("[CONFIG] Test directories created under ./TEST_INPUT and ./TEST_OUTPUT\n\n")

# ============================================================
# 1. GENERATE SYNTHETIC INPUT FILES
# ============================================================

set.seed(42)

n_genes <- 30
n_kmers <- 60

gene_ids <- paste0("ENSG", sprintf("%011d", 1:n_genes))
kmer_ids <- paste0("kmer_", sprintf("%04d", 1:n_kmers))

cat("[SYNTH] Generating synthetic reference files...\n")

# kmer → gene map (each kmer maps to 1–3 genes)
kmer_gene_map_syn <- data.table(
  kmer       = rep(kmer_ids, each=2),
  ensembl_id = sample(gene_ids, n_kmers * 2, replace=TRUE)
)
kmer_gene_map_syn <- unique(kmer_gene_map_syn)
fwrite(kmer_gene_map_syn, kmer_map_file)
cat("  Written:", kmer_map_file, "\n")

# kmer features
kmer_features_syn <- data.table(
  kmer         = kmer_ids,
  gc_content   = runif(n_kmers, 0.3, 0.7),
  kmer_entropy = runif(n_kmers, 1.5, 3.0),
  splice_motif = sample(c(TRUE, FALSE), n_kmers, replace=TRUE),
  gene_specific = sample(c(TRUE, FALSE), n_kmers, replace=TRUE),
  position     = sample(c("exon","intron","junction"), n_kmers, replace=TRUE)
)
fwrite(kmer_features_syn, kmer_features_file)
cat("  Written:", kmer_features_file, "\n")

# Generate DEG set files and matching kmer DESeq2 files
# Filename pattern: all_DEGs_<fastq_trimmer>_..._within_<deg_list_trimmer>_<counter>.csv
#   col 1 = aligner, last col = ensembl_id

test_combos <- expand.grid(
  fastq_trimmer    = trimmers,
  deg_list_trimmer = c("fastp","trimmomatic"),
  counter          = counters,
  stringsAsFactors = FALSE
)

cat("\n[SYNTH] Generating", nrow(test_combos), "DEG set + kmer file pairs...\n")

for (i in seq_len(nrow(test_combos))) {
  ftrim  <- test_combos$fastq_trimmer[i]
  dltrim <- test_combos$deg_list_trimmer[i]
  ctr    <- test_combos$counter[i]

  # Build filename matching parse_metadata expectations
  fname  <- paste0("all_DEGs_", ftrim, "_overlap_within_", dltrim, "_", ctr, ".csv")

  # DEG file: col 1 = aligner, last col = ensembl_id
  n_degs <- sample(8:15, 1)
  deg_syn <- data.table(
    aligner    = sample(aligners, n_degs, replace=TRUE),
    log2FC     = rnorm(n_degs),
    padj       = runif(n_degs, 0, 0.05),
    ensembl_id = sample(gene_ids, n_degs, replace=TRUE)
  )
  fwrite(deg_syn, file.path(deg_sets_dir, fname))

  # Matching kmer DESeq2 file (same filename, different directory)
  n_k <- sample(20:35, 1)
  kmer_syn <- data.table(
    kmer       = sample(kmer_ids, n_k, replace=TRUE),
    baseMean   = runif(n_k, 10, 1000),
    log2FC     = rnorm(n_k),
    padj       = runif(n_k, 0, 0.05),
    GC         = runif(n_k, 0.3, 0.7),
    length     = sample(20:150, n_k, replace=TRUE),
    fastq_version = ftrim   # column representing which fastq version produced this kmer
  )
  fwrite(kmer_syn, file.path(kmer_deseq_dir, fname))
}

cat("  DEG files written to:  ", deg_sets_dir, "\n")
cat("  Kmer files written to: ", kmer_deseq_dir, "\n\n")

# ============================================================
# 2. LOAD REFERENCE FILES
# ============================================================

cat("[STEP 1] Loading reference files...\n")
kmer_gene_map <- fread(kmer_map_file)
kmer_features <- fread(kmer_features_file)
cat("  kmer_gene_map rows:", nrow(kmer_gene_map), "\n")
cat("  kmer_features rows:", nrow(kmer_features), "\n\n")

# ============================================================
# 3. PARSE FILENAME METADATA — with parse verification
# ============================================================

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

  counter_match <- parts[parts %in% counters]
  counter       <- if (length(counter_match) > 0) counter_match[1] else "none"

  list(
    fastq_trimmer    = fastq_trimmer,
    grouping_type    = grouping_type,
    deg_list_trimmer = deg_list_trimmer,
    counter          = counter
  )
}

cat("[STEP 2] Verifying parse_metadata on all filenames...\n")
deg_files <- list.files(deg_sets_dir, pattern="\\.csv$", full.names=TRUE)
cat("  Found", length(deg_files), "DEG files\n")

parse_check <- rbindlist(lapply(deg_files, function(f) {
  m <- parse_metadata(f)
  data.table(
    file             = basename(f),
    fastq_trimmer    = m$fastq_trimmer,
    deg_list_trimmer = m$deg_list_trimmer %||% NA_character_,
    counter          = m$counter,
    grouping_type    = m$grouping_type
  )
}))
cat("\n  Parse results (first 10 rows):\n")
print(parse_check[1:min(10, .N)])

# Check for unexpected NAs or "none" counters
n_bad_counter <- sum(parse_check$counter == "none")
n_bad_dlt     <- sum(is.na(parse_check$deg_list_trimmer))
cat("\n  Files with counter='none'          :", n_bad_counter,
    if (n_bad_counter == 0) "[OK]" else "[WARNING]", "\n")
cat("  Files with missing deg_list_trimmer:", n_bad_dlt,
    if (n_bad_dlt == 0) "[OK]" else "[WARNING]", "\n\n")

# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

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

match_kmer_file <- function(deg_file) {
  file.path(kmer_deseq_dir, basename(deg_file))
}

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
                     total_bg, total_obs,
                     lower.tail = FALSE))
  })
  res_df <- as.data.table(do.call(rbind, res))
  res_df[, `-log10p` := -log10(as.numeric(pval) + 1e-12)]
  return(res_df[order(as.numeric(pval))])
}

compute_jaccard <- function(method_list) {
  mnames <- names(method_list)
  mat <- matrix(0, nrow=length(mnames), ncol=length(mnames),
                dimnames=list(mnames, mnames))
  for (i in seq_along(mnames))
    for (j in seq_along(mnames)) {
      inter     <- length(intersect(method_list[[i]], method_list[[j]]))
      union_sz  <- length(unique(c(method_list[[i]], method_list[[j]])))
      mat[i,j]  <- if (union_sz > 0) inter / union_sz else 0
    }
  return(mat)
}

# ============================================================
# 5. MAIN LOADING + PER-DEG-SET ANALYSIS LOOP
# ============================================================

cat("[STEP 3] Running main loading loop...\n\n")

full_dataset          <- list()
gene_tracker          <- list()
all_summaries         <- list()
kmers_by_method       <- list()
all_kmer_gene_lists   <- list()
fastq_sensitivity_all <- list()

for (deg_file in deg_files) {

  meta     <- parse_metadata(deg_file)
  deg_name <- gsub("\\.csv$", "", basename(deg_file))

  deg_dt    <- load_deg_aligners(deg_file)
  kmer_file <- match_kmer_file(deg_file)
  if (!file.exists(kmer_file)) {
    cat("  [SKIP]", deg_name, "— kmer file missing\n")
    next
  }

  kmer_dt <- fread(kmer_file)
  kmer_dt <- merge(kmer_dt, kmer_gene_map, by="kmer", all.x=TRUE)
  kmer_dt <- merge(kmer_dt, kmer_features, by="kmer", all.x=TRUE)

  merged  <- merge(kmer_dt, deg_dt, by="ensembl_id", allow.cartesian=TRUE)

  if (nrow(merged) == 0) {
    cat("  [SKIP]", deg_name, "— empty merge result\n")
    next
  }

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

  # 5a. FASTQ-TRIMMER SENSITIVITY vs DEG-LIST TRIMMER
  fastq_sens <- kmers_dt[, .(
    n_genes = uniqueN(ensembl_id),
    n_kmers = uniqueN(kmer),
    n_rows  = .N
  ), by=.(deg_list_trimmer, fastq_trimmer, aligner)]
  fastq_sens[, genes_norm := n_genes / max(n_genes, na.rm=TRUE), by=deg_list_trimmer]
  fastq_sens[, kmers_norm := n_kmers / max(n_kmers, na.rm=TRUE), by=deg_list_trimmer]
  fwrite(fastq_sens,
         file.path(deg_out,
                   paste0("fastq_trimmer_sensitivity_vs_deglist_trimmer_", today, ".csv")))
  fastq_sensitivity_all[[deg_name]] <- fastq_sens

  p_sens <- ggplot(fastq_sens, aes(x=fastq_trimmer, y=genes_norm, fill=fastq_trimmer)) +
    geom_col() +
    facet_grid(aligner ~ deg_list_trimmer, scales="free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(title=paste0("FASTQ-trimmer sensitivity\n", deg_name),
         y="Normalised gene recovery", x="FASTQ trimmer")
  ggsave(file.path(deg_out, paste0("fastq_trimmer_sensitivity_", today, ".pdf")),
         p_sens, width=10, height=7)

  # 5b. KMER & GENE STABILITY ACROSS FASTQ TRIMMERS
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

  p_stab <- ggplot(gene_stability, aes(x=stability_index, fill=deg_list_trimmer)) +
    geom_histogram(binwidth=1/length(trimmers), boundary=0, colour="white") +
    facet_wrap(~ counter, scales="free_y") +
    theme_classic() +
    labs(title=paste0("Gene stability across FASTQ trimmers\n", deg_name),
         x="Stability index", y="Number of genes")
  ggsave(file.path(deg_out, paste0("gene_stability_histogram_", today, ".pdf")),
         p_stab, width=8, height=5)

  # 5c. MOTIF ENRICHMENT
  if ("gene_specific" %in% names(kmers_dt) && sum(kmers_dt$gene_specific, na.rm=TRUE) > 0) {
    motif_res <- motif_enrichment(kmers_dt[gene_specific == TRUE], kmers_dt)
    fwrite(motif_res,
           file.path(deg_out, paste0("motif_enrichment_", today, ".csv")))
  }

  # 5d. CROSS-DEG REPRODUCIBILITY
  reproducibility <- kmers_dt[, .(N = .N), by=.(ensembl_id)]
  reproducibility[, reproducibility_score := N / length(unique(kmers_dt$method))]
  fwrite(reproducibility,
         file.path(deg_out, paste0("cross_deg_reproducibility_", today, ".csv")))

  # 5e. GENE-STRATIFIED K-MER NETWORK + PERTURBATION ROBUSTNESS
  net_edges <- kmers_dt[!is.na(ensembl_id), .(kmer, ensembl_id)]
  if (nrow(net_edges) > 1) {
    g <- graph_from_data_frame(net_edges, directed=FALSE)
    saveRDS(g, file.path(deg_out, paste0("gene_kmer_network_", today, ".rds")))
    network_robustness <- replicate(20, {   # 20 replicates in test (100 in production)
      remove <- sample(V(g)$name, max(1, floor(vcount(g) * 0.1)))
      g2     <- delete_vertices(g, remove)
      length(components(g2)$no)
    })
    fwrite(data.table(n_components = network_robustness),
           file.path(deg_out, paste0("network_robustness_", today, ".csv")))
  }

  # 5f. UPSET PLOT (skip if only 1 method so far)
  if (length(kmers_by_method) > 1) {
    pdf(file.path(deg_out, paste0("upset_methods_", today, ".pdf")),
        width=10, height=6)
    upset(fromList(kmers_by_method))
    dev.off()
  }

  # 5g. GC-CONTENT DISTRIBUTION
  if ("gc_content" %in% names(kmers_dt)) {
    p_gc <- ggplot(kmers_dt, aes(gc_content, fill=fastq_trimmer)) +
      geom_density(alpha=0.4) +
      facet_wrap(~ deg_list_trimmer) +
      theme_classic() +
      labs(title=paste0("GC-content by FASTQ trimmer\n", deg_name))
    suppressMessages(ggsave(
      file.path(deg_out, paste0("gc_density_", today, ".pdf")), p_gc))
  }

  # 5h. SUMMARY TABLE
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

cat("\n[STEP 3] Loop complete.\n")
cat("  DEG sets processed:", length(full_dataset), "/", length(deg_files), "\n\n")

analysis_dt <- rbindlist(full_dataset, fill=TRUE)
cat("  analysis_dt rows:", nrow(analysis_dt),
    "| cols:", ncol(analysis_dt), "\n\n")

# ============================================================
# 6. GLOBAL FASTQ-TRIMMER SENSITIVITY SUMMARY
# ============================================================

cat("[STEP 4] Global fastq-trimmer sensitivity summary...\n")

fastq_sensitivity_dt <- rbindlist(fastq_sensitivity_all, idcol="deg_set")

global_fastq_sens <- fastq_sensitivity_dt[, .(
  mean_genes_norm = mean(genes_norm, na.rm=TRUE),
  mean_kmers_norm = mean(kmers_norm, na.rm=TRUE),
  sd_genes_norm   = sd(genes_norm,   na.rm=TRUE)
), by=.(deg_list_trimmer, fastq_trimmer)]

fwrite(global_fastq_sens,
       file.path(base_output_dir,
                 paste0("global_fastq_sensitivity_summary_", today, ".csv")))
cat("  Rows:", nrow(global_fastq_sens), "\n")
print(global_fastq_sens)

p_global_sens <- ggplot(global_fastq_sens,
                         aes(x=fastq_trimmer, y=mean_genes_norm,
                             ymin=mean_genes_norm - sd_genes_norm,
                             ymax=mean_genes_norm + sd_genes_norm,
                             colour=fastq_trimmer)) +
  geom_pointrange() +
  facet_wrap(~ deg_list_trimmer) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title="Global FASTQ-trimmer sensitivity by DEG-list trimmer",
       y="Mean normalised gene recovery", x="FASTQ trimmer")
ggsave(file.path(base_output_dir,
                 paste0("global_fastq_sensitivity_", today, ".pdf")),
       p_global_sens, width=10, height=5)
cat("  [OK]\n\n")

# ============================================================
# 7. GENE ROBUSTNESS STABILITY CURVES
# ============================================================

cat("[STEP 5] Gene robustness stability curves...\n")

gene_freq <- table(unlist(gene_tracker))
gene_robustness_global <- data.table(
  ensembl_id = names(gene_freq),
  freq        = as.numeric(gene_freq)
)
gene_robustness_global[, robustness := freq / length(gene_tracker)]

cat("  Genes tracked:", nrow(gene_robustness_global), "\n")
cat("  Robustness range: [",
    round(min(gene_robustness_global$robustness), 3), ",",
    round(max(gene_robustness_global$robustness), 3), "]\n")

png(file.path(output_dir, paste0("GeneStabilityCurve_", today, ".png")),
    width=2400, height=2000, res=300)
print(
  ggplot(gene_robustness_global, aes(x=robustness)) +
    geom_density(fill="steelblue", alpha=0.5) +
    theme_classic() +
    labs(title="Gene robustness across all DEG sets",
         x="Robustness (fraction of DEG sets gene appeared in)")
)
dev.off()
cat("  [OK]\n\n")

# ============================================================
# 8. CROSS-DEG CONSENSUS & ROBUSTNESS RANKING
# ============================================================

cat("[STEP 6] Cross-DEG consensus and robustness ranking...\n")

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
cat("  Consensus rows:", nrow(consensus), "\n")
cat("  [OK]\n\n")

# ============================================================
# 9. FULL GENE × PIPELINE COMBINATION MATRIX
# ============================================================

cat("[STEP 7] Building full gene x pipeline matrix...\n")

all_genes    <- unique(analysis_dt$ensembl_id)
all_aligners <- unique(analysis_dt$aligner)
all_trimrs   <- unique(analysis_dt$fastq_trimmer)
all_counters <- unique(analysis_dt$counter)

cat("  Genes:", length(all_genes),
    "| Aligners:", length(all_aligners),
    "| Trimmers:", length(all_trimrs),
    "| Counters:", length(all_counters), "\n")
cat("  Expected matrix rows:", length(all_genes) * length(all_aligners) *
                               length(all_trimrs) * length(all_counters), "\n")

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

cat("  full_dt rows:", nrow(full_dt),
    "| DEG_binary=1:", sum(full_dt$DEG_binary),
    sprintf("(%.1f%%)", 100 * mean(full_dt$DEG_binary)), "\n")

reserved_cols  <- c("ensembl_id","aligner","fastq_trimmer","counter",
                    "DEG_binary","method","grouping_type","deg_list_trimmer","kmer")
gene_meta_cols <- setdiff(names(analysis_dt), reserved_cols)
gene_meta      <- unique(analysis_dt[, c("ensembl_id", gene_meta_cols), with=FALSE])
full_dt        <- merge(full_dt, gene_meta, by="ensembl_id", all.x=TRUE)
cat("  full_dt cols after gene meta merge:", ncol(full_dt), "\n")
cat("  [OK]\n\n")

# ============================================================
# 10. KMER GLOBAL STABILITY
# ============================================================

cat("[STEP 8] Kmer global stability...\n")

kmer_global_stability <- analysis_dt[, .(
  n_deg_sets          = uniqueN(method),
  n_fastq_trimmers    = uniqueN(fastq_trimmer),
  n_counters          = uniqueN(counter),
  n_aligners          = uniqueN(aligner),
  n_deg_list_trimmers = uniqueN(deg_list_trimmer),
  mean_baseMean       = mean(baseMean, na.rm=TRUE)
), by=.(kmer, ensembl_id)]

kmer_global_stability[, kmer_stability_score :=
  (n_fastq_trimmers    / length(trimmers)) *
  (n_counters          / length(counters)) *
  (n_aligners          / length(aligners)) *
  (n_deg_sets          / length(deg_files))]

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

cat("  kmer_global_stability rows:", nrow(kmer_global_stability), "\n")
cat("  Stability score range: [",
    round(min(kmer_global_stability$kmer_stability_score), 4), ",",
    round(max(kmer_global_stability$kmer_stability_score), 4), "]\n")
cat("  Kmers with score >= 0.5:",
    sum(kmer_global_stability$kmer_stability_score >= 0.5), "\n")
cat("  [OK]\n\n")

# ============================================================
# 11. FINALIZED BIOMARKER OUTPUT
# ============================================================

cat("[STEP 9] Finalizing biomarker table...\n")

biomarker_dt <- merge(gene_robustness_global, gene_kmer_stability,
                      by="ensembl_id", all=TRUE)

gene_detection <- full_dt[, .(
  detection_rate       = mean(DEG_binary),
  n_pipelines_detected = sum(DEG_binary),
  n_pipelines_total    = .N
), by=ensembl_id]

biomarker_dt <- merge(biomarker_dt, gene_detection, by="ensembl_id", all=TRUE)

gene_expr <- kmer_global_stability[,
  .(mean_baseMean = mean(mean_baseMean, na.rm=TRUE)), by=ensembl_id]
biomarker_dt <- merge(biomarker_dt, gene_expr, by="ensembl_id", all.x=TRUE)

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

cat("  Biomarker table rows:", nrow(biomarker_dt), "\n")
cat("  Tier summary:\n")
print(biomarker_dt[, .N, by=tier])
cat("  Top 5 biomarkers:\n")
print(biomarker_dt[1:5, .(ensembl_id, robustness, mean_kmer_stability,
                           detection_rate, biomarker_score, tier)])

png(file.path(output_dir, "figures",
              paste0("biomarker_score_distribution_", today, ".png")),
    width=1800, height=1400, res=200)
print(
  ggplot(biomarker_dt, aes(x=biomarker_score, fill=tier)) +
    geom_histogram(bins=40, colour="white") +
    theme_classic() +
    labs(title="Finalized biomarker score distribution",
         x="Biomarker score", y="Number of genes")
)
dev.off()

png(file.path(output_dir, "figures",
              paste0("detection_vs_kmer_stability_", today, ".png")),
    width=1800, height=1400, res=200)
print(
  ggplot(biomarker_dt,
         aes(x=mean_kmer_stability, y=detection_rate, colour=tier)) +
    geom_point(alpha=0.6, size=2) +
    theme_classic() +
    labs(title="Detection rate vs kmer stability",
         x="Mean kmer stability score", y="Pipeline detection rate")
)
dev.off()
cat("  [OK]\n\n")

# ============================================================
# 12. ROBUSTNESS REPORT HEATMAPS (no Bayesian model in test)
# ============================================================

cat("[STEP 10] Robustness report heatmaps...\n")

pdf(file.path(base_output_dir,
              paste0("Robustness_Report_with_Heatmaps_", today, ".pdf")),
    width=14, height=10)

heatmap_mat <- reshape2::acast(all_summaries_dt,
                               method ~ deg_set,
                               value.var = "unique_genes",
                               fill      = 0)
draw(Heatmap(heatmap_mat,
             name         = "Unique Genes",
             column_title = "Unique genes per method x DEG set",
             col          = colorRamp2(
               c(min(heatmap_mat, na.rm=TRUE), max(heatmap_mat, na.rm=TRUE)),
               c("white","steelblue"))))

barplot(consensus$total_genes,
        names.arg = consensus$method,
        las=2, main="Total genes per method", col="steelblue")

for (dn in names(all_kmer_gene_lists)) {
  jac_mat <- compute_jaccard(all_kmer_gene_lists[[dn]])
  draw(Heatmap(jac_mat,
               name         = "Jaccard",
               column_title = paste0("Method similarity: ", dn),
               col          = colorRamp2(c(0,1), c("white","darkgreen"))))
}

dev.off()
cat("  [OK]\n\n")

# NOTE: Bayesian model (brms), LOO, kfold, leave-one-aligner-out, and
# ROC/AUC sections are skipped in this test version. Run the production
# script for those after verifying outputs here.

# ============================================================
# FINAL OUTPUT INVENTORY
# ============================================================

cat("============================================================\n")
cat("  TEST COMPLETE — Output inventory\n")
cat("============================================================\n\n")

all_outputs <- list.files(c(base_output_dir, output_dir),
                           recursive=TRUE, full.names=TRUE)
cat("Total files written:", length(all_outputs), "\n\n")

output_summary <- data.table(
  path = all_outputs,
  size_kb = round(file.size(all_outputs) / 1024, 1)
)
output_summary[, dir := dirname(path)]
output_summary[, file := basename(path)]
print(output_summary[, .(file, size_kb)][order(size_kb, decreasing=TRUE)])

cat("\nKey outputs to inspect:\n")
key_files <- c(
  file.path(output_dir, "tables", paste0("FINALIZED_BIOMARKERS_",       today, ".tsv")),
  file.path(output_dir, "tables", paste0("kmer_global_stability_",      today, ".tsv")),
  file.path(output_dir, "tables", paste0("gene_kmer_stability_",        today, ".tsv")),
  file.path(base_output_dir,      paste0("global_fastq_sensitivity_summary_", today, ".csv")),
  file.path(base_output_dir,      paste0("ROBUSTNESS_RANKING_",         today, ".csv"))
)
for (f in key_files) {
  exists_flag <- if (file.exists(f)) "[OK]" else "[MISSING]"
  cat(" ", exists_flag, f, "\n")
}

cat("\nAll parse_metadata checks:\n")
print(parse_check)
