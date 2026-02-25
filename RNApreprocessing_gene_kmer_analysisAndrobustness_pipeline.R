#!/usr/bin/env Rscript
# RNApreprocessing_gene_kmer_analysisAndrobustness_pipeline.R
# Full Robustness Pipeline — HPC-Optimized + Parallel + Bayesian + Publication Ready

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
})

# Prevent dplyr from masking data.table functions
first <- data.table::first
last  <- data.table::last

# ============================================================
# 0. CONFIGURATION
# ============================================================

n_cores <- 12      # increase to 32 for full HPC runs
testrun <- FALSE   # set FALSE for full production run
handlers("txtprogressbar")

# Input paths
deg_sets_dir       <- "/cjc/data/Gage_align_bamFiles/count_DEG_analysis/WithConfounds/HOMER/overlap_DEGs/DEG_sets"
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
kmer_gene_map <- kmer_gene_map[, .(kmer, kmer_chrom, kmer_start, kmer_end,
                                    ensembl_id, gene_symbol)]
kmer_features <- fread(kmer_features_file)

cat("  kmer_gene_map:", nrow(kmer_gene_map), "rows\n")
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
  gene_symbol       = gene_symbol[gene_symbol != "unassigned"][1]
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
  within_idx <- which(parts == "within")
  fixed_value <- if (length(within_idx) > 0 && within_idx + 1 <= length(parts)) {
    remaining <- paste(parts[(within_idx + 1):length(parts)], collapse="_")
    # Strip trailing date (e.g. _11-11-2025)
    sub("_?\\d{2}-\\d{2}-\\d{4}$", "", remaining)
  } else {
    NA_character_
  }
  
  # Identify what is fixed by matching fixed_value against ALL known dimension lists
  # Do NOT assume — a counters file could fix an aligner OR a trimmer
  matched_trimmer <- trimmers[sapply(trimmers, function(t)
    grepl(paste0("(^|_)", t, "(_|$)"), fixed_value))]
  matched_aligner <- aligners[sapply(aligners, function(a)
    grepl(paste0("(^|_)", a, "(_|$)"), fixed_value))]
  matched_counter <- counters[sapply(counters, function(c)
    grepl(paste0("(^|_)", c, "(_|$)"), fixed_value))]
  
  # Take longest match first to handle overlapping names
  # (e.g. homer_gene_condensed_unique before homer_gene)
  matched_trimmer <- sort(matched_trimmer, decreasing=TRUE)
  matched_aligner <- sort(matched_aligner, decreasing=TRUE)
  matched_counter <- sort(matched_counter, decreasing=TRUE)
  
  deg_list_trimmer <- if (length(matched_trimmer) > 0) matched_trimmer[1] else NA_character_
  deg_list_aligner <- if (length(matched_aligner) > 0) matched_aligner[1] else NA_character_
  deg_list_counter <- if (length(matched_counter) > 0) matched_counter[1] else NA_character_
  
  # Warn if fixed_value matched nothing or matched multiple dimensions unexpectedly
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
  
  # fastq_trimmer: the trimmer applied to FASTQs before alignment
  # Only extractable when a trimmer is fixed in the filename
  fastq_trimmer <- deg_list_trimmer %||% NA_character_
  
  # Counter from filename (varying in counter files, or present in aligner/trimmer files)
  counter_match <- sort(counters[sapply(counters, grepl, x=base, fixed=TRUE)],
                        decreasing=TRUE)
  counter <- if (length(counter_match) > 0) counter_match[1] else "all"
  
  # n_varying: denominator for stability indices
  n_varying <- switch(grouping_type,
                      "aligner" = length(aligners),
                      "trimmer" = length(trimmers),
                      "counter" = length(counters),
                      1L
  )
  
  list(
    grouping_type    = grouping_type,    # dimension varying across col 1
    fixed_value      = fixed_value,      # raw within_ string
    deg_list_trimmer = deg_list_trimmer, # fixed trimmer if present, else NA
    deg_list_aligner = deg_list_aligner, # fixed aligner if present, else NA
    deg_list_counter = deg_list_counter, # fixed counter if present, else NA
    fastq_trimmer    = fastq_trimmer,    # alias for deg_list_trimmer for FASTQ context
    counter          = counter,          # counter matched from full filename
    n_varying        = n_varying         # denominator for stability indices
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

match_kmer_file <- function(deg_file, fastq_trimmer) {
  new_name <- sub("\\.csv$", "_AD_vs_CTRL.csv", basename(deg_file))
  new_name <- gsub("-", "_", new_name)
  list(
    deseq2   = file.path(kmer_deseq_dir, paste0("kmc_deseq2_",     fastq_trimmer, "_", new_name)),
    postfilter = file.path(kmer_deseq_dir, paste0("kmc_postfilter_", fastq_trimmer, "_", 
                                                  sub("_AD_vs_CTRL\\.csv$", ".csv", new_name))),
    filtered   = file.path(kmer_deseq_dir, paste0("kmc_filtered_",   fastq_trimmer, "_", 
                                                  sub("_AD_vs_CTRL\\.csv$", ".csv.zst", new_name)))
  )
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
                  ceiling(seq_along(unique(kmers)) / 50000))  # 50k kmers per chunk
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

deg_files <- list.files(deg_sets_dir, pattern="^all_DEGs_.*\\.csv$", full.names=TRUE)
cat("Found", length(deg_files), "DEG files\n")

if (testrun) {
  deg_files <- deg_files[seq(1, length(deg_files), by=5)][1:3]
  cat("TEST RUN with", length(deg_files), "DEG files\n")
}

# Accumulator lists
full_dataset           <- list()
gene_tracker           <- list()
all_summaries          <- list()
all_deg_dt             <- list()   # raw DEG gene lists — for cross-DEG analyses
kmers_by_method        <- list()
all_kmer_gene_lists    <- list()
all_kmer_aligner_lists <- list()
fastq_sensitivity_all  <- list()
motif_results_all      <- list()
subkmer_results_all    <- list()

for (deg_file in deg_files) {

  meta     <- parse_metadata(deg_file)
  deg_name <- gsub("\\.csv$", "", basename(deg_file))

  cat("\n============================================================\n")
  cat("Processing:", deg_name, "\n")
  cat("  fastq_trimmer   :", meta$fastq_trimmer,    "\n")
  cat("  deg_list_trimmer:", meta$deg_list_trimmer,  "\n")
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
  # Store raw DEG gene list BEFORE kmer merging
  # This is the correct input for cross-DEG analyses:
  # which genes were called as DEGs per aligner per DEG set,
  # independent of kmer detection.
  # ----------------------------------------------------------
  all_deg_dt[[deg_name]] <- deg_dt[, .(
    aligner          = aligner,
    ensembl_id       = ensembl_id,
    deg_name         = deg_name,
    deg_list_trimmer = meta$deg_list_trimmer,
    counter          = meta$counter,
    grouping_type    = meta$grouping_type
  )]

  # ----------------------------------------------------------
  # Load kmer DESeq2 files for each FASTQ trimmer
  # ----------------------------------------------------------
  kmer_dt_list <- list()
  for (fastq_trimmer in trimmers) {
    kmer_file       <- match_kmer_file(deg_file, fastq_trimmer)
    kmer_file_match <- Sys.glob(sub("\\.csv$", "*.csv", kmer_file$deseq2))
    if (length(kmer_file_match) == 0) {
      cat("  [SKIP] kmer file not found:", kmer_file, "\n")
      next
    }
    dt <- fread(kmer_file_match[1])
    dt[, fastq_trimmer := fastq_trimmer]
    kmer_dt_list[[fastq_trimmer]] <- dt
  }

  if (length(kmer_dt_list) == 0) {
    cat("  [SKIP] No kmer files found for", deg_name, "\n")
    next
  }

  kmer_dt <- rbindlist(kmer_dt_list, fill=TRUE)

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
  # Merge sequence features + gene annotations into kmer_dt
  # ----------------------------------------------------------
  kmer_dt <- merge(kmer_dt, all_kmers_dt, by="kmer", all.x=TRUE)

  # Gene assignments — many-to-many intentional
  kmer_dt <- merge(kmer_dt,
                   kmer_gene_map[, .(kmer, ensembl_id, gene_symbol)],
                   by="kmer", all.x=TRUE, allow.cartesian=TRUE)

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
  
  # Step 1: tag each kmer with whether its gene appears in any DEG list
  # This is a simple flag — no expansion needed
  kmer_dt[, is_deg_gene  := ensembl_id %in% deg_dt$ensembl_id]
  kmer_dt[, n_aligners_called := {
    # For each kmer's gene, how many aligners called it as a DEG
    deg_dt[match(ensembl_id, deg_dt$ensembl_id), .N]
  }]
  
  cat("  kmers in DEG genes:", sum(kmer_dt$is_deg_gene, na.rm=TRUE), "\n")
  
  # Step 2: for the subset of kmers in DEG genes, join aligner information
  # This limits the cartesian expansion to only annotated kmers
  kmer_dt_deg <- kmer_dt[is_deg_gene == TRUE]
  
  cat("  kmer_dt_deg rows before aligner join:", nrow(kmer_dt_deg), "\n")
  
  # Build merged separately for gene/aligner-level analyses only
  # Cartesian is now bounded: only DEG-annotated kmers × aligners
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
  
  # Verify no NAs in method before split
  cat("  NA methods:", sum(is.na(merged$method)), "\n")
  cat("  Unique methods:", uniqueN(merged$method), "\n")
  cat("  nrow merged:", nrow(merged), "\n")
  cat("  length unique kmer:", length(unique(merged$kmer)), "\n")
  
  # Force-fill any remaining NA methods before split
  merged[is.na(method), method := paste(meta$grouping_type,
                                        meta$fixed_value %||% "unknown",
                                        meta$counter     %||% "all",
                                        sep="_")]
  
  full_dataset[[deg_name]]           <- merged
  gene_tracker[[deg_name]]           <- unique(merged$ensembl_id)
  kmers_by_method[[deg_name]] <- lapply(
    split(merged$kmer, merged$method),
    unique
  )
  all_kmer_gene_lists[[deg_name]]    <- split(merged$ensembl_id,   merged$method)
  all_kmer_aligner_lists[[deg_name]] <- split(merged$ensembl_id,   merged$aligner)

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
    fastq_sens[, genes_norm := n_genes / max(n_genes, na.rm=TRUE)]
    fastq_sens[, kmers_norm := n_kmers / max(n_kmers, na.rm=TRUE)]
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
    trimmer_sens[, genes_norm := n_genes / max(n_genes, na.rm=TRUE)]
    trimmer_sens[, kmers_norm := n_kmers / max(n_kmers, na.rm=TRUE)]
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
    counter_sens[, genes_norm := n_genes / max(n_genes, na.rm=TRUE)]
    counter_sens[, kmers_norm := n_kmers / max(n_kmers, na.rm=TRUE)]
    fwrite(counter_sens,
           file.path(deg_out, paste0("counter_sensitivity_vs_fixed_dim_",
                                     today, ".csv")))
    fastq_sensitivity_all[[deg_name]] <- counter_sens
  }

  # -------------------------------------------------------
  # 6b. KMER STABILITY ACROSS FASTQ TRIMMERS
  # -------------------------------------------------------
  # Gene-level symbol lookup — one symbol per ensembl_id from kmers_dt
  gene_symbol_lookup <- unique(kmers_dt[!is.na(gene_symbol) & gene_symbol != "unassigned",
                                        .(ensembl_id, gene_symbol)])
  gene_symbol_lookup <- gene_symbol_lookup[, .(gene_symbol = gene_symbol[1]), by=ensembl_id]
  
  varying_col <- meta$grouping_type  # "aligner", "trimmer", or "counter"
  
  gene_stability <- kmers_dt[, .(
    n_varying_detected = uniqueN(get(varying_col))    # was n_trimmers_detected
  ), by=.(ensembl_id, fixed_value, counter)]
  gene_stability[, stability_index := n_varying_detected / meta$n_varying]
  
  gene_stability <- merge(gene_stability, gene_symbol_lookup, by="ensembl_id", all.x=TRUE)
  setcolorder(gene_stability, c("ensembl_id", "gene_symbol",
                                "fixed_value", "counter",
                                "n_varying_detected", "stability_index"))
  fwrite(gene_stability,
         file.path(deg_out, paste0("gene_stability_across_fastq_", today, ".csv")))
  
  kmer_stability <- kmers_dt[, .(
    n_varying_detected = uniqueN(get(varying_col))    # was n_trimmers_detected
  ), by=.(kmer, ensembl_id, fixed_value, counter)]
  kmer_stability[, stability_index := n_varying_detected / meta$n_varying]
  kmer_stability <- merge(kmer_stability, gene_symbol_lookup, by="ensembl_id", all.x=TRUE)
  setcolorder(kmer_stability, c("kmer", "ensembl_id", "gene_symbol",
                                "fixed_value", "counter",
                                "n_varying_detected", "stability_index"))
  fwrite(kmer_stability,
         file.path(deg_out, paste0("kmer_stability_across_fastq_", today, ".csv")))
  
  # -------------------------------------------------------
  # Combined gene x kmer stability
  # Genes that are stable at the gene level AND have stable
  # kmers within them — the intersection of both signals.
  # Aggregates kmer stability per gene then joins to gene
  # stability so each row is one gene with both metrics.
  # -------------------------------------------------------
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
  
  # Fill NAs for genes with no mapped kmers
  gene_kmer_stability_merged[is.na(n_total_kmers),       n_total_kmers      := 0L]
  gene_kmer_stability_merged[is.na(n_stable_kmers),      n_stable_kmers     := 0L]
  gene_kmer_stability_merged[is.na(n_partial_kmers),     n_partial_kmers    := 0L]
  gene_kmer_stability_merged[is.na(mean_kmer_stability), mean_kmer_stability := 0]
  gene_kmer_stability_merged[is.na(max_kmer_stability),  max_kmer_stability  := 0]
  gene_kmer_stability_merged[is.na(stable_kmers_list),   stable_kmers_list  := ""]
  
  # Compute derived columns BEFORE setcolorder
  gene_kmer_stability_merged[, joint_stable         := stability_index >= 1.0 & n_stable_kmers > 0]
  gene_kmer_stability_merged[, joint_stability_score := stability_index * mean_kmer_stability]
  
  # Now all columns exist — safe to reorder
  desired_cols <- c("ensembl_id", "gene_symbol",
                    "fixed_value", "counter",
                    "stability_index", "n_varying_detected",
                    "mean_kmer_stability", "max_kmer_stability",
                    "n_stable_kmers", "n_partial_kmers", "n_total_kmers",
                    "joint_stable", "joint_stability_score",
                    "most_stable_kmer", "stable_kmers_list")
  # Only reorder columns that actually exist to prevent future colnamesInt errors
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
  # 6c. MOTIF ENRICHMENT (kmer-level + sub-kmer + universalmotif)
  # -------------------------------------------------------
  n_sig <- sum(kmer_dt$padj < 0.05, na.rm=TRUE)
  cat(sprintf("\n  Significant kmers (padj<0.05): %d\n", n_sig))
  
  if (n_sig > 0) {
    foreground_up   <- kmer_dt[padj < 0.05 & log2FoldChange > 0]
    foreground_down <- kmer_dt[padj < 0.05 & log2FoldChange < 0]
    
    # === BACKGROUND: load DESeq2 script outputs as conditioned background ===
    # Prefer postfilter (DESeq2-normalized, same detectability threshold as foreground)
    # Fall back to kmer_dt itself if outputs not found
    
    if (file.exists(kmer_file$postfilter) && file.exists(kmer_file$deseq2)) {
      cat("  Loading DESeq2 script outputs as background controls...\n")
      deseq2_res_bg    <- fread(kmer_file$deseq2)
      background_kmers <- deseq2_res_bg[is.na(padj) | padj >= 0.05, kmer]
      postfilter_wide  <- fread(kmer_file$postfilter)
      kmer_bg          <- postfilter_wide[kmer %in% background_kmers, .(kmer)]
      cat(sprintf("  Background: %d non-sig k-mers (postfilter, padj>=0.05 or NA)\n", nrow(kmer_bg)))
    } else {
      cat("  [WARN] DESeq2 background outputs not found — falling back to full kmer_dt\n")
      kmer_bg <- kmer_dt[, .(kmer)]
    }
    
    cat(sprintf("  Motif enrichment: %s (up=%d, down=%d kmers, bg=%d kmers)\n",
                deg_name, nrow(foreground_up), nrow(foreground_down), nrow(kmer_bg)))
    
    # Kmer-level enrichment
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
    
    # Sub-kmer enrichment across all k sizes in parallel
    cat(sprintf("  Sub-kmer enrichment (k=%s) in parallel...\n",
                paste(subkmer_sizes, collapse=",")))
    
    with_progress({
      p <- progressor(steps = length(subkmer_sizes) * 2)  # up + down per k size
      
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
    
    # Sequence logos for significant sub-kmers
    top_motifs <- subkmer_results[padj < 0.05][
      order(pval), .SD[1:min(.N, 20)], by=.(k, direction)]
    for (k_size in subkmer_sizes) {
      plot_seqlogo_robust(top_motifs, k_size, deg_out, today, deg_name)
    }
    
    # universalmotif PWM-based enrichment
    cat("  Running universalmotif enrichment...\n")
    run_universalmotif_enrichment(foreground_up,   kmer_bg, deg_out, deg_name, "up",   today)
    run_universalmotif_enrichment(foreground_down, kmer_bg, deg_out, deg_name, "down", today)
    
  } else {
    cat("  [SKIP] No significant kmers for motif enrichment\n")
  }

  # -------------------------------------------------------
  # 6d. WITHIN-DEG-SET REPRODUCIBILITY (across methods)
  # -------------------------------------------------------
  reproducibility <- kmers_dt[, .(N = .N), by=.(ensembl_id)]
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
    p_gc <- ggplot(kmers_dt, aes(gc_content, fill=fastq_trimmer)) +
      geom_density(alpha=0.4) +
      facet_wrap(~ deg_list_trimmer) +
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

  cat("  [DONE]", deg_name, "\n")
}

# ============================================================
# POST-LOOP: COMBINE ALL RESULTS
# ============================================================

cat("\nCombining all DEG sets...\n")
analysis_dt <- rbindlist(full_dataset, fill=TRUE)

# ============================================================
# CROSS-DEG ANALYSES
# Uses all_deg_dt — raw DEG gene lists BEFORE kmer merging.
# Answers: which genes are reproducibly called as DEGs across
# independent DEG set comparisons, and does this vary by
# aligner, counter, or trimmer context?
# ============================================================

all_deg_dt_combined <- rbindlist(all_deg_dt, fill=TRUE)

fwrite(all_deg_dt_combined,
       file.path(output_dir, "tables",
                 paste0("all_deg_tables_combined_", today, ".csv")))

# Overall: how many DEG sets does each gene appear in
gene_cross_deg <- all_deg_dt_combined[, .(
  n_deg_sets       = uniqueN(deg_name),
  n_aligners       = uniqueN(aligner),
  deg_sets_list    = paste(unique(deg_name), collapse="|"),
  aligners_list    = paste(unique(aligner),  collapse="|")
), by=ensembl_id]
gene_cross_deg[, cross_deg_fraction := n_deg_sets / length(deg_files)]

fwrite(gene_cross_deg,
       file.path(output_dir, "tables",
                 paste0("gene_cross_deg_reproducibility_", today, ".csv")))

# Per aligner: controls for aligner-specific bias
# Genes robust within an aligner across DEG sets are aligner-stable
gene_cross_deg_aligner <- all_deg_dt_combined[, .(
  n_deg_sets         = uniqueN(deg_name),
  cross_deg_fraction = uniqueN(deg_name) / length(deg_files)
), by=.(ensembl_id, aligner)]

fwrite(gene_cross_deg_aligner,
       file.path(output_dir, "tables",
                 paste0("gene_cross_deg_per_aligner_", today, ".csv")))

# Per counter: isolates counting-method effect on DEG reproducibility
gene_cross_deg_counter <- all_deg_dt_combined[, .(
  n_deg_sets         = uniqueN(deg_name),
  cross_deg_fraction = uniqueN(deg_name) / length(deg_files)
), by=.(ensembl_id, counter)]

fwrite(gene_cross_deg_counter,
       file.path(output_dir, "tables",
                 paste0("gene_cross_deg_per_counter_", today, ".csv")))

# Per deg_list_trimmer: isolates trimmer-context effect
gene_cross_deg_trimmer <- all_deg_dt_combined[, .(
  n_deg_sets         = uniqueN(deg_name),
  cross_deg_fraction = uniqueN(deg_name) / length(deg_files)
), by=.(ensembl_id, deg_list_trimmer)]

fwrite(gene_cross_deg_trimmer,
       file.path(output_dir, "tables",
                 paste0("gene_cross_deg_per_trimmer_", today, ".csv")))

# Cross-DEG summary plots
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

p_cross_deg_counter <- ggplot(gene_cross_deg_counter,
                               aes(x=cross_deg_fraction, fill=counter)) +
  geom_density(alpha=0.4) +
  theme_classic() +
  labs(title="Cross-DEG gene reproducibility by counter",
       x="Fraction of DEG sets gene appeared in", y="Density")
ggsave(file.path(output_dir, "figures",
                 paste0("cross_deg_reproducibility_by_counter_", today, ".pdf")),
       p_cross_deg_counter, width=10, height=6)

cat("\nCross-DEG summary:\n")
cat("  Total unique genes across all DEG sets:", uniqueN(all_deg_dt_combined$ensembl_id), "\n")
cat("  Genes in all DEG sets (fraction=1)    :", sum(gene_cross_deg$cross_deg_fraction == 1), "\n")
cat("  Genes in >50% of DEG sets             :", sum(gene_cross_deg$cross_deg_fraction > 0.5), "\n")

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

# ============================================================
# 8. GLOBAL FASTQ-TRIMMER SENSITIVITY SUMMARY
# ============================================================

fastq_sensitivity_dt <- rbindlist(fastq_sensitivity_all, fill=TRUE)

# Summarise separately per grouping type since the varying dimension differs
for (gtype in c("aligner", "trimmer", "counter")) {
  sub <- fastq_sensitivity_dt[grouping_type == gtype]
  if (nrow(sub) == 0) next
  
  vary_col  <- gtype
  fixed_col <- switch(gtype,
                      "aligner" = "deg_list_trimmer",
                      "trimmer" = "deg_list_aligner",
                      "counter" = "fixed_value")
  
  if (!fixed_col %in% names(sub)) next
  
  global_sens <- sub[, .(
    mean_genes_norm = mean(genes_norm, na.rm=TRUE),
    mean_kmers_norm = mean(kmers_norm, na.rm=TRUE),
    sd_genes_norm   = sd(genes_norm,   na.rm=TRUE)
  ), by=c(fixed_col, vary_col)]
  
  fwrite(global_sens,
         file.path(output_dir, "tables",
                   paste0("global_", gtype, "_sensitivity_summary_", today, ".csv")))
}

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
ggsave(file.path(output_dir, "figures",
                 paste0("global_fastq_sensitivity_", today, ".pdf")),
       p_global_sens, width=12, height=6)

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
), by=.(method, counter, deg_list_trimmer)]

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

all_genes    <- unique(analysis_dt$ensembl_id)
all_aligners <- unique(analysis_dt$aligner)
all_trimrs   <- unique(analysis_dt$fastq_trimmer)
all_counters <- unique(analysis_dt$counter)

full_dt <- unique(analysis_dt[, .(ensembl_id, aligner, fastq_trimmer, counter,
                                  grouping_type, fixed_value, DEG_binary)])
full_dt[, DEG_binary := 0L]
full_dt[analysis_dt,
        on=.(ensembl_id, aligner, trimmer=fastq_trimmer, counter),
        DEG_binary := 1L]

gene_feature_cols <- intersect(
  c("ensembl_id","gc_content","kmer_entropy","length","has_exon_overlap"),
  names(analysis_dt))
gene_meta <- analysis_dt[, lapply(.SD, function(x) {
  if (is.numeric(x)) mean(x, na.rm=TRUE) else x[1]
}), .SDcols=setdiff(gene_feature_cols, "ensembl_id"), by=ensembl_id]
full_dt <- merge(full_dt, gene_meta, by="ensembl_id", all.x=TRUE)

cat("Full gene x pipeline matrix:", nrow(full_dt), "rows |",
    "DEG_binary=1:", sum(full_dt$DEG_binary),
    sprintf("(%.1f%%)\n", 100*mean(full_dt$DEG_binary)))

# ============================================================
# 12. KMER GLOBAL STABILITY
# ============================================================

kmer_global_stability <- analysis_dt[, .(
  n_deg_sets          = uniqueN(method),
  n_fastq_trimmers    = uniqueN(fastq_trimmer),
  n_counters          = uniqueN(counter),
  n_aligners          = uniqueN(aligner),
  n_deg_list_trimmers = uniqueN(deg_list_trimmer),
  mean_baseMean       = mean(baseMean, na.rm=TRUE)
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
# Integrates: gene robustness (kmer-based), kmer stability,
# pipeline detection rate, cross-DEG reproducibility (DEG-based),
# and mean expression into a composite biomarker score.
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

# Cross-DEG fraction: genes reproducible across independent DEG set
# comparisons score higher regardless of which pipeline produced the call
biomarker_dt <- merge(biomarker_dt,
                      gene_cross_deg[, .(ensembl_id, cross_deg_fraction, n_aligners)],
                      by="ensembl_id", all.x=TRUE)
biomarker_dt[is.na(cross_deg_fraction), cross_deg_fraction := 0]

# Composite biomarker score
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

has_aligner <- "aligner"      %in% names(full_dt) && uniqueN(full_dt$aligner)      > 1
has_trimmer <- "fastq_trimmer" %in% names(full_dt) && uniqueN(full_dt$fastq_trimmer) > 1
has_counter <- "counter"      %in% names(full_dt) && uniqueN(full_dt$counter)      > 1

formula_full <- bf(
  DEG_binary ~
    kmer_entropy + gc_content + length +
    (if (has_aligner) "aligner" else NULL) +
    (if (has_trimmer) "fastq_trimmer" else NULL) +
    (if (has_counter) "counter" else NULL) +
    (1 | ensembl_id),
  family = bernoulli()
)

fit_full <- brm(formula_full, data=full_dt,
                chains=4, cores=n_cores, iter=4000, backend="cmdstanr")
saveRDS(fit_full, file.path(output_dir, "models/full_model.rds"))

# ============================================================
# 15. MODEL COMPARISON (LOO / WAIC)
# ============================================================

loo_full  <- loo(fit_full)
waic_full <- waic(fit_full)
saveRDS(loo_full,  file.path(output_dir, "models/loo_full.rds"))
saveRDS(waic_full, file.path(output_dir, "models/waic_full.rds"))

# ============================================================
# 16. K-FOLD CROSS-VALIDATION
# ============================================================

kfold_res <- kfold(fit_full, K=10)
saveRDS(kfold_res, file.path(output_dir, "models/kfold.rds"))

# ============================================================
# 17. LEAVE-ONE-ALIGNER-OUT VALIDATION
# ============================================================

unique_aligners <- unique(full_dt$aligner)
loo_align <- future_lapply(unique_aligners, function(a) {
  brm(formula_full, data=full_dt[aligner != a],
      chains=2, cores=2, iter=2000, backend="cmdstanr")
})
names(loo_align) <- unique_aligners
saveRDS(loo_align, file.path(output_dir, "models/leave_one_aligner_out.rds"))

# ============================================================
# 18. ROC / AUC PER ALIGNER
# ============================================================

fitted_probs <- fitted(fit_full)[, 1]
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
# 19. VARIANCE DECOMPOSITION
# ============================================================

vc <- VarCorr(fit_full)
saveRDS(vc, file.path(output_dir, "tables/variance_decomposition.rds"))

# ============================================================
# 20. SHAP-STYLE VARIANCE ATTRIBUTION
# ============================================================

posterior_draws <- as_draws_df(fit_full)
var_components  <- colVars(as.matrix(posterior_draws))
fwrite(data.table(parameter=names(var_components), variance=var_components),
       file.path(output_dir, "tables/variance_attribution.tsv"))

# ============================================================
# 21. HIERARCHICAL CLUSTERING
# ============================================================

clust_mat <- dcast(full_dt, ensembl_id ~ aligner + trimmer + counter,
                   value.var="DEG_binary", fill=0)
dist_mat  <- dist(as.matrix(clust_mat[, -1]))
hc        <- hclust(dist_mat)

pdf(file.path(output_dir, "figures/hierarchical_dendrogram.pdf"))
plot(hc, labels=FALSE, main="Gene clustering by pipeline detection pattern")
dev.off()

# ============================================================
# 22. POSTERIOR PREDICTIVE CHECKS
# ============================================================

pdf(file.path(output_dir, "diagnostics/posterior_predictive.pdf"))
pp_check(fit_full)
dev.off()

# ============================================================
# 23. INTERACTION VISUALIZATION
# ============================================================

pdf(file.path(output_dir, "figures/interaction_plots.pdf"))
plot(conditional_effects(fit_full, effects="aligner:trimmer"))
dev.off()

# ============================================================
# 24. PDF ROBUSTNESS REPORT WITH METHOD SIMILARITY HEATMAPS
# ============================================================

pdf(file.path(output_dir,
              paste0("Robustness_Report_with_Heatmaps_", today, ".pdf")),
    width=14, height=10)

heatmap_mat <- reshape2::acast(all_summaries_dt,
                               method ~ deg_set, value.var="unique_genes")
if (!is.null(heatmap_mat) && nrow(heatmap_mat) > 0) {
  draw(Heatmap(heatmap_mat,
               name         = "Unique Genes",
               column_title = "Unique genes per method x DEG set",
               col          = colorRamp2(
                 c(min(heatmap_mat, na.rm=TRUE), max(heatmap_mat, na.rm=TRUE)),
                 c("white","steelblue"))))
}

barplot(consensus$total_genes,
        names.arg=consensus$method, las=2,
        main="Total genes per method", col="steelblue")

for (dn in names(all_kmer_gene_lists)) {
  jac_mat <- compute_jaccard(all_kmer_gene_lists[[dn]])
  if (nrow(jac_mat) > 0) {
    draw(Heatmap(jac_mat,
                 name         = "Jaccard",
                 column_title = paste0("Method similarity: ", dn),
                 col          = colorRamp2(c(0,1), c("white","darkgreen"))))
  }
}

dev.off()

# ============================================================
# 25. ZSTD COMPRESSION
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

cat("\n============================================================\n")
cat("Pipeline completed successfully.\n")
cat("Key outputs:\n")
cat("  Finalized biomarkers            :",
    file.path(output_dir, "tables", paste0("FINALIZED_BIOMARKERS_",               today, ".tsv")), "\n")
cat("  All DEG tables combined         :",
    file.path(output_dir, "tables", paste0("all_deg_tables_combined_",            today, ".csv")), "\n")
cat("  Gene cross-DEG reproducibility  :",
    file.path(output_dir, "tables", paste0("gene_cross_deg_reproducibility_",     today, ".csv")), "\n")
cat("  Cross-DEG per aligner           :",
    file.path(output_dir, "tables", paste0("gene_cross_deg_per_aligner_",         today, ".csv")), "\n")
cat("  Cross-DEG per counter           :",
    file.path(output_dir, "tables", paste0("gene_cross_deg_per_counter_",         today, ".csv")), "\n")
cat("  Cross-DEG per trimmer           :",
    file.path(output_dir, "tables", paste0("gene_cross_deg_per_trimmer_",         today, ".csv")), "\n")
cat("  Kmer global stability           :",
    file.path(output_dir, "tables", paste0("kmer_global_stability_",              today, ".tsv")), "\n")
cat("  Gene kmer stability             :",
    file.path(output_dir, "tables", paste0("gene_kmer_stability_",                today, ".tsv")), "\n")
cat("  Motif enrichment (all DEG sets) :",
    file.path(output_dir, "tables", paste0("motif_enrichment_all_degsets_",       today, ".csv")), "\n")
cat("  Sub-kmer enrichment (all DEG)   :",
    file.path(output_dir, "tables", paste0("subkmer_enrichment_all_degsets_",     today, ".csv")), "\n")
cat("  Global fastq sensitivity        :",
    file.path(output_dir, "tables", paste0("global_fastq_sensitivity_summary_",   today, ".csv")), "\n")
cat("  Robustness ranking              :",
    file.path(output_dir, "diagnostics", paste0("ROBUSTNESS_RANKING_",            today, ".csv")), "\n")
