### Updated DESeq2 analysis - Parallel HPC version
### Skip completed, verify 27 samples, auto-detect DEG lists
### Parallelises across trim_dir x deg_list combinations

library(DESeq2)
library(data.table)
library(dplyr)
library(future)
library(future.apply)
library(parallel)

# ============================================================
# HPC PARALLEL CONFIG
# ============================================================
N_CORES      <- 12L   # total cores to use on this node
N_WORKERS    <- 8L    # parallel DEG list combinations at once
DESEQ2_CORES <- max(1L, floor(N_CORES / N_WORKERS))  # cores per DESeq2 job

cat(sprintf("HPC config: %d total cores | %d parallel workers | %d cores per DESeq2\n",
            N_CORES, N_WORKERS, DESEQ2_CORES))

plan(multicore, workers = N_WORKERS)
options(mc.cores = DESEQ2_CORES)

cat("🚀 PARALLEL HPC MODE - Skip completed, verify 27 samples, 5M k-mer cap\n")

# ============================================================
# AUTO-DETECT ACTUAL STRUCTURE
# ============================================================
cat("🔍 Scanning kmerator_results for actual directories...\n")
sample_dirs <- list.dirs("kmerator_results", full.names=FALSE, recursive=FALSE)
allowed_samples <- intersect(sample_dirs, c(
  "X53yr_KP_F", "X56yr_840C_F", "X_yr_HR_F", "X65yr_8150_F_S84", "X85yr_14096_F",
  "X81yr_8097_F", "X80yr_3367_F", "X80yr_2800_F", "X83yr_3093_F", "X81yr_3053_F",
  "X81yr_8175_F_S82", "X83yr_3113_F_S78", "X84yr_40_F", "X85yr_3056_F",
  "X86yr_8011_F", "X86yr_27_F", "X87yr_8078_F_S83", "X89yr_2991_F",
  "X67yr_2608_F_S81", "X79yr_3121_F", "X82yr_14095_F_S85", "X81yr_8020_F",
  "X78yr_2785_F_S80", "X83yr_8149LG_F", "X80yr_3131_F_S79", "X81yr_3383_F",
  "X88yr_654LG_F"
))
cat("✅ Found", length(allowed_samples), "/27 samples\n")

# ============================================================
# AUTO-DETECT DEG LISTS
# ============================================================
if (length(allowed_samples) > 0) {
  candidate_roots <- c(
    "/netapp/snl/scratch25/AHAllen_projects/Gage2019/kmerator_results",
    "/netapp/snl/scratch25/AHAllen_Projects/Gage2019/kmerator_results"
  )

  samples_to_check <- unique(c("X53yr_KP_F", allowed_samples))
  deg_lists        <- character(0)

  for (root in candidate_roots) {
    for (samp in samples_to_check) {
      for (trimmer in c("fastp", "trimmomatic", "trimgalore", "untrimmed")) {
        scan_path <- file.path(root, samp, trimmer)
        if (!dir.exists(scan_path)) next

        deg_dirs <- list.dirs(scan_path, full.names=FALSE, recursive=FALSE)
        deg_dirs <- deg_dirs[grep("^all_DEGs_.*11-11-2025$", deg_dirs)]
        deg_dirs <- deg_dirs[vapply(deg_dirs, function(d) {
          full <- file.path(scan_path, d)
          file.exists(file.path(full, "kmc_counts.txt")) ||
          file.exists(file.path(full, "kmc_counts.txt.zst"))
        }, logical(1))]

        deg_lists <- unique(c(deg_lists, deg_dirs))
      }
      if (length(deg_lists) > 0) break
    }
    if (length(deg_lists) > 0) break
  }

  cat(sprintf("  Auto-detected %d DEG lists with kmc_counts files\n", length(deg_lists)))

} else {
  deg_lists <- character(0)
}

# Priority ordering
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

priority_available     <- intersect(priority_deg_lists, deg_lists)
non_priority_available <- setdiff(deg_lists, priority_deg_lists)
priority_missing       <- setdiff(priority_deg_lists, deg_lists)
deg_lists              <- c(priority_available, non_priority_available)

cat(sprintf("  Priority available  : %d\n", length(priority_available)))
cat(sprintf("  Priority missing    : %d\n", length(priority_missing)))
cat(sprintf("  Additional available: %d\n", length(non_priority_available)))
cat(sprintf("  Total to process    : %d\n", length(deg_lists)))

# ============================================================
# SETUP
# ============================================================
trim_dirs <- c("untrimmed", "fastp", "trimmomatic", "trimgalore")
dir.create("kmerator_results/deseq2_results", showWarnings=FALSE, recursive=TRUE)

metadata_csv <- "Gage_merged_MetaData_updated_10-23-2024_10-6-2025.csv"
meta         <- fread(metadata_csv)
meta_subset  <- meta[, .(`Sample label`, Diag, Sex, Age)][
  `Sample label` %in% allowed_samples]
setnames(meta_subset, "Sample label", "sample")
cat("✅ Metadata aligned:", nrow(meta_subset), "samples\n")

max_kmers    <- 5000000L
target_samples <- length(allowed_samples)

# ============================================================
# BUILD COMBO LIST — all trim x deg combinations to process
# ============================================================
all_combos <- data.table(
  CJ(trim_dir=trim_dirs, deg_list=deg_lists)
)

# Pre-filter already completed combos
all_combos[, safe_trim := gsub("[^A-Za-z0-9_]", "_", trim_dir)]
all_combos[, safe_deg  := gsub("[^A-Za-z0-9_]", "_", deg_list)]
all_combos[, result_path := sprintf(
  "kmerator_results/deseq2_results/kmc_deseq2_%s_%s_AD_vs_CTRL.csv",
  safe_trim, safe_deg)]
all_combos[, completed := file.exists(result_path)]

n_already_done <- sum(all_combos$completed)
combos_to_run  <- all_combos[completed == FALSE]

cat(sprintf("\n📊 Total combinations   : %d\n", nrow(all_combos)))
cat(sprintf("✅ Already completed    : %d\n", n_already_done))
cat(sprintf("🔄 Remaining to process : %d\n", nrow(combos_to_run)))

if (nrow(combos_to_run) == 0) {
  cat("🎉 All combinations already completed!\n")
  quit(save="no")
}

# ============================================================
# CORE PROCESSING FUNCTION (one trim x deg combo)
# ============================================================
process_combo <- function(combo_row, meta_subset, allowed_samples,
                          max_kmers, target_samples, deseq2_cores) {

  trim_dir    <- combo_row$trim_dir
  deg_list    <- combo_row$deg_list
  safe_trim   <- combo_row$safe_trim
  safe_deg    <- combo_row$safe_deg
  result_path <- combo_row$result_path

  # Double-check not completed (another worker may have finished it)
  if (file.exists(result_path)) {
    return(list(status="skipped", combo=paste(trim_dir, deg_list, sep="__")))
  }

  log_prefix <- sprintf("[%s__%s]", trim_dir, deg_list)
  cat(sprintf("%s Starting\n", log_prefix))

  # --- LOAD KMC FILES ---
  kmc_pattern  <- file.path("kmerator_results", "*", trim_dir, deg_list, "kmc_counts.txt*")
  kmc_files_all <- Sys.glob(kmc_pattern)
  kmc_files    <- kmc_files_all[
    basename(dirname(dirname(dirname(kmc_files_all)))) %in% allowed_samples]

  if (length(kmc_files) < target_samples * 0.8) {
    cat(sprintf("%s Only %d/%d samples — skipping\n",
                log_prefix, length(kmc_files), target_samples))
    return(list(status="insufficient_samples",
                combo=paste(trim_dir, deg_list, sep="__"),
                n_samples=length(kmc_files)))
  }

  cat(sprintf("%s Loading %d/%d samples\n",
              log_prefix, length(kmc_files), target_samples))

  # --- LOAD COUNTS ---
  counts_list <- lapply(kmc_files, function(f) {
    sample_name <- basename(dirname(dirname(dirname(f))))
    cmd <- if (grepl("\\.zst$", f)) paste("zstdcat", shQuote(f)) else f
    dt  <- tryCatch(
      fread(cmd=cmd, col.names=c("kmer", "count")),
      error = function(e) NULL
    )
    if (is.null(dt)) return(NULL)
    dt[, sample := sample_name][]
  })
  counts_list <- counts_list[!vapply(counts_list, is.null, logical(1))]

  if (length(counts_list) == 0) {
    cat(sprintf("%s No readable files — skipping\n", log_prefix))
    return(list(status="no_readable_files",
                combo=paste(trim_dir, deg_list, sep="__")))
  }

  counts_long <- rbindlist(counts_list)
  cat(sprintf("%s %d obs | %d k-mers\n",
              log_prefix, nrow(counts_long), uniqueN(counts_long$kmer)))

  # --- FILTER TOP KMERS ---
  cat(sprintf("%s Selecting top %d variable k-mers\n", log_prefix, max_kmers))
  kmer_stats  <- counts_long[, .(
    total_count = sum(count),
    n_samples   = uniqueN(sample),
    cv          = sd(count) / mean(count)
  ), by=kmer][order(-cv)][, .(kmer, score=total_count * n_samples * cv)]
  top_kmers   <- head(kmer_stats$kmer, max_kmers)
  counts_long <- counts_long[kmer %in% top_kmers]

  # --- FINAL SAMPLE FILTER ---
  avail_samples <- intersect(
    intersect(unique(counts_long$sample), meta_subset$sample),
    allowed_samples)

  if (length(avail_samples) < 10) {
    cat(sprintf("%s Only %d samples after metadata filter — skipping\n",
                log_prefix, length(avail_samples)))
    return(list(status="insufficient_samples_post_filter",
                combo=paste(trim_dir, deg_list, sep="__"),
                n_samples=length(avail_samples)))
  }

  counts_long <- counts_long[sample %in% avail_samples]
  cat(sprintf("%s FINAL: %d/%d samples\n",
              log_prefix, length(avail_samples), target_samples))

  # --- BUILD MATRIX ---
  counts_wide  <- dcast(counts_long, kmer ~ sample, value.var="count", fill=0)
  count_matrix <- as.matrix(counts_wide[, -"kmer", with=FALSE])
  rownames(count_matrix) <- counts_wide$kmer

  # --- ALIGN METADATA ---
  col_data <- meta_subset[sample %in% colnames(count_matrix)]
  col_data <- col_data[match(colnames(count_matrix), col_data$sample)]

  # Convert to data.frame with rownames for DESeq2
  col_data_df           <- as.data.frame(col_data)
  rownames(col_data_df) <- col_data_df$sample
  col_data_df$Diag      <- relevel(factor(col_data_df$Diag), ref="CTRL")
  col_data_df$Sex       <- factor(col_data_df$Sex)
  col_data_df$Age       <- as.numeric(scale(as.numeric(col_data_df$Age)))

  # --- FILTER LOW-COUNT KMERS ---
  keep         <- rowSums(count_matrix) >= 10
  count_matrix <- count_matrix[keep, ]
  col_data_df  <- col_data_df[match(colnames(count_matrix),
                                     rownames(col_data_df)), ]

  cat(sprintf("%s Matrix: %d k-mers x %d samples\n",
              log_prefix, nrow(count_matrix), ncol(count_matrix)))

  # --- SAVE FILTERED COUNTS ---
  counts_out <- as.data.table(count_matrix, keep.rownames="kmer")
  fwrite(counts_out,
         sprintf("kmerator_results/deseq2_results/kmc_filtered_%s_%s.csv.zst",
                 safe_trim, safe_deg))

  # --- DESeq2 ---
  cat(sprintf("%s Running DESeq2 with %d cores\n", log_prefix, deseq2_cores))

  tryCatch({
    # DESeq2 uses BiocParallel internally — set workers per job
    BiocParallel::register(BiocParallel::MulticoreParam(deseq2_cores))

    dds <- DESeqDataSetFromMatrix(
      countData = round(count_matrix),
      colData   = col_data_df,
      design    = ~ Age + Sex + Diag
    )
    dds      <- DESeq(dds, fitType="local", parallel=TRUE)
    res_main <- results(dds, contrast=c("Diag", "AD", "CTRL"))

    sig_count <- sum(res_main$padj < 0.05, na.rm=TRUE)
    cat(sprintf("%s %d sig k-mers (padj<0.05)\n", log_prefix, sig_count))

    # --- SAVE RESULTS ---
    res_df       <- as.data.frame(res_main)
    res_df$kmer  <- rownames(res_df)
    fwrite(res_df, result_path)

    post_counts <- assay(dds)
    post_df     <- data.table(kmer=rownames(post_counts),
                               as.data.table(as.data.frame(post_counts)))
    fwrite(post_df,
           sprintf("kmerator_results/deseq2_results/kmc_postfilter_%s_%s.csv",
                   safe_trim, safe_deg))

    writeLines(
      paste(sort(colnames(count_matrix)), collapse=","),
      sprintf("kmerator_results/deseq2_results/samples_%s_%s.txt",
              safe_trim, safe_deg)
    )

    cat(sprintf("%s ✅ COMPLETED (%d samples, %d sig k-mers)\n",
                log_prefix, ncol(count_matrix), sig_count))

    return(list(status      = "completed",
                combo       = paste(trim_dir, deg_list, sep="__"),
                n_samples   = ncol(count_matrix),
                n_sig_kmers = sig_count))

  }, error = function(e) {
    cat(sprintf("%s ❌ DESeq2 ERROR: %s\n", log_prefix, e$message))
    return(list(status = "error",
                combo  = paste(trim_dir, deg_list, sep="__"),
                error  = e$message))
  })
}

# ============================================================
# PARALLEL EXECUTION
# ============================================================
cat(sprintf("\n🚀 Launching %d combinations across %d workers...\n",
            nrow(combos_to_run), N_WORKERS))

# Convert to list of rows for future_lapply
combo_list <- split(combos_to_run, seq_len(nrow(combos_to_run)))

results_list <- future_lapply(combo_list, function(combo_row) {
  process_combo(
    combo_row      = combo_row,
    meta_subset    = meta_subset,
    allowed_samples = allowed_samples,
    max_kmers      = max_kmers,
    target_samples = target_samples,
    deseq2_cores   = DESEQ2_CORES
  )
}, future.seed=TRUE)

# ============================================================
# SUMMARY
# ============================================================
results_dt <- rbindlist(lapply(results_list, as.data.table), fill=TRUE)

cat("\n============================================================\n")
cat("🎉 PIPELINE COMPLETE\n")
cat("============================================================\n")
cat(sprintf("✅ Already done      : %d\n", n_already_done))
cat(sprintf("✅ Completed now     : %d\n", sum(results_dt$status == "completed")))
cat(sprintf("⏭️  Skipped (done)   : %d\n", sum(results_dt$status == "skipped")))
cat(sprintf("❌ Insufficient samp : %d\n", sum(grepl("insufficient", results_dt$status))))
cat(sprintf("❌ Errors            : %d\n", sum(results_dt$status == "error")))

if (any(results_dt$status == "error")) {
  cat("\nFailed combinations:\n")
  print(results_dt[status == "error", .(combo, error)])
}

if (any(results_dt$status == "completed")) {
  cat(sprintf("\nMean sig k-mers per completed combo: %.0f\n",
              mean(results_dt[status == "completed"]$n_sig_kmers, na.rm=TRUE)))
}

fwrite(results_dt,
       sprintf("kmerator_results/deseq2_results/run_summary_%s.csv",
               format(Sys.Date(), "%Y-%m-%d")))
