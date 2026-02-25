### Updated DESeq2 analysis - Skip completed, verify 27 samples, auto-detect DEG lists
library(DESeq2)
library(data.table)
library(dplyr)

cat('🚀 EFFICIENT MODE - Skip completed combos, verify 27 samples, 5M k-mer cap\n')

# === AUTO-DETECT ACTUAL STRUCTURE ===
cat('🔍 Scanning kmerator_results for actual directories...\n')
sample_dirs <- list.dirs("kmerator_results", full.names = FALSE, recursive = FALSE)
allowed_samples <- intersect(sample_dirs, c(
  "X53yr_KP_F", "X56yr_840C_F", "X_yr_HR_F", "X65yr_8150_F_S84", "X85yr_14096_F",
  "X81yr_8097_F", "X80yr_3367_F", "X80yr_2800_F", "X83yr_3093_F", "X81yr_3053_F",
  "X81yr_8175_F_S82", "X83yr_3113_F_S78", "X84yr_40_F", "X85yr_3056_F",
  "X86yr_8011_F", "X86yr_27_F", "X87yr_8078_F_S83", "X89yr_2991_F",
  "X67yr_2608_F_S81", "X79yr_3121_F", "X82yr_14095_F_S85", "X81yr_8020_F",
  "X78yr_2785_F_S80", "X83yr_8149LG_F", "X80yr_3131_F_S79", "X81yr_3383_F", "X88yr_654LG_F"
))
cat('✅ Found', length(allowed_samples), '/27 samples:', paste(allowed_samples, collapse = ", "), '\n')

# Auto-detect ACTUAL DEG lists from first sample
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

        cat(sprintf("  Scanning %s\n", scan_path))
        deg_dirs <- list.dirs(scan_path, full.names=FALSE, recursive=FALSE)
        deg_dirs <- deg_dirs[grep("^all_DEGs_.*11-11-2025$", deg_dirs)]

        # Keep only dirs that contain kmc_counts files
        deg_dirs <- deg_dirs[vapply(deg_dirs, function(d) {
          full <- file.path(scan_path, d)
          file.exists(file.path(full, "kmc_counts.txt")) ||
          file.exists(file.path(full, "kmc_counts.txt.zst"))
        }, logical(1))]

        # Keep original hyphenated date — do NOT substitute to underscores
        # unless downstream DESeq2 files confirmed to use underscore dates
        found     <- unique(deg_dirs)
        deg_lists <- unique(c(deg_lists, found))
      }
      if (length(deg_lists) > 0) break  # stop after first sample with data
    }
    if (length(deg_lists) > 0) break    # stop after first root with data
  }

  cat(sprintf("  Auto-detected %d DEG lists with kmc_counts files:\n",
              length(deg_lists)))
  for (d in deg_lists) cat(sprintf("    %s\n", d))

} else {
  deg_lists <- character(0)
  cat("  No allowed samples — deg_lists is empty\n")
}

cat(sprintf("🔍 Auto-detected %d DEG lists: %s\n",
            length(deg_lists), paste(deg_lists, collapse = ", ")))

# Priority DEG lists to process first — stripped of .csv suffix for consistency
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

# Intersect auto-detected with priority list
# Priority list entries that are available get processed first,
# remaining auto-detected entries appended after
priority_available    <- intersect(priority_deg_lists, deg_lists)
non_priority_available <- setdiff(deg_lists, priority_deg_lists)
priority_missing      <- setdiff(priority_deg_lists, deg_lists)

deg_lists <- c(priority_available, non_priority_available)

cat(sprintf("  Priority DEG lists available  : %d\n", length(priority_available)))
cat(sprintf("  Priority DEG lists missing    : %d\n", length(priority_missing)))
cat(sprintf("  Additional DEG lists available: %d\n", length(non_priority_available)))
cat(sprintf("  Total to process              : %d\n", length(deg_lists)))

if (length(priority_missing) > 0) {
  cat("  Missing priority DEG lists (no kmc_counts yet):\n")
  for (d in priority_missing) cat(sprintf("    %s\n", d))
}

trim_dirs <- c("untrimmed", "fastp", "trimmomatic", "trimgalore")
dir.create("kmerator_results/deseq2_results", showWarnings = FALSE, recursive = TRUE)

# === LOAD METADATA ===
metadata_csv <- "Gage_merged_MetaData_updated_10-23-2024_10-6-2025.csv"
meta <- fread(metadata_csv)
meta_subset <- meta[, .(`Sample label`, Diag, Sex, Age)][`Sample label` %in% allowed_samples]
setnames(meta_subset, "Sample label", "sample")
cat('✅ Metadata aligned:', nrow(meta_subset), 'samples\n')

max_kmers <- 5000000
total_combos <- length(trim_dirs) * length(deg_lists)
completed_combos <- 0
skipped_combos <- 0

# === MAIN ANALYSIS LOOP ===
for (trim_dir in trim_dirs) {
  cat(sprintf('\n🏗️  TRIM: %s (%d/%d)\n', trim_dir, which(trim_dirs == trim_dir), length(trim_dirs)))
  
  for (deg_list in deg_lists) {
    # === CHECK IF ALREADY COMPLETED ===
    safe_trim <- gsub("[^A-Za-z0-9_]", "_", trim_dir)
    safe_deg  <- gsub("[^A-Za-z0-9_]", "_", deg_list)
    result_path <- sprintf("kmerator_results/deseq2_results/kmc_deseq2_%s_%s_AD_vs_CTRL.csv", safe_trim, safe_deg)
    
    if (file.exists(result_path)) {
      cat(sprintf('    ✅ SKIP: %s__%s (already completed)\n', trim_dir, deg_list))
      completed_combos <- completed_combos + 1
      next
    }
    
    cat(sprintf('    🔄 [%d/%d] Processing %s__%s\n', 
                completed_combos + skipped_combos + 1, total_combos, trim_dir, deg_list))
    
    # === LOAD KMC FILES ===
    kmc_pattern <- file.path("kmerator_results", "*", trim_dir, deg_list, "kmc_counts.txt*")
    kmc_files_all <- Sys.glob(kmc_pattern)
    kmc_files <- kmc_files_all[basename(dirname(dirname(dirname(kmc_files_all)))) %in% allowed_samples]
    
    target_samples <- length(allowed_samples)
    if (length(kmc_files) < target_samples * 0.8) {  # Require 80% coverage
      cat(sprintf("    ❌ Only %d/%d samples (need ~%d) - skipping\n", 
                  length(kmc_files), target_samples, ceiling(target_samples * 0.8)))
      skipped_combos <- skipped_combos + 1
      next
    }
    
    cat(sprintf("    📂 Found %d/%d samples (%.0f%% coverage) ✅\n", 
                length(kmc_files), target_samples, length(kmc_files)/target_samples*100))
    
    # === LOAD COUNTS ===
    counts_list <- lapply(kmc_files, function(f) {
      sample_name <- basename(dirname(dirname(dirname(f))))
      cat(sprintf("      📂 %s\n", sample_name))
      # Handle both .zst and uncompressed
      cmd <- if (grepl("\\.zst$", f)) paste("zstdcat", shQuote(f)) else f
      dt <- tryCatch(fread(cmd = cmd, col.names = c("kmer", "count")), 
                     error = function(e) NULL)
      if (is.null(dt)) return(NULL)
      dt[, sample := sample_name][]
    })
    counts_list <- counts_list[!vapply(counts_list, is.null, logical(1))]
    
    if (length(counts_list) == 0) {
      cat("    ❌ No readable files - skipping\n")
      skipped_combos <- skipped_combos + 1
      next
    }
    
    counts_long <- rbindlist(counts_list)
    cat(sprintf("    ✅ %d obs | %d k-mers\n", nrow(counts_long), uniqueN(counts_long$kmer)))
    
    # === FILTER TOP KMERS ===
    cat(sprintf("    🔧 Top %d variable k-mers...\n", max_kmers))
    kmer_stats <- counts_long[, .(
      total_count = sum(count), n_samples = uniqueN(sample), cv = sd(count)/mean(count)
    ), by = kmer][order(-cv)][, .(kmer, score = total_count * n_samples * cv)]
    top_kmers <- head(kmer_stats$kmer, max_kmers)
    counts_long <- counts_long[kmer %in% top_kmers]
    
    # === FINAL SAMPLE FILTER ===
    avail_samples <- intersect(intersect(unique(counts_long$sample), meta_subset$sample), allowed_samples)
    if (length(avail_samples) < 10) {  # Minimum 10 samples
      cat(sprintf("    ❌ Only %d samples after metadata filter - skipping\n", length(avail_samples)))
      skipped_combos <- skipped_combos + 1
      next
    }
    
    counts_long <- counts_long[sample %in% avail_samples]
    cat(sprintf("    👥 FINAL: %d/%d samples (%.0f%%) → %s\n", 
                length(avail_samples), target_samples, 
                length(avail_samples)/target_samples*100, 
                paste(head(sort(avail_samples), 3), collapse = ", ")))
    
    # === BUILD MATRIX ===
    counts_wide <- dcast(counts_long, kmer ~ sample, value.var = "count", fill = 0)
    count_matrix <- as.matrix(counts_wide[, -"kmer", with=FALSE])
    rownames(count_matrix) <- counts_wide$kmer
    
    # === ALIGN METADATA ===
    col_data <- meta_subset[sample %in% colnames(count_matrix)]
    col_data <- col_data[match(colnames(count_matrix), col_data$sample)]
    col_data$Diag <- relevel(factor(col_data$Diag), ref = "CTRL")
    col_data$Sex <- factor(col_data$Sex)
    col_data$Age <- as.numeric(scale(as.numeric(col_data$Age)))
    
    # Filter low-count k-mers
    keep <- rowSums(count_matrix) >= 10
    count_matrix <- count_matrix[keep, ]
    col_data <- col_data[match(colnames(count_matrix), col_data$sample)]
    
    cat(sprintf("    📊 Matrix: %d k-mers × %d samples\n", nrow(count_matrix), ncol(count_matrix)))
    
    # === SAVE FILTERED COUNTS ===
    counts_out <- as.data.table(count_matrix, keep.rownames = "kmer")
    fwrite(counts_out, sprintf("kmerator_results/deseq2_results/kmc_filtered_%s_%s.csv.zst", safe_trim, safe_deg))
    
    # === DESeq2 ===
    cat("    🚀 DESeq2...\n")
    # Convert to data.frame with rownames before DESeq2
    col_data_df <- as.data.frame(col_data)
    rownames(col_data_df) <- col_data_df$sample
    col_data_df$Diag <- relevel(factor(col_data_df$Diag), ref="CTRL")
    col_data_df$Sex  <- factor(col_data_df$Sex)
    col_data_df$Age  <- as.numeric(scale(as.numeric(col_data_df$Age)))
    
    dds <- DESeqDataSetFromMatrix(round(count_matrix), col_data_df, ~ Age + Sex + Diag)
    dds <- DESeq(dds, fitType = "local")
    res_main <- results(dds, contrast = c("Diag", "AD", "CTRL"))
    
    sig_count <- sum(res_main$padj < 0.05, na.rm = TRUE)
    cat(sprintf("    🎯 %d sig k-mers (padj<0.05)!\n", sig_count))
    
    # === SAVE RESULTS ===
    res_df <- as.data.frame(res_main); res_df$kmer <- rownames(res_df)
    fwrite(res_df, result_path)
    
    post_counts <- assay(dds)
    post_df <- data.frame(kmer = rownames(post_counts), as.data.frame(post_counts))
    fwrite(post_df, sprintf("kmerator_results/deseq2_results/kmc_postfilter_%s_%s.csv", safe_trim, safe_deg))
    
    # Save sample list
    writeLines(paste(sort(colnames(count_matrix)), collapse = ","), 
               sprintf("kmerator_results/deseq2_results/samples_%s_%s.txt", safe_trim, safe_deg))
    
    cat(sprintf("    💾 COMPLETED %s__%s (%d samples, %d sig)\n\n", trim_dir, deg_list, ncol(count_matrix), sig_count))
    gc()
    completed_combos <- completed_combos + 1
  }
}

cat(sprintf('🎉 PIPELINE COMPLETE!\n'))
cat(sprintf('✅ %d/%d combos completed, %d skipped (already done), %d skipped (insufficient samples)\n', 
            completed_combos, total_combos, skipped_combos, total_combos - completed_combos - skipped_combos))

