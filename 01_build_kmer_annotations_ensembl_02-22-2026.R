#!/usr/bin/env Rscript
# 01_build_kmer_annotations_ensembl.R
#
# Produces two files using a single Ensembl GTF throughout:
#   kmer_rMATS_annotations/kmer_gene_map_from_gtf_<date>.csv
#   kmer_features/kmer_splice_motif_position_features_<date>.csv
#
# Both will have true ENSG* ensembl_id values and gene_symbol (gene_name)
# from the Ensembl annotation, not NCBI RefSeq "gene-GENENAME" style IDs.
#
# Run with:
#   nohup Rscript 01_build_kmer_annotations_ensembl_02-22-2026.R > annotation_$(date +%Y%m%d).log 2>&1 &
# process 1568484
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(stringr)
  library(stringi)
})

# ============================================================
# 0. CONFIGURATION — edit these paths
# ============================================================

gtf_file   <- "/home/thilsabeck/Documents/genomes/Reference/Homo_sapiens.GRCh38.112.gtf"
coords_dir <- "kmer_rMATS_annotations"   # contains kmer_coords_main_chromosomes_*.csv
map_outdir <- "kmer_rMATS_annotations"
feat_outdir <- "kmer_features"
today      <- format(Sys.Date(), "%Y%m%d")

dir.create(map_outdir,  recursive=TRUE, showWarnings=FALSE)
dir.create(feat_outdir, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. LOAD KMER COORDINATES
# ============================================================

cat("=== STEP 1: Loading kmer coordinates ===\n")

coords_files <- list.files(coords_dir,
                           pattern="^kmer_coords_main_chromosomes_.*\\.csv$",
                           full.names=TRUE)
if (length(coords_files) == 0) stop("No kmer_coords_main_chromosomes_*.csv found in ", coords_dir)

kmer_coords <- fread(coords_files[1])
setnames(kmer_coords,
         old = intersect(names(kmer_coords), c("kmer_start","kmer_end")),
         new = c("start","end")[seq_len(sum(c("kmer_start","kmer_end") %in% names(kmer_coords)))])

cat(sprintf("  Loaded %d kmer coordinates\n", nrow(kmer_coords)))
cat(sprintf("  Example chroms: %s\n", paste(head(unique(kmer_coords$chrom)), collapse=", ")))

# ============================================================
# 2. LOAD AND PARSE ENSEMBL GTF
#    Uses Homo_sapiens.GRCh38.112.gtf throughout — no NCBI GTF
# ============================================================

cat("\n=== STEP 2: Loading Ensembl GTF ===\n")

gtf <- fread(
  gtf_file,
  sep       = "\t",
  header    = FALSE,
  skip      = 5,          # skip GTF header comment lines
  col.names = c("chrom","source","feature","start","end",
                "score","strand","frame","attributes"),
  colClasses = "character",
  fill       = TRUE
)

# Keep exons only
gtf <- gtf[feature == "exon"]
gtf[, `:=`(start = as.integer(start), end = as.integer(end))]

# Parse ensembl_id (gene_id) and gene_symbol (gene_name) from attributes
gtf[, ensembl_id  := str_match(attributes, 'gene_id "([^"]+)"')[, 2]]
gtf[, gene_symbol := str_match(attributes, 'gene_name "([^"]+)"')[, 2]]

# Strip version suffix: ENSG00000123456.13 → ENSG00000123456
gtf[, ensembl_id := sub("\\..*$", "", ensembl_id)]

# Remove rows with missing ensembl_id
gtf <- gtf[!is.na(ensembl_id) & ensembl_id != ""]

# Convert Ensembl-style numeric chroms (1, 2, X, Y) to UCSC (chr1, chr2, chrX, chrY)
gtf[, chrom_ucsc := ifelse(grepl("^([0-9]+|X|Y|MT)$", chrom),
                           paste0("chr", chrom),
                           chrom)]
# Rename chrMT → chrM to match typical kmer coordinate files
gtf[chrom_ucsc == "chrMT", chrom_ucsc := "chrM"]
gtf <- gtf[!is.na(chrom_ucsc)]

cat(sprintf("  %d exons across %d unique genes\n",
            nrow(gtf), uniqueN(gtf$ensembl_id)))

# Verify chromosome overlap with kmer coordinates
kmer_chroms  <- unique(kmer_coords$chrom)
gtf_chroms   <- unique(gtf$chrom_ucsc)
common_chroms <- intersect(kmer_chroms, gtf_chroms)
cat(sprintf("  Kmer chroms: %d | GTF chroms: %d | Common: %d\n",
            length(kmer_chroms), length(gtf_chroms), length(common_chroms)))
if (length(common_chroms) == 0) stop("No chromosome overlap between kmer coords and GTF — check chrom naming")

# ============================================================
# 3. BUILD GRANGES
# ============================================================

cat("\n=== STEP 3: Building GRanges objects ===\n")

kmer_gr <- GRanges(
  seqnames = kmer_coords$chrom,
  ranges   = IRanges(kmer_coords$start, kmer_coords$end),
  kmer     = kmer_coords$kmer
)

exon_gr <- GRanges(
  seqnames   = gtf$chrom_ucsc,
  ranges     = IRanges(gtf$start, gtf$end),
  strand     = gtf$strand,
  ensembl_id  = gtf$ensembl_id,
  gene_symbol = gtf$gene_symbol
)

# ============================================================
# 4. EXACT EXON OVERLAPS → kmer_gene_map
#    All kmer-to-exon overlaps, one row per kmer-gene pair
#    (many-to-many intentionally preserved)
# ============================================================

cat("\n=== STEP 4: Finding exact exon overlaps (kmer_gene_map) ===\n")

suppressWarnings({
  exact_ov <- findOverlaps(kmer_gr, exon_gr, ignore.strand=TRUE)
})
cat(sprintf("  %d kmer-exon overlap pairs\n", length(exact_ov)))

kmer_gene_map <- data.table(
  kmer        = mcols(kmer_gr)$kmer[queryHits(exact_ov)],
  kmer_chrom  = as.character(seqnames(kmer_gr))[queryHits(exact_ov)],
  kmer_start  = start(kmer_gr)[queryHits(exact_ov)],
  kmer_end    = end(kmer_gr)[queryHits(exact_ov)],
  ensembl_id  = mcols(exon_gr)$ensembl_id[subjectHits(exact_ov)],
  gene_symbol = mcols(exon_gr)$gene_symbol[subjectHits(exact_ov)]
)

# Remove fully duplicate rows but preserve genuine multi-gene mappings
kmer_gene_map <- unique(kmer_gene_map)
setorder(kmer_gene_map, kmer)

map_outfile <- file.path(map_outdir,
                         sprintf("kmer_gene_map_from_gtf_%s.csv", today))
fwrite(kmer_gene_map, map_outfile)
cat(sprintf("  SAVED: %s\n", map_outfile))
cat(sprintf("  %d unique kmers | %d kmer-gene pairs | %.1f%% of all kmers covered\n",
            uniqueN(kmer_gene_map$kmer),
            nrow(kmer_gene_map),
            100 * uniqueN(kmer_gene_map$kmer) / nrow(kmer_coords)))

# ============================================================
# 5. PROXIMAL OVERLAPS (1kb window) → for kmer_features annotation
#    Take the single best (closest) gene per kmer for features table
# ============================================================

cat("\n=== STEP 5: Finding proximal overlaps for kmer_features (1kb window) ===\n")

suppressWarnings({
  prox_ov <- findOverlaps(kmer_gr, exon_gr, ignore.strand=TRUE, maxgap=1000L)
})
cat(sprintf("  %d proximal kmer-exon pairs\n", length(prox_ov)))

kmer_pos_annot <- data.table(
  kmer        = mcols(kmer_gr)$kmer[queryHits(prox_ov)],
  ensembl_id  = mcols(exon_gr)$ensembl_id[subjectHits(prox_ov)],
  gene_symbol = mcols(exon_gr)$gene_symbol[subjectHits(prox_ov)]
)

# Take first match per kmer (closest by GRanges hit order)
kmer_pos_annot <- kmer_pos_annot[, .SD[1], by=kmer]

# ============================================================
# 6. BUILD kmer_features BASE TABLE
# ============================================================

cat("\n=== STEP 6: Building kmer_features ===\n")

kmer_features <- data.table(
  kmer       = kmer_coords$kmer,
  chrom      = kmer_coords$chrom,
  start      = kmer_coords$start,
  end        = kmer_coords$end,
  length     = kmer_coords$end - kmer_coords$start + 1L,
  gc_content = (stri_count_fixed(kmer_coords$kmer, "G") +
                stri_count_fixed(kmer_coords$kmer, "C")) /
               nchar(kmer_coords$kmer)
)

# ============================================================
# 7. ATTACH EXACT OVERLAP FLAG
#    A kmer has an exact exon overlap if it appears in kmer_gene_map
# ============================================================

exact_kmers <- unique(kmer_gene_map[, .(kmer, ensembl_id, gene_symbol)])

# For the features table, collapse to one representative gene per kmer
# (exact overlap, first hit alphabetically by ensembl_id for determinism)
exact_per_kmer <- exact_kmers[order(kmer, ensembl_id)][, .SD[1], by=kmer]
setnames(exact_per_kmer,
         c("ensembl_id","gene_symbol"),
         c("gene_id_exact","gene_symbol_exact"))
exact_per_kmer[, has_exact_overlap := TRUE]

kmer_features <- merge(kmer_features, exact_per_kmer, by="kmer", all.x=TRUE)
kmer_features[is.na(has_exact_overlap), has_exact_overlap := FALSE]

# ============================================================
# 8. ATTACH PROXIMAL GENE ANNOTATION
#    Prefer exact match; fall back to proximal; fall back to "unassigned"
# ============================================================

kmer_features <- merge(kmer_features, kmer_pos_annot, by="kmer", all.x=TRUE)

kmer_features[, has_gene := !is.na(ensembl_id) | !is.na(gene_id_exact)]

# Final ensembl_id: exact first, then proximal, then unassigned
kmer_features[, ensembl_id := fifelse(
  !is.na(gene_id_exact),  gene_id_exact,
  fifelse(!is.na(ensembl_id), ensembl_id, "unassigned")
)]

# Final gene_symbol: exact first, then proximal, then unassigned
kmer_features[, gene_symbol := fifelse(
  !is.na(gene_symbol_exact), gene_symbol_exact,
  fifelse(!is.na(gene_symbol),  gene_symbol,  "unassigned")
)]

# ============================================================
# 9. POSITION FEATURES (dist to TSS / TES)
# ============================================================

cat("\n=== STEP 9: Computing position features ===\n")

kmer_features[has_gene == TRUE, `:=`(
  gene_min_start = min(start),
  gene_max_end   = max(end)
), by=gene_symbol]

kmer_features[has_gene == TRUE, `:=`(
  dist_to_tss = pmin(abs(start - gene_min_start), abs(end - gene_min_start)),
  dist_to_tes = pmin(abs(start - gene_max_end),   abs(end - gene_max_end))
)]

kmer_features[has_gene == FALSE, `:=`(
  dist_to_tss = 1000000000L,
  dist_to_tes = 1000000000L
)]

kmer_features[, `:=`(gene_min_start = NULL, gene_max_end = NULL)]

# rMATS overlap == exact exon overlap
kmer_features[, has_rMATS_overlap := has_exact_overlap]

# ============================================================
# 10. SANITY CHECKS
# ============================================================

cat("\n=== STEP 10: Sanity checks ===\n")

cat(sprintf("  Total kmers in features table   : %d\n", nrow(kmer_features)))
cat(sprintf("  Kmers with exact exon overlap   : %d (%.1f%%)\n",
            sum(kmer_features$has_exact_overlap),
            100 * mean(kmer_features$has_exact_overlap)))
cat(sprintf("  Kmers with any gene annotation  : %d (%.1f%%)\n",
            sum(kmer_features$has_gene),
            100 * mean(kmer_features$has_gene)))
cat(sprintf("  Kmers with true ENSG* ID        : %d (%.1f%%)\n",
            sum(grepl("^ENSG", kmer_features$ensembl_id)),
            100 * mean(grepl("^ENSG", kmer_features$ensembl_id))))
cat(sprintf("  Kmers still 'unassigned'        : %d (%.1f%%)\n",
            sum(kmer_features$ensembl_id == "unassigned"),
            100 * mean(kmer_features$ensembl_id == "unassigned")))

# Confirm no "gene-" style IDs remain
n_gene_prefix <- sum(grepl("^gene-", kmer_features$ensembl_id))
if (n_gene_prefix > 0) {
  cat(sprintf("  [WARNING] %d rows still have 'gene-' style IDs — check GTF source\n",
              n_gene_prefix))
} else {
  cat("  [OK] No 'gene-' style IDs found — all annotations are Ensembl\n")
}

cat(sprintf("\n  kmer_gene_map: %d unique kmers mapped to %d gene pairs\n",
            uniqueN(kmer_gene_map$kmer), nrow(kmer_gene_map)))
cat(sprintf("  Example ensembl_ids: %s\n",
            paste(head(unique(kmer_features$ensembl_id[
              kmer_features$ensembl_id != "unassigned"]), 5),
              collapse=", ")))

# ============================================================
# 11. SAVE kmer_features
# ============================================================

feat_outfile <- file.path(feat_outdir,
                          sprintf("kmer_splice_motif_position_features_%s.csv", today))
fwrite(kmer_features, feat_outfile)
cat(sprintf("\n  SAVED: %s\n", feat_outfile))

cat("\n=== COMPLETE ===\n")
cat(sprintf("  kmer_gene_map : %s\n", map_outfile))
cat(sprintf("  kmer_features : %s\n", feat_outfile))
