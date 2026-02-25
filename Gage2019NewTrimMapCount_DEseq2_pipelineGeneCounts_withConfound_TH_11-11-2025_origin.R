# TH 10-16-2024 Performing integration approach following tutorial from https://satijalab.org/seurat/articles/integration_introduction.html of Gage 2024 data in 2 batches
# TH 4-3-2025 Plotting significant DEseq2 results from each of the mapping pipelines used with Gage 2019 iN data for F32 grant preliminary data
# TH 11-6-2025 Analyzing all preprocessing pipelines for fibroblasts, including HOMER counts

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("Rsubread")
# BiocManager::install("ensembldb")
# BiocManager::install("sva")

# library(Rsubread)

#Load libraries
# source("~/dimurali/wahl_scRNA/scripts/library.R")
# 
# if (!dir.exists("./DiffExpGenes/Preprocessing/")) {
#   dir.create("./DiffExpGenes/Preprocessing/")#./seurat/D3to5Adult")
# }
# setwd("./DiffExpGenes/Preprocessing/")#./seurat/D3to5Adult")

# # install BiocManager
# install.packages("BiocManager")
# 
# # install Bioconductor core packages
# BiocManager::install()
# 
# # install additional packages:
# BiocManager::install(c("WGCNA", "igraph", "devtools", "GeneOverlap", "ggrepel", "UCell"))
# devtools::install_github("NightingaleHealth/ggforestplot", force = TRUE)
# 
# # install.packages("devtools")
# devtools::install_github("immunogenomics/presto")
# 
# # install Seurat v5 
# install.packages("Seurat", dependencies=TRUE)
# 
# devtools::install_github('smorabit/hdWGCNA', ref='dev')
# 
# BiocManager::install('multtest')
# install.packages('metap')
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
# 
# # Install sctransform from CRAN
# install.packages("sctransform")
# #Install DESeq2
# BiocManager::install('DESeq2')
# #Install RVenn
# BiocManager::install('RVenn')
# install.packages(c("tidyverse", "magrittr", "WGCNA))
library(magrittr)      # provides the %>% operator

# Load the parallel package (comes with base R)
library(parallel)

# single-cell analysis package
# library(Seurat)

# plotting and data science packages
library(cowplot)
library(patchwork)
library(ggVennDiagram)
library(RVenn)
library(ggplot2)
library(plotly)

#!/usr/bin/env Rscript

library(biomaRt)
library(DESeq2)
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(matrixStats)
library(tools)
library(VennDiagram)
library(grid)

# Today's Date
today <- format(Sys.Date(), "%m-%d-%Y")

# # Loading WGCNA processing functions
# source("/home/thilsabeck/Documents/utils_TH/WGCNA_processing_2-21-2025_8-18-2025mod_noSCTorhetcor_enrichR_irGSEA_traitcor.R")#"/home/thilsabeck/Documents/utils_TH/WGCNA_processing_2-21-2025_7-18-2025mod_noSCTorhetcor_enrichR_irGSEA_traitcor.R")#WGCNA_processing_2-21-2025_5-9-2025mod_noSCTorhetcor_enrichR_irGSEA_traitcor.R")# Using newer version that handles single vs multiple trait/gene correlation inputs 5-9-2025 TH #WGCNA_processing_2-21-2025_4-25-2025mod_noSCTorhetcor_enrichR_irGSEA_traitcor.R")#source("/home/thilsabeck/Documents/utils_TH/WGCNA_processing_2-21-2025_4-9-2025mod_noSCT_enrichR_irGSEA.R")#3-24-2025mod_noSCT_enrichR.R")#3-18-2025mod_noSCT_enrichR.R")

## ==== CONFIGURATION ====
counts_dir <- "/home/thilsabeck/Documents/Gage2019/mapping/rmats_analysis/server/counts"
metadata_file <- "/home/thilsabeck/Documents/Gage2019/RNAseq Counts/SeqFibroblast_fibroblasts_metadata_reduced_6-17-2025.csv"
output_dir <- "/home/thilsabeck/Documents/Gage2019/mapping/rmats_analysis/server/WithConfounds/HOMER" # Adding /WithConfounds for reanalysis with confounding variables 8-25-2025
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## ==== LOAD METADATA ====
Gage.metadata_fib <- read.csv(metadata_file, stringsAsFactors = TRUE)
# Ensure factors
Gage.metadata_fib$Diag <- factor(Gage.metadata_fib$Diag)
if ("Batch" %in% colnames(Gage.metadata_fib))
  Gage.metadata_fib$Batch <- factor(Gage.metadata_fib$Batch)

## ==== BIOMART Mapping ====
cat("Fetching transcript to gene mapping...\n")
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
tx2gene <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
                 mart = mart) %>% as_tibble()

## ==== HELPER FUNCTIONS ====
library(dplyr)
library(tibble)

load_and_aggregate_counts <- function(file_path) {
  # Detect delimiter by first line
  first_line <- readLines(file_path, n = 1)
  delimiter <- ifelse(grepl(",", first_line), ",", "\t")
  
  if (delimiter == ",") {
    counts <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    counts <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # Check if first cell contains HOMER command string pattern
  first_colname <- colnames(counts)[1]
  first_cell <- counts[[1]][1]
  # Fix first column if first cell is HOMER command string
  first_col_vals <- counts[[1]]
  if (grepl("^Transcript/RepeatID.*analyzeRepeats.pl", first_col_vals[1])) {
    counts[[1]] <- c(first_col_vals[-1], NA)
    counts <- counts[-nrow(counts), ]
  }
  
  # Assign rownames from first column (transcript)
  rownames(counts) <- counts[[1]]
  counts <- counts[, -1, drop = FALSE]
  
  # Replace "Y", "X", or other non-numeric with 0 and convert to numeric
  counts[] <- lapply(counts, function(col) {
    col <- ifelse(grepl("^[0-9.]+$", col), col, "0")
    as.numeric(col)
  })
  
  counts <- as_tibble(counts, rownames = "transcript")
  num_cols <- setdiff(colnames(counts), "transcript")
  
  if (any(startsWith(counts$transcript, "ENST"))) {
    counts_df <- counts %>%
      left_join(tx2gene, by = c("transcript" = "ensembl_transcript_id")) %>%
      filter(!is.na(ensembl_gene_id)) %>%
      group_by(ensembl_gene_id) %>%
      summarise(across(all_of(num_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop")
    mat <- as.data.frame(counts_df)
    rownames(mat) <- mat$ensembl_gene_id
  } else {
    counts_df <- counts
    mat <- as.data.frame(counts_df)
    rownames(mat) <- mat$transcript
  }
  
  mat <- mat[, num_cols, drop = FALSE]
  as.matrix(mat)
}

variance_filter <- function(mat, threshold = 0.01) {
  mat[rowVars(mat) > threshold, , drop = FALSE]
}

library(tools)  # for file_path_sans_ext
library(lubridate)  # for today()

library(ggplot2)
library(plotly)
library(htmlwidgets)
library(tools)  # for file_path_sans_ext

pca_plot_data <- function(pca_res, metadata_df, plot_title, f, output_dir) {
  pca_df <- as.data.frame(pca_res$x)
  pca_df$Sample_label <- rownames(pca_df)
  merged_df <- merge(pca_df, metadata_df, by.x = "Sample_label", by.y = "Sample", all.x = TRUE)
  
  # Ensure Gender is a factor for shapes
  merged_df$Gender <- factor(merged_df$Gender)
  
  pc1_cor <- cor.test(merged_df$PC1, as.numeric(merged_df$Diag))
  pc2_cor <- cor.test(merged_df$PC2, as.numeric(merged_df$Diag))
  cat(sprintf("[%s] PC1 p=%.3g, PC2 p=%.3g\n", plot_title, pc1_cor$p.value, pc2_cor$p.value))
  
  p <- ggplot(merged_df, aes(PC1, PC2, color = Diag, shape = Gender)) +
    geom_point(size = 3) +
    # scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +  # Adjust as needed based on your Gender levels
    labs(title = sprintf("%s\nPC1 p=%.2g, PC2 p=%.2g", plot_title, pc1_cor$p.value, pc2_cor$p.value)) +
    theme_minimal()
  print(p)
  
  # Prepare output paths replicating file directory structure for PNG and HTML
  today_date <- format(Sys.Date(), "%Y%m%d")
  # relative_path <- dirname(f)
  output_path <- file.path(output_dir, 'PCA_Plots')
  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  
  plot_png_file <- file.path(output_path, paste0(file_path_sans_ext(basename(f)), "_PCA_", today_date, ".png"))
  plot_html_file <- file.path(output_path, paste0(file_path_sans_ext(basename(f)), "_PCA_3D_", today_date, ".html"))
  
  # Save ggplot 2D PCA as PNG high resolution
  ggsave(plot_png_file, plot = p, width = 8, height = 6, dpi = 300)
  cat("Saved 2D PCA plot to:", plot_png_file, "\n")
  
  # Create 3D PCA plotly figure
  fig <- plot_ly(merged_df, x = ~PC1, y = ~PC2, z = ~PC3,
                 color = ~Diag, colors = c("CTRL" = "blue", "AD" = "red"),
                 symbol = ~Gender,
                 type = "scatter3d", mode = "markers") %>%
    layout(title = paste(plot_title, "(3D)"),
           scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")))
  
  # Save 3D plot as standalone HTML preserving interactivity
  saveWidget(fig, file = plot_html_file, selfcontained = TRUE)
  cat("Saved 3D PCA plot as HTML to:", plot_html_file, "\n")
  
  return(list(pca_plot_2d = p, pca_plot_3d = fig))
}

# Attempt DESeq run with error handling, adding a design_formula parameter to include confounding variables 8-25-2025
run_deseq_safe <- function(dds, design_formula) {
  tryCatch({
    design(dds) <- design_formula
    # Relevel the factor of interest if present (example given with 'Diag')
    if ("Diag" %in% colnames(colData(dds))) {
      dds$Diag <- relevel(dds$Diag, ref = "CTRL")
    }
    dds <- DESeq(dds, parallel = TRUE)
    return(dds)
  }, error = function(e) {
    message("DESeq failed with error: ", e$message)
    message("Filtering genes with zero counts across all samples and retrying...")
    
    # Filter genes with zero counts in all samples
    keep <- rowSums(counts(dds)) > 0
    dds_filtered <- dds[keep, ]
    if (nrow(dds_filtered) == 0) {
      warning("No genes remain after filtering zero-count genes, cannot run DESeq.")
      return(NULL)
    }
    design(dds_filtered) <- design_formula
    if ("Diag" %in% colnames(colData(dds_filtered))) {
      dds_filtered$Diag <- relevel(dds_filtered$Diag, ref = "CTRL")
    }
    
    # Retry DESeq
    tryCatch({
      dds_filtered <- DESeq(dds_filtered, parallel = TRUE)
      return(dds_filtered)
    }, error = function(e2) {
      warning("DESeq failed again after filtering: ", e2$message)
      return(NULL)
    })
  })
}

## ==== VENN DIAGRAM HELPERS ====
load_deg_genes <- function(filepath, pval_cutoff=0.05) {
  if(!file.exists(filepath)) return(character(0))
  degs <- read.csv(filepath)
  gene_col <- if("gene" %in% colnames(degs)) "gene" else "gene_id"
  genes <- degs %>%
    filter(!is.na(padj), padj <= pval_cutoff) %>%
    pull(!!sym(gene_col)) %>%
    unique()
  return(genes)
}

save_venn_plot <- function(gene_lists, title, filename, base_fill_color = "blue") {
  # Replace "No significant DEGs" with empty sets
  for (i in seq_along(gene_lists)) {
    if (identical(gene_lists[[i]], "No significant DEGs") ||
        (length(gene_lists[[i]]) == 1 && gene_lists[[i]][1] == "No significant DEGs")) {
      gene_lists[[i]] <- character(0)
    }
  }

  # Create ggVennDiagram plot, venn_obj
  p <- ggVennDiagram(gene_lists, label_alpha = 0) +
    scale_fill_gradient(low = "white", high = base_fill_color) +
    ggtitle(title) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 18)
    ) +
    scale_x_continuous(expand = expansion(mult = 0.2))

  ggsave(filename, plot = p, width = 6, height = 6, dpi = 300)
  message(paste("Saved Venn diagram:", filename))
}

save_upset_plot <- function(gene_lists, title, filename, base_fill_color = "blue") {
  # Clean input: replace "No significant DEGs" with empty sets
  gene_lists <- lapply(gene_lists, function(x) {
    if (identical(x, "No significant DEGs") ||
        (length(x) == 1 && x[1] == "No significant DEGs")) character(0) else x
  })

  filename <- sub("Venn", "UpSet", filename)

  # Use UpSetR::fromList to convert list of sets to binary matrix
  upset_input <- UpSetR::fromList(gene_lists)

  png(filename, width = 900, height = 600, res = 120)
  print(UpSetR::upset(
    upset_input,
    sets = colnames(upset_input),
    sets.bar.color = "#56B4E9",
    order.by = "freq",
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Set Size"
  ))
  grid::grid.force()
  grid.text(title, x = 0.5, y = 0.96, gp = gpar(fontsize = 18, fontface = "bold")) # optional title
  dev.off()
  message(paste("Saved UpSet plot:", filename))
}

# library(ComplexUpset)
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(tibble)
# 
# save_upset_plot <- function(gene_lists, title, filename, base_fill_color = "blue") {
#   # Clean input
#   gene_lists <- lapply(gene_lists, function(x) {
#     if (identical(x, "No significant DEGs") ||
#         (length(x) == 1 && x[1] == "No significant DEGs"))
#       character(0)
#     else x
#   })
#   
#   filename <- sub("Venn", "UpSet", filename)
#   
#   # Convert list into tidy dataframe
#   df <- enframe(gene_lists, name = "set") %>%
#     unnest_longer(value) %>%
#     filter(!is.na(value)) %>%
#     rename(gene = value) %>%
#     mutate(present = TRUE)
#   
#   # Convert to wide format: one row per gene, columns are membership TRUE/FALSE
#   df_wide <- df %>%
#     pivot_wider(names_from = set, values_from = present, values_fill = FALSE)
#   
#   upset_sets <- colnames(df_wide)[-1]
#   
#   p <- ggplot(df_wide) +
#     ComplexUpset::upset(
#       data = df_wide,
#       intersect = upset_sets,
#       base_annotations = list(
#         'Intersection size' = ComplexUpset::intersection_size(),
#         'Set size' = ComplexUpset::upset_set_size()
#       )
#     ) +
#     ggtitle(title) +
#     theme(plot.title = element_text(hjust = 0.5, size = 18))
#   
#   pdf(filename, width = 11, height = 8.5)
#   print(p)
#   dev.off()
#   
#   message(paste("Saved ComplexUpset plot:", filename))
# }

## ==== MAIN PIPELINE ====

files <- list.files(counts_dir, pattern = "\\.(csv|tsv)$", full.names = TRUE)

design_formula = ~ Diag + Age + Gender
reduced_formula <- ~ Age + Gender

# Parse trimmer and aligner from filenames, assuming consistent naming: trimmer_aligner_*.csv
extract_trimmer_aligner <- function(fname) {
  base <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(fname)))
  parts <- unlist(strsplit(base, "_"))
  # Adjust depending on your naming; here first two parts:
  trimmer <- parts[1]
  if (toupper(parts[2]) %in% c("STAR", "GSNAP")) {
    if (toupper(parts[2]) == "GSNAP" && parts[3] == "2024genome") {
      aligner <- paste(parts[2], parts[3], sep="_")
      # Filter out "concatenated" and "counts" from parts[4:length(parts)]
      filtered_parts <- parts[4:length(parts)][!parts[4:length(parts)] %in% c("concatenated", "counts", "merged")]
    } else if (toupper(parts[2]) == "GSNAP" && parts[3] != "2024genome") {
      aligner <- parts[2]
      # Filter out "concatenated" and "counts" from parts[4:length(parts)]
      filtered_parts <- parts[3:length(parts)][!parts[3:length(parts)] %in% c("concatenated", "counts", "merged")]
    } else if (toupper(parts[2]) == "STAR") {
      aligner <- paste(parts[2], parts[3], sep="_")
      # Filter out "concatenated" and "counts" from parts[4:length(parts)]
      filtered_parts <- parts[4:length(parts)][!parts[4:length(parts)] %in% c("concatenated", "counts", "merged")]
      # Replace "ReadsPerGene" in filtered_parts with "star"
      filtered_parts <- gsub("ReadsPerGene", "star", filtered_parts)
    }
    # Paste together with underscore separator
    counter <- paste(filtered_parts, collapse = "_")
  } else{
    aligner <- parts[2]
    # Filter out "concatenated" and "counts" from parts[4:length(parts)]
    filtered_parts <- parts[4:length(parts)][!parts[4:length(parts)] %in% c("concatenated", "counts", "merged")]
    # Paste together with underscore separator
    counter <- paste(filtered_parts, collapse = "_")
  }
  list(trimmer=trimmer, aligner=aligner, counter=counter)
}

extract_trimmer_aligner <- function(fname) {#, trimmer_list, aligner_list, counter_list) {
  trimmer_list <- c("trimmomatic", "trimgalore", "fastp") # example trimmers
  aligner_list <- c("STAR", "STAR", "GSNAP", "hisat2")
  counter_list <- c("ReadsPerGene", "star", "htseq", "homer_exon_counts", "homer_exon_counts_unique", "homer_gene_counts", "homer_gene_counts_unique", "homer_gene_condensed_counts", "homer_gene_condensed_counts_unique") # if applicable
  base <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(fname)))
  parts <- unlist(strsplit(base, "_"))
  
  # Find matching trimmer from parts, default to first part if none matched
  trimmer <- if(any(parts %in% trimmer_list)) parts[parts %in% trimmer_list][1] else parts[1]
  
  # Determine aligner combining next parts if needed
  aligner_candidate_parts <- parts[-1]
  aligner <- NULL
  if (toupper(aligner_candidate_parts[1]) %in% aligner_list) {
    # Special case for gsnap with "2024genome"
    if (toupper(aligner_candidate_parts[1]) == "GSNAP" && length(aligner_candidate_parts) > 1 && aligner_candidate_parts[2] == "2024genome") {
      aligner <- paste(aligner_candidate_parts[1:2], collapse = "_")
      filtered_parts <- aligner_candidate_parts[-(1:2)]
    } else if (toupper(aligner_candidate_parts[1]) == "STAR" && length(aligner_candidate_parts) > 1 && aligner_candidate_parts[2] %in% c("single", "twopass")) {
      aligner <- paste(aligner_candidate_parts[1:2], collapse = "_")
      filtered_parts <- aligner_candidate_parts[-(1:2)]
    } else {
      aligner <- aligner_candidate_parts[1]
      filtered_parts <- aligner_candidate_parts[-1]
    }
  } else {
    aligner <- aligner_candidate_parts[1]
    filtered_parts <- aligner_candidate_parts[-1]
  }
  
  # Filter out undesired counter tokens
  filtered_parts <- filtered_parts[!filtered_parts %in% c("concatenated", "counts", "merged")]
  
  # Replace e.g. "ReadsPerGene" with "star" in counter (optional)
  filtered_parts <- gsub("ReadsPerGene", "star", filtered_parts)
  
  # Reconstruct counter using only allowed tokens if provided
  if (!missing(counter_list)) {
    filtered_parts <- filtered_parts[
      sapply(filtered_parts, function(fp) any(grepl(fp, counter_list)))
    ]
  }
  counter <- paste(filtered_parts, collapse = "_")
  
  list(trimmer = trimmer, aligner = aligner, counter = counter)
}

# Collect DEGs info for Venn after DEG extraction
deg_summary <- tibble(
  file = character(),
  trimmer = character(),
  aligner = character(),
  counter = character(),
  degs_wald = list(),
  degs_lrt = list()
)

# Log file for count file errors
countFileErrorsLog = c()
countFileErrorsLog_path <- file.path(output_dir, paste0("CountFileErrorsLog_", today, ".txt"))

for (f in files) {#[97:length(files)]) {
  tryCatch({
    cat("Processing file:", f, "\n")
    
    # Not skipping homer files 11-6-2025 # Skipping if "homer" is in file name
    # if (grepl("homer", f, ignore.case = TRUE)) {
    #   cat("Skipping homer file:", f, "\n")
    #   next
    # }
    
    mat <- load_and_aggregate_counts(f)
    
    samples <- colnames(mat)
    
    coldata <- Gage.metadata_fib[Gage.metadata_fib$Sample %in% samples, ]
    rownames(coldata) <- coldata$Sample
    
    # HOMER counts could have X for transcripts that did not have any reads
    # This aligns with HOMER's counting output conventions, where "X" is a missing count marker rather than a numeric value.
    
    coldata <- coldata[match(colnames(mat), rownames(coldata)), ]
    if (any(is.na(rownames(coldata)))) {
      stop("Sample names not matching between metadata and counts")
    }
    
    # Remove any NaN from mat
    mat[is.na(mat)] <- 0
    
    # Create DESeqDataSet with full data and design
    dds <- DESeqDataSetFromMatrix(countData = round(mat)+1,  # Adding 1 to avoid zeros
                                  colData = coldata,
                                  design = design_formula)
    dds$Diag <- relevel(dds$Diag, ref = "CTRL")
    
    # Run DESeq2 with LRT to test global significance of condition (Diag)
    dds_lrt <- DESeq(dds, test = "LRT", reduced = reduced_formula, parallel = TRUE)
    
    # Get all coefficient names after LRT model fit
    coef_names <- resultsNames(dds_lrt)
    
    # Prepare empty list to hold per-coefficient significant gene data frames
    sig_list <- list()
    
    # Loop through each coefficient except intercept
    for (coef in coef_names[-1]) {
      # Extract Wald test results for this coefficient
      res_wald <- results(dds_lrt, name = coef)
      
      # Filter for padj < 0.05 significance
      sig_genes <- which(res_wald$padj < 0.05 & !is.na(res_wald$padj))
      if (length(sig_genes) > 0) {
        res_df <- as.data.frame(res_wald[sig_genes, ])
        res_df$Gene <- rownames(res_wald)[sig_genes]
        res_df$Coefficient <- coef
        sig_list[[coef]] <- res_df
      }
    }
    
    # Combine all significant gene data frames into one table
    res_sig_all <- do.call(rbind, sig_list)
    
    # Optionally reorder columns to have Gene and Coefficient up front
    res_sig_all <- res_sig_all[, c("Gene", "Coefficient", setdiff(colnames(res_sig_all), c("Gene", "Coefficient")))]
    
    # Write combined results to CSV
    deg_lrt_file <- file.path(output_dir, paste0(file_path_sans_ext(basename(f)), "_LRT_DEGs_AllCoeffs_", today, ".csv"))
    write.csv(res_sig_all, deg_lrt_file, row.names = FALSE)
    
    cat(sprintf("Saved combined significant DE genes for all coefficients to %s\n", deg_lrt_file))
    
    # Run Wald test DESeq safely on full data (all genes)
    dds_safe <- run_deseq_safe(dds, design_formula = design_formula)
    if (is.null(dds_safe)) {
      cat("Skipping DESeq analysis due to errors.\n")
      next
    } else {
      dds <- dds_safe
    }
    
    # Variance stabilizing transform for visualization
    vst_data <- tryCatch(
      vst(dds, blind = FALSE),
      error = function(e) {
        message("vst() failed, using varianceStabilizingTransformation(): ", e$message)
        varianceStabilizingTransformation(dds, blind = FALSE)
      }
    )
    vst_counts <- assay(vst_data)
    out_vst <- file.path(output_dir, paste0(file_path_sans_ext(basename(f)), "_VST_", today, ".csv"))
    write.csv(vst_counts, out_vst)
    
    vst_counts_filt <- variance_filter(vst_counts, threshold = 0.01)
    pca_res <- prcomp(t(vst_counts_filt), center = TRUE, scale. = TRUE)
    pca_plot_data(pca_res, coldata, paste("PCA -", file_path_sans_ext(basename(f))), f, output_dir)
    
    # Extract Wald test results from full data
    res <- results(dds)
    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene")
    
    # Labeling genes if significant in LRT
    res_df <- res_df %>%
      mutate(significant_in_LRT = ifelse(gene %in% res_sig_all$Gene, TRUE, FALSE))
    
    # Filter Wald results to LRT significant genes for downstream analysis
    res_sig_df <- res_df[order(-res_df$significant_in_LRT, res_df$padj), ]
    
    deg_file <- file.path(output_dir, paste0(file_path_sans_ext(basename(f)), "_DEGs_", today, ".csv"))
    write.csv(res_sig_df, deg_file, row.names = FALSE)
    
    # Store DEG gene list info for Venn plotting
    ta <- extract_trimmer_aligner(f)
    deg_genes <- res_df %>% filter(!is.na(padj), padj <= 0.05) %>% pull(gene) %>% unique()
    if (length(deg_genes) == 0) {
      cat("No significant DEGs found for file:", f, "\n")
      deg_genes <- c("No significant DEGs")
    }
    deg_summary <- deg_summary %>% add_row(
      file = f,
      trimmer = ta$trimmer,
      aligner = ta$aligner,
      counter = ta$counter,
      degs_wald = list(deg_genes),
      degs_lrt = list(rownames(mat)[sig_genes])
    )
  }, error = function(e) {
    message(paste("Error loading counts from file:", f, " - ", e$message))
    countFileErrorsLog <<- c(countFileErrorsLog, paste("Error loading counts from file:", f, " - ", e$message))
    return(NULL)
  })
}
# Write deg_summary to .csv
deg_summary <- deg_summary[deg_summary$counter != "",]
# Convert list-columns to character strings by collapsing elements with ";"
deg_summary_flat <- deg_summary %>%
  mutate(
    degs_wald = sapply(degs_wald, function(x) paste(x, collapse = ";")),
    degs_lrt = sapply(degs_lrt, function(x) paste(x, collapse = ";"))
  )

deg_summary_file <- file.path(output_dir, paste0("DEG_Summary_", today, ".csv"))
write.csv(deg_summary_flat, deg_summary_file, row.names = FALSE)

# Write countFileErrorsLog to file
if (length(countFileErrorsLog) > 0) {
  writeLines(countFileErrorsLog, con = countFileErrorsLog_path)
  cat("Logged count file errors to:", countFileErrorsLog_path, "\n")
}

cat("DEG analysis complete, generating Venn diagrams...\n")

# Prepare lists for venn: aligners within each trimmer
venn_by_trimmer <- list()
for(t in unique(deg_summary$trimmer)) {
  sub <- deg_summary %>% filter(trimmer == t)
  gene_lists <- list()
  for(a in unique(sub$aligner)) {
    genes <- sub %>% filter(aligner == a) %>% pull(degs_wald) %>% unlist()
    gene_lists[[toupper(a)]] <- unique(genes)
  }
  venn_by_trimmer[[t]] <- gene_lists
}

# Prepare lists for venn: trimmers within each aligner
venn_by_aligner <- list()
for(a in unique(deg_summary$aligner)) {
  sub <- deg_summary %>% filter(aligner == a)
  gene_lists <- list()
  for(t in unique(sub$trimmer)) {
    genes <- sub %>% filter(trimmer == t) %>% pull(degs_wald) %>% unlist()
    gene_lists[[toupper(t)]] <- unique(genes)
  }
  venn_by_aligner[[a]] <- gene_lists
}

# Prepare lists for venn: gene sets grouped by unique trimmer-aligner combination (counter)
venn_by_counter <- list()
padj_threshold <- 0.05

for(cntr in unique(deg_summary$counter)) {
  sub <- deg_summary %>% filter(counter == cntr)
  gene_lists <- list()
  for(t in unique(sub$trimmer)) {
    for(a in unique(sub$aligner)) {
      # Extract all DESeq2 result data frames or tables
      res_list <- sub %>% filter(trimmer == t, aligner == a) %>% pull(degs_wald)
      # Check if res_list is a list of duplicate lists and remove one
      if(length(res_list) == 2 && identical(res_list[[1]], res_list[[2]])) {
        res_list <- res_list[1]
      }
      
      # Combine and filter by padj, then get unique gene names
      sig_genes <- unique(unlist(lapply(res_list, function(df) {
        if(is.data.frame(df) && "padj" %in% colnames(df) && "gene" %in% colnames(df)) {
          genes_sig <- df$gene[df$padj < padj_threshold]
          return(genes_sig)
        } else {
          # fallback if gene names only
          return(character(0))
        }
      })))
      
      if(length(sig_genes) > 0) {
        combined_name <- paste0(toupper(t), "_", toupper(a))
        gene_lists[[combined_name]] <- sig_genes
      }
    }
  }
  venn_by_counter[[cntr]] <- gene_lists
}

# Prepare lists for venn: aligners within each trimmer and counter, named by "trimmer_counter"
venn_by_trimmer_counter <- list()
for (t in unique(deg_summary$trimmer)) {
  sub_t <- deg_summary %>% filter(trimmer == t)
  counters <- unique(sub_t$counter)  # replace 'counter' with your actual column name
  counter_lists <- list()
  for (c in counters) {
    sub_c <- sub_t %>% filter(counter == c)
    gene_lists <- list()
    for (a in unique(sub_c$aligner)) {
      genes <- sub_c %>% filter(aligner == a) %>% pull(degs_wald) %>% unlist()
      gene_lists[[toupper(a)]] <- unique(genes)
    }
    # Use concatenated name as key
    counter_name <- paste0(t, "_", c)
    counter_lists[[counter_name]] <- gene_lists
  }
  venn_by_trimmer_counter[[t]] <- counter_lists
}

# Prepare lists for venn: trimmers within each aligner and counter, named by "aligner_counter"
venn_by_aligner_counter <- list()
for (a in unique(deg_summary$aligner)) {
  sub_a <- deg_summary %>% filter(aligner == a)
  counters <- unique(sub_a$counter)
  counter_lists <- list()
  for (c in counters) {
    sub_c <- sub_a %>% filter(counter == c)
    gene_lists <- list()
    for (t in unique(sub_c$trimmer)) {
      genes <- sub_c %>% filter(trimmer == t) %>% pull(degs_wald) %>% unlist()
      gene_lists[[toupper(t)]] <- unique(genes)
    }
    counter_name <- paste0(a, "_", c)
    counter_lists[[counter_name]] <- gene_lists
  }
  venn_by_aligner_counter[[a]] <- counter_lists
}

# Prepare lists for venn: aligners within each trimmer
venn_by_trimmer <- list()
for(t in unique(deg_summary$trimmer)) {
  sub <- deg_summary %>% filter(trimmer == t)
  gene_lists <- list()
  for(a in unique(sub$aligner)) {
    genes <- sub %>% filter(aligner == a) %>% pull(degs_wald) %>% unlist() %>% unique()
    gene_lists[[toupper(a)]] <- genes
  }
  venn_by_trimmer[[t]] <- gene_lists
}

# Prepare lists for venn: trimmers within each aligner
venn_by_aligner <- list()
for(a in unique(deg_summary$aligner)) {
  sub <- deg_summary %>% filter(aligner == a)
  gene_lists <- list()
  for(t in unique(sub$trimmer)) {
    genes <- sub %>% filter(trimmer == t) %>% pull(degs_wald) %>% unlist() %>% unique()
    gene_lists[[toupper(t)]] <- genes
  }
  venn_by_aligner[[a]] <- gene_lists
}

# Prepare lists for venn: aligners within each trimmer and counter, named by "trimmer_counter"
venn_by_trimmer_counter <- list()
for (t in unique(deg_summary$trimmer)) {
  sub_t <- deg_summary %>% filter(trimmer == t)
  counters <- unique(sub_t$counter)
  counter_lists <- list()
  for (c in counters) {
    sub_c <- sub_t %>% filter(counter == c)
    gene_lists <- list()
    for (a in unique(sub_c$aligner)) {
      genes <- sub_c %>% filter(aligner == a) %>% pull(degs_wald) %>% unlist() %>% unique()
      gene_lists[[toupper(a)]] <- genes
    }
    counter_name <- paste0(t, "_", c)
    counter_lists[[counter_name]] <- gene_lists
  }
  venn_by_trimmer_counter[[t]] <- counter_lists
}

# Prepare lists for venn: aligners within each trimmer and counter, named by "trimmer_aligner"
venn_by_trimmer_aligner <- list()
for (t in unique(deg_summary$trimmer)) {
  sub_t <- deg_summary %>% filter(trimmer == t)
  aligners <- unique(sub_t$aligner)
  aligner_lists <- list()
  for (a in aligners) {
    sub_a <- sub_t %>% filter(aligner == a)
    gene_lists <- list()
    for (c in unique(sub_a$counter)) {
      genes <- sub_a %>% filter(counter == c) %>% pull(degs_wald) %>% unlist() %>% unique()
      gene_lists[[toupper(c)]] <- genes
    }
    aligner_name <- paste0(t, "_", a)
    aligner_lists[[aligner_name]] <- gene_lists
  }
  venn_by_trimmer_aligner[[t]] <- aligner_lists
}

# Prepare lists for venn: trimmers within each aligner and counter, named by "aligner_counter"
venn_by_aligner_counter <- list()
for (a in unique(deg_summary$aligner)) {
  sub_a <- deg_summary %>% filter(aligner == a)
  counters <- unique(sub_a$counter)
  counter_lists <- list()
  for (c in counters) {
    sub_c <- sub_a %>% filter(counter == c)
    gene_lists <- list()
    for (t in unique(sub_c$trimmer)) {
      genes <- sub_c %>% filter(trimmer == t) %>% pull(degs_wald) %>% unlist() %>% unique()
      gene_lists[[toupper(t)]] <- genes
    }
    counter_name <- paste0(a, "_", c)
    counter_lists[[counter_name]] <- gene_lists
  }
  venn_by_aligner_counter[[a]] <- counter_lists
}

# Create output directory for venn plots
venn_dir <- file.path(output_dir, "venn_plots")
dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)

# Setup directories for saving overlap CSVs
overlap_dir <- file.path(output_dir, "overlap_DEGs")
dir.create(overlap_dir, showWarnings = FALSE, recursive = TRUE)

# Fetch ENSEMBL-to-symbol mapping (change 'hsapiens_gene_ensembl' for other species)
cat("Fetching gene symbol mapping from biomaRt...\n")
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = mart)

# Helper: annotate Ensembl IDs with gene symbol
get_symbols <- function(overlap_df, gene_map) {
  overlap_df %>%
    left_join(gene_map, by = "ensembl_gene_id")
}

# Helper function to convert GSNAP RNA IDs to Ensembl gene IDs
convert_gsnap_to_ensembl <- function(gene_ids, mart) {
  # Remove "rna-" prefix
  accession_ids <- sub("^rna-", "", gene_ids)
  
  # Remove version suffix from RefSeq accessions (e.g., NM_001008.4 -> NM_001008)
  accession_ids_noversion <- sub("\\.\\d+$", "", accession_ids)
  
  # Query Ensembl biomart for RefSeq mRNA and ncRNA accessions
  mapping_mrna <- getBM(
    attributes = c("ensembl_gene_id", "refseq_mrna"),
    filters = "refseq_mrna",
    values = accession_ids_noversion,
    mart = mart
  )
  
  mapping_ncrna <- getBM(
    attributes = c("ensembl_gene_id", "refseq_ncrna"),
    filters = "refseq_ncrna",
    values = accession_ids_noversion,
    mart = mart
  )
  
  mapping <- dplyr::bind_rows(mapping_mrna, mapping_ncrna) %>%
    dplyr::distinct()
  
  # Create mapping of RefSeq to Ensembl gene, merging both refseq_mrna and refseq_ncrna in one vector
  unique_map <- mapping %>%
    tidyr::pivot_longer(cols = c("refseq_mrna", "refseq_ncrna"), names_to = "type", values_to = "refseq") %>%
    dplyr::filter(!is.na(refseq)) %>%
    dplyr::distinct(refseq, ensembl_gene_id)
  
  # Match input accession IDs to Ensembl gene IDs, introducing NAs for unmapped
  ensembl_ids <- unique_map$ensembl_gene_id[match(accession_ids_noversion, unique_map$refseq)]
  return(ensembl_ids)
}

# Modified save_gene_lists_and_overlap function with conversion step
save_gene_lists_and_overlap <- function(gene_lists, comparison_name, fill_cols, venn_dir, overlap_dir, mart, gene_map) {
  
  # Convert GSNAP gene IDs to Ensembl IDs within each list if necessary
  gene_lists_ensembl <- lapply(gene_lists, function(genes) {
    # Detect if genes look like GSNAP IDs by regex starting with rna-
    if(length(genes) > 0) {
      # Identify which genes are transcript-like (starting with "rna-")
      is_transcript <- grepl("^rna-", genes)
      
      if(all(is_transcript)) {
        # All are transcripts: convert all
        ensembl_ids <- convert_gsnap_to_ensembl(genes, mart)
        ensembl_ids <- ensembl_ids[!is.na(ensembl_ids)]
        return(ensembl_ids)
      } else if(any(is_transcript)) {
        # Mixed transcript and gene IDs
        # Convert only transcript-like IDs
        converted <- convert_gsnap_to_ensembl(genes[is_transcript], mart)
        converted <- converted[!is.na(converted)]
        
        # Keep existing Ensembl IDs as is
        ensembl_ids <- c(converted, genes[!is_transcript])
        # Remove possible duplicates
        ensembl_ids <- unique(ensembl_ids)
        return(ensembl_ids)
      } else {
        # No transcripts, all assumed Ensembl IDs, return as is
        return(genes)
      }
    } else {
      return(genes)
    }
  })
  
  # Prepare dataframe with gene lists using Ensembl IDs
  gene_lists_df <- gene_lists_ensembl %>%
    purrr::imap_dfr(~ tibble(sample = .y, ensembl_gene_id = .x), .id = NULL) %>%
    left_join(gene_map, by = "ensembl_gene_id")
  
  write.csv(gene_lists_df,
            file.path(overlap_dir, paste0("all_DEGs_", comparison_name, "_", today, ".csv")),
            row.names = FALSE)
  
  valid_lists <- gene_lists_ensembl[sapply(gene_lists_ensembl, function(x) length(x) > 0 && all(x != "No significant DEGs"))]
  
  # valid_lists: named list of gene vectors
  # get the list names
  list_names <- names(valid_lists)
  n <- length(valid_lists)
  
  # Initialize list to hold overlap data frames
  all_overlaps <- list()
  
  # Loop for combination size from 2 to n
  for (combo_size in 2:n) {
    combos <- combn(list_names, combo_size, simplify = FALSE)
    for (combo in combos) {
      # perform intersection across >2 sets with Reduce
      overlap_genes <- Reduce(intersect, gene_lists_ensembl[combo])
      combo_name <- paste(combo, collapse = "_and_")
      all_overlaps[[combo_name]] <- overlap_genes
    }
  }
  
  # Combine results in a single data frame
  overlap_df <- do.call(rbind, lapply(names(all_overlaps), function(nm) {
    if(length(all_overlaps[[nm]]) > 0)
      data.frame(combination = nm, ensembl_gene_id = all_overlaps[[nm]], stringsAsFactors = FALSE)
  }))
  
  # Save to CSV
  overlap_annotated <- get_symbols(overlap_df, gene_map)
  csv_file <- file.path(overlap_dir, paste0("overlap_DEGs_", comparison_name, "_", today, ".csv"))
  write.csv(overlap_annotated, csv_file, row.names = FALSE)
  
  overlap_row <- tibble(
    comparison = comparison_name,
    overlap_genes = if(length(overlap_genes) > 0) paste(overlap_genes, collapse = ";") else "No overlap"
  )
  write.csv(overlap_row,
            file.path(overlap_dir, paste0("overlap_DEGs_", comparison_name, "_summary.csv")),
            row.names = FALSE)
  
  png_file <- file.path(venn_dir, paste0("Venn_DEGs_", comparison_name, "_", today, ".png"))
  n_sets <- length(gene_lists_ensembl)
  if (n_sets <= 6) {
    save_venn_plot(gene_lists_ensembl, paste("DEGs by", gsub("_", " ", comparison_name)), png_file)
  } else if (n_sets == 7) {
    save_venn_plot(gene_lists_ensembl, paste("DEGs by", gsub("_", " ", comparison_name)), png_file)
    save_upset_plot(gene_lists_ensembl, paste("DEGs by", gsub("_", " ", comparison_name)), png_file)
  } else {
    save_upset_plot(gene_lists_ensembl, paste("DEGs by", gsub("_", " ", comparison_name)), png_file)
  }
}

# Venn and overlap for aligners within each trimmer
for (t in names(venn_by_trimmer)) {
  gene_lists <- venn_by_trimmer[[t]]
  if (length(gene_lists) > 1 && sum(sapply(gene_lists, length) > 0) >= 2) {
    # fill_cols <- c("red", "green", "blue")[seq_along(gene_lists)]
    save_gene_lists_and_overlap(gene_lists,
                                paste0("aligners_within_", t),
                                fill_cols,
                                venn_dir,
                                overlap_dir,
                                mart,
                                gene_map)
  }
}

# Venn and overlap for trimmers within each aligner
for (a in names(venn_by_aligner)) {
  gene_lists <- venn_by_aligner[[a]]
  if (length(gene_lists) > 1 && sum(sapply(gene_lists, length) > 0) >= 2) {
    # fill_cols <- c("orange", "purple", "cyan")[seq_along(gene_lists)]
    save_gene_lists_and_overlap(gene_lists,
                                paste0("trimmers_within_", a),
                                fill_cols,
                                venn_dir,
                                overlap_dir,
                                mart,
                                gene_map)
  }
}

# Venn and overlap for trimmer-aligner combinations split by counter
for (counter_name in names(venn_by_counter)) {
  gene_lists <- venn_by_counter[[counter_name]]
  if (length(gene_lists) > 1 && sum(sapply(gene_lists, length) > 0) >= 2) {
    # Define fill colors or use default
    # fill_cols <- c("orange", "purple", "cyan")[seq_along(gene_lists)]
    save_gene_lists_and_overlap(
      gene_lists = gene_lists,
      comparison_name = paste0("Trimmer_Aligner_within_", counter_name),
      fill_cols = fill_cols,
      venn_dir = venn_dir,
      overlap_dir = overlap_dir,
      mart,
      gene_map
    )
  }
}

# Venn and overlap for aligners within each trimmer_counter combination
for (t in names(venn_by_trimmer_counter)) {
  combo_lists <- venn_by_trimmer_counter[[t]]
  if (length(combo_lists) > 1 && sum(sapply(combo_lists, length) > 0) >= 2) {
    for (c in names(combo_lists)) {
      gene_lists <- combo_lists[[c]]
      if (length(gene_lists) > 1 && sum(sapply(gene_lists, length) > 0) >= 2) {
        # fill_cols <- c("red", "green", "blue")[seq_along(gene_lists)]
        save_gene_lists_and_overlap(gene_lists,
                                    paste0("aligners_within_", c),
                                    fill_cols,
                                    venn_dir,
                                    overlap_dir,
                                    mart,
                                    gene_map)
      }
    }
  }
}

# Venn and overlap for trimmers within each aligner_counter combination
for (a in names(venn_by_aligner_counter)) {
  combo_lists <- venn_by_aligner_counter[[a]]
  if (length(combo_lists) > 1 && sum(sapply(combo_lists, length) > 0) >= 2) {
    for (c in names(combo_lists)) {
      gene_lists <- combo_lists[[c]]
      if (length(gene_lists) > 1 && sum(sapply(gene_lists, length) > 0) >= 2) {
        # fill_cols <- c("red", "green", "blue")[seq_along(gene_lists)]
        save_gene_lists_and_overlap(gene_lists,
                                    paste0("trimmers_within_", c),
                                    fill_cols,
                                    venn_dir,
                                    overlap_dir,
                                    mart,
                                    gene_map)
      }
    }
  }
}

# Venn and overlap for trimmers within each aligner_counter combination
for (t in names(venn_by_trimmer_aligner)) {
  combo_lists <- venn_by_trimmer_aligner[[t]]
  if (length(combo_lists) > 1 && sum(sapply(combo_lists, length) > 0) >= 2) {
    for (c in names(combo_lists)) {
      gene_lists <- combo_lists[[c]]
      if (length(gene_lists) > 1 && sum(sapply(gene_lists, length) > 0) >= 2) {
        # fill_cols <- c("red", "green", "blue")[seq_along(gene_lists)]
        save_gene_lists_and_overlap(gene_lists,
                                    paste0("counters_within_", c),
                                    fill_cols,
                                    venn_dir,
                                    overlap_dir,
                                    mart,
                                    gene_map)
      }
    }
  }
}

# 
# # Plot and save venn diagrams for aligners within trimmers
# for(t in names(venn_by_trimmer)) {
#   gene_lists <- venn_by_trimmer[[t]]
#   if(length(gene_lists) > 1 && sum(sapply(gene_lists, length) > 0) >= 2) {
#     png_file <- file.path(venn_dir, paste0("Venn_DEGs_aligners_within_", t, ".png"))
#     fill_cols <- c("red", "green", "blue")[seq_along(gene_lists)]
#     save_venn_plot(gene_lists, paste("DEGs by Aligners within", t), png_file, fill_cols)
#   }
# }
# 
# # Plot and save venn diagrams for trimmers within aligners
# for(a in names(venn_by_aligner)) {
#   gene_lists <- venn_by_aligner[[a]]
#   if(length(gene_lists) > 1 && sum(sapply(gene_lists, length) > 0) >= 2) {
#     png_file <- file.path(venn_dir, paste0("Venn_DEGs_trimmers_within_", a, ".png"))
#     fill_cols <- c("orange", "purple", "cyan")[seq_along(gene_lists)]
#     save_venn_plot(gene_lists, paste("DEGs by Trimmers within", toupper(a)), png_file, fill_cols)
#   }
# }

cat("Venn diagram generation complete. Pipeline finished.\n")
