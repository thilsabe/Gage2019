library(vegan)
library(stringr)

# Function to extract clean sample basename without suffixes like _R1_paired.fq
extract_clean_sample_name <- function(path_str) {
  base <- basename(path_str)
  # Remove suffix like _R1_paired.fq or _R1_paired.fastq.gz if present
  clean_name <- str_replace(base, "_R1_paired.*$", "")
  return(clean_name)
}

# Load mash distances long table
mash_tsv <- "mash_analysis_3000sketchsize_15kmersize/mash_distances.tsv"
df <- read.delim(mash_tsv, header=FALSE, stringsAsFactors=FALSE)
colnames(df) <- c("query", "reference", "distance", "pval", "shared")

# Extract clean sample names from query and reference
df$query_sample <- sapply(df$query, extract_clean_sample_name)
df$reference_sample <- sapply(df$reference, extract_clean_sample_name)

# Create unique sorted sample list
samples <- sort(unique(c(df$query_sample, df$reference_sample)))
n <- length(samples)

# Initialize empty matrix
mat <- matrix(NA, nrow=n, ncol=n, dimnames=list(samples, samples))

# Fill matrix: assign distances to correct positions
for(i in seq_len(nrow(df))) {
  mat[df$query_sample[i], df$reference_sample[i]] <- df$distance[i]
}

# Make symmetric and fill diagonal zeros
for(i in seq_len(n)) {
  for(j in seq_len(n)) {
    if(is.na(mat[i,j])) mat[i,j] <- mat[j,i]
  }
}
diag(mat)[is.na(diag(mat))] <- 0

# Convert to dist object
dist_matrix <- as.dist(mat)

# Load metadata
metadata_csv <- "Gage2019_samples_groups.csv"
metadata <- read.csv(metadata_csv, stringsAsFactors=FALSE)
metadata$Sample <- as.character(metadata$Sample)
rownames(metadata) <- metadata$Sample

# Filter metadata to samples present in distance matrix
metadata <- metadata[rownames(metadata) %in% samples, ]

# Convert group variables to factors
metadata$Group <- factor(metadata$Group)
metadata$Sex <- factor(metadata$Sex)

# Run adonis2 with Group, Sex, and their interaction
result <- adonis2(dist_matrix ~ Group * Sex, data=metadata, permutations=999)

print(result)

# Save to output file
output_file <- "Gage2019_mashDistance_3000sketchsize_15kmersize_GroupSexInter_permanova_results_10_7_2025.txt"
capture.output(result, file=output_file)
