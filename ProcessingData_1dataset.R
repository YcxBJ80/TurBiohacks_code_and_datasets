# ==========================================================
# R Script: GSE25097 Series Matrix → CSV with Probe → Gene Mapping + Train/Validation Split
# Process GSE25097 dataset, mark tumor samples as 1, other samples (healthy, cirrhotic, non_tumor) as 0
# ==========================================================

library(GEOquery)

# ------------------------------
# 1. Set input and output paths
# ------------------------------
input_file <- "/Users/yangchengxuan/Downloads/Oncology Dataset/GSE25097 Series Matrix.txt"
output_file <- "GSE25097.csv"         # Full dataset CSV
train_file  <- "GSE25097_train.csv"   # Training set CSV
valid_file  <- "GSE25097_valid.csv"   # Validation set CSV
gpl_id <- "GPL10687"                   # GEO platform ID

# ------------------------------
# 2. Safe function to read Series Matrix
# ------------------------------
read_series_matrix_safe <- function(file_path) {
  lines <- readLines(file_path)
  start <- which(grepl("^!series_matrix_table_begin", lines)) + 1
  end   <- which(grepl("^!series_matrix_table_end", lines)) - 1
  
  expr_lines <- lines[start:end]
  expr_lines <- expr_lines[!grepl("^!", expr_lines)]  # remove comment lines
  expr_lines <- expr_lines[nchar(expr_lines) > 0]     # remove empty lines
  
  expr_con <- textConnection(expr_lines)
  exprs_raw <- read.table(expr_con, sep = "\t", header = TRUE,
                          check.names = FALSE, stringsAsFactors = FALSE)
  close(expr_con)
  
  rownames(exprs_raw) <- exprs_raw[,1]  # first column = probe ID
  exprs_raw <- exprs_raw[,-1]
  return(exprs_raw)
}

# ------------------------------
# 3. Safe function to map probe → gene
# ------------------------------
map_probe_to_gene_safe <- function(exprs_matrix, gpl_id) {
  cat("🔗 Mapping probe IDs to gene symbols for", gpl_id, "...\n")
  
  gpl <- getGEO(gpl_id, destdir = tempdir())
  probe_table <- Table(gpl)
  
  if(!"GeneSymbol" %in% colnames(probe_table)){
    stop("❌ GPL platform does not have 'GeneSymbol' column")
  }
  
  probe2gene <- probe_table[, c("ID", "GeneSymbol")]
  probe2gene$GeneSymbol <- as.character(probe2gene$GeneSymbol)
  
  probe2gene <- probe2gene[probe2gene$ID %in% rownames(exprs_matrix), ]
  
  probe2gene$Gene <- sapply(probe2gene$GeneSymbol, function(x){
    x <- trimws(x)
    if(is.na(x) || x == "" || x == " ") return(NA)
    strsplit(x, "///")[[1]][1]
  })
  
  probe2gene <- probe2gene[!is.na(probe2gene$Gene) & probe2gene$Gene != "", ]
  
  if(nrow(probe2gene) == 0){
    stop("❌ No valid GeneSymbol found for this platform")
  }
  
  exprs_matrix <- exprs_matrix[probe2gene$ID, , drop=FALSE]
  exprs_matrix <- cbind(exprs_matrix, Gene = probe2gene$Gene)
  
  df <- as.data.frame(exprs_matrix)
  sample_cols <- setdiff(colnames(df), "Gene")
  aggregated <- aggregate(. ~ Gene, data = df, FUN = function(x) mean(as.numeric(x), na.rm = TRUE))
  
  rownames(aggregated) <- aggregated$Gene
  merged_exprs <- as.matrix(aggregated[, sample_cols])
  
  # transpose: row = sample, col = gene
  merged_exprs_t <- t(merged_exprs)
  return(merged_exprs_t)
}

# ------------------------------
# 4. Read Series Matrix
# ------------------------------
exprs_raw <- read_series_matrix_safe(input_file)

# ------------------------------
# 5. Map probe ID → gene & transpose
# ------------------------------
exprs_mapped <- map_probe_to_gene_safe(exprs_raw, gpl_id)

# ------------------------------
# 6. Add Label column (only tumor samples as 1, others as 0)
# ------------------------------
sample_names <- rownames(exprs_mapped)

# Extract sample title information from original file
cat("🏷️ Extracting sample title information...\n")
lines <- readLines(input_file)
sample_title_line <- lines[grepl("^!Sample_title", lines)]

if(length(sample_title_line) > 0) {
  # 解析样本标题
  sample_titles <- unlist(strsplit(sample_title_line, "\t"))
  sample_titles <- sample_titles[-1]  # Remove first element ("!Sample_title")
  sample_titles <- gsub('"', '', sample_titles)  # Remove quotes
  
  # 获取样本ID行
  sample_id_line <- lines[grepl("^!Sample_geo_accession", lines)]
  sample_ids <- unlist(strsplit(sample_id_line, "\t"))
  sample_ids <- sample_ids[-1]  # Remove first element
  sample_ids <- gsub('"', '', sample_ids)  # Remove quotes
  
  # Create ID to title mapping
  id_to_title <- setNames(sample_titles, sample_ids)
  
  # Assign labels based on sample titles: only tumor samples as 1, others as 0
  labels <- sapply(sample_names, function(x) {
    title <- id_to_title[x]
    if(is.na(title)) {
      cat("Warning: Sample title not found for", x, "\n")
      return(0)
    }
    # Exact match: only "tumor sample" starting samples are tumor, exclude "non_tumor sample"
    if(grepl("^tumor sample", title, ignore.case = TRUE)) return(1)    # Only tumor as 1
    else return(0)                                                     # All other types as 0
  })
} else {
  cat("Warning: Sample title information not found, using sample names for matching\n")
  labels <- sapply(sample_names, function(x) {
    if(grepl("^tumor", x, ignore.case = TRUE)) return(1)    # Only tumor starting samples as 1
    else return(0)                                          # All other types as 0
  })
}

# Convert to numeric to avoid becoming object type in Python
labels <- as.numeric(labels)

# Confirm count for each category
cat("Label distribution:\n")
print(table(labels))
cat("Label 0 (Non-tumor): ", sum(labels == 0), "\n")
cat("Label 1 (Tumor): ", sum(labels == 1), "\n")

# Display some sample titles and labels
cat("\nFirst 10 sample titles and labels:\n")
for(i in 1:min(10, length(sample_names))) {
  sample_id <- sample_names[i]
  title <- if(exists("id_to_title")) id_to_title[sample_id] else sample_id
  cat(sprintf("%s: %s -> Label: %d\n", sample_id, title, labels[i]))
}

exprs_mapped <- cbind(exprs_mapped, Label = labels)

# ------------------------------
# 7. Split 8:2 train/validation (maintain class distribution)
# ------------------------------
set.seed(123)
train_idx <- unlist(lapply(unique(labels), function(lbl) {
  idx <- which(labels == lbl)
  sample(idx, size = floor(0.8 * length(idx)))
}))

train_df <- exprs_mapped[train_idx, , drop=FALSE]
valid_df <- exprs_mapped[-train_idx, , drop=FALSE]

# ------------------------------
# 8. Save CSV
# ------------------------------
write.csv(exprs_mapped, output_file, row.names = TRUE)
write.csv(train_df, train_file, row.names = TRUE)
write.csv(valid_df, valid_file, row.names = TRUE)

cat("✅ GSE25097 data conversion and splitting completed.\n")
cat("   - Complete dataset:", output_file, "\n")
cat("   - Training set (80%):", train_file, "\n")
cat("   - Validation set (20%):", valid_file, "\n")
cat("   - Label rule: tumor=1, others=0\n")
