# ==========================================================
# R Script: Process GSE14520 + GSE121248
# Automatically take intersected genes, merge expression matrices,
# clean Labels, map probe IDs to gene symbols,
# and split 80/20 train/validation dataset
# ==========================================================

# ------------------------------
# 1. Install dependencies
# ------------------------------
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  install.packages("GEOquery")
}
library(GEOquery)

# ------------------------------
# 2. Set local file paths
# ------------------------------
file_gse14520_a <- "/Users/yangchengxuan/Downloads/Oncology Dataset/GSE14520 GPL571 Series Matrix.txt"
file_gse14520_b <- "/Users/yangchengxuan/Downloads/Oncology Dataset/GSE14520 GPL3921 Series Matrix.txt"
file_gse121248  <- "/Users/yangchengxuan/Downloads/Oncology Dataset/GSE121248 Series Matrix.txt"

# ------------------------------
# 3. Safe function to read Series Matrix
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
  
  rownames(exprs_raw) <- exprs_raw[,1]  # first column is probe ID
  exprs_raw <- exprs_raw[,-1]           # remove probe ID column
  return(exprs_raw)
}

# ------------------------------
# 4. Extract and clean Labels function
# ------------------------------
extract_labels <- function(file_path, exprs_matrix) {
  lines <- readLines(file_path)
  
  # Try to get from !Sample_title first
  title_line <- grep("^!Sample_title", lines, value = TRUE)
  if (length(title_line) > 0) {
    parts <- unlist(strsplit(title_line, "\t"))
    labels <- parts[-1]  # remove "!Sample_title"
  } else {
    # fallback: use !Sample_characteristics_ch1
    sample_lines <- lines[grepl("^!Sample_characteristics_ch1", lines)]
    labels <- sapply(sample_lines, function(x) sub("!Sample_characteristics_ch1 = ", "", x))
  }
  
  # Align with GSM order
  gsm_ids <- colnames(exprs_matrix)
  names(labels) <- gsm_ids[seq_along(labels)]
  
  # ---- Cleaning rules ----
  clean_labels <- tolower(labels)
  clean_labels[grepl("tumor", clean_labels) & !grepl("non", clean_labels) & !grepl("adjacent", clean_labels)] <- "HBV-HCC"
  clean_labels[grepl("hcc", clean_labels)] <- "HBV-HCC"
  clean_labels[grepl("non", clean_labels) | grepl("adjacent", clean_labels)] <- "Normal"
  clean_labels[grepl("normal", clean_labels) & !grepl("adjacent", clean_labels)] <- "Normal"
  
  return(clean_labels)
}

# ------------------------------
# 5. Function to map probe IDs to gene symbols
# ------------------------------
map_probe_to_gene <- function(exprs_matrix, gpl_id) {
  cat("🔗 Mapping probe IDs to gene symbols for", gpl_id, "...\n")
  
  gpl <- getGEO(gpl_id, destdir = tempdir())
  probe2gene <- Table(gpl)[, c("ID", "Gene Symbol")]
  
  # keep only probes present in expression matrix
  probe2gene <- probe2gene[probe2gene$ID %in% rownames(exprs_matrix), ]
  
  # clean Gene Symbol: take first gene if multiple separated by '///'
  probe2gene$Gene <- sapply(probe2gene$`Gene Symbol`, function(x) {
    if (is.na(x) || x == "") return(NA)
    strsplit(x, "///")[[1]][1]
  })
  
  # remove probes with missing gene
  valid_idx <- !is.na(probe2gene$Gene) & probe2gene$Gene != ""
  probe2gene <- probe2gene[valid_idx, ]
  exprs_matrix <- exprs_matrix[probe2gene$ID, , drop = FALSE]
  
  # add Gene column
  exprs_matrix <- cbind(exprs_matrix, Gene = probe2gene$Gene)
  
  # convert to data.frame to aggregate by Gene
  df <- as.data.frame(exprs_matrix)
  
  # all sample columns
  sample_cols <- setdiff(colnames(df), "Gene")
  
  # aggregate by Gene: take mean of duplicate genes
  aggregated <- aggregate(. ~ Gene, data = df, FUN = function(x) mean(as.numeric(x), na.rm = TRUE))
  
  # convert back to matrix
  rownames(aggregated) <- aggregated$Gene
  merged_exprs <- as.matrix(aggregated[, sample_cols])
  
  return(merged_exprs)
}





# ------------------------------
# 6. Read three datasets
# ------------------------------
cat("📂 Reading GSE14520 GPL571 ...\n")
exprs14520_a <- read_series_matrix_safe(file_gse14520_a)
labels14520_a <- extract_labels(file_gse14520_a, exprs14520_a)
exprs14520_a <- map_probe_to_gene(exprs14520_a, "GPL571")
cat("✅ Completed GSE14520 GPL571\n")

cat("📂 Reading GSE14520 GPL3921 ...\n")
exprs14520_b <- read_series_matrix_safe(file_gse14520_b)
labels14520_b <- extract_labels(file_gse14520_b, exprs14520_b)
exprs14520_b <- map_probe_to_gene(exprs14520_b, "GPL3921")
cat("✅ Completed GSE14520 GPL3921\n")

cat("📂 Reading GSE121248 ...\n")
exprs121248 <- read_series_matrix_safe(file_gse121248)
labels121248 <- extract_labels(file_gse121248, exprs121248)
exprs121248 <- map_probe_to_gene(exprs121248, "GPL570")  # example GPL, replace with real
cat("✅ Completed GSE121248\n")

# ------------------------------
# 7. Take common genes
# ------------------------------
common_genes <- Reduce(intersect, list(rownames(exprs14520_a),
                                       rownames(exprs14520_b),
                                       rownames(exprs121248)))
cat("🔗 Number of common genes: ", length(common_genes), "\n")

exprs14520_a <- exprs14520_a[common_genes, ]
exprs14520_b <- exprs14520_b[common_genes, ]
exprs121248  <- exprs121248[common_genes, ]

# ------------------------------
# 8. Merge expression matrices + Labels
# ------------------------------
merged_exprs <- cbind(exprs14520_a, exprs14520_b, exprs121248)
merged_labels <- c(labels14520_a, labels14520_b, labels121248)

if (length(merged_labels) != ncol(merged_exprs)) {
  stop("❌ Number of labels does not match number of samples!")
}

merged_df <- as.data.frame(t(merged_exprs))
merged_df$Label <- merged_labels

# ------------------------------
# 9. Split train (80%) / validation (20%)
# ------------------------------
set.seed(123)
n <- nrow(merged_df)
train_idx <- sample(seq_len(n), size = 0.8 * n)

train_df <- merged_df[train_idx, ]
valid_df <- merged_df[-train_idx, ]

# ------------------------------
# 10. Save CSV files
# ------------------------------
write.csv(train_df, "train_dataset.csv", row.names = TRUE)
write.csv(valid_df, "valid_dataset.csv", row.names = TRUE)
write.csv(merged_df, "merged_features_labels.csv", row.names = TRUE)

cat("✅ Data preparation completed! Files generated:\n")
cat("   - train_dataset.csv (80%)\n")
cat("   - valid_dataset.csv (20%)\n")
cat("   - merged_features_labels.csv (full)\n")
