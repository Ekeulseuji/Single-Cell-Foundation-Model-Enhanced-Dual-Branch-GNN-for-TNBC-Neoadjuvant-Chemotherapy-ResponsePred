library(tidyverse)
library(GEOquery)
library(org.Hs.eg.db)
library(AnnotationDbi)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

tpm_gz <- "GSE123845_exp_tpm_matrix.csv.gz"
matrix_gz <- "GSE123845_series_matrix.txt.gz"

if (!file.exists(tpm_gz)) stop("TPM file not found: ", tpm_gz)
if (!file.exists(matrix_gz)) stop("series_matrix.gz file not found: ", matrix_gz)

message("Reading TPM matrix...")
tpm_data <- read.csv(tpm_gz, row.names = 1, check.names = FALSE)

message("Starting gene symbol standardization...")
valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
genes <- rownames(tpm_data)
keep <- genes %in% valid_symbols
message("Original genes: ", length(genes), ", matched to official symbols: ", sum(keep))
if (sum(keep) == 0) stop("No matching symbols found")
tpm_data <- tpm_data[keep, , drop = FALSE]
if (any(duplicated(rownames(tpm_data)))) {
  message("Duplicate symbols detected, aggregating by mean.")
  tpm_data <- aggregate(tpm_data, by = list(Symbol = rownames(tpm_data)), FUN = mean)
  rownames(tpm_data) <- tpm_data$Symbol
  tpm_data$Symbol <- NULL
}
message("Standardized matrix dimensions: ", nrow(tpm_data), " genes x ", ncol(tpm_data), " samples")

lines <- readLines(matrix_gz)
sample_id_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
sample_ids <- strsplit(sample_id_line, "\t")[[1]][-1]
sample_ids <- gsub('^"|"$', '', sample_ids)

title_line <- grep("^!Sample_title", lines, value = TRUE)
titles <- strsplit(title_line, "\t")[[1]][-1]
titles <- gsub('^"|"$', '', titles)
gsm_to_title <- setNames(titles, sample_ids)

char_lines_idx <- grep("^!Sample_characteristics_ch", lines)
feat_lines <- list()
for (idx in char_lines_idx) {
  line <- lines[idx]
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) < 2) next
  first_val <- gsub('^"|"$', '', parts[2]) |> trimws()
  if (!grepl(":", first_val)) {
    feat_name <- paste0("UNKNOWN_", idx)
  } else {
    feat_name <- trimws(strsplit(first_val, ":")[[1]][1])
  }
  feat_lines[[feat_name]] <- line
}

extract_from_line <- function(line) {
  parts <- strsplit(line, "\t")[[1]]
  values <- parts[-1]
  values <- sapply(values, function(v) {
    v <- gsub('^"|"$', '', v) |> trimws()
    if (grepl(":", v)) v <- trimws(sub("^[^:]+:\\s*", "", v))
    if (v %in% c("NA", "", "null", "NULL", "na")) return(NA_character_)
    return(v)
  }, USE.NAMES = FALSE)
  return(values)
}

clin_df <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
for (fname in names(feat_lines)) {
  clin_df[[fname]] <- extract_from_line(feat_lines[[fname]])
}
colnames(clin_df) <- gsub(" ", "_", colnames(clin_df))
colnames(clin_df) <- gsub("\\(|\\)", "", colnames(clin_df))

subtype_col <- grep("subtype_consensus", colnames(clin_df), ignore.case = TRUE)[1]
time_col <- grep("^time$", colnames(clin_df), ignore.case = TRUE)[1]
if (is.na(time_col)) time_col <- grep("timepoint", colnames(clin_df), ignore.case = TRUE)[1]
pcr_col <- grep("pcr_status", colnames(clin_df), ignore.case = TRUE)[1]
if (is.na(subtype_col) || is.na(time_col) || is.na(pcr_col)) 
  stop("Critical clinical columns not found")

clin_df$is_tnbc <- clin_df[[subtype_col]] == "TN"
clin_df$is_t1 <- clin_df[[time_col]] == "T1" | clin_df[[time_col]] == "1"
clin_tnbc_t1 <- clin_df[clin_df$is_tnbc & clin_df$is_t1 & 
                          !is.na(clin_df$is_tnbc) & !is.na(clin_df$is_t1), ]
if (nrow(clin_tnbc_t1) == 0) {
  warning("No TNBC+T1 samples found, using all TNBC samples.")
  clin_tnbc_t1 <- clin_df[clin_df$is_tnbc, ]
}

pcr_raw <- clin_tnbc_t1[[pcr_col]]
clin_tnbc_t1$pcr_rd <- NA_character_
clin_tnbc_t1$pcr_rd[pcr_raw == "1"] <- "pCR"
clin_tnbc_t1$pcr_rd[pcr_raw == "0"] <- "RD"
clin_tnbc_t1$pcr_rd[is.na(clin_tnbc_t1$pcr_rd) & !is.na(pcr_raw)] <- "RD"
clin_clean <- clin_tnbc_t1[!is.na(clin_tnbc_t1$pcr_rd), ]
if (nrow(clin_clean) == 0) stop("No samples with valid labels")

clin_clean$title <- gsm_to_title[clin_clean$sample_id]
common_titles <- intersect(colnames(tpm_data), clin_clean$title)
if (length(common_titles) == 0) stop("No overlap between expression columns and titles")
clin_clean <- clin_clean[match(common_titles, clin_clean$title), ]
final_expr <- tpm_data[, common_titles, drop = FALSE]
colnames(final_expr) <- clin_clean$sample_id
stopifnot(all(colnames(final_expr) == clin_clean$sample_id))

write.csv(final_expr, "GSE123845_TNBC_expression.csv", quote = FALSE)
write.csv(clin_clean[, c("sample_id", "pcr_rd")], 
          "GSE123845_TNBC_labels.csv", row.names = FALSE)

clin_full <- data.frame(sample_id = sample_ids, 
                        title = gsm_to_title[sample_ids], 
                        stringsAsFactors = FALSE)
for (fname in names(feat_lines)) {
  clin_full[[fname]] <- extract_from_line(feat_lines[[fname]])
}
write.csv(clin_full, "GSE123845_all_clinical_data.csv", row.names = FALSE)

message("Preprocessing complete. Output files:")
message("  - GSE123845_TNBC_expression.csv")
message("  - GSE123845_TNBC_labels.csv")
message("  - GSE123845_all_clinical_data.csv")