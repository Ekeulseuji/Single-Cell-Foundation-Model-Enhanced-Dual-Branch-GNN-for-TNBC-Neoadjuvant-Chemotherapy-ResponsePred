library(tidyverse)
library(GEOquery)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

tpm_gz  <- "GSE163882_all.data.tpms_222Samples.csv.gz"
matrix_gz <- "GSE163882_series_matrix.txt.gz"

if (!file.exists(tpm_gz)) stop("TPM file not found: ", tpm_gz)
if (!file.exists(matrix_gz)) stop("series_matrix.gz file not found: ", matrix_gz)

message("Reading TPM matrix...")
tpm_data <- read.csv(tpm_gz, row.names = 1, check.names = FALSE)
message("TPM matrix dimensions: ", nrow(tpm_data), " genes x ", ncol(tpm_data), " samples")
message("TPM column examples: ", paste(head(colnames(tpm_data)), collapse = ", "))

lines <- readLines(matrix_gz)

sample_id_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
sample_ids <- strsplit(sample_id_line, "\t")[[1]][-1]
sample_ids <- gsub('^"|"$', '', sample_ids)
message("Clinical samples: ", length(sample_ids))

title_line <- grep("^!Sample_title", lines, value = TRUE)
if (length(title_line) == 0) {
  title_line <- grep("^!Sample_description", lines, value = TRUE)
}
titles <- strsplit(title_line, "\t")[[1]][-1]
titles <- gsub('^"|"$', '', titles)
message("Title examples: ", paste(head(titles), collapse = ", "))

extract_ba <- function(t) {
  m <- regexpr("BA[0-9]+", t)
  if (m > 0) {
    substr(t, m, m + attr(m, "match.length") - 1)
  } else {
    NA_character_
  }
}
ba_ids <- sapply(titles, extract_ba)

mapping <- data.frame(GSM = sample_ids, BA = ba_ids, stringsAsFactors = FALSE)
message("Mapping preview:")
print(head(mapping))

if (any(is.na(mapping$BA))) {
  warning("Some samples failed to extract BA ID. Please check Title format.")
}

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
    if (grepl(":", v)) {
      v <- trimws(sub("^[^:]+:\\s*", "", v))
    }
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

resp_col <- grep("response", colnames(clin_df), ignore.case = TRUE)[1]
er_col <- grep("estrogen_receptor", colnames(clin_df), ignore.case = TRUE)[1]
pr_col <- grep("progesterone_receptor", colnames(clin_df), ignore.case = TRUE)[1]
her2_col <- grep("her2_receptor", colnames(clin_df), ignore.case = TRUE)[1]

message("Response column: ", resp_col, " (", colnames(clin_df)[resp_col], ")")
message("ER column: ", er_col, " (", colnames(clin_df)[er_col], ")")
message("PR column: ", pr_col, " (", colnames(clin_df)[pr_col], ")")
message("HER2 column: ", her2_col, " (", colnames(clin_df)[her2_col], ")")

if (!is.na(resp_col)) {
  message("Response values: ", paste(unique(clin_df[[resp_col]]), collapse = ", "))
}
if (!is.na(er_col)) {
  message("ER values: ", paste(unique(clin_df[[er_col]]), collapse = ", "))
}
if (!is.na(pr_col)) {
  message("PR values: ", paste(unique(clin_df[[pr_col]]), collapse = ", "))
}
if (!is.na(her2_col)) {
  message("HER2 values: ", paste(unique(clin_df[[her2_col]]), collapse = ", "))
}

clean_receptor <- function(x) {
  x <- toupper(trimws(x))
  x[x %in% c("N", "NEG", "NEGATIVE")] <- "NEG"
  x[x %in% c("P", "POS", "POSITIVE")] <- "POS"
  x[grepl("^NEG", x)] <- "NEG"
  x[grepl("^POS", x)] <- "POS"
  return(x)
}

if (!is.na(er_col) && !is.na(pr_col) && !is.na(her2_col)) {
  clin_df$er_clean <- clean_receptor(clin_df[[er_col]])
  clin_df$pr_clean <- clean_receptor(clin_df[[pr_col]])
  clin_df$her2_clean <- clean_receptor(clin_df[[her2_col]])
  
  tnbc_idx <- clin_df$er_clean == "NEG" &
    clin_df$pr_clean == "NEG" &
    clin_df$her2_clean == "NEG" &
    !is.na(clin_df$er_clean) &
    !is.na(clin_df$pr_clean) &
    !is.na(clin_df$her2_clean)
  
  clin_tnbc <- clin_df[tnbc_idx, ]
  message("TNBC samples: ", nrow(clin_tnbc))
  
  if (nrow(clin_tnbc) == 0) {
    warning("No TNBC samples found. Please check receptor status columns.")
  }
} else {
  stop("Complete ER/PR/HER2 feature columns not found")
}

if (!is.na(resp_col)) {
  resp_raw <- clin_tnbc[[resp_col]]
  clin_tnbc$pcr_rd <- NA_character_
  
  unique_resp <- unique(resp_raw[!is.na(resp_raw)])
  message("Response unique values: ", paste(unique_resp, collapse = ", "))
  
  if (any(grepl("pCR|PCR", unique_resp, ignore.case = TRUE))) {
    clin_tnbc$pcr_rd[grepl("pCR|PCR", resp_raw, ignore.case = TRUE)] <- "pCR"
    clin_tnbc$pcr_rd[grepl("RD|R$", resp_raw, ignore.case = TRUE) & !grepl("pCR|PCR", resp_raw, ignore.case = TRUE)] <- "RD"
  } else if (any(unique_resp %in% c("R", "NR"))) {
    clin_tnbc$pcr_rd[resp_raw == "R"] <- "RD"
    clin_tnbc$pcr_rd[resp_raw == "NR"] <- "pCR"
  } else {
    clin_tnbc$pcr_rd[grepl("pCR|complete|responder", resp_raw, ignore.case = TRUE)] <- "pCR"
    clin_tnbc$pcr_rd[grepl("RD|residual|non-responder", resp_raw, ignore.case = TRUE)] <- "RD"
  }
  
  clin_tnbc$pcr_rd[is.na(clin_tnbc$pcr_rd) & !is.na(resp_raw)] <- "RD"
  
  clin_clean <- clin_tnbc[!is.na(clin_tnbc$pcr_rd), ]
  message("TNBC samples with valid labels: ", nrow(clin_clean))
  print(table(clin_clean$pcr_rd, useNA = "ifany"))
  
  if (nrow(clin_clean) == 0) {
    stop("No samples with valid labels")
  }
} else {
  stop("Response column not found")
}

clin_clean$ba_id <- mapping$BA[match(clin_clean$sample_id, mapping$GSM)]

if (any(is.na(clin_clean$ba_id))) {
  warning("Some samples missing BA ID mapping")
}

common_ba <- intersect(colnames(tpm_data), clin_clean$ba_id)
message("Matched BA IDs: ", length(common_ba))

if (length(common_ba) == 0) {
  message("TPM column examples: ", paste(head(colnames(tpm_data)), collapse = ", "))
  message("Clinical BA ID examples: ", paste(head(clin_clean$ba_id), collapse = ", "))
  stop("No overlap between expression columns and BA IDs")
}

final_expr <- tpm_data[, common_ba, drop = FALSE]
final_clin <- clin_clean[match(common_ba, clin_clean$ba_id), ]

colnames(final_expr) <- final_clin$sample_id

stopifnot(all(colnames(final_expr) == final_clin$sample_id))

message("Final expression matrix: ", nrow(final_expr), " genes x ", ncol(final_expr), " samples")
message("Final label distribution:")
print(table(final_clin$pcr_rd))

write.csv(final_expr, "GSE163882_TNBC_expression_TPM.csv", quote = FALSE)
write.csv(final_clin[, c("sample_id", "pcr_rd")], 
          "GSE163882_TNBC_labels.csv", row.names = FALSE)

clin_full <- data.frame(sample_id = sample_ids, GSM = sample_ids, BA = mapping$BA, stringsAsFactors = FALSE)
for (fname in names(feat_lines)) {
  clin_full[[fname]] <- extract_from_line(feat_lines[[fname]])
}
write.csv(clin_full, "GSE163882_all_clinical_data.csv", row.names = FALSE)

message("Preprocessing complete. Output files:")
message("  - GSE163882_TNBC_expression_TPM.csv")
message("  - GSE163882_TNBC_labels.csv")
message("  - GSE163882_all_clinical_data.csv")