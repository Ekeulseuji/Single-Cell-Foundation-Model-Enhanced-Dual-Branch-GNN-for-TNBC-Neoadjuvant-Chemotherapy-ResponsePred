library(tidyverse)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

expr_file <- "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz"
matrix_files <- c("GSE194040-GPL20078_series_matrix.txt.gz",
                  "GSE194040-GPL30493_series_matrix.txt.gz")

if (!file.exists(expr_file)) stop("Expression matrix file not found")
for (f in matrix_files) if (!file.exists(f)) stop("Series matrix file not found: ", f)

message("Reading expression matrix...")
expr_data <- read.delim(gzfile(expr_file), header = TRUE, row.names = 1,
                        check.names = FALSE, stringsAsFactors = FALSE)
message("Matrix dimensions: ", nrow(expr_data), " genes x ", ncol(expr_data), " samples")
colnames(expr_data) <- as.character(colnames(expr_data))

extract_clin_from_series <- function(series_file) {
  message("Reading series matrix: ", series_file)
  lines <- readLines(gzfile(series_file))
  
  sample_id_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
  sample_ids <- strsplit(sample_id_line, "\t")[[1]][-1]
  sample_ids <- gsub('^"|"$', '', sample_ids)
  
  title_line <- grep("^!Sample_title", lines, value = TRUE)
  titles <- strsplit(title_line, "\t")[[1]][-1]
  titles <- gsub('^"|"$', '', titles)
  
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
      return(v)
    }, USE.NAMES = FALSE)
    return(values)
  }
  
  clin_df <- data.frame(sample_id = sample_ids, title = titles, stringsAsFactors = FALSE)
  for (fname in names(feat_lines)) {
    clin_df[[fname]] <- extract_from_line(feat_lines[[fname]])
  }
  colnames(clin_df) <- gsub(" ", "_", colnames(clin_df))
  return(clin_df)
}

clin_list <- lapply(matrix_files, extract_clin_from_series)
clin_df <- bind_rows(clin_list)
message("Total samples after merge: ", nrow(clin_df))

message("Clinical data columns: ", paste(colnames(clin_df), collapse = ", "))

pcr_col <- grep("^pcr$", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
hr_col <- grep("^hr$", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
her2_col <- grep("^her2$", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]
patient_col <- grep("patient_id", colnames(clin_df), ignore.case = TRUE, value = TRUE)[1]

if (is.na(pcr_col) || is.na(hr_col) || is.na(her2_col)) {
  stop("Required columns (pcr, hr, her2) not found")
}

clin_df$pcr_status <- as.numeric(clin_df[[pcr_col]])
clin_df$hr_status <- as.numeric(clin_df[[hr_col]])
clin_df$her2_status <- as.numeric(clin_df[[her2_col]])
clin_df$patient_id <- as.character(clin_df[[patient_col]])

clin_tnbc <- clin_df[clin_df$hr_status == 0 & clin_df$her2_status == 0 &
                       !is.na(clin_df$hr_status) & !is.na(clin_df$her2_status), ]
message("TNBC samples: ", nrow(clin_tnbc))

clin_clean <- clin_tnbc[!is.na(clin_tnbc$pcr_status), ]
message("TNBC samples with pCR labels: ", nrow(clin_clean))

if (nrow(clin_clean) == 0) stop("No usable samples")

clin_clean$pcr_rd <- ifelse(clin_clean$pcr_status == 1, "pCR", "RD")
message("Label distribution:")
print(table(clin_clean$pcr_rd))

common_patients <- intersect(colnames(expr_data), clin_clean$patient_id)
message("Common patients: ", length(common_patients))

if (length(common_patients) == 0) {
  stop("No overlap between expression columns and patient IDs")
}

expr_final <- expr_data[, common_patients, drop = FALSE]
clin_final <- clin_clean[match(common_patients, clin_clean$patient_id), ]

colnames(expr_final) <- clin_final$sample_id

stopifnot(all(colnames(expr_final) == clin_final$sample_id))

message("Final expression matrix: ", nrow(expr_final), " genes x ", ncol(expr_final), " samples")

write.csv(expr_final, "GSE194040_TNBC_expression.csv", quote = FALSE)
write.csv(clin_final[, c("sample_id", "pcr_rd")],
          "GSE194040_TNBC_labels.csv", row.names = FALSE)

message("Preprocessing complete. Output files:")
message("  - GSE194040_TNBC_expression.csv")
message("  - GSE194040_TNBC_labels.csv")