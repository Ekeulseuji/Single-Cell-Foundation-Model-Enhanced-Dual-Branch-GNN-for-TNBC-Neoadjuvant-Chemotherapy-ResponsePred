library(affy)
library(hgu133plus2.db)
library(AnnotationDbi)
library(tidyverse)
library(GEOquery)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

raw_folder <- "GSE23988_RAW"
matrix_gz  <- "GSE23988_series_matrix.txt.gz"

if (!dir.exists(raw_folder)) stop("RAW folder not found: ", raw_folder)
if (!file.exists(matrix_gz)) stop("series_matrix.gz file not found: ", matrix_gz)

cel_gz_files <- list.files(raw_folder, pattern = "\\.CEL\\.gz$", full.names = TRUE)
if (length(cel_gz_files) > 0) {
  message("Uncompressing .CEL.gz files...")
  for (f in cel_gz_files) system(paste("gunzip", shQuote(f)))
}

cel_files <- list.files(raw_folder, pattern = "\\.CEL$", ignore.case = TRUE, full.names = TRUE)
if (length(cel_files) == 0) stop("No CEL files found")

message("Reading CEL files and performing RMA normalization...")
raw_data <- ReadAffy(filenames = cel_files)
eset <- rma(raw_data)
expr_probe <- exprs(eset)

probe_ids <- rownames(expr_probe)
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

keep <- !is.na(gene_symbols)
expr_probe <- expr_probe[keep, ]
gene_symbols <- gene_symbols[keep]

expr_gene <- aggregate(expr_probe, by = list(gene = gene_symbols), FUN = max)
rownames(expr_gene) <- expr_gene$gene
expr_gene <- expr_gene[, -1]

colnames(expr_gene) <- sub("\\.CEL$", "", colnames(expr_gene))
colnames(expr_gene) <- basename(colnames(expr_gene))
colnames(expr_gene) <- sub("\\.CEL$", "", colnames(expr_gene))

message("Gene expression matrix: ", nrow(expr_gene), " genes x ", ncol(expr_gene), " samples")

lines <- readLines(matrix_gz)

sample_id_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
sample_ids <- strsplit(sample_id_line, "\t")[[1]][-1]
sample_ids <- gsub('^"|"$', '', sample_ids)

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

find_feature <- function(keywords) {
  for (kw in keywords) {
    exact <- names(feat_lines)[tolower(names(feat_lines)) == tolower(kw)]
    if (length(exact) > 0) return(feat_lines[[exact[1]]])
    partial <- names(feat_lines)[grepl(kw, names(feat_lines), ignore.case = TRUE)]
    if (length(partial) > 0) return(feat_lines[[partial[1]]])
  }
  return(NULL)
}

feat_er   <- find_feature(c("er positive vs negative", "er_status", "er"))
feat_resp <- find_feature(c("pcr.v.rd", "pcr_rd", "response", "pcr"))

if (is.null(feat_er)) {
  message("Available features: ", paste(names(feat_lines), collapse = ", "))
  stop("ER feature not found. Please check feat_er keywords.")
}

if (is.null(feat_resp)) {
  message("Available features: ", paste(names(feat_lines), collapse = ", "))
  stop("Response feature not found. Please check feat_resp keywords.")
}

clin_df <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
clin_df$er_status <- extract_from_line(feat_er)
clin_df$response <- extract_from_line(feat_resp)

message("ER unique values: ", paste(unique(clin_df$er_status), collapse = ", "))
message("Response unique values: ", paste(unique(clin_df$response), collapse = ", "))

tnbc_idx <- clin_df$er_status == "ERneg" & !is.na(clin_df$er_status)
clin_tnbc <- clin_df[tnbc_idx, ]
message("TNBC samples (ERneg): ", nrow(clin_tnbc))

if (nrow(clin_tnbc) == 0) {
  stop("No samples passed TNBC filtering. Please check ER column values.")
}

clin_tnbc$pcr_rd <- NA_character_
clin_tnbc$pcr_rd[clin_tnbc$response == "pCR"] <- "pCR"
clin_tnbc$pcr_rd[clin_tnbc$response == "RD"] <- "RD"

clin_clean <- clin_tnbc[!is.na(clin_tnbc$pcr_rd), ]
message("TNBC samples with valid labels: ", nrow(clin_clean))
print(table(clin_clean$pcr_rd))

if (nrow(clin_clean) == 0) {
  stop("No samples with valid labels. Please check response column values.")
}

common_ids <- intersect(colnames(expr_gene), clin_clean$sample_id)
if (length(common_ids) == 0) {
  message("Expression column examples: ", paste(head(colnames(expr_gene)), collapse = ", "))
  message("Clinical sample ID examples: ", paste(head(clin_clean$sample_id), collapse = ", "))
  stop("No overlap between expression columns and sample IDs")
}

final_expr <- expr_gene[, common_ids, drop = FALSE]
final_clin <- clin_clean[match(common_ids, clin_clean$sample_id), ]
stopifnot(all(colnames(final_expr) == final_clin$sample_id))

message("Final expression matrix: ", nrow(final_expr), " genes x ", ncol(final_expr), " samples")
message("Final label distribution:")
print(table(final_clin$pcr_rd))

write.csv(final_expr, "GSE23988_TNBC_expression.csv", quote = FALSE)
write.csv(final_clin[, c("sample_id", "pcr_rd")], 
          "GSE23988_TNBC_labels.csv", row.names = FALSE)

clin_full <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
for (fname in names(feat_lines)) {
  clin_full[[fname]] <- extract_from_line(feat_lines[[fname]])
}
write.csv(clin_full, "GSE23988_all_clinical_data.csv", row.names = FALSE)

message("Preprocessing complete. Output files:")
message("  - GSE23988_TNBC_expression.csv")
message("  - GSE23988_TNBC_labels.csv")
message("  - GSE23988_all_clinical_data.csv")