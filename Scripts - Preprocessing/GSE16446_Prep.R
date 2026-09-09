library(affy)
library(hgu133plus2.db)
library(AnnotationDbi)
library(tidyverse)
library(GEOquery)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

raw_folder <- "GSE16446_RAW"
matrix_gz  <- "GSE16446_series_matrix.txt.gz"

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
message("Expression column names example: ", paste(head(colnames(expr_gene)), collapse = ", "))

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

message("Available feature names:")
print(names(feat_lines))

feat_her2 <- find_feature(c("her2fishbin", "her2fish", "HER2fishbin"))
feat_pcr  <- find_feature(c("pcr", "response", "pCR"))

if (is.null(feat_her2)) stop("'her2fishbin' feature not found")
if (is.null(feat_pcr))   stop("'pcr' response feature not found")

clin_df <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
clin_df$her2fishbin <- extract_from_line(feat_her2)
clin_df$pcr_raw <- extract_from_line(feat_pcr)

message("her2fishbin unique values: ", paste(unique(clin_df$her2fishbin), collapse = ", "))
message("pcr_raw unique values: ", paste(unique(clin_df$pcr_raw), collapse = ", "))

tnbc_idx <- clin_df$her2fishbin == "0" & !is.na(clin_df$her2fishbin)
clin_tnbc <- clin_df[tnbc_idx, ]
message("TNBC samples (her2fishbin=0): ", nrow(clin_tnbc))

clin_tnbc$pcr_rd <- NA_character_
clin_tnbc$pcr_rd[clin_tnbc$pcr_raw == "1"] <- "pCR"
clin_tnbc$pcr_rd[clin_tnbc$pcr_raw == "0"] <- "RD"

clin_clean <- clin_tnbc[!is.na(clin_tnbc$pcr_rd), ]
message("TNBC samples with valid labels: ", nrow(clin_clean))
print(table(clin_clean$pcr_rd))

if (nrow(clin_clean) == 0) {
  stop("No samples with valid labels. Please check the pcr column values.")
}

common_ids <- intersect(colnames(expr_gene), clin_clean$sample_id)
if (length(common_ids) == 0) {
  message("Expression column name examples: ", paste(head(colnames(expr_gene)), collapse = ", "))
  message("Clinical sample ID examples: ", paste(head(clin_clean$sample_id), collapse = ", "))
  stop("No overlap between expression columns and sample IDs")
}

final_expr <- expr_gene[, common_ids, drop = FALSE]
final_clin <- clin_clean[match(common_ids, clin_clean$sample_id), ]
stopifnot(all(colnames(final_expr) == final_clin$sample_id))

message("Final expression matrix: ", nrow(final_expr), " genes x ", ncol(final_expr), " samples")
message("Final label distribution:")
print(table(final_clin$pcr_rd))

write.csv(final_expr, "GSE16446_TNBC_expression.csv", quote = FALSE)
write.csv(final_clin[, c("sample_id", "pcr_rd")], 
          "GSE16446_TNBC_labels.csv", row.names = FALSE)

message("Preprocessing complete. Output files:")
message("  - GSE16446_TNBC_expression.csv")
message("  - GSE16446_TNBC_labels.csv")