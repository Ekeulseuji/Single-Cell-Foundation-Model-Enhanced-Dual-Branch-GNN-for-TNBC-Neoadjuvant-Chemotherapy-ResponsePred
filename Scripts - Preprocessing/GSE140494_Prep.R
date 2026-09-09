library(affy)
library(hgu133plus2.db)
library(AnnotationDbi)
library(tidyverse)
library(GEOquery)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

raw_folder <- "GSE140494_RAW"
matrix_gz  <- "GSE140494_series_matrix.txt.gz"

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
colnames(expr_gene) <- sub("_.*", "", colnames(expr_gene))
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

feat_er   <- find_feature(c("er_status", "er_status_ihc", "er", "er-ihc"))
feat_pr   <- find_feature(c("pr_status", "pr_status_ihc", "pr", "pr-ihc"))
feat_her2 <- find_feature(c("her2_status", "her2", "her2-ihc"))
feat_resp <- find_feature(c("pathological response", "pathologic_response_pcr_rd", "pcr_rd", "response"))

if (is.null(feat_resp)) stop("'pathological response' feature not found")

clin_df <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
if (!is.null(feat_er))   clin_df$er_status <- extract_from_line(feat_er)
if (!is.null(feat_pr))   clin_df$pr_status <- extract_from_line(feat_pr)
if (!is.null(feat_her2)) clin_df$her2_status <- extract_from_line(feat_her2)
clin_df$response <- extract_from_line(feat_resp)

tnbc_fix_ids <- c("GSM4171496", "GSM4171517")
idx_fix <- clin_df$sample_id %in% tnbc_fix_ids
if (any(idx_fix)) {
  clin_df$er_status[idx_fix]   <- "0"
  clin_df$pr_status[idx_fix]   <- "0"
  clin_df$her2_status[idx_fix] <- "0"
  message("Corrected receptor status for ", sum(idx_fix), " samples to negative")
}

clean_status <- function(x) {
  x <- toupper(trimws(x))
  x[x %in% c("N", "NEG", "NEGATIVE", "0", "NEGATIVE (0-1+)")] <- "NEG"
  x[x %in% c("P", "POS", "POSITIVE", "1+", "2+", "3+")] <- "POS"
  x[x == "I"] <- "IND"
  return(x)
}

if (!is.null(feat_er) && !is.null(feat_pr) && !is.null(feat_her2)) {
  clin_df$er_clean  <- clean_status(clin_df$er_status)
  clin_df$pr_clean  <- clean_status(clin_df$pr_status)
  clin_df$her2_clean <- clean_status(clin_df$her2_status)
  
  tnbc_idx <- clin_df$er_clean == "NEG" &
    clin_df$pr_clean == "NEG" &
    clin_df$her2_clean == "NEG" &
    !is.na(clin_df$er_clean) &
    !is.na(clin_df$pr_clean) &
    !is.na(clin_df$her2_clean)
  
  clin_df <- clin_df[tnbc_idx, ]
  message("TNBC samples after filtering (including corrected samples): ", nrow(clin_df))
} else {
  stop("ER/PR/HER2 features missing")
}

resp_raw <- clin_df$response
clin_df$pcr_rd <- NA_character_
clin_df$pcr_rd[resp_raw == "pCR"] <- "pCR"
clin_df$pcr_rd[!is.na(resp_raw) & resp_raw != "pCR"] <- "RD"

clin_clean <- clin_df[!is.na(clin_df$pcr_rd), ]
message("TNBC samples with valid labels: ", nrow(clin_clean))
print(table(clin_clean$pcr_rd))

if (nrow(clin_clean) == 0) stop("No samples with valid labels")

common_ids <- intersect(colnames(expr_gene), clin_clean$sample_id)
if (length(common_ids) == 0) stop("No common sample IDs between expression matrix and clinical data")

final_expr <- expr_gene[, common_ids, drop = FALSE]
final_clin <- clin_clean[match(common_ids, clin_clean$sample_id), ]
stopifnot(all(colnames(final_expr) == final_clin$sample_id))

message("Final expression matrix: ", nrow(final_expr), " genes x ", ncol(final_expr), " samples")
message("Final label distribution:")
print(table(final_clin$pcr_rd))

write.csv(final_expr, "GSE140494_TNBC_expression.csv", quote = FALSE)
write.csv(final_clin[, c("sample_id", "pcr_rd")], 
          "GSE140494_TNBC_labels.csv", row.names = FALSE)

message("Preprocessing complete. Output files:")
message("  - GSE140494_TNBC_expression.csv")
message("  - GSE140494_TNBC_labels.csv")