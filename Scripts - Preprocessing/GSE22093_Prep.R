library(affy)
library(hgu133a.db)
library(AnnotationDbi)
library(tidyverse)
library(GEOquery)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

raw_folder <- "GSE22093_RAW"
matrix_gz  <- "GSE22093_series_matrix.txt.gz"

if (!dir.exists(raw_folder)) stop("RAW folder not found: ", raw_folder)
if (!file.exists(matrix_gz)) stop("series_matrix.gz file not found: ", matrix_gz)

lines <- readLines(matrix_gz)

gpl_line <- grep("^!Series_platform_id", lines, value = TRUE)
gpl_id <- if (length(gpl_line) > 0) {
  strsplit(gpl_line, "\t")[[1]][2] |> gsub('"', '', x = _) |> trimws()
} else NA_character_
message("Detected platform GPL: ", gpl_id)

annotation_pkg <- "hgu133a.db"
if (!requireNamespace(annotation_pkg, quietly = TRUE)) {
  BiocManager::install(annotation_pkg)
}
library(annotation_pkg, character.only = TRUE)
message("Using annotation package: ", annotation_pkg)

cel_gz_files <- list.files(raw_folder, pattern = "\\.CEL\\.gz$", full.names = TRUE)
if (length(cel_gz_files) > 0) {
  message("Uncompressing .CEL.gz files...")
  for (f in cel_gz_files) system(paste("gunzip", shQuote(f)))
}

cel_files <- list.files(raw_folder, pattern = "\\.CEL$", ignore.case = TRUE, full.names = TRUE)
if (length(cel_files) == 0) stop("No CEL files found")
message("Found ", length(cel_files), " CEL files.")

message("Reading CEL files and performing RMA normalization...")
raw_data <- ReadAffy(filenames = cel_files)
eset <- rma(raw_data)
expr_probe <- exprs(eset)

probe_ids <- rownames(expr_probe)
gene_symbols <- mapIds(hgu133a.db,
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

feat_er_ihc <- find_feature(c("er positive vs negative by immunohistochemistry", 
                              "er_ihc", "er_status_ihc", "ER IHC"))
feat_er_mrna <- find_feature(c("er positive vs negative by esr1 mrna gene expression (probe 205225_at)",
                               "er_mrna", "ESR1 mRNA", "er by esr1"))
feat_resp  <- find_feature(c("pcr.v.rd", "pcr_rd", "response", "PCR.V.RD", "pCR"))

clin_df <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)

if (!is.null(feat_er_ihc))   clin_df$er_ihc <- extract_from_line(feat_er_ihc)
if (!is.null(feat_er_mrna))  clin_df$er_mrna <- extract_from_line(feat_er_mrna)
if (!is.null(feat_resp))     clin_df$response <- extract_from_line(feat_resp)

if (!is.null(feat_er_ihc))   message("ER IHC unique values: ", paste(unique(clin_df$er_ihc), collapse = ", "))
if (!is.null(feat_er_mrna))  message("ER mRNA unique values: ", paste(unique(clin_df$er_mrna), collapse = ", "))
if (!is.null(feat_resp))     message("Response unique values: ", paste(unique(clin_df$response), collapse = ", "))

if (is.null(feat_resp)) {
  stop("Response feature not found. Please check feat_resp keywords.")
}

clean_er <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("^.*?:\\s*", "", x)
  return(x)
}

if (!is.null(feat_er_ihc) && !is.null(feat_er_mrna)) {
  clin_df$er_ihc_clean <- clean_er(clin_df$er_ihc)
  clin_df$er_mrna_clean <- clean_er(clin_df$er_mrna)
  
  clin_df$er_neg <- clin_df$er_ihc_clean == "ERNEG" & 
    clin_df$er_mrna_clean == "ERNEG" &
    !is.na(clin_df$er_ihc_clean) &
    !is.na(clin_df$er_mrna_clean)
  
  message("ER negative classification (both columns ERneg):")
  print(table(clin_df$er_neg, useNA = "ifany"))
  
  ihc_neg <- clin_df$er_ihc_clean == "ERNEG"
  mrna_neg <- clin_df$er_mrna_clean == "ERNEG"
  discordant <- (ihc_neg & !mrna_neg) | (!ihc_neg & mrna_neg)
  if (any(discordant, na.rm = TRUE)) {
    message("Discordant samples: ", sum(discordant, na.rm = TRUE))
  }
  
} else {
  stop("Complete ER IHC and ER mRNA columns not found")
}

tnbc_idx <- clin_df$er_neg == TRUE & !is.na(clin_df$er_neg)
clin_tnbc <- clin_df[tnbc_idx, ]
message("TNBC samples: ", nrow(clin_tnbc))

if (nrow(clin_tnbc) == 0) {
  message("ER values cross-tabulation:")
  print(table(clin_df$er_ihc_clean, clin_df$er_mrna_clean, useNA = "ifany"))
  stop("No TNBC samples found")
}

resp_raw <- clin_tnbc$response
clin_tnbc$pcr_rd <- NA_character_

unique_resp <- unique(resp_raw[!is.na(resp_raw)])
message("Response unique values: ", paste(unique_resp, collapse = ", "))

if (any(grepl("pCR|PCR", unique_resp, ignore.case = TRUE))) {
  clin_tnbc$pcr_rd[toupper(resp_raw) == "PCR"] <- "pCR"
  clin_tnbc$pcr_rd[toupper(resp_raw) == "RD"] <- "RD"
  other_idx <- !is.na(resp_raw) & !(toupper(resp_raw) %in% c("PCR", "RD"))
  clin_tnbc$pcr_rd[other_idx] <- "RD"
} else if (any(resp_raw %in% c("1", "0"))) {
  clin_tnbc$pcr_rd[resp_raw == "1"] <- "pCR"
  clin_tnbc$pcr_rd[resp_raw == "0"] <- "RD"
} else if (any(grepl("pCR", resp_raw, ignore.case = TRUE))) {
  clin_tnbc$pcr_rd[grepl("pCR", resp_raw, ignore.case = TRUE)] <- "pCR"
  clin_tnbc$pcr_rd[!is.na(resp_raw) & !grepl("pCR", resp_raw, ignore.case = TRUE)] <- "RD"
} else {
  message("Response format not recognized. Attempting manual mapping...")
  pcr_idx <- resp_raw %in% c("pCR", "complete", "responder")
  clin_tnbc$pcr_rd[pcr_idx] <- "pCR"
  clin_tnbc$pcr_rd[!is.na(resp_raw) & !pcr_idx] <- "RD"
}

clin_clean <- clin_tnbc[!is.na(clin_tnbc$pcr_rd), ]
message("TNBC samples with valid labels: ", nrow(clin_clean))
print(table(clin_clean$pcr_rd, useNA = "ifany"))

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

write.csv(final_expr, "GSE22093_TNBC_expression.csv", quote = FALSE)
write.csv(final_clin[, c("sample_id", "pcr_rd")], 
          "GSE22093_TNBC_labels.csv", row.names = FALSE)

clin_full <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
for (fname in names(feat_lines)) {
  clin_full[[fname]] <- extract_from_line(feat_lines[[fname]])
}
write.csv(clin_full, "GSE22093_all_clinical_data.csv", row.names = FALSE)

message("Preprocessing complete. Output files:")
message("  - GSE22093_TNBC_expression.csv")
message("  - GSE22093_TNBC_labels.csv")
message("  - GSE22093_all_clinical_data.csv")