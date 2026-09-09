library(affy)
library(hgu133a.db)
library(AnnotationDbi)
library(tidyverse)
library(GEOquery)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

raw_folder <- "GSE20194_RAW"
matrix_gz  <- "GSE20194_series_matrix.txt.gz"

if (!dir.exists(raw_folder)) stop("RAW folder not found: ", raw_folder)
if (!file.exists(matrix_gz)) stop("series_matrix.gz file not found: ", matrix_gz)

lines <- readLines(matrix_gz)

gpl_line <- grep("^!Series_platform_id", lines, value = TRUE)
gpl_id <- if (length(gpl_line) > 0) {
  strsplit(gpl_line, "\t")[[1]][2] |> gsub('"', '', x = _) |> trimws()
} else NA_character_
message("Detected platform GPL: ", gpl_id)

annotation_pkg <- switch(
  gpl_id,
  "GPL96"  = "hgu133a.db",
  "GPL570" = "hgu133plus2.db",
  if (grepl("96", gpl_id, ignore.case = TRUE)) "hgu133a.db" else "hgu133plus2.db"
)

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
gene_symbols <- mapIds(get(annotation_pkg),
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
colnames(expr_gene) <- sub("\\.CEL$", "", colnames(expr_gene))

message("Gene expression matrix: ", nrow(expr_gene), " genes x ", ncol(expr_gene), " samples")
message("Expression column examples: ", paste(head(colnames(expr_gene)), collapse = ", "))

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

feat_resp <- find_feature(c("pcr", "response", "pCR", "RD", "pathologic_response", "prc or rd"))

clin_df <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
if (!is.null(feat_resp)) {
  clin_df$response <- extract_from_line(feat_resp)
} else {
  message("Response column not found. Searching for columns containing pCR/RD...")
  for (fname in names(feat_lines)) {
    vals <- extract_from_line(feat_lines[[fname]])
    if (any(grepl("pCR|RD", vals, ignore.case = TRUE), na.rm = TRUE)) {
      message("Found potential response column: ", fname)
      feat_resp <- fname
      clin_df$response <- extract_from_line(feat_lines[[fname]])
      break
    }
  }
  
  if (is.null(feat_resp)) {
    stop("Response column not found")
  }
}

if (!is.null(feat_resp)) {
  message("Response unique values: ", paste(unique(clin_df$response), collapse = ", "))
}

tnbc_sample_ids <- c(
  "GSM505331", "GSM505333", "GSM505344", "GSM505345", "GSM505350", "GSM505351",
  "GSM505357", "GSM505359", "GSM505360", "GSM505370", "GSM505372", "GSM505411",
  "GSM505413", "GSM505416", "GSM505421", "GSM505423", "GSM505428", "GSM505432",
  "GSM505438", "GSM505440", "GSM505450", "GSM505454", "GSM505459", "GSM505462",
  "GSM505466", "GSM505467", "GSM505470", "GSM505475", "GSM505476", "GSM505478",
  "GSM505487", "GSM505489", "GSM505490", "GSM505491", "GSM505496", "GSM505497",
  "GSM505498", "GSM505499", "GSM505500", "GSM505503", "GSM505514", "GSM505516",
  "GSM505527", "GSM505528", "GSM505533", "GSM505535", "GSM505536", "GSM505537",
  "GSM505540", "GSM505541", "GSM505543", "GSM505545", "GSM505552", "GSM505553",
  "GSM505554", "GSM505555", "GSM505558", "GSM505562", "GSM505563", "GSM505564",
  "GSM505566", "GSM505570", "GSM505572", "GSM505573", "GSM505574", "GSM505575",
  "GSM505580", "GSM505583", "GSM505584", "GSM505588", "GSM505592"
)

message("Manually defined TNBC samples: ", length(tnbc_sample_ids))

tnbc_ids_available <- intersect(tnbc_sample_ids, colnames(expr_gene))
message("TNBC samples available in expression matrix: ", length(tnbc_ids_available))

clin_tnbc <- clin_df[clin_df$sample_id %in% tnbc_ids_available, ]
message("TNBC samples in clinical data: ", nrow(clin_tnbc))

if (nrow(clin_tnbc) == 0) {
  stop("No TNBC samples found in clinical data")
}

if (!is.null(feat_resp)) {
  resp_raw <- clin_tnbc$response
  clin_tnbc$pcr_rd <- NA_character_
  
  unique_resp <- unique(resp_raw[!is.na(resp_raw)])
  message("Response unique values in TNBC samples: ", paste(unique_resp, collapse = ", "))
  
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
    pcr_idx <- resp_raw %in% c("pCR", "complete", "responder", "1")
    clin_tnbc$pcr_rd[pcr_idx] <- "pCR"
    clin_tnbc$pcr_rd[!is.na(resp_raw) & !pcr_idx] <- "RD"
  }
  
  clin_clean <- clin_tnbc[!is.na(clin_tnbc$pcr_rd), ]
  message("TNBC samples with valid labels: ", nrow(clin_clean))
  print(table(clin_clean$pcr_rd, useNA = "ifany"))
  
  if (nrow(clin_clean) == 0) {
    stop("No samples with valid labels")
  }
} else {
  stop("Response column not found")
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

if ("GSM505496" %in% final_clin$sample_id) {
  final_clin$pcr_rd[final_clin$sample_id == "GSM505496"] <- "pCR"
  message("Manual correction: GSM505496 label set to pCR")
} else {
  message("GSM505496 not in final data")
}

message("Final expression matrix: ", nrow(final_expr), " genes x ", ncol(final_expr), " samples")
message("Final label distribution (after correction):")
print(table(final_clin$pcr_rd))

write.csv(final_expr, "GSE20194_TNBC_expression.csv", quote = FALSE)
write.csv(final_clin[, c("sample_id", "pcr_rd")], 
          "GSE20194_TNBC_labels.csv", row.names = FALSE)

clin_full <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
for (fname in names(feat_lines)) {
  clin_full[[fname]] <- extract_from_line(feat_lines[[fname]])
}
write.csv(clin_full, "GSE20194_all_clinical_data.csv", row.names = FALSE)

message("Preprocessing complete. Output files:")
message("  - GSE20194_TNBC_expression.csv")
message("  - GSE20194_TNBC_labels.csv")
message("  - GSE20194_all_clinical_data.csv")