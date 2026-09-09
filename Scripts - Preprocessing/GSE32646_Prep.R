library(affy)
library(hgu133plus2.db)
library(AnnotationDbi)
library(tidyverse)
library(GEOquery)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

raw_folder <- "GSE32646_RAW"
matrix_gz  <- "GSE32646_series_matrix.txt.gz"

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
  "GPL570" = "hgu133plus2.db",
  "GPL96"  = "hgu133a.db",
  if (grepl("570", gpl_id, ignore.case = TRUE)) "hgu133plus2.db" else "hgu133plus2.db"
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

feat_er    <- find_feature(c("er_status", "er", "ER", "ER STATUS", "er_ihc", "ER IHC", "er by ihc"))
feat_pr    <- find_feature(c("pr_status", "pr", "PR", "PR STATUS", "pr_ihc", "PR IHC", "pr by ihc"))
feat_her2  <- find_feature(c("her2_status", "her2", "HER2", "HER2 STATUS", "her2_ihc", "HER2 IHC", "her2 by ihc"))
feat_resp  <- find_feature(c("pcr", "response", "pCR", "RD", "pathologic_response", "pathological_response", "pCR_status"))

clin_df <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
if (!is.null(feat_er))    clin_df$er_status <- extract_from_line(feat_er)
if (!is.null(feat_pr))    clin_df$pr_status <- extract_from_line(feat_pr)
if (!is.null(feat_her2))  clin_df$her2_status <- extract_from_line(feat_her2)
if (!is.null(feat_resp))  clin_df$response <- extract_from_line(feat_resp)

if (!is.null(feat_er))    message("ER unique values: ", paste(unique(clin_df$er_status), collapse = ", "))
if (!is.null(feat_pr))    message("PR unique values: ", paste(unique(clin_df$pr_status), collapse = ", "))
if (!is.null(feat_her2))  message("HER2 unique values: ", paste(unique(clin_df$her2_status), collapse = ", "))
if (!is.null(feat_resp))  message("Response unique values: ", paste(unique(clin_df$response), collapse = ", "))

if (is.null(feat_resp)) {
  message("Response column not found. Searching for columns containing pCR/RD...")
  for (fname in names(feat_lines)) {
    vals <- extract_from_line(feat_lines[[fname]])
    if (any(grepl("pCR|RD|PCR", vals, ignore.case = TRUE), na.rm = TRUE)) {
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

clean_status <- function(x) {
  x <- toupper(trimws(x))
  x[x == "N"] <- "NEG"
  x[x == "P"] <- "POS"
  x[x %in% c("NEG", "NEGATIVE", "0", "NEGATIVE (0-1+)", "NEGATIVE 0")] <- "NEG"
  x[x %in% c("POS", "POSITIVE", "1", "1+", "2+", "3+", "POSITIVE 3+")] <- "POS"
  x[grepl("^NEG", x)] <- "NEG"
  x[grepl("^POS", x)] <- "POS"
  x[tolower(x) == "positive"] <- "POS"
  x[tolower(x) == "negative"] <- "NEG"
  x[grepl("ER?NEG", x)] <- "NEG"
  return(x)
}

if (!is.null(feat_er) && !is.null(feat_pr) && !is.null(feat_her2)) {
  clin_df$er_clean  <- clean_status(clin_df$er_status)
  clin_df$pr_clean  <- clean_status(clin_df$pr_status)
  clin_df$her2_clean <- clean_status(clin_df$her2_status)
  
  message("Cleaned status distributions:")
  print(table(clin_df$er_clean, useNA = "ifany"))
  print(table(clin_df$pr_clean, useNA = "ifany"))
  print(table(clin_df$her2_clean, useNA = "ifany"))
  
  tnbc_idx <- clin_df$er_clean == "NEG" &
    clin_df$pr_clean == "NEG" &
    clin_df$her2_clean == "NEG" &
    !is.na(clin_df$er_clean) &
    !is.na(clin_df$pr_clean) &
    !is.na(clin_df$her2_clean)
  
  clin_tnbc <- clin_df[tnbc_idx, ]
  message("TNBC samples: ", nrow(clin_tnbc))
  
  if (nrow(clin_tnbc) == 0) {
    if (!is.null(feat_er)) message("ER original values: ", paste(unique(clin_df$er_status), collapse = ", "))
    if (!is.null(feat_pr)) message("PR original values: ", paste(unique(clin_df$pr_status), collapse = ", "))
    if (!is.null(feat_her2)) message("HER2 original values: ", paste(unique(clin_df$her2_status), collapse = ", "))
    message("Using all samples.")
    clin_tnbc <- clin_df
  }
} else {
  message("Complete ER/PR/HER2 features not detected. Using all samples.")
  clin_tnbc <- clin_df
}

if (!is.null(feat_resp)) {
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

message("Final expression matrix: ", nrow(final_expr), " genes x ", ncol(final_expr), " samples")
message("Final label distribution:")
print(table(final_clin$pcr_rd))

write.csv(final_expr, "GSE32646_TNBC_expression.csv", quote = FALSE)
write.csv(final_clin[, c("sample_id", "pcr_rd")], 
          "GSE32646_TNBC_labels.csv", row.names = FALSE)

clin_full <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
for (fname in names(feat_lines)) {
  clin_full[[fname]] <- extract_from_line(feat_lines[[fname]])
}
write.csv(clin_full, "GSE32646_all_clinical_data.csv", row.names = FALSE)

message("Preprocessing complete. Output files:")
message("  - GSE32646_TNBC_expression.csv")
message("  - GSE32646_TNBC_labels.csv")
message("  - GSE32646_all_clinical_data.csv")