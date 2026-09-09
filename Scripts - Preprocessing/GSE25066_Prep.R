library(affy)
library(hgu133a.db)
library(AnnotationDbi)
library(tidyverse)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

tar_file <- "GSE25066_RAW.tar"
cel_dir  <- "GSE25066_RAW"

if (!dir.exists(cel_dir)) {
  message("Extracting CEL files...")
  untar(tar_file, exdir = cel_dir)
}

cel_gz_files <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE)
if (length(cel_gz_files) > 0) {
  for (f in cel_gz_files) system(paste("gunzip", f))
}

cel_files <- list.files(cel_dir, pattern = "\\.CEL$", ignore.case = TRUE, full.names = TRUE)
if (length(cel_files) == 0) stop("No CEL files found")

raw_data <- ReadAffy(filenames = cel_files)
eset     <- rma(raw_data)
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

colnames(expr_gene) <- sub("_.*", "", colnames(expr_gene))
message("Expression matrix: ", nrow(expr_gene), " genes × ", ncol(expr_gene), " samples")

# -----------------------------------------------------------------------------
# Parse series_matrix.txt
# -----------------------------------------------------------------------------
if (!file.exists("GSE25066_series_matrix.txt")) {
  system("gunzip GSE25066_series_matrix.txt.gz")
}

lines <- readLines("GSE25066_series_matrix.txt")

sample_id_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
sample_ids <- strsplit(sample_id_line, "\t")[[1]][-1]
sample_ids <- gsub('^"|"$', '', sample_ids)

char_lines_idx <- grep("^!Sample_characteristics_ch", lines)
feat_lines <- list()
for (idx in char_lines_idx) {
  line <- lines[idx]
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) < 2) next
  first_val <- gsub('^"|"$', '', parts[2])
  first_val <- trimws(first_val)
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
    v <- gsub('^"|"$', '', v)
    v <- trimws(v)
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

feat_er   <- find_feature(c("er_status_ihc", "er_status"))
feat_pr   <- find_feature(c("pr_status_ihc", "pr_status"))
feat_her2 <- find_feature(c("her2_status"))
feat_resp <- find_feature(c("pathologic_response_pcr_rd", "pcr_rd", "response"))

if (is.null(feat_er))   stop("ER feature not found")
if (is.null(feat_pr))   stop("PR feature not found")
if (is.null(feat_her2)) stop("HER2 feature not found")
if (is.null(feat_resp)) stop("Response feature not found")

clin_df <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
clin_df$er_status_ihc <- extract_from_line(feat_er)
clin_df$pr_status_ihc <- extract_from_line(feat_pr)
clin_df$her2_status   <- extract_from_line(feat_her2)
clin_df$pathologic_response_pcr_rd <- extract_from_line(feat_resp)

# -----------------------------------------------------------------------------
# Select TNBC (ER-, PR-, HER2-)
# -----------------------------------------------------------------------------
clean_status <- function(x) {
  x <- toupper(trimws(x))
  x[x %in% c("N", "NEG", "NEGATIVE", "0", "NEGATIVE (0-1+)")] <- "NEG"
  x[x %in% c("P", "POS", "POSITIVE", "1+", "2+", "3+")] <- "POS"
  x[x == "I"] <- "IND"
  return(x)
}

clin_df$er_clean  <- clean_status(clin_df$er_status_ihc)
clin_df$pr_clean  <- clean_status(clin_df$pr_status_ihc)
clin_df$her2_clean <- clean_status(clin_df$her2_status)

tnbc_idx <- clin_df$er_clean == "NEG" &
  clin_df$pr_clean == "NEG" &
  clin_df$her2_clean == "NEG" &
  !is.na(clin_df$er_clean) &
  !is.na(clin_df$pr_clean) &
  !is.na(clin_df$her2_clean)

tnbc_clin <- clin_df[tnbc_idx, ]
message("Initial TNBC samples: ", nrow(tnbc_clin))

# -----------------------------------------------------------------------------
# Manual update of pathologic_response_pcr_rd for specific samples
# -----------------------------------------------------------------------------
update_map <- data.frame(
  sample_id = c(
    "GSM615805", "GSM615637", "GSM615639", "GSM615640", "GSM615644", "GSM615649",
    "GSM615650", "GSM615651", "GSM615657", "GSM615658", "GSM615660",
    "GSM615661", "GSM615666", "GSM615667", "GSM615668", "GSM615671",
    "GSM615672", "GSM615674", "GSM615677", "GSM615680", "GSM615681",
    "GSM615687", "GSM615689", "GSM615757", "GSM615762", "GSM615763",
    "GSM615764", "GSM615766", "GSM615769", "GSM615118", "GSM615119",
    "GSM615696", "GSM615805", "GSM615707", "GSM615714"
  ),
  new_label = c(
    "RD", "pCR", "pCR", "RD", "pCR", "RD", 
    "RD",  "RD",  "RD",  "pCR", "pCR",
    "RD",  "RD",  "RD",  "pCR", "RD", 
    "pCR", "RD",  "pCR", "pCR", "pCR",
    "RD",  "RD",  NA,    NA,    NA,
    NA,    NA,    NA,    NA,    NA,
    "RD", "RD", "RD", "RD"
  ),
  stringsAsFactors = FALSE
)

idx <- match(update_map$sample_id, tnbc_clin$sample_id)
if (any(is.na(idx))) {
  warning("Some sample IDs not found; removing them.")
  update_map <- update_map[!is.na(idx), ]
  idx <- idx[!is.na(idx)]
}
tnbc_clin$pathologic_response_pcr_rd[idx] <- update_map$new_label

# -----------------------------------------------------------------------------
# Generate pcr_rd column (combine pathologic_response_pcr_rd and RCB text)
# -----------------------------------------------------------------------------
tnbc_clin$pcr_rd <- NA_character_

# Direct match from pathologic_response_pcr_rd
tnbc_clin$pcr_rd[tnbc_clin$pathologic_response_pcr_rd == "pCR"] <- "pCR"
tnbc_clin$pcr_rd[tnbc_clin$pathologic_response_pcr_rd == "RD"]  <- "RD"

# Map from RCB text within the same column
idx_rcb_text <- grepl("RCB", tnbc_clin$pathologic_response_pcr_rd)
tnbc_clin$pcr_rd[idx_rcb_text & grepl("RCB-0/I", tnbc_clin$pathologic_response_pcr_rd)] <- "pCR"
tnbc_clin$pcr_rd[idx_rcb_text & grepl("RCB-II|RCB-III", tnbc_clin$pathologic_response_pcr_rd)] <- "RD"

# Map from separate rcb_class column (if present)
if ("pathologic_response_rcb_class" %in% colnames(tnbc_clin)) {
  idx_miss <- is.na(tnbc_clin$pcr_rd) & !is.na(tnbc_clin$pathologic_response_rcb_class)
  if (any(idx_miss)) {
    rcb2 <- tnbc_clin$pathologic_response_rcb_class[idx_miss]
    rcb2_upper <- toupper(rcb2)
    is_pcr_rcb <- grepl("RCB[-\\s]*0[/\\-]?I?", rcb2_upper) | grepl("RCB[-\\s]*0", rcb2_upper)
    is_rd_rcb  <- grepl("RCB[-\\s]*II|RCB[-\\s]*III", rcb2_upper)
    tnbc_clin$pcr_rd[idx_miss][is_pcr_rcb] <- "pCR"
    tnbc_clin$pcr_rd[idx_miss][is_rd_rcb]  <- "RD"
  }
}

tnbc_clean <- tnbc_clin[!is.na(tnbc_clin$pcr_rd), ]
message("Final TNBC samples with valid labels: ", nrow(tnbc_clean))
print(table(tnbc_clean$pcr_rd))

# -----------------------------------------------------------------------------
# Align expression matrix with labels
# -----------------------------------------------------------------------------
common_ids <- intersect(colnames(expr_gene), tnbc_clean$sample_id)
if (length(common_ids) == 0) stop("No common sample IDs between expression and clinical data")
tnbc_expr <- expr_gene[, common_ids, drop = FALSE]
tnbc_clean <- tnbc_clean[match(common_ids, tnbc_clean$sample_id), ]
stopifnot(all(colnames(tnbc_expr) == tnbc_clean$sample_id))

# -----------------------------------------------------------------------------
# Save output
# -----------------------------------------------------------------------------
write.csv(tnbc_expr, "GSE25066_TNBC_expression.csv", quote = FALSE)
write.csv(tnbc_clean[, c("sample_id", "pcr_rd")], "GSE25066_TNBC_labels.csv", row.names = FALSE)

message("Done. Output files:")
message("  - GSE25066_TNBC_expression.csv")
message("  - GSE25066_TNBC_labels.csv")