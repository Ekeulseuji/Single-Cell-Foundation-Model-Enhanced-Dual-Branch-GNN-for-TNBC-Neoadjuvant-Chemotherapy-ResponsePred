library(limma)
library(GEOquery)
library(tidyverse)

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

raw_folder <- "GSE32603_RAW"
matrix_gz  <- "GSE32603_series_matrix.txt.gz"

gpr_gz_files <- list.files(raw_folder, pattern = "\\.gpr\\.gz$", full.names = TRUE)
if (length(gpr_gz_files) > 0) {
  message("Uncompressing .gpr.gz files...")
  for (f in gpr_gz_files) system(paste("gunzip", shQuote(f)))
}

gpr_files <- list.files(raw_folder, pattern = "\\.gpr$", ignore.case = TRUE, full.names = TRUE)
if (length(gpr_files) == 0) stop("No .gpr files found")
message("Found ", length(gpr_files), " GPR files.")

col_names <- c(
  "Block", "Column", "Row", "Name", "ID", "X", "Y", "Dia.",
  "F635 Median", "F635 Mean", "F635 SD", "B635 Median", "B635 Mean", "B635 SD",
  "% > B635+1SD", "% > B635+2SD", "F635 % Sat.",
  "F532 Median", "F532 Mean", "F532 SD", "B532 Median", "B532 Mean", "B532 SD",
  "% > B532+1SD", "% > B532+2SD", "F532 % Sat.",
  "Ratio of Medians", "Ratio of Means", "Median of Ratios", "Mean of Ratios",
  "Ratios SD", "Rgn Ratio", "Rgn R�", "F Pixels", "B Pixels",
  "Sum of Medians", "Sum of Means", "Log Ratio",
  "F635 Median - B635", "F532 Median - B532", "F635 Mean - B635", "F532 Mean - B532",
  "Flags"
)

test_data <- read.table(gpr_files[1], skip = 27, header = FALSE, sep = "\t",
                        quote = "", fill = TRUE, stringsAsFactors = FALSE,
                        col.names = col_names, check.names = FALSE)

all_probe_ids <- test_data[["ID"]]
num_probes <- length(all_probe_ids)
num_files <- length(gpr_files)

expr_matrix <- matrix(NA, nrow = num_probes, ncol = num_files)
sample_names <- character(num_files)

for (i in seq_along(gpr_files)) {
  f <- gpr_files[i]
  message("Reading ", basename(f), " (", i, "/", num_files, ")")
  dt <- read.table(f, skip = 27, header = FALSE, sep = "\t",
                   quote = "", fill = TRUE, stringsAsFactors = FALSE,
                   col.names = col_names, check.names = FALSE)
  
  log_ratio <- as.numeric(dt[["Log Ratio"]])
  
  probe_ids_file <- dt[["ID"]]
  common_probes <- intersect(all_probe_ids, probe_ids_file)
  idx_base <- match(common_probes, all_probe_ids)
  idx_file <- match(common_probes, probe_ids_file)
  
  expr_matrix[idx_base, i] <- log_ratio[idx_file]
  sample_names[i] <- gsub("\\.gpr$", "", basename(f))
}

keep_rows <- apply(expr_matrix, 1, function(x) !all(is.na(x)))
expr_matrix <- expr_matrix[keep_rows, ]
probe_ids <- all_probe_ids[keep_rows]

rownames(expr_matrix) <- probe_ids
colnames(expr_matrix) <- sample_names

message("Probe-level expression matrix: ", nrow(expr_matrix), " probes x ", ncol(expr_matrix), " samples")

message("\nRetrieving platform annotation...")
gpl <- getGEO("GPL14668", destdir = ".")
tbl <- Table(gpl)

map_df <- tbl[, c("ID", "GENE SYMBOL")]
colnames(map_df) <- c("probe_id", "gene_symbol")
map_df$probe_id <- trimws(as.character(map_df$probe_id))
map_df$gene_symbol <- trimws(as.character(map_df$gene_symbol))
map_df$gene_symbol <- gsub('^"|"$', '', map_df$gene_symbol)
map_df$gene_symbol <- sapply(strsplit(map_df$gene_symbol, split = "[/,;]"), `[`, 1)
map_df$gene_symbol <- trimws(map_df$gene_symbol)
map_df <- map_df[!is.na(map_df$gene_symbol) & map_df$gene_symbol != "", ]
map_df <- map_df[!duplicated(map_df$probe_id), ]
rownames(map_df) <- map_df$probe_id
message("Obtained ", nrow(map_df), " probe-gene symbol mappings.")

expr_probes <- rownames(expr_matrix)
expr_probes_clean <- gsub('"', '', expr_probes)
expr_probes_clean <- gsub("'", '', expr_probes_clean)
expr_probes_clean <- trimws(expr_probes_clean)
valid_idx <- expr_probes_clean != "" & !is.na(expr_probes_clean)
expr_probes_clean <- expr_probes_clean[valid_idx]
expr_matrix_clean <- expr_matrix[valid_idx, , drop = FALSE]
rownames(expr_matrix_clean) <- expr_probes_clean

common <- intersect(expr_probes_clean, map_df$probe_id)
message("Matched probes: ", length(common))
if (length(common) == 0) stop("No probes matched to gene symbols. Please check platform annotation.")

expr_mapped <- expr_matrix_clean[common, , drop = FALSE]
gene_symbols <- map_df[common, "gene_symbol"]
rownames(expr_mapped) <- gene_symbols
message("Mapped expression matrix: ", nrow(expr_mapped), " probes x ", ncol(expr_mapped), " samples")

message("\nAggregating duplicate gene symbols using data.table...")

expr_dt <- as.data.table(expr_mapped, keep.rownames = "gene_symbol")
sample_cols <- setdiff(names(expr_dt), "gene_symbol")
expr_dt[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]
expr_agg <- expr_dt[, lapply(.SD, mean, na.rm = TRUE), 
                    by = gene_symbol, 
                    .SDcols = sample_cols]

final_expr <- as.matrix(expr_agg[, ..sample_cols])
rownames(final_expr) <- expr_agg[["gene_symbol"]]

message("Final gene symbol expression matrix: ", nrow(final_expr), " genes x ", ncol(final_expr), " samples")
message("First 10 gene symbols: ", paste(head(rownames(final_expr), 10), collapse = ", "))

if (any(is.na(final_expr))) {
  message("Warning: NA values found in expression matrix. Replacing with 0.")
  final_expr[is.na(final_expr)] <- 0
}

lines <- readLines(matrix_gz)

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

feat_hr <- find_feature(c("hr-positive_yes is 1", "hr"))
feat_her2 <- find_feature(c("her2-positive_yes is 1", "her2"))
feat_resp <- find_feature(c("pcr_yes is 1", "pcr"))

if (is.null(feat_hr) || is.null(feat_her2) || is.null(feat_resp)) {
  stop("Required features not found")
}

sample_id_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
sample_ids <- strsplit(sample_id_line, "\t")[[1]][-1]
sample_ids <- gsub('^"|"$', '', sample_ids)

hr_vals <- extract_from_line(feat_hr)
her2_vals <- extract_from_line(feat_her2)
resp_vals <- extract_from_line(feat_resp)

clin_full <- data.frame(
  sample_id = sample_ids,
  hr = hr_vals,
  her2 = her2_vals,
  pcr = resp_vals,
  stringsAsFactors = FALSE
)

tnbc_idx <- clin_full$hr == "0" & clin_full$her2 == "0" & 
  !is.na(clin_full$hr) & !is.na(clin_full$her2)
clin_tnbc <- clin_full[tnbc_idx, ]
message("TNBC samples: ", nrow(clin_tnbc))

clin_tnbc$pcr_rd <- NA_character_
clin_tnbc$pcr_rd[clin_tnbc$pcr == "1"] <- "pCR"
clin_tnbc$pcr_rd[clin_tnbc$pcr == "0"] <- "RD"

clin_clean <- clin_tnbc[!is.na(clin_tnbc$pcr_rd), ]
message("TNBC samples with valid labels: ", nrow(clin_clean))

common_ids <- intersect(colnames(final_expr), clin_clean$sample_id)
if (length(common_ids) == 0) stop("No common sample IDs between expression matrix and clinical data")

final_expr_tnbc <- final_expr[, common_ids, drop = FALSE]
final_clin <- clin_clean[match(common_ids, clin_clean$sample_id), ]

message("Final expression matrix: ", nrow(final_expr_tnbc), " genes x ", ncol(final_expr_tnbc), " samples")
message("Gene symbol examples: ", paste(head(rownames(final_expr_tnbc), 10), collapse = ", "))

write.csv(final_expr_tnbc, "GSE32603_TNBC_expression.csv", quote = FALSE)
write.csv(final_clin[, c("sample_id", "pcr_rd")], 
          "GSE32603_TNBC_labels.csv", row.names = FALSE)

write.csv(expr_matrix, "GSE32603_TNBC_probe_level_expression.csv", quote = FALSE)

message("\nDone. Output files:")
message("  - GSE32603_TNBC_expression.csv")
message("  - GSE32603_TNBC_labels.csv")
message("  - GSE32603_TNBC_probe_level_expression.csv")