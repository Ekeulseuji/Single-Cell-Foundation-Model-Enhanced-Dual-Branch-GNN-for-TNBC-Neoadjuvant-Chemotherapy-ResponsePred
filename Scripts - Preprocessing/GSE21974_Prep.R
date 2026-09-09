library(tidyverse)
library(limma)
library(GEOquery)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

raw_folder <- "GSE21974_RAW"

data_files <- list.files(raw_folder, pattern = "\\.txt$", full.names = TRUE)
message("Found ", length(data_files), " txt files")

parse_agilent_file <- function(file_path) {
  all_lines <- readLines(file_path, warn = FALSE)
  
  header_line <- NULL
  data_start <- -1
  
  for (i in seq_along(all_lines)) {
    line <- all_lines[i]
    if (grepl("SystematicName", line, ignore.case = TRUE) || 
        grepl("ProbeName", line, ignore.case = TRUE)) {
      header_line <- line
      data_start <- i + 1
      break
    }
  }
  
  if (is.null(header_line)) {
    stop("Cannot find header line")
  }
  
  col_names <- strsplit(header_line, "\t")[[1]]
  col_names <- trimws(col_names)
  if (length(col_names) > 0 && col_names[1] == "FEATURES") {
    col_names <- col_names[-1]
  }
  
  data_rows <- all_lines[data_start:length(all_lines)]
  data_rows <- data_rows[grepl("^DATA", data_rows)]
  
  if (length(data_rows) == 0) {
    stop("No data rows found")
  }
  
  data_list <- list()
  for (row in data_rows) {
    parts <- strsplit(row, "\t")[[1]]
    if (length(parts) > 0 && parts[1] == "DATA") {
      parts <- parts[-1]
    }
    if (length(parts) == length(col_names)) {
      data_list <- c(data_list, list(parts))
    }
  }
  
  if (length(data_list) == 0) {
    stop("No data rows parsed")
  }
  
  df <- do.call(rbind, data_list)
  colnames(df) <- col_names
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  return(df)
}

matrix_gz <- "GSE21974_series_matrix.txt.gz"
lines <- readLines(matrix_gz)

sample_id_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
sample_ids <- strsplit(sample_id_line, "\t")[[1]][-1]
sample_ids <- gsub('^"|"$', '', sample_ids)

char_lines <- grep("^!Sample_characteristics_ch", lines)
timepoint_values <- list()
for (idx in char_lines) {
  line <- lines[idx]
  if (grepl("before/after chemotherapy", line, ignore.case = TRUE)) {
    parts <- strsplit(line, "\t")[[1]]
    values <- parts[-1]
    values <- gsub('^"|"$', '', values)
    timepoint_values <- values
    break
  }
}

before_samples <- sample_ids[grepl("before", timepoint_values, ignore.case = TRUE)]
message("Before samples: ", length(before_samples))

get_gsm <- function(filename) {
  basename <- basename(filename)
  m <- regexpr("GSM[0-9]+", basename)
  if (m > 0) {
    return(substr(basename, m, m + attr(m, "match.length") - 1))
  }
  return(NA)
}

before_files <- data_files[sapply(data_files, function(f) get_gsm(f) %in% before_samples)]
message("Found ", length(before_files), " before sample files")

all_expr <- list()

for (i in seq_along(before_files)) {
  file_path <- before_files[i]
  gsm_id <- get_gsm(file_path)
  
  message("Processing sample ", i, "/", length(before_files), ": ", gsm_id)
  
  tryCatch({
    dat <- parse_agilent_file(file_path)
    
    probe_col <- grep("SystematicName", colnames(dat), ignore.case = TRUE)[1]
    if (is.na(probe_col)) probe_col <- 1
    
    green_col <- grep("gMedianSignal", colnames(dat), ignore.case = TRUE)[1]
    red_col <- grep("rMedianSignal", colnames(dat), ignore.case = TRUE)[1]
    
    if (is.na(green_col)) {
      signal_cols <- grep("Signal", colnames(dat), ignore.case = TRUE)
      if (length(signal_cols) >= 2) {
        green_col <- signal_cols[1]
        red_col <- signal_cols[2]
      } else {
        stop("Cannot identify signal columns")
      }
    }
    
    probe_ids <- as.character(dat[, probe_col])
    
    if (!is.na(red_col) && red_col != green_col) {
      green_signal <- as.numeric(dat[, green_col])
      red_signal <- as.numeric(dat[, red_col])
      expr_values <- log2(red_signal / green_signal)
    } else {
      green_signal <- as.numeric(dat[, green_col])
      expr_values <- log2(green_signal)
    }
    
    expr_values[is.infinite(expr_values) | is.nan(expr_values) | is.na(expr_values)] <- NA
    
    temp_df <- data.frame(
      probe_id = probe_ids,
      expr = expr_values,
      stringsAsFactors = FALSE
    )
    
    temp_df <- temp_df[!is.na(temp_df$probe_id) & temp_df$probe_id != "" & 
                         temp_df$probe_id != "NA" & !is.na(temp_df$expr), ]
    
    if (any(duplicated(temp_df$probe_id))) {
      temp_df <- aggregate(expr ~ probe_id, data = temp_df, FUN = median, na.rm = TRUE)
    }
    
    all_expr[[gsm_id]] <- temp_df
    message("  Retained ", nrow(temp_df), " probes")
    
  }, error = function(e) {
    message("  Processing failed: ", e$message)
  })
}

if (length(all_expr) == 0) stop("No before samples processed successfully")

message("\nMerging all before samples...")

all_probes <- unique(unlist(lapply(all_expr, function(x) x$probe_id)))
message("Total probes: ", length(all_probes))

sample_names <- names(all_expr)
expr_matrix <- matrix(NA, nrow = length(all_probes), ncol = length(sample_names))
rownames(expr_matrix) <- all_probes
colnames(expr_matrix) <- sample_names

for (i in seq_along(sample_names)) {
  sample_id <- sample_names[i]
  temp_df <- all_expr[[sample_id]]
  idx <- match(temp_df$probe_id, all_probes)
  expr_matrix[idx, i] <- temp_df$expr
}

na_prop <- rowSums(is.na(expr_matrix)) / ncol(expr_matrix)
keep_probes <- na_prop < 0.5
expr_matrix <- expr_matrix[keep_probes, ]
message("Probes after NA filtering: ", nrow(expr_matrix))

for (i in 1:nrow(expr_matrix)) {
  row_median <- median(expr_matrix[i, ], na.rm = TRUE)
  if (!is.na(row_median) && is.finite(row_median)) {
    expr_matrix[i, is.na(expr_matrix[i, ])] <- row_median
  }
}

message("\nMapping probes to gene symbols...")

gpl <- getGEO("GPL6480", destdir = getwd())
gpl_table <- Table(gpl)

probe_col <- "ID"
gene_col <- "GeneName"

if (!probe_col %in% colnames(gpl_table)) {
  probe_col <- grep("ID|Probe", colnames(gpl_table), ignore.case = TRUE)[1]
}
if (!gene_col %in% colnames(gpl_table)) {
  gene_col <- grep("GeneName|Symbol|Gene", colnames(gpl_table), ignore.case = TRUE)[1]
}
if (is.na(probe_col) || is.na(gene_col)) {
  stop("Cannot find probe ID or gene symbol columns. Available columns: ", paste(colnames(gpl_table), collapse = ", "))
}

message("Using probe column: ", probe_col, ", gene symbol column: ", gene_col)

mapping <- data.frame(
  probe_id = as.character(gpl_table[[probe_col]]),
  gene = as.character(gpl_table[[gene_col]]),
  stringsAsFactors = FALSE
)

mapping <- mapping[!is.na(mapping$gene) & mapping$gene != "", ]
mapping$gene <- trimws(mapping$gene)
mapping$gene <- sapply(strsplit(mapping$gene, " /// "), `[`, 1)
mapping$gene <- sapply(strsplit(mapping$gene, " / "), `[`, 1)
mapping$gene <- sapply(strsplit(mapping$gene, "; "), `[`, 1)
mapping <- mapping[!is.na(mapping$gene) & mapping$gene != "", ]

expr_df <- data.frame(expr_matrix, row.names = NULL, stringsAsFactors = FALSE)
expr_df$probe <- rownames(expr_matrix)

expr_df <- merge(expr_df, mapping, by.x = "probe", by.y = "probe_id", all.x = FALSE)
message("Probes retained after mapping: ", nrow(expr_df))

dup_genes <- expr_df$gene[duplicated(expr_df$gene)]
if (length(dup_genes) > 0) {
  message("Found ", length(unique(dup_genes)), " genes with multiple probes, aggregating by mean")
}

expr_df <- expr_df[, !(colnames(expr_df) %in% "probe")]
expr_gene <- aggregate(. ~ gene, data = expr_df, FUN = mean, na.rm = TRUE)
rownames(expr_gene) <- expr_gene$gene
expr_gene <- expr_gene[, -1]

message("Gene expression matrix: ", nrow(expr_gene), " genes x ", ncol(expr_gene), " samples")
message("Sample columns: ", paste(colnames(expr_gene), collapse = ", "))

if (any(is.na(rownames(expr_gene)))) {
  warning("NA gene symbols found")
}

write.csv(expr_gene, "GSE21974_TNBC_expression.csv", quote = FALSE)
message("\nBefore sample expression matrix saved with ", nrow(expr_gene), " gene symbols.")