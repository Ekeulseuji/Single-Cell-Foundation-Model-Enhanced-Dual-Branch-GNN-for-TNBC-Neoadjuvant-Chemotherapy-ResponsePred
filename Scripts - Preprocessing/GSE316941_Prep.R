library(tidyverse)

base_dir <- "/Users/ekeulseuji/Downloads"
setwd(base_dir)

tpm_file <- "GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz"
label_file <- "NeoTRIP_baseline_D1C2_pheno_and_clinical.txt"

if (!file.exists(tpm_file)) stop("TPM file not found: ", tpm_file)
if (!file.exists(label_file)) stop("Label file not found: ", label_file)

message("Reading TPM matrix...")
tpm_raw <- read.table(tpm_file, header = TRUE, check.names = FALSE, sep = "\t", stringsAsFactors = FALSE)
gene_symbols <- tpm_raw[, 1]
sample_cols <- colnames(tpm_raw)[-(1:2)]
expr_matrix <- tpm_raw[, -(1:2), drop = FALSE]
rownames(expr_matrix) <- gene_symbols
colnames(expr_matrix) <- sample_cols
message("TPM matrix dimensions: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " samples")
message("Sample column examples: ", paste(head(sample_cols), collapse = ", "))

message("Reading label file...")
labels_raw <- read.delim(label_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

labels_baseline <- labels_raw[labels_raw$Timepoint == "Baseline", ]
message("Baseline samples: ", nrow(labels_baseline))

labels_clean <- labels_baseline[, c("Sample_ID", "pCR")]
colnames(labels_clean) <- c("sample_id", "pcr_rd")

labels_clean$pcr_rd <- ifelse(labels_clean$pcr_rd == "pCR", "pCR", "RD")
message("Label distribution:")
print(table(labels_clean$pcr_rd))

common_samples <- intersect(colnames(expr_matrix), labels_clean$sample_id)
message("Common samples: ", length(common_samples))

if (length(common_samples) == 0) {
  stop("No overlap between expression columns and Sample_IDs")
}

expr_final <- expr_matrix[, common_samples, drop = FALSE]
labels_final <- labels_clean[match(common_samples, labels_clean$sample_id), ]
stopifnot(all(colnames(expr_final) == labels_final$sample_id))

message("Final expression matrix: ", nrow(expr_final), " genes x ", ncol(expr_final), " samples")
message("Final label distribution:")
print(table(labels_final$pcr_rd))

write.csv(expr_final, "GSE319641_TNBC_expression_TPM.csv", quote = FALSE)
write.csv(labels_final, "GSE319641_TNBC_labels.csv", row.names = FALSE)

message("Preprocessing complete. Output files:")
message("  - GSE319641_TNBC_expression_TPM.csv")
message("  - GSE319641_TNBC_labels.csv")