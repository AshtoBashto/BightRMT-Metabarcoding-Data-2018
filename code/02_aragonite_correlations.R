# 02_aragonite_correlations.R
# R port of 02_aragonite_correlations.m. Spearman correlation of each ASV's
# relative abundance against avgArag (15-100 m integrated aragonite saturation
# state), joined with ASV occurrence and the taxonomy-traits table.
#
# Inputs  (data/):  RelabundCO1_aligned_clean.csv, Metadata_aligned.csv,
#                   TaxonomyTraits_aligned_clean.csv
# Output  (data/):  ASVcorrelations_avgArag_taxonomytraits_clean.csv
#
# rho and p-value come from cor.test(x, y, method = "spearman"), which reproduces
# MATLAB corr(..., 'Type', 'Spearman') to floating-point precision on this data.
#
# Run from the repository root (the folder containing data/ and figures/).
# Required packages: dplyr

library(dplyr)
stopifnot(dir.exists("data"))

otu      <- read.csv("data/RelabundCO1_aligned_clean.csv", check.names = FALSE)
metadata <- read.csv("data/Metadata_aligned.csv", check.names = FALSE)
taxtrait <- read.csv("data/TaxonomyTraits_aligned_clean.csv", check.names = FALSE)

metadata$SampleID <- as.character(metadata$SampleID)
sample_cols <- setdiff(names(otu), "ASV")

# ----- Align metadata to OTU sample order -----
if (!all(sample_cols %in% metadata$SampleID)) {
  stop("Some OTU sample names are not found in metadata. Check for mismatches.")
}
metadata <- metadata[match(sample_cols, metadata$SampleID), , drop = FALSE]
message("Sample names aligned between OTU and metadata tables.")

if (!("avgArag" %in% names(metadata))) stop("avgArag column is missing from metadata.")
arag <- metadata$avgArag
message(sprintf("avgArag range: %.4f to %.4f", min(arag), max(arag)))

otu_mat <- as.matrix(otu[, sample_cols])
rownames(otu_mat) <- as.character(otu$ASV)
n <- ncol(otu_mat)

# ----- Spearman correlation per ASV -----
spearman_with_p <- function(x, y) {
  if (sd(x) == 0 || anyNA(x)) return(c(rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(tryCatch(cor.test(x, y, method = "spearman"),
                                  error = function(e) NULL))
  if (is.null(ct)) return(c(rho = NA_real_, p = NA_real_))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

res <- t(apply(otu_mat, 1, spearman_with_p, y = arag))
corr_tbl <- data.frame(
  ASV = rownames(otu_mat),
  SpearmanRho = res[, "rho"],
  pValue = res[, "p"],
  row.names = NULL,
  stringsAsFactors = FALSE
)

n_nan <- sum(is.na(corr_tbl$SpearmanRho))
message("ASVs with undefined correlation (dropped): ", n_nan)
corr_tbl <- corr_tbl[!is.na(corr_tbl$SpearmanRho), , drop = FALSE]
message("ASVs significantly correlated with avgArag (p < 0.05): ",
        sum(corr_tbl$pValue < 0.05))

# ----- Occurrence -----
occ_count <- rowSums(otu_mat > 0)
occ_tbl <- data.frame(
  ASV = names(occ_count),
  Occurrence = as.integer(occ_count),
  OccurrencePct = 100 * occ_count / n,
  row.names = NULL,
  stringsAsFactors = FALSE
)

# ----- Join and sort -----
merged <- corr_tbl %>%
  inner_join(occ_tbl, by = "ASV") %>%
  inner_join(taxtrait %>% mutate(ASV = as.character(ASV)), by = "ASV") %>%
  arrange(desc(SpearmanRho))

write.csv(merged, "data/ASVcorrelations_avgArag_taxonomytraits_clean.csv",
          row.names = FALSE, na = "NA")
message("Saved data/ASVcorrelations_avgArag_taxonomytraits_clean.csv  (",
        nrow(merged), " ASVs)")
