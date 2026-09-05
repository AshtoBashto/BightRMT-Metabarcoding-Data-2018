# 01_relative_abundance.R
# R port of 01_relative_abundance.m. Produces the same cleaned, sample-aligned
# CO1 ASV relative abundance table and taxonomy-traits table.
#
# Workflow:
#   1. Load aligned ASV counts, metadata, taxonomy, and Order-level traits.
#   2. Align metadata sample order to the counts table.
#   3. Merge taxonomy with Order-level traits by Order.
#   4. Normalize the full counts table to relative abundance, by sample.
#   5. Remove ASVs flagged in OrderTraits.csv (Removal = "Yes terrestrial" or
#      "Yes non-target"). Relative abundances are NOT recomputed after removal.
#
# Inputs  (data/):  ASVCountsCO1_aligned.csv, Metadata_aligned.csv,
#                   Taxonomy_aligned.csv, OrderTraits.csv
# Outputs (data/):  TaxonomyTraits_aligned.csv, TaxonomyTraits_aligned_clean.csv,
#                   RelabundCO1_clean.csv, RelabundCO1_aligned_clean.csv
#
# Run from the repository root (the folder containing data/ and figures/).
# Required packages: dplyr

library(dplyr)
stopifnot(dir.exists("data"))

# ----- Load -----
counts   <- read.csv("data/ASVCountsCO1_aligned.csv", check.names = FALSE)
metadata <- read.csv("data/Metadata_aligned.csv", check.names = FALSE)
taxonomy <- read.csv("data/Taxonomy_aligned.csv", check.names = FALSE)
traits   <- read.csv("data/OrderTraits.csv", check.names = FALSE, fileEncoding = "UTF-8-BOM")

metadata$SampleID <- as.character(metadata$SampleID)
sample_cols <- setdiff(names(counts), "ASV")

# ----- Align metadata sample order with the counts columns -----
if (!all(sample_cols %in% metadata$SampleID)) {
  stop("Some count sample names are not found in metadata. Check for mismatches.")
}
metadata <- metadata[match(sample_cols, metadata$SampleID), , drop = FALSE]
message("Sample names aligned between counts and metadata tables.")

# ----- Merge taxonomy with Order-level traits -----
taxonomy$Order <- as.character(taxonomy$Order)
traits$Order   <- as.character(traits$Order)
taxonomytraits <- left_join(taxonomy, traits, by = "Order")
message("Taxonomy merged with traits by Order.")

write.csv(taxonomytraits, "data/TaxonomyTraits_aligned.csv", row.names = FALSE, na = "NA")
message("Saved data/TaxonomyTraits_aligned.csv")

# ----- Normalize counts to relative abundance (per sample, before filtering) -----
count_mat <- as.matrix(counts[, sample_cols])
rownames(count_mat) <- as.character(counts$ASV)

col_sums <- colSums(count_mat)
if (any(col_sums == 0)) {
  warning("Some samples have zero total counts and are dropped: ",
          paste(names(col_sums)[col_sums == 0], collapse = ", "))
  keep <- col_sums != 0
  count_mat <- count_mat[, keep, drop = FALSE]
  col_sums <- col_sums[keep]
  sample_cols <- sample_cols[keep]
  metadata <- metadata[metadata$SampleID %in% sample_cols, , drop = FALSE]
}

relab <- sweep(count_mat, 2, col_sums, "/")
stopifnot(all(abs(colSums(relab) - 1) < 1e-10))
message("All sample columns normalized to relative abundance (sums ~ 1).")

# ----- Remove flagged ASVs (after normalization, no re-normalization) -----
remove_asv <- taxonomytraits$ASV[taxonomytraits$Removal %in% c("Yes terrestrial", "Yes non-target")]
message("ASVs flagged for removal: ", length(remove_asv))

relab_clean <- relab[!(rownames(relab) %in% as.character(remove_asv)), , drop = FALSE]
taxonomytraits_clean <- taxonomytraits[!(taxonomytraits$ASV %in% remove_asv), , drop = FALSE]
message("OTU table: ", nrow(relab), " ASVs before, ", nrow(relab_clean), " ASVs after filtering.")

write.csv(taxonomytraits_clean, "data/TaxonomyTraits_aligned_clean.csv", row.names = FALSE, na = "NA")
message("Saved data/TaxonomyTraits_aligned_clean.csv")

# ----- Save cleaned relative abundance table -----
relab_clean_df <- data.frame(ASV = rownames(relab_clean), relab_clean, check.names = FALSE)
write.csv(relab_clean_df, "data/RelabundCO1_clean.csv", row.names = FALSE)
message("Saved data/RelabundCO1_clean.csv")

# ----- Save cleaned, metadata-aligned relative abundance table -----
if (!all(metadata$SampleID %in% colnames(relab_clean))) {
  stop("Some SampleIDs in metadata are not found in the cleaned OTU table.")
}
relab_aligned <- relab_clean[, metadata$SampleID, drop = FALSE]
relab_aligned_df <- data.frame(ASV = rownames(relab_aligned), relab_aligned, check.names = FALSE)
write.csv(relab_aligned_df, "data/RelabundCO1_aligned_clean.csv", row.names = FALSE)
message("Saved data/RelabundCO1_aligned_clean.csv")
