# Analysis code

Author-generated analysis code and processed data for Bandy et al., "Zooplankton
DNA metabarcoding supports regional ocean acidification monitoring and identifies
novel candidate bioindicators in the Southern California Bight" (PLOS ONE,
PONE-D-26-23798).

Raw sequence processing was done with REVAMP (McAllister et al. 2023; Cutadapt ->
DADA2 -> BLAST vs. NCBI nucleotide -> LCA taxonomy) and is not reproduced here.
These scripts start from the resulting ASV count table and take the analysis
through to the manuscript figures, tables, and statistics.

## Running

Run every script **from the repository root** (the folder that contains `data/`
and `figures/`). Each script checks for `data/` and stops with a message if the
working directory is wrong. Scripts write regenerated CSVs back into `data/` and
all plots and derived tables into `figures/`.

Run them in numeric order to regenerate everything from `data/ASVCountsCO1_aligned.csv`,
or run any single script on its own (`data/` ships with every intermediate file
already present).

**All five steps are R.** Steps 01-03 were originally run in MATLAB for the
manuscript's initial submission; the R port was verified against that MATLAB
analysis on the originally submitted data before being adopted as the sole
pipeline (not shipped here): the regenerated `RelabundCO1_aligned_clean.csv`
and `TaxonomyTraits_aligned_clean.csv` were identical to floating-point
precision, correlation coefficients and p-values in
`ASVcorrelations_avgArag_taxonomytraits_clean.csv` matched to ~1e-15, and the
PCA gave the same variance explained either way.

With the `data/` shipped in this repository (taxonomy reassigned against a
current NCBI nucleotide database and restricted to family-or-deeper
assignments, see "Note on taxonomy" below), running `01`-`03` gives 4505 ASVs
retained after filtering and PCA variance explained of PC1 = 71.78%,
PC2 = 11.46%, PC3 = 3.06%.

Requirements: R (4.2 or newer). Package lists are in the header of each script;
the base pipeline (01, 02) needs only `dplyr`.

## Pipeline

| Step | Script | Reads | Writes | Manuscript output |
|---|---|---|---|---|
| 01 | `01_relative_abundance.R` | `ASVCountsCO1_aligned.csv`, `Metadata_aligned.csv`, `Taxonomy_aligned.csv`, `OrderTraits.csv` | `RelabundCO1_aligned_clean.csv`, `TaxonomyTraits_aligned_clean.csv` (+ unaligned variants) | community-composition inputs |
| 02 | `02_aragonite_correlations.R` | `RelabundCO1_aligned_clean.csv`, `Metadata_aligned.csv`, `TaxonomyTraits_aligned_clean.csv` | `ASVcorrelations_avgArag_taxonomytraits_clean.csv` | source for S2/S3 Table and Fig 5 |
| 03 | `03_pca.R` | `RelabundCO1_aligned_clean.csv`, `Metadata_aligned.csv`, `TaxonomyTraits_aligned_clean.csv` | `figures/PCA_OTU_ColorBy_Season2*.png`, `figures/PCA_OTU_ColorBy_CalanoidaRatio2*.png` | **S3 Fig**, **Fig 4** |
| 04 | `04_trait_group_stats.R` | `ASVcorrelations_avgArag_taxonomytraits_clean.csv` | `figures/DunnWithLetters_*.csv`, `figures/GroupLetters_*.csv` | **Fig 6** statistics (Kruskal-Wallis / Dunn's / BH letters) |
| 05 | `05_forest_plot.R` | `ASVcorrelations_avgArag_taxonomytraits_clean.csv` | `figures/forest_posneg_onepanel_occColor.png` (+ per-direction panels) | **Fig 5**; the top-15 sets are **S2 Table** and **S3 Table** |

Carbonate chemistry estimation (the `avgArag`, `avgTA`, `avgDIC`, `avgpH`, and
`avgCTD*` columns in `Metadata_aligned.csv`, and Fig 2 / S1 Fig / S2 Fig) was done
by co-author C.A. Frieder following Alin et al. (2012), with uncertainty
propagation by N. Lombardo. That code and the underlying CTD and
bottle-chemistry data are handled separately; see the root README.

## Methods notes

- **Relative abundance.** `01` normalizes the full ASV count table to relative
  abundance per sample, then removes ASVs flagged in `OrderTraits.csv`
  (`Removal` = "Yes terrestrial" or "Yes non-target"). Relative abundances are
  **not** recomputed after removal, so each value is an ASV's fraction of the
  total reads in that sample.
- **Correlations.** `02` computes a Spearman correlation of every ASV's relative
  abundance against `avgArag` (water-column integrated aragonite saturation
  state, 15-100 m). ASVs with an undefined correlation (zero variance) are
  dropped. The R version uses `cor.test(..., method = "spearman")`, which matches
  MATLAB's `corr(..., 'Type', 'Spearman')` on this data.
- **S2 / S3 Table selection.** `05` selects ASVs with occurrence > 50%, sorts by
  absolute Spearman rho, and takes the top 15 positive and top 15 negative. S2/S3
  Table are those two sets; Fig 5 shows them together.
- **Trait groups.** `04` compares the Spearman coefficients across Order-level
  trait groups (calcification status, mineral form, larval structure) with
  Kruskal-Wallis followed by Dunn's post-hoc, using both Benjamini-Hochberg and
  Holm adjustment; group letters use a 0.05 threshold.

## Data files (`data/`)

| File | Contents |
|---|---|
| `ASVCountsCO1_aligned.csv` | ASV x sample read counts, column order aligned to the metadata |
| `Taxonomy_aligned.csv` | per-ASV taxonomy from REVAMP (Kingdom -> Species) |
| `OrderTraits.csv` | Order-level calcification / structure / mineral-form assignments and the `Removal` flag |
| `Metadata_aligned.csv` | one row per sample: station, agency, season, coordinates, and 15-100 m integrated carbonate chemistry and CTD variables |
| `RelabundCO1_aligned_clean.csv` | output of `01`: relative abundances after trait-based ASV removal |
| `TaxonomyTraits_aligned_clean.csv` | output of `01`: taxonomy joined to Order-level traits, same ASV set as above |
| `ASVcorrelations_avgArag_taxonomytraits_clean.csv` | output of `02`: per-ASV Spearman rho vs avgArag, occurrence, taxonomy, and traits |

## Note on taxonomy

`Taxonomy_aligned.csv` was regenerated against a current NCBI nucleotide
database in response to reviewer comments (original submission used a nt
snapshot from December 2022; this repository ships the rerun against nt June
2026). Per Reviewer 2's recommendation, ASVs whose deepest confident rank
assignment is Order or shallower are treated as unassigned here (all ranks
blanked, so they drop out via the existing non-target/terrestrial removal
rule in step `01`), restricting the analysis to family-or-deeper taxonomic
assignments. `OrderTraits.csv` was updated alongside the taxonomy rerun to add
trait assignments for Order-level taxa that only appear under the new
database.
