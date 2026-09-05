# 03_pca.R
# R port of 03_pca.m (manuscript figures only). PCA of the cleaned CO1 ASV
# relative abundance table.
#
# Figures produced:
#   figures/PCA_OTU_ColorBy_Season2_R.png          -> S3 Fig (PCA by season)
#   figures/PCA_OTU_ColorBy_CalanoidaRatio2_R.png  -> Fig 4  (PCA by Calanoida:Euphausiacea ratio)
#
# PCA matches MATLAB pca(): columns centered, not scaled (prcomp center = TRUE,
# scale. = FALSE). Principal component signs are arbitrary and may be mirrored
# relative to the MATLAB figure; the variance explained and the sample grouping
# are identical. (03_pca.m also produces additional exploratory PCA plots that
# are not in the manuscript.)
#
# Inputs (data/): RelabundCO1_aligned_clean.csv, Metadata_aligned.csv,
#                 TaxonomyTraits_aligned_clean.csv
#
# Run from the repository root (the folder containing data/ and figures/).
# Required packages: ggplot2

library(ggplot2)
stopifnot(dir.exists("data"))
dir.create("figures", showWarnings = FALSE)

otu      <- read.csv("data/RelabundCO1_aligned_clean.csv", check.names = FALSE)
metadata <- read.csv("data/Metadata_aligned.csv", check.names = FALSE)
taxtrait <- read.csv("data/TaxonomyTraits_aligned_clean.csv", check.names = FALSE)

metadata$SampleID <- as.character(metadata$SampleID)
sample_cols <- setdiff(names(otu), "ASV")
metadata <- metadata[match(sample_cols, metadata$SampleID), , drop = FALSE]

otu_mat <- as.matrix(otu[, sample_cols])          # ASVs x samples
rownames(otu_mat) <- as.character(otu$ASV)

# ----- PCA (samples as observations) -----
pc <- prcomp(t(otu_mat), center = TRUE, scale. = FALSE)
explained <- 100 * pc$sdev^2 / sum(pc$sdev^2)
message(sprintf("PC1 %.2f%%  PC2 %.2f%%  PC3 %.2f%%", explained[1], explained[2], explained[3]))

scores <- data.frame(
  SampleID = sample_cols,
  PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC3 = pc$x[, 3],
  Season = factor(metadata$Season,
                  levels = c("Spring 2019", "Summer 2019", "Fall 2019", "Winter 2020")),
  row.names = NULL
)

# ----- Calanoida : Euphausiacea ratio -----
asv_order <- setNames(as.character(taxtrait$Order), as.character(taxtrait$ASV))
cal_asv <- names(asv_order)[asv_order == "Calanoida"]
eup_asv <- names(asv_order)[asv_order == "Euphausiacea"]
cal_asv <- intersect(cal_asv, rownames(otu_mat))
eup_asv <- intersect(eup_asv, rownames(otu_mat))
message(sprintf("Matched %d Calanoida and %d Euphausiacea ASVs", length(cal_asv), length(eup_asv)))

sum_cal <- colSums(otu_mat[cal_asv, , drop = FALSE])
sum_eup <- colSums(otu_mat[eup_asv, , drop = FALSE])
denom <- sum_cal + sum_eup
scores$CalanoidaRatio <- ifelse(denom == 0, NA_real_, sum_cal / denom)

season_cols <- c("Spring 2019" = rgb(107, 142, 35, maxColorValue = 255),
                 "Summer 2019" = rgb(255, 215,  0, maxColorValue = 255),
                 "Fall 2019"   = rgb(255, 140,  0, maxColorValue = 255),
                 "Winter 2020" = rgb( 64, 224, 208, maxColorValue = 255))

parula <- c("#352a87", "#0363e1", "#1485d4", "#06a7c6", "#38b99e",
            "#92bf73", "#d9ba56", "#fcce2e", "#f9fb0e")

xlab <- sprintf("PC1 (%.2f%% Variance Explained)", explained[1])
ylab <- sprintf("PC2 (%.2f%% Variance Explained)", explained[2])

# ----- S3 Fig: PCA colored by season -----
g_season <- ggplot(scores, aes(PC1, PC2, fill = Season)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(shape = 21, size = 3, colour = "grey20") +
  scale_fill_manual(values = season_cols) +
  labs(x = xlab, y = ylab) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.position = c(0.88, 0.16))
ggsave("figures/PCA_OTU_ColorBy_Season2_R.png", g_season, width = 8, height = 6, dpi = 300)

# ----- Fig 4: PCA colored by Calanoida:Euphausiacea ratio -----
g_ratio <- ggplot(scores, aes(PC1, PC2, colour = CalanoidaRatio)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(size = 3) +
  scale_colour_gradientn(colours = parula, limits = c(0, 1),
                         name = "Calanoida/(Calanoida + Euphausiacea)",
                         guide = guide_colourbar(title.position = "top", barwidth = 12)) +
  labs(x = xlab, y = ylab) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom")
ggsave("figures/PCA_OTU_ColorBy_CalanoidaRatio2_R.png", g_ratio, width = 8, height = 6.5, dpi = 300)

message("Saved figures/PCA_OTU_ColorBy_Season2_R.png and figures/PCA_OTU_ColorBy_CalanoidaRatio2_R.png")
