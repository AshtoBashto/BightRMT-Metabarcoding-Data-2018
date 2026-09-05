# 05_forest_plot.R
# Forest plot of the ASVs most strongly correlated with avgArag (Fig 5), and the
# per-direction top-15 panels. The top-15 sets are the source of S2 Table
# (positive) and S3 Table (negative).
#
# Selection rule: ASVs with occurrence > 50%, sorted by absolute Spearman rho,
# top 15 positive and top 15 negative.
#
# Input   (data/):    ASVcorrelations_avgArag_taxonomytraits_clean.csv
# Outputs (figures/):  forest_15pos.png, forest_15neg.png,
#                      forest_posneg_onepanel.png, forest_posneg_onepanel_occColor.png (Fig 5)
#
# Run from the repository root (the folder containing data/ and figures/).
# Required packages: dplyr, ggplot2, readr, stringr, ggtext, patchwork, viridis
#   install.packages(c("dplyr","ggplot2","readr","stringr","ggtext","patchwork","viridis"))

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(ggtext)
library(patchwork)
library(viridis)

stopifnot(dir.exists("data"))
dir.create("figures", showWarnings = FALSE)

# ----- Build the top-15 positive/negative sets from the full correlation table -----
corr <- read_csv("data/ASVcorrelations_avgArag_taxonomytraits_clean.csv")

sel <- corr %>%
  filter(OccurrencePct > 50) %>%
  transmute(
    ASV, Phylum, Class, Order, Family, Genus, Species,
    `avgArag Corr` = SpearmanRho,
    pValue,
    `Occ (%)` = OccurrencePct
  )

pos <- sel %>% filter(`avgArag Corr` > 0) %>% arrange(desc(`avgArag Corr`)) %>% slice_head(n = 15) %>% mutate(n = 80L)
neg <- sel %>% filter(`avgArag Corr` < 0) %>% arrange(`avgArag Corr`)       %>% slice_head(n = 15) %>% mutate(n = 80L)

make_forest <- function(df, outfile){
  df1 <- df %>%
    rename(
      corr = `avgArag Corr`, p = pValue, occ = `Occ (%)`,
      phylum = Phylum, cls = Class, ord = Order, family = Family,
      genus = Genus, species = Species, asv = ASV
    ) %>%
    mutate(
      species_lab = case_when(
        !is.na(species) & species != "" ~ species,
        !is.na(genus)   & genus   != "" ~ paste0(genus, " sp."),
        !is.na(ord)     & ord     != "" ~ ord,
        TRUE ~ "Unknown taxon"
      ),
      label_md = paste0("ASV ", asv, " – <i>", species_lab, "</i>")
    )

  if ("n" %in% names(df1)) {
    z  <- atanh(df1$corr)
    se <- 1 / sqrt(pmax(df1$n - 3, 1))
    df1$lo <- tanh(z - 1.96 * se)
    df1$hi <- tanh(z + 1.96 * se)
  }

  df1 <- df1 %>% arrange(desc(abs(corr))) %>% mutate(idx = row_number())
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

  g <- ggplot(df1, aes(x = corr, y = idx)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0) +
    geom_point(shape = 16, size = 2.8) +
    coord_cartesian(xlim = c(-1, 1), clip = "off") +
    scale_y_continuous(breaks = df1$idx, labels = df1$label_md) +
    geom_text(aes(x = 1.02, y = idx,
                  label = ifelse(is.na(occ), "", paste0("Occ = ", formatC(occ, format = "f", digits = 0), "%"))),
              hjust = 0, size = 3) +
    labs(x = "Correlation (Spearman ρ)", y = NULL, title = NULL,
         caption = paste0("Points = ρ; bars = 95% CIs via Fisher z (approx). N = ",
                          if ("n" %in% names(df1)) unique(df1$n) else "?")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(), legend.position = "none",
          plot.margin = margin(5.5, 78, 5.5, 5.5, "pt"), axis.text.y = element_markdown())

  ggsave(outfile, g, width = 4.6, height = 5.0, units = "in")
}

make_forest(pos, "figures/forest_15pos.png")
make_forest(neg, "figures/forest_15neg.png")

# ----- One-panel versions (positives top -> negatives bottom) -----
prep <- function(df){
  df %>%
    rename(
      corr = `avgArag Corr`, p = pValue, occ = `Occ (%)`,
      phylum = Phylum, cls = Class, ord = Order, family = Family,
      genus = Genus, species = Species, asv = ASV
    ) %>%
    mutate(
      species_lab = dplyr::case_when(
        !is.na(species) & species != "" ~ species,
        !is.na(genus)   & genus   != "" ~ paste0(genus, " sp."),
        !is.na(ord)     & ord     != "" ~ ord,
        TRUE ~ "Unknown taxon"
      ),
      label_md = paste0("ASV ", asv, " – <i>", species_lab, "</i>")
    )
}

make_forest_onepanel <- function(pos, neg, outfile, xlim = c(-1, 1), width = 6.8, height = 7.2){
  both <- dplyr::bind_rows(prep(pos), prep(neg))
  if ("n" %in% names(both)) {
    z  <- atanh(both$corr)
    se <- 1 / sqrt(pmax(both$n - 3, 1))
    both$lo <- tanh(z - 1.96 * se)
    both$hi <- tanh(z + 1.96 * se)
  }
  both <- both %>% arrange(desc(corr)) %>% mutate(idx = row_number())
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

  g <- ggplot(both, aes(x = corr, y = idx)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0) +
    geom_point(shape = 16, size = 2.8) +
    coord_cartesian(xlim = xlim, clip = "off") +
    scale_y_reverse(breaks = both$idx, labels = both$label_md) +
    geom_text(aes(x = xlim[2] + 0.02,
                  label = ifelse(is.na(occ), "", paste0("Occ = ", formatC(occ, format = "f", digits = 0), "%"))),
              hjust = 0, size = 3) +
    labs(x = "Correlation (Spearman ρ)", y = NULL, title = NULL,
         caption = paste0("Points = ρ; bars = 95% CIs via Fisher z (approx). N = ",
                          if ("n" %in% names(both)) unique(both$n) else "?")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(), legend.position = "none",
          plot.margin = margin(5.5, 90, 5.5, 5.5, "pt"), axis.text.y = element_markdown())

  ggsave(outfile, g, width = width, height = height, units = "in")
}

make_forest_onepanel(pos, neg, "figures/forest_posneg_onepanel.png")

# ----- Fig 5: one panel, points colored by occurrence -----
make_forest_onepanel_occColor <- function(pos, neg, outfile, xlim = c(-1, 1), width = 6.8, height = 7.2){
  both <- dplyr::bind_rows(prep(pos), prep(neg))
  if ("n" %in% names(both)) {
    z  <- atanh(both$corr)
    se <- 1 / sqrt(pmax(both$n - 3, 1))
    both$lo <- tanh(z - 1.96 * se)
    both$hi <- tanh(z + 1.96 * se)
  }
  both <- both %>% arrange(desc(corr)) %>% mutate(idx = dplyr::row_number())
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

  g <- ggplot(both, aes(x = corr, y = idx)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0) +
    geom_point(aes(color = occ), shape = 16, size = 2.9, na.rm = TRUE) +
    coord_cartesian(xlim = xlim, clip = "off") +
    scale_y_reverse(breaks = both$idx, labels = both$label_md) +
    scale_color_viridis_c(name = "Occurrence (%)", option = "D", end = 0.95, na.value = "grey80") +
    labs(x = "Correlation (Spearman ρ)", y = NULL, title = NULL,
         caption = paste0("Points = ρ; bars = 95% CIs via Fisher z (approx). N = ",
                          if ("n" %in% names(both)) unique(both$n) else "?")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(), legend.position = "right",
          plot.margin = margin(5.5, 18, 5.5, 5.5, "pt"), axis.text.y = element_markdown())

  ggsave(outfile, g, width = width, height = height, units = "in")
}

make_forest_onepanel_occColor(pos, neg, "figures/forest_posneg_onepanel_occColor.png")
