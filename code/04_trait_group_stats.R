# 04_trait_group_stats.R
# Kruskal-Wallis + Dunn's post-hoc (Benjamini-Hochberg and Holm) tests comparing
# ASV-avgArag Spearman correlation coefficients across trait groups
# (calcification status, mineral form, larval structure). Produces the group
# letters and pairwise comparison tables that annotate Fig 6.
#
# Input  (data/):    ASVcorrelations_avgArag_taxonomytraits_clean.csv
# Outputs (figures/): DunnWithLetters_<group>.csv, GroupLetters_<group>.csv
#
# Run from the repository root (the folder containing data/ and figures/).
# Required packages: dunn.test, FSA, readr, dplyr, tidyr, rcompanion, tibble, multcompView
#   install.packages(c("dunn.test","FSA","readr","dplyr","tidyr","rcompanion","tibble","multcompView"))

library(dunn.test)
library(FSA)
library(readr)
library(dplyr)
library(tidyr)
library(rcompanion)
library(tibble)
library(multcompView)

stopifnot(dir.exists("data"))
dir.create("figures", showWarnings = FALSE)

# ----- Load Data -----
corrtraits <- read_csv("data/ASVcorrelations_avgArag_taxonomytraits_clean.csv")

# ----- Filter for Significant Correlations -----
sigOnly <- corrtraits %>% filter(pValue < 0.05)

run_kruskal_dunn_export <- function(df, value_col, group_col, title, out_prefix) {
  message("\n------ ", title, " ------")

  df <- df %>% filter(!is.na(.data[[group_col]]))
  df[[group_col]] <- as.factor(df[[group_col]])

  # BH-adjusted Dunn test
  bh_raw <- dunnTest(as.formula(paste(value_col, "~", group_col)), data = df, method = "bh")$res
  bh_raw <- bh_raw %>% filter(!is.na(Comparison) & grepl(" - ", Comparison))
  bh_df <- bh_raw %>%
    separate(Comparison, into = c("Group1", "Group2"), sep = " - ") %>%
    mutate(
      Method = "BH",
      p.adj = as.numeric(P.adj),
      Signif = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        p.adj < 0.05  ~ "*",
        TRUE          ~ ""
      ),
      Comparison = paste(Group1, "-", Group2)
    )

  # Holm-adjusted Dunn test
  holm_raw <- dunnTest(as.formula(paste(value_col, "~", group_col)), data = df, method = "holm")$res
  holm_raw <- holm_raw %>% filter(!is.na(Comparison) & grepl(" - ", Comparison))
  holm_df <- holm_raw %>%
    separate(Comparison, into = c("Group1", "Group2"), sep = " - ") %>%
    mutate(
      Method = "Holm",
      p.adj = as.numeric(P.adj),
      Comparison = paste(Group1, "-", Group2)
    )

  # Assign compact-letter-display group letters based on adjusted p-values
  to_matrix <- function(df, all_groups) {
    mat <- matrix(1, nrow = length(all_groups), ncol = length(all_groups),
                  dimnames = list(all_groups, all_groups))
    for (i in seq_len(nrow(df))) {
      mat[df$Group1[i], df$Group2[i]] <- df$p.adj[i]
      mat[df$Group2[i], df$Group1[i]] <- df$p.adj[i]
    }
    mat
  }

  all_groups <- sort(unique(c(bh_df$Group1, bh_df$Group2)))
  bh_letters <- multcompLetters(to_matrix(bh_df, all_groups), threshold = 0.05)$Letters %>%
    enframe(name = "Group", value = "Letter_BH")
  holm_letters <- multcompLetters(to_matrix(holm_df, all_groups), threshold = 0.05)$Letters %>%
    enframe(name = "Group", value = "Letter_Holm")

  letters_df <- full_join(bh_letters, holm_letters, by = "Group") %>% arrange(Group)

  export_df <- bh_df %>%
    select(Group1, Group2, Z, P.unadj, p.adj, Signif) %>%
    left_join(letters_df, by = c("Group1" = "Group")) %>%
    rename(Letter_BH_1 = Letter_BH, Letter_Holm_1 = Letter_Holm) %>%
    left_join(letters_df, by = c("Group2" = "Group")) %>%
    rename(Letter_BH_2 = Letter_BH, Letter_Holm_2 = Letter_Holm) %>%
    select(Group1, Group2, Z, P.unadj, p.adj, Signif,
           Letter_BH_1, Letter_Holm_1, Letter_BH_2, Letter_Holm_2)

  write_csv(export_df, file.path("figures", paste0("DunnWithLetters_", out_prefix, ".csv")))
  write_csv(letters_df, file.path("figures", paste0("GroupLetters_", out_prefix, ".csv")))

  group_counts <- df %>%
    count(.data[[group_col]]) %>%
    rename(Group = 1, N = n) %>%
    arrange(Group)
  print(group_counts)

  list(comparisons = export_df, letters = letters_df)
}

# ----- Run Tests -----
results <- list()
results$CalcifierGeneral_sigOnly <- run_kruskal_dunn_export(sigOnly,   "SpearmanRho", "CalcifierGeneral", "CalcifierGeneral (sigOnly)", "CalcifierGeneral_sigOnly")
results$CalcifierGeneral_all     <- run_kruskal_dunn_export(corrtraits, "SpearmanRho", "CalcifierGeneral", "CalcifierGeneral (all)",     "CalcifierGeneral_all")
results$MineralForm_sigOnly      <- run_kruskal_dunn_export(sigOnly,   "SpearmanRho", "MineralForm",      "MineralForm (sigOnly)",      "MineralForm_sigOnly")
results$MineralForm_all          <- run_kruskal_dunn_export(corrtraits, "SpearmanRho", "MineralForm",      "MineralForm (all)",          "MineralForm_all")
results$StructureLarval_sigOnly  <- run_kruskal_dunn_export(sigOnly,   "SpearmanRho", "StructureLarval",  "StructureLarval (sigOnly)",  "StructureLarval_sigOnly")
results$StructureLarval_all      <- run_kruskal_dunn_export(corrtraits, "SpearmanRho", "StructureLarval",  "StructureLarval (all)",      "StructureLarval_all")
