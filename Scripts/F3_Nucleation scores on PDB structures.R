#--------------------------
# Load required packages
#--------------------------
library(tidyverse)  # includes dplyr, ggplot2, tidyr, purrr, stringr
library(ggpubr)     # for stat_cor, geom_signif
library(reshape2)   # for melt
library(readxl)     # if Excel files are used

#--------------------------
# Load datasets
#--------------------------
load("All_mutants.RData")  # contains 'all_ds'

#--------------------------
# Panel A: Violin plot of mean nucleation scores per position grouped by residue exposure
# Residues with relative ASA > 0.25 are considered exposed
#--------------------------

# Load and process relative solvent accessibility (RASA) data
RASA <- read.csv("merged_ss_sasa.csv")
rasa_cols <- grep("SASA", colnames(RASA), value = TRUE)

df_rasa <- RASA %>%
  select(AAPos, all_of(rasa_cols)) %>%
  mutate(Pos = as.numeric(str_extract(AAPos, "\\d+"))) %>%
  select(-AAPos) %>%
  pivot_longer(cols = -Pos, names_to = "structure", values_to = "ASA") %>%
  mutate(
    structure = str_remove(structure, "SASA_"),
    ASA = as.numeric(ASA)
  ) %>%
  group_by(structure) %>%
  complete(Pos = 6:37) %>%  # ensure all positions included
  ungroup()

# Filter single amino acid substitutions from all_ds
singles_IAPP <- all_ds %>%
  filter(Mutation_type == "Singles") %>%
  mutate(Pos = as.numeric(str_extract(HGVS_nomenclature, "\\d+")))

# Compute mean nucleation score per position
singles_mean <- singles_IAPP %>%
  group_by(Pos) %>%
  summarize(Mean = mean(fitness_merged, na.rm = TRUE)) %>%
  ungroup()

# Merge mean NS and ASA for plotting
df_not_av <- left_join(singles_mean, df_rasa, by = "Pos") %>%
  filter(!is.na(ASA)) %>%
  mutate(
    buried = if_else(ASA > 0.25, "Exposed", "Buried"),
    structure = factor(structure, levels = c("7M61", "7M62", "7M64", "7M65", "6Y1A", "6ZRF", "8R4I", "6VW2"))
  )

# Generate violin plot per structure
p_strand <- ggplot(df_not_av, aes(x = buried, y = Mean)) +
  geom_violin(fill = "white") +
  geom_boxplot(width = 0.3) +
  geom_point(color = "maroon") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~structure, nrow = 1) +
  geom_signif(
    comparisons = list(c("Buried", "Exposed")),
    test = "t.test",
    map_signif_level = TRUE,
    col = "black",
    size = 0.4,
    y_position = 0.3
  ) +
  labs(
    x = NULL,
    y = "Mean nucleation score per position"
  ) +
  ylim(-4, 1) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.y = element_line(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45)
  )

ggsave("Mean_nucleation_scores_buried_exposed.pdf", p_strand, width = 10, height = 3)

#--------------------------
# Panel C: Correlation between ASA and mean nucleation score
#--------------------------

# Reshape df_rasa to wide format for correlation
df_sasa_wide <- df_rasa %>%
  pivot_wider(names_from = structure, values_from = ASA)

# Merge mean nucleation scores with ASA
NSmean_asa <- left_join(singles_mean, df_sasa_wide, by = "Pos")

# Define PDB IDs of interest
pdb_ids <- c("7M61", "7M62", "7M64", "7M65", "6Y1A", "6ZRF", "8R4I", "6VW2")

# Compute correlation for each structure
cor_results <- map_dfr(pdb_ids, function(pdb) {
  test <- cor.test(NSmean_asa$Mean, NSmean_asa[[pdb]], method = "pearson")
  tibble(PDB = pdb, R = test$estimate, p_val = test$p.value)
}) %>%
  mutate(
    p_sig = p_val < 0.05,
    PDB = factor(PDB, levels = rev(pdb_ids))
  )

# Barplot of correlation coefficients
p_corr <- ggplot(cor_results, aes(x = R, y = PDB)) +
  geom_col(fill = "maroon", width = 0.5) +
  geom_text(data = filter(cor_results, p_sig), aes(label = "*"), hjust = -0.5, vjust = 0.5) +
  labs(x = "ASA and NS corr. coef. (R)", y = "PDB ID") +
  coord_cartesian(xlim = c(0, 0.6)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill = alpha('white', 0.5))
  )

ggsave("Barplot_ASA_vs_NS_corr.pdf", p_corr, width = 5, height = 3)

#--------------------------
# Point-wise correlation plots per PDB - These pltos are also used in Supplementary Figure 3. 
#--------------------------
merged_df_melt <- left_join(singles_mean, df_rasa, by = "Pos") %>%
  filter(!is.na(ASA))

p_point_corr <- ggplot(merged_df_melt, aes(x = Mean, y = ASA)) +
  geom_point() +
  facet_wrap(~structure, nrow = 8, scales = "free") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  xlab("Mean nucleation score") +
  ylab("Relative ASA (from PDB structure)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.1),
    axis.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 14)
  )

ggsave("Pointwise_ASA_vs_NS_correlation.pdf", p_point_corr, width = 2.2, height = 16)
