# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(stringr)

# --- Load SASA data ---
sasa_files <- list.files(".", pattern = "merged")
df_sasa_raw <- read.csv(sasa_files[1])

# Extract SASA columns and positions
sasa_cols <- grep("SASA", colnames(df_sasa_raw), value = TRUE)
df_sasa <- df_sasa_raw[, sasa_cols]
df_sasa$Pos <- as.numeric(str_extract(df_sasa_raw$AAPos, "[0-9]+"))

# Reshape SASA data to long format
df_sasa_long <- melt(df_sasa, id.vars = "Pos", variable.name = "structure", value.name = "ASA")
df_sasa_long$structure <- str_remove(df_sasa_long$structure, "SASA_")
df_sasa_long$Pos <- as.numeric(df_sasa_long$Pos)

# Fill missing position data (6 to 37) per structure
df_sasa_long <- df_sasa_long %>% group_by(structure) %>% complete(Pos = 6:37)

# --- Load nucleation scores from INDEL data ---
load("INDEL_datasets.RData")

# Compute mean and SD of nucleation scores per position
ns_summary <- Singles.df %>%
  group_by(Pos) %>%
  summarize(
    mean = mean(nscore_c, na.rm = TRUE),
    sd = sd(nscore_c, na.rm = TRUE)
  ) %>%
  mutate(nscore_label = "NS mean")

# Merge SASA and nucleation score data
df_combined <- left_join(df_sasa_long, ns_summary, by = "Pos") %>%
  filter(!is.na(structure)) %>%
  mutate(structure = factor(structure, levels = c("7m61", "7m62", "7m64", "7m65", "6y1a", "6zrf", "8r4i", "6vw2")))

# Facet labels for structure
structure_labels <- c(
  "7m61" = "PDB: 7m61", "7m62" = "PDB: 7m62", "7m64" = "PDB: 7m64", "7m65" = "PDB: 7m65",
  "6y1a" = "PDB: 6y1a", "6zrf" = "PDB: 6zrf", "8r4i" = "PDB: 8r4i", "6vw2" = "PDB: 6vw2"
)

# --- Plot correlation between ASA and nucleation score ---
p <- ggplot(df_combined, aes(x = mean, y = ASA)) +
  geom_point() +
  facet_wrap(~structure, nrow = 8, scales = "free", labeller = labeller(structure = structure_labels)) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  xlab("Mean nucleation score") +
  ylab("Relative ASA (from PDB structure)") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.1),
    axis.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 14))

ggsave(filename = "Correlations_of_relative_ASA_and_mean_nucleation_score_IAPP_substitutions.pdf",
  plot = p,
  width = 2.5,
  height = 16)

