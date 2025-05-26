
library(tidyverse)
library(ggpubr)
library(reshape2)

### Panel a: Violin plot of mean nucleation scores per position grouped by residue exposure
#Residues with relative ASA > 0.25 are considered exposed.

# Load and process SASA data
sasa_file <- list.files(".", pattern = "merged")[1]
df_sasa_raw <- read.csv(sasa_file)

# Extract SASA columns and reshape to long format
sasa_cols <- grep("SASA", colnames(df_sasa_raw), value = TRUE)
df_sasa <- df_sasa_raw[, sasa_cols]
df_sasa$Pos <- as.numeric(str_extract(df_sasa_raw$AAPos, "[0-9]+"))

df_sasa_long <- melt(df_sasa, id.vars = "Pos", variable.name = "structure", value.name = "ASA") %>%
  mutate(structure = str_remove(structure, "SASA_"),
    Pos = as.numeric(Pos)) %>%
  group_by(structure) %>%
  complete(Pos = 6:37)

# Convert to wide format for correlation analysis
df_sasa_wide <- spread(df_sasa_long, key = structure, value = ASA)

# Load nucleation scores and compute per-position mean
load("../INDEL_datasets.RData")
ns_summary <- Singles.df %>%
  group_by(Pos) %>%
  summarize(Mean = mean(nscore_c, na.rm = TRUE), .groups = "drop")

NSmean_asa <- left_join(ns_summary, df_sasa_wide, by = "Pos")

# Define PDBs of interest
pdb_ids <- c("7m61", "7m62", "7m64", "7m65", "6y1a", "6zrf", "8r4i", "6vw2")

# Merge and label exposure status for non-averaged plot
df_not_av <- left_join(ns_summary, df_sasa_long, by = "Pos") %>%
  filter(!is.na(ASA)) %>%
  mutate(
    buried = if_else(ASA > 0.25, "Exposed", "Buried"),
    structure = factor(structure, levels = pdb_ids)
  )

# Generate Panel A plot
p_strand <- ggplot(df_not_av, aes(x = buried, y = Mean)) +
  geom_violin(fill = "white") +
  geom_point(color = "maroon") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~structure, nrow = 1) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  geom_signif(
    comparisons = list(c("Buried", "Exposed")),
    test = "wilcox.test",
    map_signif_level = TRUE,
    col = "black",
    size = 0.4,
    y_position = 0
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

ggsave("Mean nucleation scores in buried and exposed residues.pdf", p_strand, width = 10, height = 3)

### Panel b: 6y1a painted structure appears in this figure as an example. Scripts to reproduce it: S4.1-S4.4

### Panel c: Barplot of correlation coefficients (ASA vs NS)

# Compute correlation between ASA and NS mean for each structure
cor_results <- map_dfr(pdb_ids, function(pdb) {
  cor_test <- cor.test(NSmean_asa$Mean, NSmean_asa[[pdb]])
  tibble(
    PDB = pdb,
    R = cor_test$estimate,
    p_val = cor_test$p.value
  )
}) %>%
  mutate(
    p_sig = p_val < 0.05,
    PDB = factor(PDB, levels = rev(pdb_ids))
  )

# Generate Panel C plot
p_corr <- ggplot(cor_results, aes(x = R, y = PDB)) +
  geom_col(fill = "maroon", width = 0.5) +
  geom_text(data = filter(cor_results, p_sig), aes(label = "*"), hjust = -0.5, vjust = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dotted") +
  labs(
    x = "ASA and NS corr. coef. (R)",
    y = "PDB ID"
  ) +
  coord_cartesian(xlim = c(0, 0.6)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank(), axis.line.y = element_line(), axis.line.x = element_line(), legend.position = c(0.8, 0.2), legend.background = element_rect(fill = alpha('white', 0.5))) 
ggsave("Barplot of correlation of ASA and nucleation score", p_corr, width = 5, height = 3)
