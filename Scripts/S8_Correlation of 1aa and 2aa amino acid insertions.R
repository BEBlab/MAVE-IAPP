
# Load required packages
library(ggplot2)
library(dplyr)
library(ggpubr)

# Load data
load("../INDEL_datasets.RData")

# Parse and preprocess PolSlip.df
# Extract insertion position from name
PolSlip.df$ins_pos <- as.numeric(sapply(PolSlip.df$name, function(x) strsplit(x, "_")[[1]][3]))

# Extract inserted 2 amino acids
PolSlip.df$ins2aa <- apply(PolSlip.df, 1, function(row) {
  substring(text = row["aa_seq"], first = as.numeric(row["ins_pos"]) + 1, last = as.numeric(row["ins_pos"]) + 2)})

# Extract insertion type
PolSlip.df$type <- sapply(PolSlip.df$name, function(x) strsplit(x, "_")[[1]][2])

# Filter direct and tandem repeat insertions
polslips_direct <- PolSlip.df %>% filter(type == "2i")
polslips_tandem <- PolSlip.df %>% filter(type == "2c")

# Extract first aa from the 2-aa insertion
polslips_direct$ins_aa <- sapply(polslips_direct$ins2aa, function(x) strsplit(x, "")[[1]][1])

df_pol_direct <- polslips_direct %>%
  select(ins_aa, nscore_c, sigma, ins_pos) %>%
  mutate(id_comb = paste(ins_pos, ins_aa, sep = ""))

polslips_tandem <- polslips_tandem %>%
  mutate(ins_firstaa = sapply(ins2aa, function(x) strsplit(x, "")[[1]][1]),
         ins_secondaa = sapply(ins2aa, function(x) strsplit(x, "")[[1]][2]))

df_pol_tandem <- polslips_tandem %>%
  select(ins_firstaa, ins_secondaa, nscore_c, sigma, ins_pos) %>%
  mutate(id_comb_first = paste(ins_pos, ins_firstaa, sep = ""),
         id_comb_second = paste(ins_pos, ins_secondaa, sep = ""))


# Join with single insertions data

# Prepare IDs for matching
SingleInsertions.df <- SingleInsertions.df %>%
  mutate(
    id_comb_after = paste(ins_pos, ins_aa, sep = ""),
    id_comb_before = paste(as.numeric(ins_pos) + 1, ins_aa, sep = ""),
    ins_pos = as.numeric(ins_pos))

df_ins <- SingleInsertions.df %>%
  select(nscore_c, sigma, ins_pos, id_comb_after, id_comb_before, ins_aa) %>%
  rename(
    nscore_c_insertions = nscore_c,
    sigma_insertions = sigma,
    Mut = ins_aa)

# Merge datasets on matching insertion site
df_pol_ins_direct <- inner_join(df_ins, df_pol_direct, by = c("id_comb_after" = "id_comb", "ins_pos"))
df_pol_tandemfirst_ins <- inner_join(df_ins, df_pol_tandem, by = c("id_comb_after" = "id_comb_first"))
df_pol_tandemsecond_ins <- inner_join(df_ins, df_pol_tandem, by = c("id_comb_after" = "id_comb_second"))

# Correlation aa of Direct Repeat with 1 aa insertions

p_direct <- ggplot(df_pol_ins_direct, aes(x = nscore_c, y = nscore_c_insertions)) +
  geom_errorbar(aes(xmin = nscore_c - sigma, xmax = nscore_c + sigma), width = 0, linewidth = 0.1) +
  geom_errorbar(aes(ymin = nscore_c_insertions - sigma_insertions, ymax = nscore_c_insertions + sigma_insertions), width = 0, linewidth = 0.1) +
  geom_point(alpha = 0.5) +
  xlab("2 aa insertion (direct repeat)") +
  ylab("1 aa insertion") +
  xlim(-8, 3.5) +
  ylim(-8, 3.5) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  ) +
  stat_cor(label.x = Inf, label.y = Inf, vjust = 1, hjust = 1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted")

ggsave(p_direct, filename = "Corr_insertions_polslips_all_direct.pdf", width = 3, height = 3)


# Correlation first aa of Tandem Repeat with 1 aa insertions

p_first <- ggplot(df_pol_tandemfirst_ins, aes(x = nscore_c, y = nscore_c_insertions)) +
  geom_errorbar(aes(xmin = nscore_c - sigma, xmax = nscore_c + sigma), width = 0, linewidth = 0.1) +
  geom_errorbar(aes(ymin = nscore_c_insertions - sigma_insertions, ymax = nscore_c_insertions + sigma_insertions), width = 0, linewidth = 0.1) +
  geom_point(alpha = 0.5) +
  xlab("2 aa insertion\n(first aa of tandem repeat)") +
  ylab("1 aa insertion") +
  xlim(-8, 3.5) +
  ylim(-8, 3.5) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  ) +
  stat_cor(label.x = Inf, label.y = Inf, vjust = 1, hjust = 1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted")

ggsave(p_first, filename = "Corr_insertions_polslips_all_tandem_first.pdf", width = 3, height = 3)

# Correlation second aa of Tandem Repeat with 1 aa insertions

p_second <- ggplot(df_pol_tandemsecond_ins, aes(x = nscore_c, y = nscore_c_insertions)) +
  geom_errorbar(aes(xmin = nscore_c - sigma, xmax = nscore_c + sigma), width = 0, linewidth = 0.1) +
  geom_errorbar(aes(ymin = nscore_c_insertions - sigma_insertions, ymax = nscore_c_insertions + sigma_insertions), width = 0, linewidth = 0.1) +
  geom_point(alpha = 0.5) +
  xlab("2 aa insertion\n(second aa of tandem repeat)") +
  ylab("1 aa insertion") +
  xlim(-8, 3.5) + ylim(-8, 3.5) + theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()) +
  stat_cor(label.x = Inf, label.y = Inf, vjust = 1, hjust = 1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted")

ggsave(p_second, filename = "Corr_insertions_polslips_all_tandem_second.pdf", width = 3.1, height = 3)
