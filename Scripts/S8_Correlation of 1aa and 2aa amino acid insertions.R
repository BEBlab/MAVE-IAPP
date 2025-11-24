# Load required packages
library(ggplot2)
library(dplyr)
library(ggpubr)
library(stringr)
library(patchwork) # for combining plots

# Load data
load("All_mutants.RData")

# ---------------------------
# Filter polymerase slip data
# ---------------------------
Polslip.df <- all_ds[all_ds$Mutation_type == "Polymerase slip", ]

# Extract inserted 2 amino acids from HGVS notation
Polslip.df$ins_2aa <- sapply(Polslip.df$HGVS_nomenclature, function(x) {
  parts <- strsplit(x, "(?i)ins")[[1]]
  if(length(parts) < 2) return(NA)
  after_ins <- strsplit(parts[2], ";")[[1]][1]
  after_ins
})

# Extract insertion position from HGVS
Polslip.df$ins_pos <- sapply(Polslip.df$HGVS_nomenclature, function(x) {
  match <- regexpr("\\d+(?=_)", x, perl = TRUE)
  if(match == -1) return(NA)
  regmatches(x, match)
})

# Categorize insertion type: same AA repeated vs different AAs
Polslip.df$ins_category <- sapply(Polslip.df$ins_2aa, function(x) {
  if(is.na(x) || x == "") return(NA)
  aas <- unlist(regmatches(x, gregexpr("[A-Z][a-z]{2}", x)))
  if(length(aas) != 2) return(NA)
  if(aas[1] == aas[2]) "Same AA repeated" else "Different AAs"
})

# Split into direct repeats and tandem repeats
polslips_direct <- Polslip.df %>% filter(ins_category == "Same AA repeated")
polslips_tandem <- Polslip.df %>% filter(ins_category == "Different AAs")

# Extract first and second AA for tandem repeats
polslips_direct$ins_aa_first <- substr(polslips_direct$ins_2aa, 1, 3)
polslips_tandem$ins_aa_first <- substr(polslips_tandem$ins_2aa, 1, 3)
polslips_tandem$ins_aa_second <- substr(polslips_tandem$ins_2aa, 4, 6)

# Prepare IDs for joining with single insertions
df_pol_direct <- polslips_direct %>%
  select(ins_aa_first, fitness_merged, sigma_scaled, ins_pos) %>%
  mutate(id_comb = paste(ins_pos, ins_aa_first, sep = "")) %>%
  mutate(ins_pos = as.numeric(ins_pos))

df_pol_tandem <- polslips_tandem %>%
  select(ins_aa_first, ins_aa_second, fitness_merged, sigma_scaled, ins_pos) %>%
  mutate(id_comb_first = paste(ins_pos, ins_aa_first, sep = ""),
         id_comb_second = paste(ins_pos, ins_aa_second, sep = ""))

# ---------------------------
# Prepare single aa insertions
# ---------------------------
SingleInsertions.df <- all_ds[all_ds$Mutation_type == "Insertions", ]

SingleInsertions.df$ins_pos <- sapply(SingleInsertions.df$HGVS_nomenclature, function(i) as.numeric(str_extract(i, "\\d+")))
SingleInsertions.df$ins_2aa <- sapply(SingleInsertions.df$HGVS_nomenclature, function(x) {
  parts <- strsplit(x, "(?i)ins")[[1]]
  if(length(parts) < 2) return(NA)
  strsplit(parts[2], ";")[[1]][1]
})

# Prepare IDs for joining
SingleInsertions.df <- SingleInsertions.df %>%
  mutate(
    id_comb_after = paste(ins_pos, ins_2aa, sep = ""),
    id_comb_before = paste(ins_pos + 1, ins_2aa, sep = ""),
    ins_pos = as.numeric(ins_pos)
  )

df_ins <- SingleInsertions.df %>%
  select(fitness_merged, sigma_scaled, ins_pos, id_comb_after, id_comb_before, ins_2aa)

# ---------------------------
# Merge datasets on matching insertion site
# ---------------------------
df_pol_ins_direct <- inner_join(df_ins, df_pol_direct, by = c("id_comb_after" = "id_comb", "ins_pos"))

df_pol_tandemfirst_ins <- inner_join(df_ins, df_pol_tandem, by = c("id_comb_after" = "id_comb_first"))
df_pol_tandemsecond_ins <- inner_join(df_ins, df_pol_tandem, by = c("id_comb_after" = "id_comb_second"))

# ---------------------------
# Plotting function
# ---------------------------
plot_correlation <- function(df, x_col, y_col, x_label, y_label, xlim_vals = c(-6,4), ylim_vals = c(-6,4)) {
  ggplot(df, aes_string(x = x_col, y = y_col)) +
    geom_errorbar(aes_string(xmin = paste0(x_col, " - sigma_scaled.x"),
                             xmax = paste0(x_col, " + sigma_scaled.x")), width = 0, linewidth = 0.1) +
    geom_errorbar(aes_string(ymin = paste0(y_col, " - sigma_scaled.y"),
                             ymax = paste0(y_col, " + sigma_scaled.y")), width = 0, linewidth = 0.1) +
    geom_point(alpha = 0.5) +
    xlab(x_label) +
    ylab(y_label) +
    xlim(xlim_vals) +
    ylim(ylim_vals) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line()) +
    stat_cor(label.x = Inf, label.y = Inf, vjust = 1, hjust = 1) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted")
}

# Generate plots
p_direct <- plot_correlation(df_pol_ins_direct, "fitness_merged.x", "fitness_merged.y",
                             "2 aa insertion (direct repeat)", "1 aa insertion")

p_first <- plot_correlation(df_pol_tandemfirst_ins, "fitness_merged.x", "fitness_merged.y",
                            "2 aa insertion\n(first aa of tandem repeat)", "1 aa insertion")

p_second <- plot_correlation(df_pol_tandemsecond_ins, "fitness_merged.x", "fitness_merged.y",
                             "2 aa insertion\n(second aa of tandem repeat)", "1 aa insertion")

# ---------------------------
# Combine plots into one figure
# ---------------------------
combined_plot <- p_direct + p_first + p_second + plot_layout(ncol = 3, guides = "collect")
combined_plot

# Save combined figure
ggsave("Corr_insertions_polslips_all_combined.pdf", combined_plot, width = 9, height = 3)

