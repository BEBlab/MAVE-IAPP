
# Load required packages
require(ggplot2)
require(ggpubr)
require(stringr)

# Load dataset
load("../INDEL_datasets.RData")

# Assign nucleation score for synonymous mutations
Synonymous.df$nscore_c <- Synonymous.df$fitness

# Combine relevant mutation datasets
df_to_write <- rbind(Synonymous.df[, c("aa_seq", "nscore_c", "mut_type", "name")],
  Deletions.df[, c("aa_seq", "nscore_c", "mut_type", "name")],
  PolSlip.df[, c("aa_seq", "nscore_c", "mut_type", "name")],
  SingleDeletions.df[, c("aa_seq", "nscore_c", "mut_type", "name")],
  SingleInsertions.df[, c("aa_seq", "nscore_c", "mut_type", "name")],
  Singles.df[, c("aa_seq", "nscore_c", "mut_type", "name")],
  Truncations.df[, c("aa_seq", "nscore_c", "mut_type", "name")],
  Stops.df[, c("aa_seq", "nscore_c", "mut_type", "name")])

# Count cysteines per sequence (subtract 2 for relative to WT)
df_to_write$cys_count <- as.numeric(lapply(df_to_write$aa_seq, function(seq) str_count(seq, "C"))) - 2

# Remove duplicated amino acid sequences
df_to_write <- df_to_write[!duplicated(df_to_write$aa_seq), ]

# Panel A: Cysteine count plot

# Define group comparisons for significance testing
my_comparisons <- list(c("0", "1"), c("-2", "0"), c("-1", "0"), c("-2", "-1"))

# Convert cysteine count to factor
df_to_write$cys_count <- factor(df_to_write$cys_count)

# Plot violin plot for cysteine count
p_cys_count <- ggplot(df_to_write, aes(x = cys_count, y = nscore_c, group = cys_count)) +
  geom_violin() +
  geom_point(color = "gray60", shape = 21) +
  xlab(expression(Delta * " cysteines w.r.t WT")) +
  ylab("Nucleation score") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
    axis.ticks.x.top = element_line(),
    axis.line = element_line(linewidth = 0.2),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, colour = c("black", "black", "red", "black", "black")),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 14),
    legend.key.size = unit(1.3, "line"),
    legend.justification = 1,
    strip.text = element_text(size = 14)) +
  ylim(-7, 8) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0, step.increase = 0.15, method = "t.test")

ggsave(p_cys_count, filename = "Cysteine_number.png", height = 3, width = 3)

# Panel B: Cysteine spacing when the total cysteine count in the sequence = 2. 

# Subset for exactly 2 cysteines (i.e., WT-like)
cys2df <- df_to_write[df_to_write$cys_count == 0 & df_to_write$mut_type %in% c("Deletion", "Single deletion"), ]

# Compute spacing between cysteines
cys2df$cys_distance <- factor(sapply(cys2df$aa_seq, function(seq) {
  pos <- gregexpr("C", seq)[[1]]
  if (length(pos) >= 2) pos[2] - pos[1] - 1 else NA}))

# Define groups for comparisons
my_comparisons <- list(c("1", "4"), c("2", "4"), c("3", "4"))

# Plot violin plot for cysteine spacing
p_cys_dis <- ggplot(cys2df, aes(x = cys_distance, y = nscore_c, group = cys_distance)) +
  geom_violin() +
  geom_point(color = "gray60", shape = 21) +
  xlab("Residues in between cysteines") +
  ylab("Nucleation score") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.x.top = element_line(),
    axis.line = element_line(linewidth = 0.2),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, colour = c("black", "black", "black", "black", "red")),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 14),
    legend.key.size = unit(1.3, "line"),
    legend.justification = 1,
    strip.text = element_text(size = 14)
  ) +
  ylim(-7, 8) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0, step.increase = 0.15, method = "t.test")

ggsave(p_cys_dis, filename = "Cysteine_distance_n=2.pdf", height = 3, width = 3)

# Panel C: Adjacent vs. non-adjacent cysteines

# Identify sequences with adjacent cysteines ("CC")
df_to_write$cc_dupla <- unlist(lapply(df_to_write$aa_seq, function(seq) gregexpr("CC", seq)[[1]]))
df_to_write$adj_cys <- "Non-adjacent\ncysteines"
df_to_write[df_to_write$cc_dupla != -1, "adj_cys"] <- "Adjacent\ncysteines"

# Plot comparison
my_comparisons <- list(c("Adjacent\ncysteines", "Non-adjacent\ncysteines"))

p_cys_dis <- ggplot(df_to_write, aes(x = adj_cys, y = nscore_c, group = adj_cys)) +
  geom_violin() +
  geom_point(color = "gray60", shape = 21) +
  xlab("") +
  ylab("Nucleation score") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
    axis.ticks.x.top = element_line(),
    axis.line = element_line(linewidth = 0.2),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.key.size = unit(1.3, "line"),
    legend.justification = 1,
    strip.text = element_text(size = 14)) + ylim(-7, 8) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0, step.increase = 0.15, method = "t.test")

ggsave(p_cys_dis, filename = "Cysteine_distance_adj.pdf", height = 3, width = 3)
