# Cleaned, complete script producing Panels A-E (keeps all aesthetics & labels)

# ---------------------- Libraries ---------------------- #
library(ggplot2)
library(ggpubr)    # stat_compare_means
library(stringr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(Biostrings)

# ---------------------- Load data & constants ---------------------- #
load("All_mutants.RData")         # expects objects used below (all_ds, myColors, levels, etc.)
colors_fdr <- rev(c("gray90", "darkblue", "darkred"))  # for 3-level FDR


fitness_vals <- all_ds$fitness_merged
min_val <- abs(min(fitness_vals, na.rm = TRUE))
max_val <- abs(max(fitness_vals, na.rm = TRUE))

cols <- c(
  colorRampPalette(c("darkred", "grey95"))((min_val/(min_val+max_val)*100)-0.5),
  "grey95",
  colorRampPalette(c("grey95","darkblue"))((max_val/(min_val+max_val)*100)-0.5)
)

myColors_gray <- c("white","white","black")
names(myColors_gray) <- levels

all_aa <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H","P","*")
IAPP_wt <- "-CNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY-"
IAPPseq <- strsplit(IAPP_wt, "")[[1]]
IAPPseq_pos <- paste0(IAPPseq, "\n", 1:38)

WT <- "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"  # reference used in sequence reconstruction

aa3_to_1 <- c(
  Ala="A", Arg="R", Asn="N", Asp="D",
  Cys="C", Gln="Q", Glu="E", Gly="G",
  His="H", Ile="I", Leu="L", Lys="K",
  Met="M", Phe="F", Pro="P", Ser="S",
  Thr="T", Trp="W", Tyr="Y", Val="V"
)

common_theme <- theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.x.top = element_line(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13)
  )

# ---------------------- Panel A: single amino-acid insertions ---------------------- #
SingleInsertions.df <- all_ds %>%
  filter(Mutation_type == "Insertions") %>%
  mutate(
    ins_pos = as.numeric(str_extract(HGVS_nomenclature, "\\d+")) + 1,
    ins_aa = sapply(HGVS_nomenclature, function(x) {
      raw <- strsplit(x, "(?i)ins")[[1]][2]
      raw <- strsplit(raw, ";")[[1]][1]
      aa3_to_1[raw]
    })
  )

# Extract single-letter insertions from sequences (expand to rows)
extract_insertions_rows <- function(ref, seq) {
  ref_chars <- strsplit(ref, "")[[1]]
  seq_chars <- strsplit(seq, "")[[1]]
  i <- j <- 1
  insertions <- positions <- c()
  while (j <= length(seq_chars)) {
    if (i <= length(ref_chars) && seq_chars[j] == ref_chars[i]) { i <- i + 1; j <- j + 1 }
    else { insertions <- c(insertions, seq_chars[j]); positions <- c(positions, i); j <- j + 1 }
  }
  if (length(insertions) == 0) return(data.frame(aa_insertion = NA, aa_position = NA))
  data.frame(aa_insertion = insertions, aa_position = positions)
}

inserts_NNK <- data.frame(aa_seq = SingleInsertions.df$aa_seq, stringsAsFactors = FALSE)
final_results <- bind_rows(lapply(inserts_NNK$aa_seq, function(s) extract_insertions_rows(WT, s)), .id = "seq_id")
NNK_insertions <- cbind(final_results, inserts_NNK)

loc <- NNK_insertions %>%
  group_by(aa_insertion) %>%
  complete(aa_position = 1:38) %>%
  rowwise() %>%
  mutate(aa_sequence = paste0(substr(WT, 1, aa_position - 1), aa_insertion, substr(WT, aa_position, nchar(WT)))) %>%
  ungroup()

SingleInsertions.df <- left_join(loc, SingleInsertions.df, by = c("aa_sequence" = "aa_seq"))

heat_ins_df <- SingleInsertions.df %>% group_by(aa_insertion) %>% complete(aa_position = 1:38)
heat_ins_df$Pos <- as.integer(heat_ins_df$aa_position)
heat_ins_df <- heat_ins_df[!(is.na(heat_ins_df$aa_insertion)), ]
aa_list <- strsplit(as.character(Biostrings::AAString("GAVLMIFYWKRDESTCNQHP")), "")[[1]]
heat_ins_df$ins_aa <- factor(heat_ins_df$aa_insertion, levels = rev(aa_list))

# Heatmap: fitness
p_heatmap_ins <- ggplot(heat_ins_df) +
  geom_tile(aes(factor(Pos, levels = 1:38), factor(aa_insertion, levels = rev(all_aa)), fill = fitness_merged)) +
  scale_fill_gradientn(colours = cols, limits = c(-7, 3.5), na.value = "grey60", name = "Nucleation\nscore") +
  labs(x = "Position of inserted amino acid", y = "Inserted amino acid") +
  common_theme
ggsave(p_heatmap_ins, file = "Heatmap_single_aa_insertions.pdf", width = 12, height = 6.5)

# Violin (distribution per position) — Panel B
ins_median_df <- SingleInsertions.df %>% group_by(aa_position) %>% summarise(median = median(fitness_merged, na.rm = TRUE))
ins_plot_df <- left_join(SingleInsertions.df, ins_median_df, by = "aa_position")

p_ins_violin <- ggplot(ins_plot_df, aes(x = factor(aa_position, levels = 1:38), y = fitness_merged, group = factor(aa_position))) +
  geom_hline(yintercept = 0, size = 0.4) +
  geom_violin(scale = "width", aes(fill = median), size = 0.2) +
  geom_boxplot(width = 0.15, outlier.shape = NA, size = 0.2) +
  annotate("rect", xmin = 0.5,  xmax = 10.5,  ymin = -Inf, ymax = Inf, fill = "#dec8c8", alpha = 0.4) +
  annotate("rect", xmin = 10.5, xmax = 21.5,  ymin = -Inf, ymax = Inf, fill = "#9895c2", alpha = 0.4) +
  annotate("rect", xmin = 21.5, xmax = 38.5,  ymin = -Inf, ymax = Inf, fill = "gray70", alpha = 0.4) +
  theme_bw() +
  labs(x = "Position of inserted amino acid", y = "Nucleation score", fill = "Median\nnucleation\nscores") +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13)
  ) +
  scale_fill_gradientn(colours = cols, limits = c(-7, 3.5), na.value = "grey60")
ggsave(p_ins_violin, file = "Distribution of nucleation scores per position - insertions.pdf", width = 12, height = 2.2)

# ---------------------- Panel C: distributions by region ---------------------- #
SingleInsertions.df$region <- ""
SingleInsertions.df$region[SingleInsertions.df$ins_pos %in% seq(1, 10, 1)] <- "1-10"
SingleInsertions.df$region[SingleInsertions.df$ins_pos %in% seq(11, 22, 1)] <- "11-21"
SingleInsertions.df$region[SingleInsertions.df$ins_pos %in% seq(23, 38, 1)] <- "22-38"
SingleInsertions.df$region <- factor(SingleInsertions.df$region, levels = c("1-10", "11-21", "22-38"))

p_region <- ggplot(SingleInsertions.df, aes(y = fitness_merged, x = region, colour = region)) +
  geom_violin(colour = "black", aes(fill = region), alpha = 0.4) +
  geom_jitter(aes(colour = region), alpha = 0.3) +
  theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.x.top = element_line(),
    plot.title = element_blank(),
    axis.line = element_line(linewidth = 0.2),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 13),
    legend.key.size = unit(1.3, "line"),
    legend.position = "none"
  ) +
  scale_colour_manual(values = c("#dec8c8", "#9895c2", "gray70")) +
  scale_fill_manual(values = c("#dec8c8", "#9895c2", "gray70")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Nucleation score") + xlab("") +
  ylim(c(-7, 8)) +
  stat_compare_means(comparisons = list(c("1-10", "11-21"), c("1-10", "22-38"), c("11-21", "22-38")),
                     label = "p.signif", method = "t.test", tip.length = 0, p.adjust.method = "bonferroni", step.increase = 0.15)
ggsave(p_region, filename = "Distributions of single amino acid insertions grouped by region.pdf", height = 3.5, width = 2.5)

# ---------------------- Panel D: distributions by region & amino acid type ---------------------- #
SingleInsertions.df$aa_type <- ""
SingleInsertions.df$aa_type[SingleInsertions.df$ins_aa %in% c("W", "Y", "F")] <- "Aromatic"
SingleInsertions.df$aa_type[SingleInsertions.df$ins_aa %in% c("D", "E", "R", "K")] <- "Charged"
SingleInsertions.df$aa_type[SingleInsertions.df$ins_aa %in% c("P", "G")] <- "P/G"
SingleInsertions.df$aa_type[SingleInsertions.df$ins_aa %in% c("S", "T", "Q", "N", "H")] <- "Polar"
SingleInsertions.df$aa_type[SingleInsertions.df$ins_aa %in% c("A", "I", "L", "M", "V")] <- "Aliphatic"
SingleInsertions.df$aa_type[SingleInsertions.df$ins_aa == "C"] <- "Cysteine"

p_type <- ggplot(SingleInsertions.df, aes(y = fitness_merged, x = region, colour = region)) +
  geom_violin(colour = "black", aes(fill = region), alpha = 0.4) +
  geom_jitter(aes(colour = region), alpha = 0.3) +
  theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    axis.ticks.x.top = element_line(),
    plot.title = element_blank(),
    axis.line = element_line(linewidth = 0.2),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 13),
    legend.position = "none"
  ) +
  scale_colour_manual(values = c("#dec8c8", "#9895c2", "gray70")) +
  scale_fill_manual(values = c("#dec8c8", "#9895c2", "gray70")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Nucleation score") + xlab("") +
  ylim(c(-7, 8)) +
  stat_compare_means(comparisons = list(c("1-10", "11-21"), c("1-10", "22-38"), c("11-21", "22-38")),
                     label = "p.signif", method = "t.test", tip.length = 0, p.adjust.method = "bonferroni", step.increase = 0.15) +
  facet_wrap(~aa_type, nrow = 1)
ggsave(p_type, filename = "Distributions of single amino acid insertions grouped by region and amino acid type.pdf", height = 3, width = 12.5)

# ---------------------- Panel E: polymerase slippage (direct & tandem) ---------------------- #
                                  
all_aa <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H","P","*")
IAPP_wt <- "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"
IAPPseq <- strsplit(IAPP_wt, "")[[1]]
IAPPseq_pos <- paste0(IAPPseq, "\n", 1:37)

WT <- "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"  # reference used in sequence reconstruction
PolSlip.df <- all_ds %>% filter(Mutation_type == "Polymerase slip")

split_to_triplets <- function(x) str_extract_all(x, ".{3}")[[1]]

PolSlip.df <- PolSlip.df %>%
  mutate(
    Pos = as.numeric(str_extract(HGVS_nomenclature, "\\d+")),
    ins_2aa = sapply(HGVS_nomenclature, function(i) {
      raw <- strsplit(i, "(?i)ins")[[1]][2]
      raw <- strsplit(raw, ";")[[1]][1]
      paste(aa3_to_1[split_to_triplets(raw)], collapse = "")
    })
  )

extract_2aa_insertions <- function(ref, seq) {
  ref_chars <- strsplit(ref, "")[[1]]
  seq_chars <- strsplit(seq, "")[[1]]
  i <- j <- 1
  insertion <- NULL
  position <- NULL
  while (j <= length(seq_chars)) {
    if (i <= length(ref_chars) && seq_chars[j] == ref_chars[i]) { i <- i + 1; j <- j + 1 }
    else {
      start_pos <- i
      block <- c()
      while (j <= length(seq_chars) && (i > length(ref_chars) || seq_chars[j] != ref_chars[i])) {
        block <- c(block, seq_chars[j]); j <- j + 1
      }
      if (length(block) == 2) { insertion <- paste(block, collapse = ""); position <- start_pos }
    }
  }
  if (is.null(insertion)) return(data.frame(aa_insertion = NA, aa_position = NA))
  data.frame(aa_insertion = insertion, aa_position = position)
}

# classify direct vs tandem by comparing letters
PolSlip_direct <- PolSlip.df %>% filter(substr(ins_2aa, 1, 1) == substr(ins_2aa, 2, 2)) %>% complete(Pos = 1:37)
PolSlip_tandem <- PolSlip.df %>% filter(substr(ins_2aa, 1, 1) != substr(ins_2aa, 2, 2)) %>% complete(Pos = 1:37)

# Corrections for mapping positions around the double-site (keep original behavior)
# (these replicate the original code's replacements to align positions 19/20 & 21/22)
if ("Pos" %in% names(PolSlip_direct)) {
  PolSlip_direct$y <- 1
  if (any(PolSlip_direct$Pos == 19, na.rm = TRUE) && any(PolSlip_direct$Pos == 20, na.rm = TRUE)) {
    PolSlip_direct[PolSlip_direct$Pos == 19, c("category_fdr", "fitness_merged", "ins_2aa")] <-
      PolSlip_direct[PolSlip_direct$Pos == 20, c("category_fdr", "fitness_merged", "ins_2aa")]
  }
  if (any(PolSlip_direct$Pos == 21, na.rm = TRUE) && any(PolSlip_direct$Pos == 22, na.rm = TRUE)) {
    PolSlip_direct[PolSlip_direct$Pos == 21, c("category_fdr", "fitness_merged", "ins_2aa")] <-
      PolSlip_direct[PolSlip_direct$Pos == 22, c("category_fdr", "fitness_merged", "ins_2aa")]
  }
}

PolSlip_tandem$y <- 1
# map particular positions to direct repeats' values as original code did
for (p in c(20, 22)) {
  if (p %in% PolSlip_tandem$Pos && p %in% PolSlip_direct$Pos) {
    PolSlip_tandem[PolSlip_tandem$Pos == p, c("category_fdr", "fitness_merged", "ins_2aa")] <-
      PolSlip_direct[PolSlip_direct$Pos == p, c("category_fdr", "fitness_merged", "ins_2aa")]
  }
}
PolSlip_tandem[PolSlip_tandem$Pos %in% seq(2, 37, 2), ]$y <- 3

# Plot direct repeats: fitness heatmap with labels
p_direct <- ggplot(PolSlip_direct, aes(x = Pos, fill = fitness_merged, y = y)) +
  geom_tile() +
  geom_tile(data = PolSlip_direct[PolSlip_direct$category_fdr != "WT-like", ], colour = "black", linejoin = "round", size = 0.4) +
  geom_text(aes(label = ins_2aa), size = 5) +
  theme_bw() +
  scale_fill_gradientn(colors = cols, limits = c(-min_val, max_val), name = "Nucleation\nscore", na.value = "white") +
  theme(
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_continuous(breaks = 1:37, labels = 1:37, expand = c(0, 0))
ggsave(p_direct, file = "Polymerase slippage variants (direct repeats).pdf", height = 1, width = 15)


# Plot tandem repeats: fitness
p_tandem <- ggplot(PolSlip_tandem, aes(xmin = Pos - 0.45, xmax = Pos + 1.45, ymin = y - 1, ymax = y + 1, fill = fitness_merged)) +
  geom_rect(size = 1) +
  geom_rect(data = PolSlip_tandem[PolSlip_tandem$category_fdr != "WT-like", ], colour = "black") +
  geom_text(aes(x = Pos + 0.5, y = y, label = ins_2aa, colour = category_fdr), size = 5, show.legend = FALSE) +
  theme_bw() +
  scale_fill_gradientn(colors = cols, limits = c(-min_val, max_val), name = "Nucleation\nscore", na.value = "white") +
  theme(
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_continuous(breaks = 1:37, labels = 1:37, expand = c(0, 0)) +
  scale_colour_manual(values = myColors_gray, na.value = "white")
ggsave(p_tandem, file = "Polymerase slippage variants (tandem repeats).pdf", height = 2, width = 13)


# ---------------------- Panel F: single amino-acid deletions ---------------------- #
SingleDeletions.df <- all_ds %>% filter(Mutation_type == "Single aa deletions") %>%
  mutate(del_pos = as.numeric(str_extract(HGVS_nomenclature, "\\d+")))

# Handle truncations that were encoded as full-length minus first/last AA
SingleDeletions.df$del_pos[SingleDeletions.df$aa_seq == "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNT"] <- 37
SingleDeletions.df$del_pos[SingleDeletions.df$aa_seq == "CNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"] <- 1

SingleDeletions.df <- SingleDeletions.df %>% complete(del_pos = 1:37)

# Mimic original code's replacement of missing position values with neighbor positions
replace_rows <- function(df, old_pos, new_pos) {
  df[df$del_pos == old_pos, c("category_fdr", "fitness_merged", "sigma_scaled")] <-
    df[df$del_pos == new_pos, c("category_fdr", "fitness_merged", "sigma_scaled")]
  df
}
SingleDeletions.df <- replace_rows(SingleDeletions.df, 21, 22)
SingleDeletions.df <- replace_rows(SingleDeletions.df, 28, 29)

p_singledeletions <- ggplot(SingleDeletions.df, aes(x = as.numeric(del_pos), y = fitness_merged)) +
  geom_hline(yintercept = 0, size = 0.1, color = "black") +
  geom_errorbar(aes(ymin = fitness_merged - sigma_scaled, ymax = fitness_merged + sigma_scaled), width = 0, size = 0.2) +
  geom_point(data = SingleDeletions.df %>% filter(category_fdr %in% c("NS+", "NS-")),
             aes(fill = fitness_merged), size = 5, shape = 21, stroke = 1.2) +
  geom_point(data = SingleDeletions.df %>% filter(!category_fdr %in% c("NS+", "NS-")),
             aes(fill = fitness_merged), size = 5, shape = 21, stroke = 0.2) +
  labs(x = "Position of deleted amino acid", y = "Nucleation score", fill = "Nucleation\nscore") +
  scale_y_continuous(limits = c(-7, 3.5), breaks = c(-6, -4, -2, 0, 2)) +
  scale_fill_gradientn(colors = cols, limits = c(-min_val, max_val), name = "Nucleation\nscore") +
  theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.x.top = element_line(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 15),
    axis.line.y = element_line()
  ) +
  scale_x_continuous(breaks = seq(1, 37, 1), labels = IAPPseq_pos, expand = c(0.02, 0.02)) + ylab("Nucleation score") + xlab("Position of deleted amino acid")
ggsave(p_singledeletions, file = "p_singledeletions.pdf", width = 12, height = 2.2)


