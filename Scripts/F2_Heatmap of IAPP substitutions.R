#--------------------------
# Load required packages
#--------------------------
library(ggtext)     # For markdown/HTML in ggplot axis labels
library(dplyr)      # Data manipulation
library(ggplot2)    # Plotting
library(tidyr)      # Data tidying
library(purrr)      # Functional programming tools (map/lapply)
library(stringr)    # String manipulation
library(Biostrings) # Handling amino acid sequences
library(readxl)     # Reading Excel files

#--------------------------
# Load datasets
#--------------------------
# This should contain the 'all_ds' dataframe with nucleation scores and variant info
load("All_mutants.RData")

#--------------------------
# Define color palettes and reference sequence
#--------------------------
colors_fdr <- rev(c("white", "darkblue", "brown3"))   # Colors for FDR categories
IAPP_wt <- "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"   # Wild-type sequence
IAPPseq <- unlist(strsplit(IAPP_wt, ""))            # Split sequence into single letters
IAPPseq_pos <- paste0(IAPPseq, "\n", 1:37)          # Labels with amino acid and position

#--------------------------
# Preprocess mutation types
#--------------------------
all_ds <- all_ds %>%
  filter(Mutation_type != "WT") %>%  # Remove wild-type rows
  mutate(
    # Recode mutation types to human-readable labels
    Mutation_type = recode(Mutation_type,
                           "Singles"                 = "Single aa\nsubstitutions",
                           "Insertions"              = "Single aa\ninsertions",
                           "Single aa deletions"     = "Single aa\ndeletions",
                           "Multi aa deletions"      = "Multiple aa\ndeletions",
                           "N-terminal truncations"  = "N-terminal\ntruncations",
                           "C-terminal truncations"  = "C-terminal\ntruncations",
                           "Polymerase slip"         = "Polymerase\nslippage\nvariants"),
    # Set factor levels to control facet order in plots
    Mutation_type = factor(Mutation_type,
                           levels = c(
                             "Single aa\nsubstitutions", "Single aa\ninsertions", "Single aa\ndeletions",
                             "Multiple aa\ndeletions", "N-terminal\ntruncations", "C-terminal\ntruncations",
                             "Polymerase\nslippage\nvariants"
                           ))
  )

#--------------------------
# Plot histograms per mutation type
#--------------------------
# Count the number of variants per mutation type for annotation
n_text <- all_ds %>%
  count(Mutation_type)

# Create histogram of nucleation scores for each mutation type
hist_plot <- ggplot(all_ds, aes(x = fitness_merged)) +
  geom_histogram(aes(fill = category_fdr), binwidth = 0.4) +
  facet_wrap(~Mutation_type, nrow = 1, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dotted") +  # Add vertical line at 0
  geom_text(data = n_text, aes(x = -4.5, y = 75, label = paste0("n = ", n)), inherit.aes = FALSE) +
  scale_fill_manual(values = colors_fdr, name = "Nucleation score\n(FDR = 0.1)") +
  scale_x_continuous(limits = c(-9, 5)) +
  scale_y_continuous(limits = c(0, 90)) +
  theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.x.top = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.line = element_line(size = 0.2),
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  ylab("Counts") +
  xlab("Nucleation score")

# Save histogram plot
ggsave(hist_plot, filename = "Nucleation_histograms_per_mut_type.pdf", width = 10, height = 3.2)

#--------------------------
# Prepare NNK substitutions data
#--------------------------
# Filter single amino acid substitutions
substitutions_NNK <- all_ds %>% filter(Mutation_type == "Single aa\nsubstitutions")
ref_seq <- "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"

# Function to extract substitution rows for a sequence
extract_substitutions_rows <- function(ref, seq) {
  ref_chars <- strsplit(ref, "")[[1]]
  seq_chars <- strsplit(seq, "")[[1]]
  len <- min(length(ref_chars), length(seq_chars))
  
  subs <- tibble()
  
  for (i in seq_len(len)) {
    if (ref_chars[i] != seq_chars[i]) {
      subs <- bind_rows(subs, tibble(
        aa_ref = ref_chars[i],
        aa_mut = seq_chars[i],
        aa_position = i
      ))
    }
  }
  
  # If no substitutions, return NA row for consistency
  if (nrow(subs) == 0) {
    subs <- tibble(aa_ref = NA, aa_mut = NA, aa_position = NA)
  }
  return(subs)
}

# Apply function to all sequences
subs_list <- lapply(substitutions_NNK$aa_seq, function(s) extract_substitutions_rows(ref_seq, s))
NNK_substitutions <- bind_rows(subs_list, .id = "seq_id") %>%
  bind_cols(substitutions_NNK)

# Complete all positions for plotting and convert factor levels for ordered y-axis
aa_levels <- rev(strsplit(as.character(Biostrings::AAString("GAVLMIFYWKRDESTCNQHP")), "")[[1]])

# Convert to character/integer for safe joining
NNK_substitutions <- NNK_substitutions %>%
  mutate(aa_mut = as.character(aa_mut),
         aa_position = as.integer(aa_position))

#--------------------------
# Mark WT-like residues with asterisks
#--------------------------
df_asterisk <- data.frame(
  aa_position = 1:37,
  aa_ref = unlist(strsplit(IAPP_wt, split = "")),
  aa_mut = unlist(strsplit(IAPP_wt, split = "")),
  fitness_merged = 0,
  category = NA
) %>%
  mutate(aa_mut = as.character(aa_mut), aa_position = as.integer(aa_position))

# Join to identify WT-like residues and mark them
keys <- c("aa_position", "aa_mut")
matches <- inner_join(df_asterisk[, keys], NNK_substitutions[, keys], by = keys)

NNK_substitutions <- NNK_substitutions %>%
  left_join(matches %>% mutate(mark = TRUE), by = keys) %>%
  mutate(
    fitness_merged = ifelse(!is.na(mark), 0, fitness_merged),
    category_fdr = ifelse(!is.na(mark), "WT-like", category_fdr)
  ) %>%
  select(-mark)

#--------------------------
# Heatmap of nucleation scores
#--------------------------
min_val <- abs(min(all_ds$fitness_merged, na.rm = TRUE))
max_val <- abs(max(all_ds$fitness_merged, na.rm = TRUE))

cols <- c(
  colorRampPalette(c("darkred", "grey95"))((min_val / (min_val + max_val) * 100) - 0.5),
  colorRampPalette("grey95")(1),
  colorRampPalette(c("grey95", "darkblue"), bias = 1)((max_val / (min_val + max_val) * 100) - 0.5)
)

# Define the desired order of amino acids (reversed so that plot goes top-to-bottom)
aa_levels <- rev(strsplit(as.character(Biostrings::AAString("GAVLMIFYWKRDESTCNQHP")), "")[[1]])

# Convert aa_mut to factor with specified levels
NNK_substitutions <- NNK_substitutions %>%
  mutate(aa_mut = factor(aa_mut, levels = aa_levels))

# Also for the asterisk dataframe (to align symbols)
df_asterisk <- df_asterisk %>%
  mutate(aa_mut = factor(aa_mut, levels = aa_levels))

heatmap_plot <- ggplot(NNK_substitutions, aes(x = aa_position, y = aa_mut, fill = as.numeric(fitness_merged))) +
  geom_tile() +
  scale_fill_gradientn(colours = cols, limits = c(-7, 4), na.value = "grey60", name = "Nucleation\nscore") +
  scale_x_continuous(breaks = 1:37, labels = IAPPseq_pos, expand = c(0,0)) +
  geom_text(data = df_asterisk, aes(x = aa_position, y = aa_mut, label = "*")) +
  ylab("Mutant amino acid") +
  xlab("") +
  theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),
    axis.ticks.x.top = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

ggsave(heatmap_plot, filename = "heatmap_IAPP_substitutions.pdf", width = 12.5, height = 6.5)

#--------------------------
# Violin plot of substitutions per position
#--------------------------
subs_median_df <- NNK_substitutions %>%
  group_by(aa_position) %>%
  summarise(median = median(fitness_merged, na.rm = TRUE))

subs_plot_df <- left_join(NNK_substitutions, subs_median_df, by = "aa_position")

p_subs_violin <- ggplot(subs_plot_df,
                        aes(x = factor(aa_position, levels = 1:37, labels = IAPPseq_pos),
                            y = fitness_merged,
                            group = factor(aa_position, levels = 1:37, labels = IAPPseq_pos))) +
  geom_hline(yintercept = 0, size = 0.2) +
  geom_violin(scale = "width", aes(fill = median), size = 0.2) +
  geom_boxplot(width = 0.15, outlier.shape = NA, size = 0.2) +
  theme_bw() +
  labs(x = "IAPP WT amino acid and position",
       y = "Nucleation score",
       fill = "Median\nnucleation\nscores") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 13)) +
  scale_fill_gradientn(colours = cols, limits = c(-7, 3.5), na.value = "grey60")

ggsave(p_subs_violin, file = "p_subs_violin_substitutions.pdf", width = 12, height = 2.5)

#--------------------------
# Glycine-specific heatmap
#--------------------------
glycine_subs <- NNK_substitutions %>%
  filter(aa_mut %in% c("P", "G")) %>%
  group_by(aa_mut) %>%
  complete(aa_position = 1:37)

glycine_plot <- ggplot(glycine_subs, aes(x = aa_position, y = aa_mut, fill = category_fdr)) +
  geom_tile(color = "gray60") +
  scale_x_continuous(breaks = 1:37, labels = IAPPseq_pos, expand = c(0,0)) +
  scale_fill_manual(values = colors_fdr, na.value = "grey60", name = "Nucleation effect") +
  geom_vline(xintercept = c(14.5, 32.5), size = 2, colour = "forestgreen") +
  ylab("Mutant\namino\nacid") +
  xlab("") +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x.top = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

ggsave(glycine_plot, filename = "glycine_substitutions_heatmap.pdf", width = 12, height = 1.4)

#--------------------------
# Panel C: Relative growth of double mutants
#--------------------------
doubles_ind <- read_excel("Low-thoughput validation.xlsx", sheet = 3)

# Set factor levels for plotting order
doubles_ind$variant <- factor(doubles_ind$variant, levels = c(
  "A25T", "K1I", "H18R", "K1I, A25T", "H18R, A25T", "F23L", "V17I", "S20G",
  "S29P", "V17I, F23L", "H18R, F23L", "S20G, F23L", "F23L, S29P"
))

# Custom labels with HTML formatting for core mutations
custom_labels <- c(
  "<span style='color:darkblue;'>A25T<br>(core mutation)</span>", "<span style='color:black'>K1I</span>",
  "<span style='color:black;'>H18R</span>", "<span style='color:black;'>K1I, A25T</span>",
  "<span style='color:black'>H18R, A25T</span>", "<span style='color:darkred;'>F23L<br>(core mutation)</span>",
  "<span style='color:black'>V17I</span>", "<span style='color:black;'>S20G</span>",
  "<span style='color:black'>S29P</span>", "<span style='color:black;'>V17I, F23L</span>",
  "<span style='color:black'>H18R, F23L</span>", "<span style='color:black'>S20G, F23L</span>",
  "<span style='color:black'>F23L, S29P</span>"
)

# Plot relative growth of single and double mutants
doubles_plot <- ggplot(doubles_ind, aes(y = Percentage, x = variant)) +
  geom_rect(aes(xmin = 0.6, xmax = 5.4, ymin = -1, ymax = 22), fill = "#b5b9d6", alpha = 0.5) +
  geom_rect(aes(xmin = 5.6, xmax = 13.4, ymin = -1, ymax = 22), fill = "#e2b0a5", alpha = 0.2) +
  geom_point() +
  scale_x_discrete(labels = custom_labels) +
  xlab("") +
  ylab("Relative growth") +
  theme_bw() +
  theme(
    axis.text.x = element_markdown(angle = 30, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = 'black', size = 0.25),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, 'cm')
  )

ggsave(doubles_plot, filename = "Double_mut_relative_growth.pdf", width = 6, height = 3)
