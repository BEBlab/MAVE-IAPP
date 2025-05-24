# Load required libraries
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(seqinr)

# Load and filter gnomAD data for the specific transcript
iapp_gnomad <- read.csv("gnomAD_v4.1.0_ENSG00000121351.csv")
iapp_gnomad <- iapp_gnomad[iapp_gnomad$Transcript == "ENST00000240652.8", ]

# Extract protein position from consequence string (e.g., "p.Ser34Phe" → 34)
iapp_gnomad$Protein_position <- unlist(lapply(
  as.character(iapp_gnomad$Protein.Consequence),
  function(i) parse_number(strsplit(i, split = 'p.')[[1]][2])))

# Focus on residues 34–70 (relative to the IAPP mature form)
iapp_gnomad <- iapp_gnomad[iapp_gnomad$Protein_position > 33 & iapp_gnomad$Protein_position < 71, ]

# Extract reference and alternative amino acids from the consequence string
iapp_gnomad <- iapp_gnomad %>% mutate(
  REF_aa = str_sub(Protein.Consequence, 3, 5),
  ALT_aa = str_sub(Protein.Consequence, -3, -1))

# Adjust position relative to the amyloid core
iapp_gnomad$Protein_position_corrected <- iapp_gnomad$Protein_position - 33

# Remove unused columns and rows with missing rsIDs
iapp_gnomad <- iapp_gnomad %>%
  select(-Chromosome, -Position, -Flags) %>%
  filter(!is.na(rsIDs))

# Create mutation IDs (e.g., A1G or A1*)
# Handle missense and stop variants separately
iapp_gnomad_notstop <- iapp_gnomad[iapp_gnomad$ALT_aa != "Ter", ]
iapp_gnomad_notstop$mut_ID <- paste0(
  seqinr::a(iapp_gnomad_notstop$REF_aa),
  iapp_gnomad_notstop$Protein_position_corrected,
  seqinr::a(iapp_gnomad_notstop$ALT_aa)
)

iapp_gnomadstop <- iapp_gnomad[iapp_gnomad$ALT_aa == "Ter", ]
iapp_gnomadstop$mut_ID <- paste0(
  seqinr::a(iapp_gnomadstop$REF_aa),
  iapp_gnomadstop$Protein_position_corrected,
  "*"
)

# Merge stop and non-stop mutations back together
iapp_gnomad <- bind_rows(iapp_gnomad_notstop, iapp_gnomadstop)

# Load nucleation score datasets
load("Not_noisy_IAPP.RData")
load("INDEL_datasets.RData")

# Keep only non-noisy single variants
Singles.df <- Singles.df[Singles.df$aa_seq %in% not_noisy$aa_seq, ]
Singles.df$mut_ID <- paste0(Singles.df$WT_AA, Singles.df$Pos, Singles.df$Mut)

# Merge gnomAD data with nucleation scores
merge <- left_join(iapp_gnomad, Singles.df, by = "mut_ID") %>%
  filter(!is.na(nscore_c)) %>% mutate(label_mut = "")

# Annotate variants with appreciable allele frequency
merge$label_mut[merge$Allele.Frequency > 0.000006571] <- 
  merge$mut_ID[merge$Allele.Frequency > 0.000006571]

# Simplify and reformat category labels
merge$category_fdr <- recode(
  merge$category_fdr,
  "NS_dec" = "NS-",
  "NS_inc" = "NS+",
  "WT-like" = "WT-like")

# Define color scheme for plot
colors_fdr <- rev(c("gray60", "darkblue", "darkred"))

# Create and save plot
p_gnomad <- ggplot(merge, aes(x = nscore_c, y = Allele.Frequency, colour = category_fdr)) +
  geom_point() + geom_pointrange(aes(xmin = nscore_c-sigma, xmax = nscore_c+sigma)) + guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom", legend.title = element_blank()) +
  ylab("Allele frequency") +
  xlab("Nucleation score") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = colors_fdr, name = "Nucleation\nFDR = 0.1") +
  geom_label_repel(aes(label = label_mut), nudge_x = 0.1, max.overlaps = 100, size = 2) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  guides(color = guide_legend(ncol = 2))

ggsave(plot = p_gnomad, filename = "Supplementary Figure 15. Allele frequency and nucleation scores.pdf", width = 3, height = 4)

