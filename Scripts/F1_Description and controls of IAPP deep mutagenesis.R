# Load required packages
require(ggplot2)     # Plotting
require(dplyr)       # Data manipulation
require(reshape2)    # Reshaping data
require(scales)      # Scaling utilities
require(ggpubr)      # Statistical plotting helpers
require(ggpointdensity) # Density-colored scatterplots
require(readxl)      # Excel import
require(ggrepel)     # Non-overlapping labels
require(tidyr)       # Tidying data
library(stringr)     # String processing

# Load unified mutant dataset
load("All_mutants.RData")

# Extract mutation position from HGVS (first number in p.XnnY)
all_ds$Pos <- lapply(all_ds$HGVS_nomenclature, function(i) as.numeric(str_extract(i, "\\d+")))

# Extract WT and mutant 3-letter amino-acid codes
all_ds$WT_AA <- sapply(all_ds$HGVS_nomenclature, function(i) str_extract(i, "(?<=p\\.)[A-Z][a-z]{2}"))
all_ds$Mut   <- sapply(all_ds$HGVS_nomenclature, function(i) str_extract(i, "[A-Z][a-z]{2}$"))

# Map 3-letter AA codes to 1-letter codes
aa3_to_1 <- c(
  Ala="A", Arg="R", Asn="N", Asp="D",
  Cys="C", Gln="Q", Glu="E", Gly="G",
  His="H", Ile="I", Leu="L", Lys="K",
  Met="M", Phe="F", Pro="P", Ser="S",
  Thr="T", Trp="W", Tyr="Y", Val="V"
)

all_ds$Mut   <- aa3_to_1[all_ds$Mut]
all_ds$WT_AA <- aa3_to_1[all_ds$WT_AA]

# Construct compact mutant names (e.g., A12V)
all_ds$Mutant <- paste(all_ds$WT_AA, all_ds$Pos, all_ds$Mut, sep="")

# Keep only single-AA substitutions
singles_iapp <- all_ds[all_ds$Mutation_type == "Singles", ]

# Load in vitro validation measurements extracted from the literature
invitro <- read_excel("Figures/Supplementary Table S1. ThT measurements of IAPP variants.xlsx", sheet=1)

# Merge nucleation-scores with in vitro data
merged <- left_join(invitro, singles_iapp, by=c("Mutation ID" = "Mutant"))

# Special-case: manually assign values for K1Δ
merged[merged$`Mutation ID` == "K1Δ", ]$fitness_merged <- all_ds[all_ds$aa_seq == "CNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY", ]$fitness_merged
merged[merged$`Mutation ID` == "K1Δ", ]$category_fdr   <- all_ds[all_ds$aa_seq == "CNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY", ]$category_fdr

# Filter variants with missing category
merged <- merged[!is.na(merged$category_fdr), ]

# Color palette for categories
colors_fdr <- rev(c("grey", "darkblue", "brown3"))

# Improve display names for in vitro categories
merged[merged$`Aggregation with respect to WT` == "Increased aggregation rate", ]$`Aggregation with respect to WT` <- "Increased\n aggr. rate"
merged[merged$`Aggregation with respect to WT` == "Decreased aggregation rate", ]$`Aggregation with respect to WT` <- "Decreased\n aggr. rate"
merged[merged$`Aggregation with respect to WT` == "As WT aggregation rate", ]$`Aggregation with respect to WT` <- "As WT\n aggr. rate"

# Remove duplicates
merged <- merged[!(duplicated(merged$`Mutation ID`)), ]

# Number of total data points
n_total <- nrow(merged)

# Rename for convenience
merged$Aggregation <- merged$`Aggregation with respect to WT`

# Pairwise comparisons (Wilcoxon test + BH correction)
stat.test <- compare_means(fitness_merged ~ Aggregation,
                           data=merged,
                           method="wilcox.test",
                           p.adjust.method="BH")

# Keep only significant comparisons
stat.test.ns <- stat.test %>% filter(p.adj < 0.05)

# Assign bracket positions above data
stat.test.ns <- stat.test.ns %>%
  mutate(y.position = max(merged$fitness_merged) * 1.05 + (1:n()) * 0.1)

# Plot: nucleation-score vs in vitro phenotype
p <- ggplot(merged, aes(y=fitness_merged, x=Aggregation)) +
  geom_boxplot(width=0.3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_jitter(alpha=0.4, width=0.1) +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x.top = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "bottom",
        axis.line = element_line(linewidth=0.2)) +
  ylab("Nucleation scores") +
  xlab("In vitro\nbehavior") +
  annotate("text", x=mean(as.numeric(as.factor(merged$Aggregation))), y=2.1,
           label=paste0("n = ", n_total), vjust=1.5, size=4.2) +
  stat_pvalue_manual(stat.test.ns, label="p.signif", hide.ns=FALSE, tip.length=0.01)

p
ggsave(p, filename="Invitro_evidences.pdf", width=3.5, height=4)

###Correlation of low-throughput growth rates and high-throughput nucleation scores
lowthorughput_scores <- read_excel("Low-thoughput validation.xlsx", sheet=1)
merge_ind_ns <- left_join(lowthorughput_scores, all_ds, by="Mutant")

# Manually set WT score to zero
merge_ind_ns[merge_ind_ns$Mutant == "WT", ]$fitness_merged <- 0

# Summarize by variant
summary_stats <- merge_ind_ns %>%
  group_by(Mutant) %>%
  summarise(mean_valuens = mean(fitness_merged),
            sd_valuens = sigma_scaled,
            mean_valueper = mean(Percentage_growth),
            sd_valueper = sd(Percentage_growth)) %>%
  distinct(Mutant, .keep_all=TRUE)

# Scatter with error bars and correlation
p <- ggplot(summary_stats, aes(x=mean_valueper, y=mean_valuens)) +
  geom_smooth(method="lm", se=FALSE, color="gray", linetype="dashed") +
  geom_pointrange(aes(ymin = mean_valuens - sd_valuens,
                      ymax = mean_valuens + sd_valuens)) +
  geom_pointrange(aes(xmin = mean_valueper - sd_valueper,
                      xmax = mean_valueper + sd_valueper)) +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x.top = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=14),
        axis.text  = element_text(size=12),
        axis.line  = element_line(linewidth=0.2)) +
  xlab("Relative growth (%)\n(-URA-ADE/-URA)") +
  ylab("NS from selection\nexperiments") +
  stat_cor(label.sep="\n", label.x=16, label.y=0) +
  geom_label_repel(aes(label=Mutant), nudge_x=5) +
  annotate("text", x=20, y=-1,
           label=paste("n =", length(unique(summary_stats$Mutant))))

p
ggsave(p, filename="Indvidual validation.pdf", width=3.5, height=3.5)

### Animal-sequence panel
animal_seq <- read_excel("Low-thoughput validation.xlsx", sheet=2)
levels = c("Human IAPP", "Baboon IAPP", "Bear IAPP", "Rat IAPP")
animal_seq$Variant <- factor(animal_seq$Variant, levels=levels)
animal_seq <- animal_seq[animal_seq$Variant != "Cat IAPP", ]

animal_seq$Growth <- animal_seq$`Percentage of growth`

# Pairwise comparisons
stat.test <- compare_means(Growth ~ Variant,
                           data=animal_seq,
                           method="wilcox.test",
                           p.adjust.method="BH")

# Keep significant comparisons
stat.test.ns <- stat.test %>% filter(p.adj < 0.05) %>%
  mutate(y.position = max(animal_seq$Growth, na.rm=TRUE) * 1.05 + (1:n()) * 0.05,
         p.signif = "*") %>%
  filter(group1 != "Baboon IAPP")

# Base plot
p_animal <- ggplot(animal_seq, aes(x=Variant, y=Growth)) +
  geom_boxplot() +
  geom_point(alpha=0.7) +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x.top = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=14),
        axis.text  = element_text(size=12),
        axis.line  = element_line(linewidth=0.2)) +
  xlab("") +
  ylab("Relative growth (%)\n(-URA-ADE/-URA)") +
  coord_flip()

y_pos <- max(animal_seq$Growth, na.rm=TRUE) * 1.05

# Manual significance labels
sig_labels <- data.frame(
  Variant = unique(animal_seq$Variant),
  label = c("*", "*", "", "*"),
  y = y_pos
)

p_animal <- p_animal + geom_text(data=sig_labels, aes(x=Variant, y=y, label=label), size=6)

ggsave(p_animal, filename="Animal_sequences.pdf", width=4, height=2.5)
