
#Load required packages
require(ggplot2)
require(dplyr)
require(reshape2)
require(scales)
require(ggpubr)
require(ggpointdensity)
require(readxl)
require(ggrepel)

#Load required datasets
load("INDEL_datasets.RData")
load("INDEL.df.RData")

#I don't want to show synonymous and WT in this plot, and I need to remove variants that are duplicated (same aa sequence different ID)
Indel_fig1_plot <- INDEL.df[!(INDEL.df$mut_type %in% c("WT", "Synonymous")), ]
short_aaseq <- lapply(Indel_fig1_plot$aa_seq, function(i) strsplit(i, split = "\\*", )[[1]][1])
Indel_fig1_plot$aa_seq <- unlist(short_aaseq)
Indel_fig1_plot <- Indel_fig1_plot[!(duplicated(Indel_fig1_plot$aa_seq)), ]

Indel_fig1_plot <- Indel_fig1_plot %>%
  mutate(mut_type = case_when(
    mut_type == "Singles" ~ "Single\nsubstitutions",
    mut_type == "Single insertion" ~ "Single\ninsertions",
    mut_type == "Single deletion" ~ "Single\ndeletions",
    mut_type == "Deletion" ~ "Multiple aa\ndeletions",
    mut_type == "PolSlip" ~ "Polymerase\nslip",
    TRUE ~ mut_type  # keep unchanged if not matched
  ))

Indel_fig1_plot <- Indel_fig1_plot[!(is.na(Indel_fig1_plot$nscore_c)), ]
Indel_fig1_plot[Indel_fig1_plot$mut_type == "STOP", ]$mut_type <- "C-terminal\ntruncations"
Indel_fig1_plot[Indel_fig1_plot$mut_type == "Truncation" & endsWith(suffix = "-37", Indel_fig1_plot$name), ]$mut_type <- "N-terminal\ntruncations"
Indel_fig1_plot[Indel_fig1_plot$mut_type == "Truncation", ]$mut_type <- "C-terminal\ntruncations"

levels = c("Single\nsubstitutions", "Single\ninsertions", "Single\ndeletions", "C-terminal\ntruncations", "N-terminal\ntruncations", "Multiple aa\ndeletions", "Polymerase\nslip")
Indel_fig1_plot$mut_type <- factor(Indel_fig1_plot$mut_type, levels = levels)

###Panel b - Boxplots of relative growth of IAPP animal sequences

animal_seq <- read_excel("Low-thoughput validation.xlsx", sheet = 2)
levels = c("Human IAPP", "Baboon IAPP", "Bear IAPP", "Rat IAPP")
animal_seq$Variant <- factor(animal_seq$Variant, levels = levels)
animal_seq <- animal_seq[animal_seq$Variant != "Cat IAPP", ]

p_animal <- ggplot(animal_seq, aes(y = Variant, x = `Percentage of growth`)) + geom_boxplot() + geom_point(alpha = 0.7)
p_animal <- p_animal + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                               axis.ticks.x.top = element_line(),       
                                               plot.title = element_blank(),
                                               panel.grid.major = element_blank(), 
                                               panel.grid.minor = element_blank(),
                                               legend.title = element_text(size=15),
                                               legend.text = element_text(size=15), 
                                               axis.text = element_text(size=12), 
                                               axis.line = element_line(linewidth = 0.2)) + xlab("Relative growth") + ylab("")

p_animal
ggsave(p_animal, filename = "Animal_sequences.pdf", width = 4, height = 2.5)

###Panel c - Correlations of growth rates assessed indivvidually and nucleation scores obtained from the high-throughput experiment. 
lowthorughput_scores <- read_excel("Low-thoughput validation.xlsx", sheet = 1)
merge_ind_ns <- left_join(lowthorughput_scores, Singles.df, by = "Mutant")
merge_ind_ns[merge_ind_ns$Mutant == "WT", ]$nscore_c <- 0

summary_stats <- merge_ind_ns %>%
  group_by(Mutant) %>%
  summarise(mean_valuens = mean(nscore_c),
            sd_valuens = sigma, 
            mean_valueper = mean(Percentage_growth),
            sd_valueper = sd(Percentage_growth))

summary_stats <- summary_stats[!(duplicated(summary_stats$Mutant)),]

p <- ggplot(summary_stats, aes(x = mean_valueper, y = mean_valuens)) + geom_smooth(method = "lm", se=F, color="gray", linetype="dashed") + geom_pointrange(aes(ymin = mean_valuens - sd_valuens, ymax = mean_valuens + sd_valuens)) 
p <- p + geom_pointrange(aes(xmin = mean_valueper - sd_valueper, xmax = mean_valueper + sd_valueper)) 
p <- p + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                 axis.ticks.x.top = element_line(),       
                                 plot.title = element_blank(),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 legend.title = element_text(size=15),
                                 legend.text = element_text(size=15), 
                                 axis.text = element_text(size=12), 
                                 axis.title = element_text(size=14), 
                                 axis.line = element_line(linewidth = 0.2)) + xlab("Relative growth") + ylab("NS from selection experiments")
p <- p + stat_cor(label.sep = "\n", label.x = 16, label.y = 0) 
p <- p + geom_label_repel(aes(label = Mutant), nudge_x = 6) + annotate(geom = "text", x = 20, y = -0.8, label = paste("n = ", length(unique(summary_stats$Mutant))))
p
ggsave(p, filename = "Indvidual validation.pdf", width = 3.2, height = 3.2)

###Panel e - Distribution of nucleation scores grouped by their ThT fluorescence behaviour. 

invitro <- read_excel("Table S1. ThT measurements of IAPP variants.xlsx", sheet = 1)
merged <- left_join(invitro, Singles.df, by = c("Mutation ID" = "Mutant"))
merged[merged$`Mutation ID` == "K1Δ", ]$nscore_c <- Truncation_reps[Truncation_reps$aa_seq == "CNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY", ]$nscore_c
merged[merged$`Mutation ID` == "K1Δ", ]$category_fdr <- Truncation_reps[Truncation_reps$aa_seq == "CNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY", ]$category_fdr
merged <- merged[!is.na(merged$category_fdr), ]
merged <- merged[!is.na(merged$category_fdr), ]
colors_fdr <- rev(c("grey", "darkblue", "brown3"))
merged[merged$`Aggregation with respect to WT` == "Increased aggregation rate", ]$`Aggregation with respect to WT` <- "Increased\n aggr. rate"
merged[merged$`Aggregation with respect to WT` == "Decreased aggregation rate", ]$`Aggregation with respect to WT` <- "Decreased\n aggr. rate"
merged[merged$`Aggregation with respect to WT` == "As WT aggregation rate", ]$`Aggregation with respect to WT` <- "As WT\n aggr. rate"

#Remove or conflicting measurements
merged <- merged[merged$Reference != "10.1038/s41467-022-28660-7", ]
merged <- merged[merged$`Mutation ID` != "N21G", ]
merged <- merged[merged$`Mutation ID` != "N21A", ]
merged <- merged[!(merged$`Mutation ID` == "H18R" & merged$Reference == "10.1016/j.biochi.2017.07.015"), ]
merged <- merged[!(merged$`Mutation ID` == "S20G" & merged$Method == "Thioflavin-T fluorescence assay in ammonium acetate"), ]

p <- ggplot(merged, aes(y = nscore_c, x = `Aggregation with respect to WT`)) + geom_boxplot(width = 0.3) + geom_hline(yintercept = 0, linetype = "dotted")
p <- p + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                 axis.ticks.x.top = element_line(),       
                                 plot.title = element_blank(),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 legend.title = element_blank(),
                                 legend.text = element_text(size=12), 
                                 axis.text = element_text(size=12),
                                 legend.position = "bottom",
                                 axis.line = element_line(linewidth = 0.2), 
                                 axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Nucleation scores") + xlab("In vitro\nbehavior") + geom_point(alpha = 0.5)
ggsave(p, filename = "Invitro_evidences.pdf", width = 2.8, height = 4)


###Panel f - Barplot of all mutants of the library grouped by their nucleation class (FDR = 0.1)

red_table <- INDEL.df[, c("category_fdr", "mut_type")]
require(reshape2)
red_table <- melt(red_table)
total_rows <- sapply(unique(factor(red_table$category_fdr, levels = c("NS+", "WT-like", "NS-"))), function(i) nrow(red_table[red_table$category_fdr == i, ]))
df_small_freq <- data.frame(unique(factor(red_table$category_fdr, levels = c("NS+", "WT-like", "NS-"))), total_rows)
df_small_freq$percentage <- df_small_freq$total_rows/sum(df_small_freq$total_rows) * 100
df_small_freq$label <- ""
colnames(df_small_freq) <- c("category_fdr", "total_rows", "percentage", "label")

p <- ggplot(df_small_freq, aes(x = label, y = percentage, fill = category_fdr, group = category_fdr)) + geom_bar(position="stack", stat = "identity") 
p <- p + theme_minimal() + theme(axis.ticks.x = element_blank(),
                                 axis.ticks.y.left = element_line(),
                                 axis.line = element_line(size = 0.2),
                                 plot.title = element_blank(),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 axis.text = element_text(size=12), 
                                 strip.text = element_text(size = 12),
                                 axis.title = element_text(size = 12),
                                 legend.title = element_blank(),
                                 legend.position = "right",
                                 legend.box.margin=margin(0, 0, 0, 0))
p <- p + ylab("Frequency (%)") + xlab("") + scale_fill_manual(values = colors_fdr) 
p
ggsave(p, filename = "FDR_categories_barplot_all_mut_type.pdf", width = 3, height = 4)

###Panel g - Distributions and barplots of each mutation type grouped by their nucleation class (FDR = 0.1)

n_text<-as.data.frame(table(Indel_fig1_plot$mut_type))
colnames(n_text)<-c("mut_type", "n")
p <- ggplot(Indel_fig1_plot, aes(x = nscore_c)) + geom_histogram(fill = "#4D8346", binwidth = 0.4) + facet_wrap(~mut_type, nrow = 1, scales = "free")
p <- p + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                 axis.ticks.x.top = element_line(),       
                                 plot.title = element_blank(),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 axis.text = element_text(size=12), 
                                 axis.line = element_line(size = 0.2),
                                 strip.text = element_text(size = 12),
                                 axis.title = element_text(size = 12))
p <- p + ylab("Counts") + xlab("Nucleation score") + geom_vline(xintercept = 0, linetype = "dotted") + 
  scale_x_continuous(limits=c(-9, 5)) + scale_y_continuous(limits=c(0, 90)) 
p <- p + geom_text(data = n_text, aes(x = -4.5, y = 75, label = paste("n = ", n, sep = "")), inherit.aes = FALSE)
p
ggsave(p, filename = "Nucleation_histograms_per_mut_type.pdf", width = 10, height = 2.5)

colors_fdr <- c("brown3", "darkblue", "grey")

Indel_fig1_plot$label <- ""
red_table <- Indel_fig1_plot[, c("category_fdr", "mut_type")]
red_table <- table(red_table)
red_table <- melt(red_table)
total_rows <- rep(sapply(unique(red_table$mut_type), function(i) sum(red_table[red_table$mut_type == i, ]$value)), each = 3)
red_table$total_rows <- total_rows
red_table$percentage <- red_table$value/red_table$total_rows * 100
red_table$label <- ""
p <- ggplot(red_table, aes(y = label, x = percentage, fill = category_fdr, group = category_fdr)) + geom_bar(position="stack", stat = "identity") + facet_wrap(~mut_type, nrow = 1, scales = "free")
p <- p + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                 axis.ticks.x.top = element_line(),
                                 axis.line.x = element_line(size = 0.2),
                                 plot.title = element_blank(),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 axis.text = element_text(size=8), 
                                 strip.text = element_text(size = 12),
                                 axis.title = element_text(size = 12),
                                 legend.title = element_blank(),
                                 legend.position = "bottom")
p <- p + xlab("Frequency (%)") + ylab("") + scale_fill_manual(values = colors_fdr) 
p
ggsave(p, filename = "FDR_categories_barplot.pdf", width = 10, height = 2)

