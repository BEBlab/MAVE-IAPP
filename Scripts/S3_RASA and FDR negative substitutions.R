
load("../INDEL_datasets.RData")
# Setup
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)
library(cowplot)

# Original sequence and x-axis labels
IAPP_wt <- "TCATQRLANFLVHSSNNFGAILSSTNVGSNTY"
IAPPseq <- unlist(strsplit(IAPP_wt, ""))
IAPPseq_pos <- paste0(IAPPseq, "\n", 6:37)

# Load and process RASA
RASA <- read.csv("merged_ss_sasa.csv")
RASA_names <- grep("SASA", colnames(RASA), value = TRUE)
df_rasa <- RASA[, c("AAPos", RASA_names)]

df_rasa$Pos <- sapply(df_rasa$AAPos, function(i) as.numeric(str_extract(i, "\\d+")))
df_rasa_melt <- melt(df_rasa, id.vars = "Pos")
df_rasa_melt$Pos <- as.numeric(df_rasa_melt$Pos)
df_rasa_melt$variable <- sub("SASA_", "", df_rasa_melt$variable)
df_rasa_melt$value <- as.numeric(as.character(df_rasa_melt$value))  

# Mean nucleation score per position
singles_mean <- Singles.df %>%
  group_by(Pos) %>%
  summarize(Mean = mean(nscore_c, na.rm = TRUE))

# Merge ASA and mean NS data
df <- left_join(singles_mean, df_rasa_melt, by = "Pos")
df <- df[!is.na(df$variable), ]
df$variable <- factor(df$variable, levels = c("7m61", "7m62", "7m64", "7m65", "6y1a", "6zrf", "8r4i", "6vw2"))

# Calculate percentage of NS- substitutions per position
table_FDR <- lapply(6:37, function(i) table(Singles.df[Singles.df$Pos == i, ]$category_fdr))
all_pos <- sapply(table_FDR, sum)
NS_neg_per <- sapply(1:length(table_FDR), function(i) (table_FDR[[i]] / all_pos[[i]])[1])
ns_per <- data.frame(Pos = 6:37, NS_per = unlist(NS_neg_per))

# Merge all into a single df for plotting
dfASA_NSper <- left_join(df, ns_per, by = "Pos")
dfASA_NSper$NS_per <- as.numeric(as.character(dfASA_NSper$NS_per))

#Plot percentage of NS- substitutions and relative ASA on all available IAPP PDB structures
prof_plot <- ggplot(dfASA_NSper, aes(x = as.numeric(Pos))) + geom_col(data = ns_per, aes(y = as.numeric(NS_per)), colour = "darkred", fill = "darkred", alpha = 0.5) 
prof_plot <- prof_plot + geom_point(aes(y = value), colour = "darkgreen")
prof_plot <- prof_plot + theme_minimal() + theme(panel.border = element_blank(),
                                                 panel.grid = element_blank(),
                                                 axis.line.y = element_line(),
                                                 axis.line.x = element_line(),
                                                 axis.text.y = element_text(size=10),
                                                 axis.text.x = element_text(size=10),
                                                 axis.title.y.left = element_text(colour = "darkred"), 
                                                 axis.title.y.right = element_text(colour = "darkgreen"))
prof_plot <- prof_plot + scale_x_continuous(breaks=seq(6, 37, 1), labels = IAPPseq_pos, expand = c(0.01, 0.01)) + geom_hline(yintercept = 0.25, linetype = "dotted", colour = "darkgreen")
prof_plot <- prof_plot + scale_y_continuous(name = "Percentage of NS- substitutions", limits = c(0, 1), sec.axis = sec_axis(~ . * 1, name = "Relative ASA (from PDB structure)"))
prof_plot <- prof_plot + xlab("IAPP WT position") + geom_hline(yintercept = 0.5, colour = "darkred", linetype = "dotted")
prof_plot
ggsave(prof_plot, filename = "Profile_ASA_NSmean_per_position.pdf", width = 10, height = 3)


