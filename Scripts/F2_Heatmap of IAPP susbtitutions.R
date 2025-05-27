
## Load required packages
require(ggplot2)
require(ggpubr)
require(stringr)
require(dplyr)
require(tidyr)
require(RColorBrewer)
require(gtools)
require(readxl)
require(ggrepel)
require(ggtext)

## Import required data 
load("INDEL_datasets.RData")
load("INDEL.df.RData")

##Define constants and aesthetics for plotting
padding <- 0.2  
custom_labels <- c("<span style='color:darkblue;'>A25T<br>(core mutation)</span>", "<span style='color:black'>K1I</span>", "<span style='color:black;'>H18R</span>", "<span style='color:black;'>K1I, A25T</span>", "<span style='color:black'>H18R, A25T</span>", "<span style='color:darkred;'>F23L<br>(core mutation)</span>", "<span style='color:black'>V17I</span>", "<span style='color:black;'>S20G</span>", "<span style='color:black'>S29P</span>", "<span style='color:black;'>V17I, F23L</span>",  "<span style='color:black'>H18R, F23L</span>", "<span style='color:black'>S20G, F23L</span>", "<span style='color:black'>F23L, S29P</span>")

##Colour legend for FDR
levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")
colors <- c("darkred", "#D45B5B", "#DC8484", "#E8C2C2","grey90", "#C3C3DE","#A1A1CF","#7E7EC0","darkblue")
myColors <- colors #Use when having all FDR levels
names(myColors) <- levels
colors_fdr <- rev(c("white", "darkblue", "darkred")) #Use when having only 3 levels

##Colour gradient
min_val <- abs(min(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))
max_val <- abs(max(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))
cols <- c(colorRampPalette(c("darkred", "grey95"))((min_val/(min_val+max_val)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max_val/(min_val+max_val)*100)-0.5))

##Colour code by amino acid type
all_aa = c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P","*")
color_axis_y <- c("black","#FEC10BFF", rep("#15983DFF",6), rep("#EE0011FF", 2), rep("#0C5BB0FF", 2),rep("#9A703EFF", 3),  rep("darkgrey", 5), "#FA6B09FF")

IAPP_wt <- "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"
IAPPseq <- unlist(strsplit(IAPP_wt, ""))
IAPPseq_pos <- paste0(IAPPseq, "\n", c(1:37))

#Parse synonymous to merge with the singles substitutions dataframe for visualisation purposes
Synonymous.df$category_fdr <- "WT-like"
Synonymous.df$sig_fdr <- FALSE
Synonymous.df$zscore <- ""
Synonymous.df$p.adjust <- ""
Synonymous.df$WT_AA <- sapply(Synonymous.df$Pos, function(i) IAPPseq[i])
Synonymous.df$Mut <- sapply(Synonymous.df$Pos, function(i) IAPPseq[i])

#Merge synonymous and single substitutions data frame. 
listdf <- list(Singles.df, Synonymous.df)
SinglesSynonymous.df <- do.call(smartbind, listdf)

#Panel A (top)

# Prepare and format mutation data and synonymous annotation
heatmap_df <- SinglesSynonymous.df %>% group_by(Mut) %>% complete(Pos = 1:37)
heatmap_df$Pos <- as.integer(as.character(heatmap_df$Pos))
heatmap_df <- heatmap_df[!(is.na(heatmap_df$Mut)), ]
aa_obj <- Biostrings::AAString("GAVLMIFYWKRDESTCNQHP")
aa_list <- unlist(strsplit(as.character(aa_obj), split = ""))
heatmap_df$Mut <- factor(heatmap_df$Mut, levels = rev(aa_list))
df_asterisk <- data.frame(as.numeric(as.character(seq(1, 37, 1))), unlist(strsplit(IAPP_wt, split = "")), unlist(strsplit(IAPP_wt, split = "")), NA, NA)
colnames(df_asterisk) <- c("Pos", "WT_AA", "Mut", "nscore_c", "category")

#Plot heatmap
pheat_snp <- ggplot(heatmap_df, aes(x = as.numeric(as.character(Pos)), y = Mut, fill = as.numeric(nscore_c))) + geom_tile() 
pheat_snp <- pheat_snp + scale_fill_gradientn(colours=cols, limits=c(-7, 3.5), na.value = "grey60", name = "Nucleation\nscore") + scale_x_continuous(breaks=seq(1:37), labels = IAPPseq_pos, expand = c(0,0))
pheat_snp <- pheat_snp + ylab("Mutant amino acid")+ xlab("")
pheat_snp <- pheat_snp + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                   axis.ticks.x.top = element_line(),       
                                   plot.title = element_blank(),
                                   panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank(),
                                   axis.text = element_text(size= 12), 
                                   axis.title = element_text(size = 14),
                                   axis.text.y = element_text(size = 14), 
                                   legend.key.size = unit(1.3,"line"),
                                   legend.justification = 1)

pheat_snp <- pheat_snp + geom_text(data = df_asterisk, aes(x = as.numeric(as.character(Pos)), y = Mut, label = "*"))
pheat_snp
ggsave(pheat_snp, filename = "Heatmap_single_aa_substitutions.pdf", height = 6.5, width = 12)

#Panel A (bottom)

#Summarise positions: calculate median 
subs_median_df <- as.data.frame(Singles.df %>% group_by(Pos) %>% dplyr::summarise(median=median(nscore_c)))
subs_plot_df <- left_join(Singles.df, subs_median_df)

#Plot distribution of nucleation scores per position 

p_subs_violin <- ggplot(subs_plot_df, aes(x=factor(Pos, levels=c(1:37), labels=IAPPseq_pos), y=nscore_c, group=factor(Pos, levels=c(1:37), labels=IAPPseq_pos)))+
  geom_hline(yintercept = 0, size=0.2)+
  geom_violin(scale = "width", aes(fill=median), size=0.2)+
  geom_boxplot(width=0.15, outlier.shape = NA, size=0.2)+
  theme_bw()+
  labs(x="IAPP WT amino acid and position", y="Nucleation score", fill="Median\nnucleation\nscores")+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size = 13))+
  scale_fill_gradientn(colours=cols, limits=c(-7, 3.5), na.value = "grey60") 
p_subs_violin
ggsave(p_subs_violin, file="Distribution of nucleation scores per position - subsitutions.pdf", width = 12, height = 2.5)


#Panel B

#Select only mutations to proline and glycine 
GP_subs <- heatmap_df[heatmap_df$Mut %in% c("P", "G"), ]
GP_subs$Mut <- factor(GP_subs$Mut, levels = c("G", "P"))

#Plot FDR categories of substitutions to proline and glycine

p_GP_subs <- ggplot(GP_subs, aes(x = Pos, y = Mut, fill = category_fdr)) + geom_tile(color = "gray60") + ylab("Mutant\namino acid") + xlab("")
p_GP_subs <- p_GP_subs + scale_x_continuous(breaks=seq(1:37), labels = IAPPseq_pos, expand = c(0,0)) 
p_GP_subs <- p_GP_subs + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                                   axis.ticks.x.top = element_line(),       
                                                   plot.title = element_blank(),
                                                   panel.grid.major = element_blank(), 
                                                   panel.grid.minor = element_blank(),
                                                   axis.text = element_text(size= 12), 
                                                   axis.title = element_text(size = 14),
                                                   axis.text.y = element_text(size = 14), 
                                                   legend.text = element_text(size = 14),
                                                   legend.title = element_text(size = 14),
                                                   legend.margin=margin(t = 0, unit='cm'),
                                                   legend.key.size = unit(0.8, 'cm'))
p_GP_subs <- p_GP_subs + scale_fill_manual(values = colors_fdr, na.value = "grey60", name = "Nucleation effect", na.translate = FALSE)
p_GP_subs <- p_GP_subs + geom_vline(xintercept = 14.5, size = 2, colour = "forestgreen")
p_GP_subs <- p_GP_subs + geom_vline(xintercept = 32.5, size = 2, colour = "forestgreen")
p_GP_subs
ggsave(p_GP_subs, file="p_glycine_line.pdf", width = 12, height = 1.6)


#Panel C - relative growth of double mutants in the core
 
#Load data table with the relative growth of double mutants 

doubles_ind <- read_excel(path = "Low-thoughput validation.xlsx", sheet = 3)
doubles_ind$variant <- factor(doubles_ind$variant, levels = c("A25T", "K1I", "H18R", "K1I, A25T", "H18R, A25T", "F23L", "V17I", "S20G", "S29P", "V17I, F23L", "H18R, F23L", "S20G, F23L", "F23L, S29P"))

#Plot the percentage of growth of each single and double mutant assayed
doubles_plot <- ggplot(doubles_ind, aes(y = Percentage, x = variant)) + 
  scale_x_discrete(labels = custom_labels) + 
  # Adjust geom_rect() so it doesn't overlap with x-axis
  geom_rect(aes(xmin = 0.6, xmax = 5.4, ymin = -1, ymax = 22), 
            fill = "#b5b9d6", alpha = 0.5) + 
  geom_rect(aes(xmin = 5.6, xmax = 13.4, ymin = -1, ymax = 22), 
            fill = "#e2b0a5", alpha = 0.2) + geom_point() + 
  xlab("") + geom_point() +
  ylab("Relative growth") + 
  theme_bw() + 
  theme(axis.text.x = element_markdown(), 
        axis.text.x.bottom = element_text(angle = 45, hjust = 1)) + 
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, 'cm'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black', size = 0.25),
        axis.title = element_text(size = 13))

ggsave(doubles_plot, filename = "Double_mut_relative_growth.pdf", width = 5, height = 3)

