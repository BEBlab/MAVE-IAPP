

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

## Import required data 
load("INDEL_datasets.RData")
load("INDEL.df.RData")

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

#Colour of the text - polymerase slippage variants 
color_gray <- c("white", "black", "black", "black","black", "black","black","black","white")
myColors_gray <- color_gray
names(myColors_gray) <- levels

##Define amino acid order and IAPP sequence for insertions

all_aa = c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P","*")
IAPP_wt <- "CNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"
IAPPseq <- unlist(strsplit(IAPP_wt, ""))
IAPPseq_pos <- paste0(IAPPseq, "\n", c(2:37))

#Panel A (top)
SingleInsertions.df$ins_pos <- as.numeric(SingleInsertions.df$ins_pos)

# Prepare and format mutation data and synonymous annotation
heat_ins_df <- SingleInsertions.df %>% group_by(ins_aa) %>% complete(ins_pos = 1:36)
heat_ins_df$Pos <- as.integer(as.character(heat_ins_df$ins_pos))
heat_ins_df <- heat_ins_df[!(is.na(heat_ins_df$ins_aa)), ]
aa_obj <- Biostrings::AAString("GAVLMIFYWKRDESTCNQHP")
aa_list <- unlist(strsplit(as.character(aa_obj), split = ""))
heat_ins_df$ins_aa <- factor(heat_ins_df$ins_aa, levels = rev(aa_list))

p_heatmap_ins <- ggplot(heat_ins_df) + 
  geom_tile(aes(factor(ins_pos, levels=c(1:36), labels=c(2:37)), factor(ins_aa, levels=rev(all_aa)),fill=nscore_c))+
  theme_minimal() + labs(x="Position of inserted amino acid", y="Inserted amino acid", fill="Nucleation\nscore")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x.top = element_line(),       
        plot.title =element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=13, face="bold"),
        legend.text = element_text(size=13), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 15)) + 
        scale_fill_gradientn(colours=cols,  limits=c(-7, 3.5), na.value = "grey60", name = "Nucleation\nscore")

ggsave(p_heatmap_ins, file="Heatmap_single_aa_insertions.pdf",width=12, height=6.5)

#Panel A (middle)

#Summarise positions: calculate median nucleation scores of insertions per position

ins_median_df <- as.data.frame(SingleInsertions.df %>% group_by(ins_pos) %>% dplyr::summarise(median=median(nscore_c)))
ins_plot_df <- left_join(SingleInsertions.df, ins_median_df)

#Plot distribution of nucleation scores per position 

p_ins_violin <- ggplot(ins_plot_df, aes(x=factor(ins_pos, levels=c(2:37), labels=IAPPseq_pos), y=nscore_c, group=factor(ins_pos, levels=c(2:37), labels=IAPPseq_pos)))+
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

ggsave(p_ins_violin, file="Distribution of nucleation scores per position - insertions.pdf", width = 12, height = 2.5)

#Panel A (bottom)

#Select only mutations to proline and glycine 
GP_ins <- heat_ins_df[heat_ins_df$ins_aa %in% c("P", "G"), ]
GP_ins$ins_aa <- factor(GP_ins$ins_aa, levels = c("G", "P"))

#Plot FDR categories of insertions to proline and glycine

p_GP_ins <- ggplot(GP_ins, aes(x = ins_pos, y = ins_aa, fill = category_fdr)) + geom_tile(color = "gray60") + ylab("Mutant\namino acid") + xlab("")
p_GP_ins <- p_GP_ins + scale_x_continuous(breaks=seq(2:37), labels = IAPPseq_pos, expand = c(0,0)) 
p_GP_ins <- p_GP_ins + theme_minimal() + theme(axis.ticks.y=element_blank(),
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
p_GP_ins <- p_GP_ins + scale_fill_manual(values = colors_fdr, na.value = "grey60", name = "Nucleation effect", na.translate = FALSE)
p_GP_ins <- p_GP_ins + geom_vline(xintercept = 13.5, size = 2, colour = "forestgreen")
p_GP_ins <- p_GP_ins + geom_vline(xintercept = 31.5, size = 2, colour = "forestgreen")
p_GP_ins
ggsave(p_GP_ins, file="Proline and glycine FDR plot insertions.pdf", width = 12, height = 1.6)

#Panel B: polymerase slippage plot 

##Define amino acid order and IAPP sequence for polymerase slippage variants and insertions. 
IAPP_wt <- "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"
IAPPseq <- unlist(strsplit(IAPP_wt, ""))
IAPPseq_pos <- paste0(IAPPseq, "\n", c(1:37))

#Parse polymerase slipagge dataframe and obtain positions of the double insertions
PolSlip.df$Pos <- sapply(PolSlip.df$name, function(i) as.numeric(strsplit(i, "_")[[1]][3]))
PolSlip.df$ins_2aa <- apply(PolSlip.df, 1, function(i) substring(text = i["aa_seq"], first = as.numeric(i["Pos"]) + 1, last = as.numeric(i["Pos"]) + 2))

#Separate polymerase slippage variants in direct and tandem repeats and prepare dataframe for plotting
polslip_direct <- PolSlip.df[grep(x = PolSlip.df$name, pattern = "PolSlip_2i"), ]
polslip_direct <- polslip_direct %>% complete(Pos = 1:37) 
polslip_direct$y <- 1

polslip_tandem <- PolSlip.df[grep(x = PolSlip.df$name, pattern = "PolSlip_2c"), ]
polslip_tandem$y <- 1
polslip_tandem[polslip_tandem$Pos %in% seq(2, 37, 2), ]$y <- 3
polslip_tandem <- polslip_tandem %>% complete(Pos = 1:37) 

#Plot direct repeats
p_direct <- ggplot(polslip_direct, aes(x = Pos, fill = nscore_c, y = y)) + geom_tile() + geom_tile(data = polslip_direct[polslip_direct$sig_fdr == TRUE, ], colour = "black", linejoin = "round", size = 0.4) + geom_text(aes(label = ins_2aa), size = 5)
p_direct <- p_direct + theme_bw() + scale_fill_gradientn(colors=cols, limits=c(-min_val,max_val), name = "Nucleation\nscore", na.value = "white")
p_direct <- p_direct + theme(axis.ticks.y=element_blank(),
                   plot.title = element_blank(),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   legend.title = element_text(size=13, face="bold"),
                   legend.text = element_text(size=13), 
                   axis.text.x = element_text(size = 14),
                   axis.text.y = element_blank(),
                   axis.title = element_blank(), 
                   panel.border = element_blank())
p_direct <- p_direct + scale_x_continuous(breaks=c(1:37), labels = c(1:37), expand = c(0,0))
p_direct
ggsave(p_direct, file= "Polymerase slippage variants (direct repeats).pdf", height = 1, width = 15)

#Plot tandem repeats
p_tandem <- ggplot(polslip_tandem, aes(xmin = Pos - 0.45, xmax = Pos + 1.45, ymin = y - 1, ymax = y + 1, fill = nscore_c)) + geom_rect(size = 1) + geom_rect(data = polslip_tandem[polslip_tandem$sig_fdr == TRUE, ], colour = "black")
p_tandem <- p_tandem + geom_text(aes(x = Pos + 0.5, y = y, label = ins_2aa, colour = category_fdr), size = 5, show.legend = FALSE) + theme_bw() + scale_fill_gradientn(colors=cols, limits=c(-min_val,max_val), name = "Nucleation\nscore", na.value = "white")
p_tandem <- p_tandem + theme(axis.ticks.y=element_blank(),
                   plot.title = element_blank(),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   legend.title = element_text(size=13, face="bold"),
                   legend.text = element_text(size=13), 
                   axis.text.x = element_text(size = 14),
                   axis.text.y = element_blank(),
                   axis.title = element_blank(), 
                   panel.border = element_blank(),
                   legend.position = "bottom")
p_tandem <- p_tandem + scale_x_continuous(breaks=c(1:37), labels = c(1:37), expand = c(0,0))
p_tandem <- p_tandem + scale_colour_manual(values = myColors_gray, na.value = "white")
p_tandem
ggsave(p_tandem, file= "Polymerase slippage variants (tandem repeats).pdf", height= 2, width = 13)

#Panel C: single amino acid deletions 

#Parse single deletions 
SingleDeletions.df$del_pos <- ""
SingleDeletions.df$del_pos <- sapply(SingleDeletions.df$name, function(i) strsplit(i, split = ";")[[1]][1])
SingleDeletions.df$del_pos <- sapply(SingleDeletions.df$del_pos, function(i) strsplit(i, split = "-")[[1]][2])

#Annotate the single aa deletions that can be a result of different deletions - like deletion of residue in position 28 and 29. 
df <- c()
for (name in SingleDeletions.df$name) {
  if(grepl(pattern = ";", x = name)){
    row_to_bind <- SingleDeletions.df[SingleDeletions.df$name == name, ]
    row_to_bind$del_pos <- as.numeric(row_to_bind$del_pos) + 1
    df <- rbind(row_to_bind, df)
  } else {
    print("OK")
  }
}
df$del_pos <- as.numeric(df$del_pos)
SingleDeletions.df$del_pos <- as.numeric(SingleDeletions.df$del_pos)
SingleDeletions.df <- rbind(SingleDeletions.df, df)

#Include deletions of the first and the last residue, that are labelled as truncations

SingleDeletions.df <- SingleDeletions.df %>% complete(del_pos = 1:37)
SingleDeletions.df[SingleDeletions.df$del_pos == 1, ]$nscore_c <- Truncation_reps[Truncation_reps$ID == "kmer36_2-37", ]$nscore_c
SingleDeletions.df[SingleDeletions.df$del_pos == 1, ]$sigma <- Truncation_reps[Truncation_reps$ID == "kmer36_2-37", ]$sigma
SingleDeletions.df[SingleDeletions.df$del_pos == 1, ]$sig_fdr <- Truncation_reps[Truncation_reps$ID == "kmer36_2-37", ]$sig_fdr
SingleDeletions.df[SingleDeletions.df$del_pos == 1, ]$category_fdr <- Truncation_reps[Truncation_reps$ID == "kmer36_2-37", ]$category_fdr

SingleDeletions.df[SingleDeletions.df$del_pos == 37, ]$nscore_c <- Truncation_reps[Truncation_reps$ID == "kmer36_1-36", ]$nscore_c
SingleDeletions.df[SingleDeletions.df$del_pos == 37, ]$sigma <- Truncation_reps[Truncation_reps$ID == "kmer36_1-36", ]$sigma
SingleDeletions.df[SingleDeletions.df$del_pos == 37, ]$sig_fdr <- Truncation_reps[Truncation_reps$ID == "kmer36_1-36", ]$sig_fdr
SingleDeletions.df[SingleDeletions.df$del_pos == 37, ]$category_fdr <- Truncation_reps[Truncation_reps$ID == "kmer36_1-36", ]$category_fdr

#Plot nucleation scores of single amino acid deletions 

p_singledeletions <- ggplot(SingleDeletions.df, aes(x= as.numeric(x = del_pos), y=nscore_c)) +
  geom_hline(yintercept = 0, size=0.1, color="black") + geom_errorbar(aes(ymin=nscore_c-sigma, ymax=nscore_c+sigma),width=0, size=0.2) +
  geom_point(data=SingleDeletions.df[SingleDeletions.df$category_fdr %in% c("NS+", "NS-"),],
             aes(fill=nscore_c),size=5, shape=21, stroke=1.2) +
  geom_point(data=SingleDeletions.df[!SingleDeletions.df$category_fdr %in% c("NS+", "NS-"),],
             aes(fill=nscore_c),size=5, shape=21, stroke=0.2) +
  labs(x="Position of deleted amino acid", y="Nucleation score", fill="Nucleation\nscore")+
  scale_y_continuous(limits = c(-7, 3.5), breaks = c(-6, -4,-2,0,2))+
  scale_fill_gradientn(colors=cols,limits=c(-min_val, max_val), name = "Nucleation\nscore") + theme_minimal()+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x.top = element_line(),       
        plot.title =element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=13, face="bold"),
        legend.text = element_text(size=13), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 15), 
        axis.line.y = element_line()) + scale_x_continuous(breaks=seq(1, 37, 1), labels = IAPPseq_pos, expand = c(0.02,0.02))

ggsave(p_singledeletions, file="p_singledeletions.pdf", width = 12, height = 2.5)

#Annotate comparison pairs 
my_comparisons <- list(c("1-9", "10-21"), c("1-9", "22-37"), c("10-21", "22-37"))


#Panel D: distributions of the effect of single amino acid insertions grouped by region

#Annotate single amino acid insertions grouped by region
SingleInsertions.df$region <- ""
SingleInsertions.df[SingleInsertions.df$ins_pos %in% seq(1, 9, 1), ]$region <- "1-9"
SingleInsertions.df[SingleInsertions.df$ins_pos %in% seq(10, 21, 1), ]$region <- "10-21"
SingleInsertions.df[SingleInsertions.df$ins_pos %in% seq(22, 37, 1), ]$region <- "22-37"
SingleInsertions.df$region <- factor(SingleInsertions.df$region , levels=c("1-9", "10-21", "22-37"))

#Plot distributions of single amino acid insertions 
p <- ggplot(SingleInsertions.df, aes(y = nscore_c, x = region)) + geom_violin() + geom_point(colour = "gray", alpha = 0.7) 
p <- p + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                 axis.ticks.x.top = element_line(),       
                                 plot.title = element_blank(),
                                 axis.line = element_line(linewidth = 0.2),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 axis.text.x = element_text(size= 11, angle = 45, hjust = 1), 
                                 axis.title = element_text(size = 11),
                                 axis.text.y = element_text(size = 11), 
                                 legend.key.size = unit(1.3,"line"),
                                 legend.justification = 1)
p <- p + geom_hline(yintercept = 0, linetype = "dotted") + ylab("Nucleation score") + xlab("") + ylim(c(-7, 7)) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", tip.length = 0, step.increase = 0.15)
p
ggsave(p, filename = "Distributions of single amino acid insertions grouped by region.pdf", height = 3.5, width = 2.5)

#Panel E: distributions of the effect of single amino acid insertions grouped by region and amino acid type

#Annotate single amino acid insertions grouped by amino acid type
SingleInsertions.df$aa_type <- ""
SingleInsertions.df[SingleInsertions.df$ins_aa %in% c("W", "Y", "F"), ]$aa_type <- "Aromatic"
SingleInsertions.df[SingleInsertions.df$ins_aa %in% c("D", "E", "R", "K"), ]$aa_type <- "Charged"
SingleInsertions.df[SingleInsertions.df$ins_aa %in% c("P", "G"), ]$aa_type <- "P/G"
SingleInsertions.df[SingleInsertions.df$ins_aa %in% c("S", "T", "Q", "N", "H"), ]$aa_type <- "Polar"
SingleInsertions.df[SingleInsertions.df$ins_aa %in% c("A", "I", "L", "M", "V"), ]$aa_type <- "Aliphatic"
SingleInsertions.df[SingleInsertions.df$ins_aa == "C", ]$aa_type <- "Cysteine"

#Plot distributions of single amino acid insertions grouped by amino acid type

p <- ggplot(SingleInsertions.df, aes(y = nscore_c, x = region)) + geom_violin() + geom_point(colour = "gray", alpha = 0.7) 
p <- p + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                 axis.ticks.x.top = element_line(),       
                                 plot.title = element_blank(),
                                 axis.line = element_line(linewidth = 0.2),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 axis.text.x = element_text(size= 11, angle = 45, hjust = 1), 
                                 axis.title = element_text(size = 11),
                                 axis.text.y = element_text(size = 11), 
                                 legend.key.size = unit(1.3,"line"),
                                 legend.justification = 1) 
p <- p + geom_hline(yintercept = 0, linetype = "dotted") + ylab("Nucleation score") + xlab("") + ylim(c(-7, 7)) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", tip.length = 0, step.increase = 0.15)
p <- p + facet_wrap(~aa_type)

ggsave(p, filename = "Distributions of single amino acid insertions grouped by region and amino acid type.pdf", height = 3.5, width = 2.5)

#Panel F: distributions of the effect of polymerase slippage variants grouped by region

#Annotate polymerase slippage variants
PolSlip.df$region <- ""
PolSlip.df[PolSlip.df$Pos %in% seq(1, 9, 1), ]$region <- "1-9"
PolSlip.df[PolSlip.df$Pos %in% seq(10, 21, 1), ]$region <- "10-21"
PolSlip.df[PolSlip.df$Pos %in% seq(22, 37, 1), ]$region <- "22-37"
PolSlip.df$region <- factor(PolSlip.df$region , levels=c("1-9", "10-21", "22-37"))

#Plot distributions of polymerase slippage variants
p <- ggplot(PolSlip.df, aes(y = nscore_c, x = region)) + geom_violin() + geom_point(colour = "gray", alpha = 0.7) 
p <- p + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                 axis.ticks.x.top = element_line(),       
                                 plot.title = element_blank(),
                                 axis.line = element_line(linewidth = 0.2),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 axis.text.x = element_text(size= 11, angle = 45, hjust = 1), 
                                 axis.title = element_text(size = 11),
                                 axis.text.y = element_text(size = 11), 
                                 legend.key.size = unit(1.3,"line"),
                                 legend.justification = 1)
p <- p + geom_hline(yintercept = 0, linetype = "dotted") + ylab("Nucleation score") + xlab("") + ylim(c(-7, 7)) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", tip.length = 0, step.increase = 0.15)
p
ggsave(p, filename = "Distributions of polymerase slippage variants grouped by region.pdf", height = 3.5, width = 2.5)


#Panel G: distributions of the effect of single amino acid deletions grouped by region

#Annotate single amino acid deletions 
SingleDeletions.df$region <- ""
SingleDeletions.df[SingleDeletions.df$del_pos %in% seq(1, 9, 1), ]$region <- "1-9"
SingleDeletions.df[SingleDeletions.df$del_pos %in% seq(10, 21, 1), ]$region <- "10-21"
SingleDeletions.df[SingleDeletions.df$del_pos %in% seq(22, 37, 1), ]$region <- "22-37"
SingleDeletions.df$region <- factor(SingleDeletions.df$region , levels=c("1-9", "10-21", "22-37"))

#Plot distributions of single amino acid deletions 

p <- ggplot(SingleDeletions.df, aes(y = nscore_c, x = region)) + geom_violin() + geom_point(colour = "gray", alpha = 0.7) 
p <- p + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                 axis.ticks.x.top = element_line(),       
                                 plot.title = element_blank(),
                                 axis.line = element_line(linewidth = 0.2),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(),
                                 axis.text.x = element_text(size= 11, angle = 45, hjust = 1), 
                                 axis.title = element_text(size = 11),
                                 axis.text.y = element_text(size = 11), 
                                 legend.key.size = unit(1.3,"line"),
                                 legend.justification = 1)
p <- p + geom_hline(yintercept = 0, linetype = "dotted") + ylab("Nucleation score") + xlab("") + ylim(c(-7, 7)) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", tip.length = 0, step.increase = 0.15)
p
ggsave(p, filename = "Distributions of single amino acid deletions grouped by cluster.pdf", height = 3.5, width = 2.5)




