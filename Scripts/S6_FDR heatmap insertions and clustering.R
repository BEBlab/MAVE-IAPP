
#Load required packages
require(ggplot2)
require(ggpubr)
require(stringr)
require(dplyr)
require(tidyr)
require(ggdendro)
require(reshape2)

## Import required data 
load("INDEL_datasets.RData")

#Define aesthetics and constants for plotting
IAPP_wt <- "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"
IAPPseq <- unlist(strsplit(IAPP_wt, ""))
IAPPseq_pos <- paste0(IAPPseq, "\n", c(1:37))
aa_obj <- Biostrings::AAString("GAVLMIFYWKRDESTCNQHP")
aa_list <- unlist(strsplit(as.character(aa_obj), split = ""))
all_aa = c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P","*")

levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like","NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")
colors <- c("darkred", "#D45B5B", "#DC8484", "#E8C2C2","grey90", "#C3C3DE","#A1A1CF","#7E7EC0","darkblue")
myColors <- colors #Use when having all FDR levels
names(myColors) <- levels

min_val <- abs(min(Singles.df[!is.na(Singles.df$nscore_c),]$nscore_c))
max_val <- abs(max(Singles.df[!is.na(Singles.df$nscore_c),]$nscore_c))
cols <- c(colorRampPalette(c("darkred", "grey95"))((min_val/(min_val+max_val)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max_val/(min_val+max_val)*100)-0.5))

#Panel A (top) - heatmap of the FDR categories of each amino acid insertion. 

#FDR categories testing and assignment

fdr_ins <- SingleInsertions.df[c("aa_seq", "ins_pos", "ins_aa", "nscore_c", "sigma", "p.adjust")]
fdr_ins$p.adjust <- as.numeric(fdr_ins$p.adjust)
fdr_ins$FDR_cutoff <- "WT-like"
fdr_ins$ins_pos <- as.numeric(fdr_ins$ins_pos)
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.25 & fdr_ins$nscore_c<0),]$FDR_cutoff<- "NS- 25%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.1 & fdr_ins$nscore_c<0),]$FDR_cutoff<- "NS- 10%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.05 & fdr_ins$nscore_c<0),]$FDR_cutoff<- "NS- 5%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.01 & fdr_ins$nscore_c<0),]$FDR_cutoff<- "NS- 1%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.25 & fdr_ins$nscore_c>0),]$FDR_cutoff<- "NS+ 25%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.1 & fdr_ins$nscore_c>0),]$FDR_cutoff<- "NS+ 10%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.05 & fdr_ins$nscore_c>0),]$FDR_cutoff<- "NS+ 5%"
fdr_ins[(!is.na(fdr_ins$nscore_c) & fdr_ins$p.adjust<0.01 & fdr_ins$nscore_c>0),]$FDR_cutoff<- "NS+ 1%"
fdr_ins <- fdr_ins %>% group_by(ins_aa) %>% complete(ins_pos = 1:36)
fdr_ins$ins_aa <- factor(fdr_ins$ins_aa, levels = rev(aa_list))

#Parse dataframe to plot the heatmap 

fdr_ins$FDR_cutoff <- factor(fdr_ins$FDR_cutoff, levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like", "NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%"))
fdr_ins <- fdr_ins[!(is.na(fdr_ins$FDR_cutoff)), ]
df_asterisk <- data.frame(as.numeric(as.character(seq(1, 37, 1))), unlist(strsplit(IAPP_wt, split = "")), unlist(strsplit(IAPP_wt, split = "")), NA, "WT-like")
colnames(df_asterisk) <- c("Pos", "WT_AA", "Mut", "nscore_c", "FDR_cutoff")

p1 <- ggplot(fdr_ins, aes(x = as.numeric(as.character(ins_pos)), y = ins_aa, fill = FDR_cutoff)) + geom_tile() 
p1 <- p1 + scale_fill_manual(values = myColors, na.value = "grey60", name = "Nucleation\nscore") + scale_x_continuous(breaks=seq(2:37), labels = seq(2, 37, 1), expand = c(0,0))
p1 <- p1 + ylab("Mutant amino acid")+ xlab("")
p1 <- p1 + theme(panel.border = element_blank(),
                 panel.grid = element_blank(),
                 axis.line.y = element_line(),
                 axis.text.x = element_text(size=12),
                 axis.title = element_text(size = 13),
                 axis.text.y = element_text(size = 14),
                 legend.title= element_text(face="bold")) 
p1
ggsave(p1, filename = "Heatmap of single amino acid insertions coloured by FDR.pdf", height = 8, width = 15)

#Panel A (bottom) - barplot summarizing the FDR categories of single insertions per position. 

#Group mutations per position and compute frequency of mutations in each class

categories <- fdr_ins %>% group_by(ins_pos,FDR_cutoff) %>% dplyr::summarise(Freq=n()) 
categories <- categories[!(is.na(categories$FDR_cutoff)), ]

p_categories_ins <- ggplot(categories, aes(fill=factor(FDR_cutoff, levels=levels), x=(as.numeric(ins_pos)), y=Freq)) + 
  theme_bw()+
  theme(legend.title=element_text(size=10, face="bold"),
        legend.text = element_text(size=9),
        legend.key.size = unit(0.4, 'cm'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=11),
        axis.text = element_text(size=9))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation\n(%FDR)", values=colors)+
  labs( x= "Position of inserted amino acid",y="Counts") +
  scale_x_continuous(breaks=seq(2:37), labels = seq(2, 37, 1), expand = c(0,0)) + scale_y_continuous(expand = c(0,0))

ggsave(p_categories_ins, file="Barplot of IAPP substitutions FDR categories.pdf", width = 12, height = 3)


#Panel B - Hierarchical clustering of single amino acid substitutions. 

#Matrix generation with positions as rows and mutations introduced as columns. 

clust_ins_df <- reshape(SingleInsertions.df[,c("ins_pos", "ins_aa", "nscore_c")], idvar=c('ins_pos'), timevar='ins_aa', direction='wide')
colnames(clust_ins_df)[2:length(clust_ins_df)]<- paste0(unique(SingleInsertions.df$ins_aa))
clust_ins_df$ins_pos <- as.numeric(clust_ins_df$ins_pos)
merged_clust <- clust_ins_df

#Euclidean distance calculation and plotting of the distances in a dendrogram. 

rownames(merged_clust) <- merged_clust$Pos
dendro_df <- as.dendrogram(hclust(d = dist(x = merged_clust)))
dendro_plot <- ggdendrogram(data = dendro_df)
neworder = labels(dendro_df)
ggsave(dendro_plot, filename = "Dendroorder_insertions.pdf", width = 12, height = 3)

#Parse dataframe to plot it on a heatmap
melt_df <- melt(merged_clust, id="ins_pos")

#Reorder x-axis according euclidean distances and plot the heatmap
p_heatmap <- ggplot(melt_df) +
  geom_tile(aes(x=factor(ins_pos, levels=neworder),y=factor(variable, levels=rev(all_aa)),fill=value), size=0.1, color="white")+
  theme_minimal() + theme(panel.border = element_blank(),
                          panel.grid = element_blank(),
                          axis.line.y = element_line(),
                          axis.text.x = element_text(size=12),
                          axis.title = element_text(size = 13),
                          legend.title= element_text(face="bold")) +
  scale_fill_gradientn(colours=cols, limits=c(-7,3.5), na.value = "grey60", name = "Nucleation\nscore") + xlab("IAPP amino acid position") + ylab("Mutant amino acid")

ggsave(p_heatmap, filename = "Clustered_heatmap_insertions.pdf", height = 8, width = 15)

