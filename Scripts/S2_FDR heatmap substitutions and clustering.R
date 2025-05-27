
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


#Panel A (top) - heatmap of the FDR categories of each amino acid substitution. 

#FDR categories testing and assignment

fdr_cat_df <- Singles.df[, c("Pos", "WT_AA", "Mut", "nscore_c", "p.adjust")]
fdr_cat_df$p.adjust <- as.numeric(fdr_cat_df$p.adjust)
fdr_cat_df$FDR_cutoff <- "WT-like"
fdr_cat_df <- as.data.frame(fdr_cat_df)
fdr_cat_df[(fdr_cat_df$p.adjust<0.25 & fdr_cat_df$nscore_c<0),]$FDR_cutoff <- "NS- 25%"
fdr_cat_df[(fdr_cat_df$p.adjust<0.1 & fdr_cat_df$nscore_c<0),]$FDR_cutoff<- "NS- 10%"
fdr_cat_df[(fdr_cat_df$p.adjust<0.05 & fdr_cat_df$nscore_c<0),]$FDR_cutoff<- "NS- 5%"
fdr_cat_df[(fdr_cat_df$p.adjust<0.01 & fdr_cat_df$nscore_c<0),]$FDR_cutoff<- "NS- 1%"

fdr_cat_df[(fdr_cat_df$p.adjust<0.25 & fdr_cat_df$nscore_c>0),]$FDR_cutoff<- "NS+ 25%"
fdr_cat_df[(fdr_cat_df$p.adjust<0.1 & fdr_cat_df$nscore_c>0),]$FDR_cutoff<- "NS+ 10%"
fdr_cat_df[(fdr_cat_df$p.adjust<0.05 & fdr_cat_df$nscore_c>0),]$FDR_cutoff<- "NS+ 5%"
fdr_cat_df[(fdr_cat_df$p.adjust<0.01 & fdr_cat_df$nscore_c>0),]$FDR_cutoff<- "NS+ 1%"
fdr_cat_df <- fdr_cat_df %>% group_by(Mut) %>% complete(Pos = 1:37)
fdr_cat_df[(fdr_cat_df$WT_AA == fdr_cat_df$Mut & !(is.na(fdr_cat_df$nscore_c))),]$FDR_cutoff <- "WT-like"
fdr_cat_df$Mut <- factor(fdr_cat_df$Mut, levels = rev(aa_list))

#Parse dataframe to plot the heatmap 

fdr_cat_df$FDR_cutoff <- factor(fdr_cat_df$FDR_cutoff, levels = c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like", "NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%"))
fdr_cat_df <- fdr_cat_df[!(is.na(fdr_cat_df$FDR_cutoff)), ]
df_asterisk <- data.frame(as.numeric(as.character(seq(1, 37, 1))), unlist(strsplit(IAPP_wt, split = "")), unlist(strsplit(IAPP_wt, split = "")), NA, "WT-like")
colnames(df_asterisk) <- c("Pos", "WT_AA", "Mut", "nscore_c", "FDR_cutoff")

p1 <- ggplot(fdr_cat_df, aes(x = as.numeric(as.character(Pos)), y = Mut, fill = FDR_cutoff)) + geom_tile() 
p1 <- p1 + scale_fill_manual(values = myColors, na.value = "grey60", name = "FDR_cutoff (FDR)") + scale_x_continuous(breaks=seq(1:37), labels = IAPPseq_pos, expand = c(0,0))
p1 <- p1 + ylab("Mutant amino acid")+ xlab("")
p1 <- p1 + theme(panel.border = element_blank(),
                 panel.grid = element_blank(),
                 axis.line.y = element_line(),
                 axis.text.y = element_text(size=12),
                 axis.text.x = element_text(size=12),
                 axis.title = element_text(size = 13),
                 legend.title= element_text(face="bold")) 
p1 <- p1 + geom_text(data = df_asterisk, aes(x = as.numeric(as.character(Pos)), y = Mut, label = "*"))
p1

ggsave(p1, filename = "Heatmap of single amino acid substitutions coloured by FDR.pdf", height = 8, width = 15)

#Panel A (bottom) - barplot summarizing the FDR categories of single substitutions per position. 

#Group mutations per position and compute frequency of mutations in each class

categories <- fdr_cat_df %>% group_by(Pos,FDR_cutoff) %>% dplyr::summarise(Freq=n()) 
categories <- categories[!(is.na(categories$FDR_cutoff)), ]

p_categories_subs <- ggplot(categories, aes(x = Pos, fill=factor(FDR_cutoff, levels=levels), y=Freq)) + 
  theme_bw() + scale_x_continuous(breaks=seq(1:37), labels = IAPPseq_pos, expand = c(0,0)) +
  theme(legend.title=element_text(size=10, face="bold"),
        
        legend.text = element_text(size=9),
        legend.key.size = unit(0.4, 'cm'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color='black', size=0.25),
        axis.title = element_text(size=13),
        axis.text = element_text(size= 12))+
  geom_bar(position="stack", stat="identity", alpha=0.8, width=0.8)+
  scale_fill_manual("Nucleation\n(%FDR)", values=colors)+
  labs( x= "IAPP amino acid position",y="Counts")+
  scale_y_continuous(expand = c(0,0))
p_categories_subs

ggsave(p_categories_subs, file="Barplot of IAPP substitutions FDR categories.pdf", width = 12, height = 3)


#Panel B - Hierarchical clustering of single amino acid substitutions. 

#Matrix generation with positions as rows and mutations introduced as columns. 
clust_subs_df <- reshape(Singles.df[,c("Pos", "Mut", "nscore_c")], idvar=c("Pos"), timevar="Mut", direction="wide")
colnames(clust_subs_df)[2:length(clust_subs_df)] <- paste0(unique(Singles.df$Mut))
clust_subs_df$Pos <- as.numeric(clust_subs_df$Pos)
merged_clust <- clust_subs_df

#Euclidean distance calculation and plotting of the distances in a dendrogram. 

rownames(merged_clust) <- merged_clust$Pos
dendro_df <- as.dendrogram(hclust(d = dist(x = merged_clust)))
dendro_plot <- ggdendrogram(data = dendro_df)
neworder = labels(dendro_df)
ggsave(dendro_plot, filename = "Dendroorder_subs.pdf", width = 12, height = 3)

#Parse dataframe to plot it on a heatmap
melt_df <- melt(merged_clust, id="Pos")

#Reorder x-axis according euclidean distances and plot the heatmap

p_heatmap <- ggplot(melt_df) +
      geom_tile(aes(x=factor(Pos, levels=neworder),y=factor(variable, levels=rev(all_aa)),fill=value), size=0.1, color="white")+
      theme_minimal() + theme(panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line.y = element_line(),
            axis.text.x = element_text(size=12),
            axis.title = element_text(size = 13),
            legend.title= element_text(face="bold")) +
      scale_fill_gradientn(colours=cols, limits=c(-7,3.5), na.value = "grey60", name = "Nucleation\nscore") + xlab("IAPP amino acid position") + ylab("Mutant amino acid")

ggsave(p_heatmap, filename = "Clustered_heatmap_subsitutions.pdf", height = 8, width = 15)

