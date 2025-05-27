
## Load required packages
require(ggplot2)
require(stringr)
require(ggpubr)
require(reshape2)
require(dplyr)

## Load required data
load("../INDEL.df.RData")
load("../INDEL_datasets.RData")

#Define aesthetics and constants for plotting
IAPP_wt<-"KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"
IAPPseq<-unlist(strsplit(IAPP_wt, ""))
IAPPseq_pos<-paste0(IAPPseq, "\n", c(1:37))
IAPPseq_pos_y <-paste0(IAPPseq, c(1:37))
all_aa <- c("G","A","V","L","M","I","F","Y","W","K","R","D","E","S","T","C","N","Q","H", "P", "-")

## NS color gradient 
min_val <-abs(min(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))
max_val <-abs(max(INDEL.df[!is.na(INDEL.df$nscore_c),]$nscore_c))

cols <- c(colorRampPalette(c( "brown3", "grey95"))((min_val/(min_val+max_val)*100)-0.5), colorRampPalette("grey95")(1),
          colorRampPalette(c("grey95",  "darkblue"), bias=1)((max_val/(min_val+max_val)*100)-0.5))

colors_fdr <- (c( "darkred", "darkblue", "white")) #Use when having only 3 levels

#Annotate starting and ending positions of deletions - positions are adjusted for different deletions that result in the same amino acid sequence. 
splitted <- lapply(Deletions_reps$name, function(i) strsplit(i, split = ";")[[1]])
splitted <- lapply(splitted, function(i) strsplit(i, split = "_")[[1]][3])
Deletions_reps$del_start <- lapply(splitted, function(i) strsplit(i, split = "-")[[1]][1])
Deletions_reps$del_start <- mapply(function(id, current_val) {
  if (!is.na(id) && id != "") {
    l <- strsplit(id, split = "_")[[1]][3]
    strsplit(l, split = "-")[[1]][1]
  } else {
    current_val  # keep existing del_start value
  }}, Deletions_reps$ID, Deletions_reps$del_start, SIMPLIFY = FALSE)


Deletions_reps$del_end <- lapply(splitted, function(i) strsplit(i, split = "-")[[1]][2])
Deletions_reps$del_end <- mapply(function(id, current_val) {
  if (!is.na(id) && id != "") {
    l <- strsplit(id, split = "_")[[1]][3]
    strsplit(l, split = "-")[[1]][2]
  } else {
    current_val  # keep existing del_start value
  }}, Deletions_reps$ID, Deletions_reps$del_end, SIMPLIFY = FALSE)
Deletions_reps$del_length <- as.numeric(Deletions_reps$del_end) - as.numeric(Deletions_reps$del_start)

Deletions <- Deletions_reps[,c("aa_seq","nscore_c", "sigma","p.adjust","category_fdr", "ID", "del_length", "del_start", "del_end")]

#Add truncations and single amino acid deletions 

truncations_to_add <- Truncation_reps[Truncation_reps$k_end==37 | Truncation_reps$k_start==1,
                                   c("aa_seq","nscore_c", "sigma","p.adjust","category_fdr", "ID", "k_length", "k_start", "k_end")]

truncations_to_add$del_length <- 37 - (as.numeric(truncations_to_add$k_length))
truncations_to_add$del_start <- 0
truncations_to_add[truncations_to_add$k_end==37,]$del_start <- 1
truncations_to_add[truncations_to_add$k_start==1,]$del_start<-(truncations_to_add[truncations_to_add$k_start==1,]$k_end)+1

truncations_to_add$del_end<-0
truncations_to_add[truncations_to_add$k_start==1,]$del_end<-37
truncations_to_add[truncations_to_add$k_end==37,]$del_end<-(truncations_to_add[truncations_to_add$k_end==37,]$k_start)-1

single_dels_to_add <- Single_deletions_reps
single_dels_to_add$del_length<-1
single_dels_to_add$del_start<-single_dels_to_add$del_pos
single_dels_to_add$del_end<-single_dels_to_add$del_pos

#Merge multi-aa, single amino acid deletions and truncations
deletions_map <- rbind(Deletions, 
                     truncations_to_add[,c("aa_seq","nscore_c", "sigma","p.adjust","category_fdr", "ID", "del_length", "del_start", "del_end")],
                     single_dels_to_add[,c("aa_seq","nscore_c", "sigma", "p.adjust", "category_fdr","ID", "del_length", "del_start", "del_end")])

#Build a matrix with all possible deletions that will be filled with the deletions quantified in this assay. Missing deletions will be depicted as NA.
all_dels <- expand.grid(c(1:37), c(1:37))
colnames(all_dels)<-c("del_start", "del_end")
all_dels<-all_dels[all_dels$del_start<=all_dels$del_end,]
all_dels$del_length<-all_dels$del_end-all_dels$del_start+1
all_dels<-all_dels[all_dels$del_length < 36,]
all_dels$ID<-paste0("Del_k", all_dels$del_length, "_", all_dels$del_start, "-", all_dels$del_end)

all_dels[all_dels$del_start==1 & all_dels$del_end!=1,]$ID<-paste0("kmer", 37-all_dels[all_dels$del_start==1& all_dels$del_end!=1,]$del_end, "_",
                                                                  all_dels[all_dels$del_start==1& all_dels$del_end!=1,]$del_end+1, "-37")

all_dels[all_dels$del_end==37 & all_dels$del_start!=37,]$ID<-paste0("kmer", 37-all_dels[all_dels$del_end==37 & all_dels$del_start!=37,]$del_length, 
                                                                    "_1-",37-all_dels[all_dels$del_end== 37 & all_dels$del_start!=37,]$del_length)

all_dels$paste_ID <- paste(all_dels$del_start, "-", all_dels$del_end)
deletions_map$paste_ID <- paste(deletions_map$del_start, "-", deletions_map$del_end)

missing_dels <- all_dels[!all_dels$paste_ID %in% deletions_map$paste_ID, ]
missing_dels[, c("aa_seq", "nscore_c", "sigma", "p.adjust", "category_fdr")] <- NA


heatmap_deletions_map <- rbind(deletions_map, missing_dels)


#Plot multi amino acid deletions heatmap 

p_deletions <- ggplot(heatmap_deletions_map, aes(x=factor(del_start, levels=c(1:37)), 
                                               y=factor(del_end, levels = c(1:37)), fill=nscore_c))+
  theme_bw()+
  geom_tile(color="white")+
  #geom_tile(data=heatmap_deletions_map[heatmap_deletions_map$box=="fAD",],color="black", fill=NA, size=1)+
  scale_fill_gradientn(colors=cols, limits=c(-min_val, max_val), na.value = "grey60")+
  labs(x="First deleted position", y="Last deleted position", fill="Nucleation\nscore")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_discrete(labels=IAPPseq_pos)+
  scale_y_discrete(labels=IAPPseq_pos_y, position = "left")
ggsave(p_deletions, filename = "Multi amino acid deletions heatmap.pdf", width = 9, height = 9)

#Plot multi amino acid deletions heatmap coloured by FDR category

p_deletions_fdr <- ggplot(heatmap_deletions_map, aes(x=factor(del_start, levels=c(1:37)), 
                                                 y=factor(del_end, levels = c(1:37)), fill=category_fdr))+
  theme_bw()+
  geom_tile(color="white")+
  scale_fill_manual(values = colors_fdr, na.value = "grey60", name = "Nucleation\nscore") +
  labs(x="First deleted position", y="Last deleted position", fill="Nucleation\nscore")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_discrete(labels=IAPPseq_pos)+
  scale_y_discrete(labels=IAPPseq_pos_y, position = "left")
ggsave(p_deletions_fdr, filename = "Multi amino acid deletions heatmap by FDR.pdf", width = 9, height = 9)
