
#Load required packages 
require(ggplot2)
require(dplyr)
require(ggpubr)
require(tidyr)

#Define labels for the alignment of IAPP and AB42. Alignment of the 2 sequences is performed with T-COFFEE. 
#For subsitutions
equiv_dict <- data.frame(Pos_IAPP = c(1, 2, "-", "-", "-", "-", "-", seq(3, 36, 1), "-", 37), 
                         Pos_AB = c(seq(1, 34, 1), "-", seq(35, 42, 1)), plotting_number = seq(1, 43, 1))

complete_pos <- as.character(c(seq(1, 34, 1), "-", seq(35, 42, 1)))
complete_pos_for_plotting <- seq(1, 43, 1)
#For insertions
equiv_dict_ins <- data.frame(Pos_IAPP = c(1, 2, "-", "-", "-", "-", "-", seq(3, 35, 1), "-", 36), 
                             Pos_AB = c(seq(1, 34, 1), "-", seq(35, 41, 1)), plotting_number = seq(1, 42, 1))

#Load required data and rename
#IAPP dataset
load("INDEL_datasets.RData")
deletions_reps_IAPP <- Deletions_reps
deletions.df_IAPP <- Deletions.df
insertions.df_IAPP <-SingleInsertions.df
insertions_reps_IAPP <- Insertions_reps
single_deletions.df_IAPP <- SingleDeletions.df
single_deletions_reps_IAPP <- Single_deletions_reps
singles.df_IAPP <- Singles.df

#Parse single aa substitutions and single aa insertions IAPP dataframes 
singles_IAPP <- as.data.frame(singles.df_IAPP[,c("nscore_c", "sigma", "WT_AA", "Mut", "Pos", "category_fdr")])
colnames(singles_IAPP) <- c("nscore_c_IAPP", "sigma_IAPP", "WT_AA_IAPP", "Mut", "Pos_IAPP", "category_fdr_IAPP")
singles_IAPP$SEQ <- "IAPP"
singles_IAPP$pos_plotting <- equiv_dict$plotting_number[match(singles_IAPP$Pos, equiv_dict$Pos_IAPP)]

insertions_IAPP <- as.data.frame(SingleInsertions.df[,c("nscore_c", "sigma", "ins_aa", "ins_pos", "category_fdr")])
colnames(insertions_IAPP) <- c("nscore_c_IAPP", "sigma_IAPP", "ins_aa", "ins_pos_IAPP", "category_fdr_IAPP")
insertions_IAPP$SEQ<-"IAPP"
insertions_IAPP$pos_plotting <- equiv_dict_ins$plotting_number[match(insertions_IAPP$ins_pos, equiv_dict_ins$Pos_IAPP)]

#AB42 dataset from Seuma et al. 
load("INDEL_datasets_AB.Rdata")
deletions_reps_AB <- deletions_reps
deletions.df_AB <- deletions.df
insertions.df_AB <- insertions.df
insertions_reps_AB <- insertions_reps
single_deletions.df_AB <- single_deletions.df
single_deletions_reps_AB <- single_deletions_reps
singles.df_AB <- singles.df

#Parse single aa substitutions and single aa insertions AB42 dataframes 
singles_AB <- as.data.frame(singles.df_AB[,c("nscore_c", "sigma", "WT_AA", "Mut", "Pos", "category_fdr")])
colnames(singles_AB) <- c("nscore_c_AB", "sigma_AB", "WT_AA_AB", "Mut", "Pos_AB", "category_fdr_AB")
singles_AB$Pos_WT_AB <- singles_AB$Pos
singles_AB$SEQ <- "AB"
singles_AB$pos_plotting <- equiv_dict$plotting_number[match(singles_AB$Pos_AB, equiv_dict$Pos_AB)]
singles_AB$ID_Equiv <- as.character(singles_AB$Pos_AB)

insertions_reps_AB <- as.data.frame(insertions_reps_AB[,c("nscore_c", "sigma", "ins_aa", "ins_pos", "category_fdr")])
colnames(insertions_reps_AB) <- c("nscore_c_AB", "sigma_AB", "ins_aa", "ins_pos_AB", "category_fdr_AB")
insertions_reps_AB$SEQ <- "AB"
insertions_reps_AB$pos_plotting <- equiv_dict_ins$plotting_number[match(insertions_reps_AB$ins_pos, equiv_dict_ins$Pos_AB)]

#Merge AB42 and IAPP susbtitutions dataframe
singles_IAPP$ID_Equiv <- equiv_dict$Pos_AB[match(singles_IAPP$Pos, equiv_dict$Pos_IAPP)]
singles_AB_IAPP <- merge(singles_AB, singles_IAPP, by = c("pos_plotting", "Mut", "ID_Equiv"))
singles_AB_IAPP$labelAB <- paste("AB42 = ", singles_AB_IAPP$WT_AA_AB, singles_AB_IAPP$Pos_WT_AB, sep = "")
singles_AB_IAPP$labelIAPP <- paste("IAPP = ", singles_AB_IAPP$WT_AA_IAPP, singles_AB_IAPP$Pos_IAPP, sep = "")

#Merge AB42 and IAPP insertions dataframe
insertions_AB_IAPP <- merge(insertions_reps_AB, insertions_IAPP, by = c("pos_plotting", "ins_aa"))
insertions_AB_IAPP$labelAB <- paste("AB42 = ", insertions_AB_IAPP$ins_pos_AB, sep = "")
insertions_AB_IAPP$labelIAPP <- paste("IAPP = ", insertions_AB_IAPP$ins_pos_AB, sep = "")

#Panel B: plot correlation of nucleation scores of AB42 and IAPP per aligned position - substitutions
p_corr_all <- ggplot(singles_AB_IAPP, aes(y = nscore_c_AB, x = nscore_c_IAPP)) + geom_point(color = "gray60")
p_corr_all <- p_corr_all + 
  geom_hline(yintercept = 0, linetype="dotted", color="black")+
  geom_vline(xintercept = 0,linetype="dotted", color="black")+
  stat_cor(size = 4) 
p_corr_all <- p_corr_all + ylab("Nucleation score AB42") + xlab("Nucleation score IAPP") 
p_corr_all <- p_corr_all + theme_bw() + theme(panel.grid = element_blank(), legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 12))
p_corr_all <- p_corr_all + xlim(c(-7, 3.5)) + ylim(c(-7, 3.5))
ggsave(p_corr_all, filename = "Correlation substitutions between AB42 and IAPP.pdf", width = 3, height = 3)

#Panel C: plot correlation of nucleation scores of AB42 and IAPP per aligned position - insertions
p_corr_all <- ggplot(insertions_AB_IAPP, aes(y = nscore_c_AB, x = nscore_c_IAPP)) + geom_point(color = "gray60")
p_corr_all <- p_corr_all + 
  geom_hline(yintercept = 0, linetype="dotted", color="black")+
  geom_vline(xintercept = 0,linetype="dotted", color="black")+
  stat_cor(size = 4) 
p_corr_all <- p_corr_all + ylab("Nucleation score AB42") + xlab("Nucleation score IAPP") 
p_corr_all <- p_corr_all + theme_bw() + theme(panel.grid = element_blank(), legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 12))
p_corr_all <- p_corr_all + xlim(c(-7, 3.5)) + ylim(c(-7, 3.5))
ggsave(p_corr_all, filename = "Correlation insertions between AB42 and IAPP.pdf", width = 3, height = 3)

#Panel D: Plot correlation of nucleation scores of AB42 and IAPP per aligned position - one plot per position
p_corr_pos <- ggplot(singles_AB_IAPP, aes(y = nscore_c_AB, x = nscore_c_IAPP)) + geom_point() + facet_wrap(~as.numeric(Pos_IAPP), nrow = 4)
p_corr_pos <- p_corr_pos + 
  geom_hline(yintercept = 0, linetype="dotted", color="black")+
  geom_vline(xintercept = 0,linetype="dotted", color="black")+
  stat_cor(size = 3) 
p_corr_pos <- p_corr_pos + ylab("Nucleation score AB42") + xlab("Nucleation score IAPP") 
p_corr_pos <- p_corr_pos + geom_text(aes(label=labelIAPP, x=-1, y=-2.5), color = "#F69415")
p_corr_pos <- p_corr_pos + geom_text(aes(label=labelAB, x=-1, y=-4.5), color = "#8127BB")
p_corr_pos <- p_corr_pos + theme_bw() + theme(panel.grid = element_blank(), legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 8))
ggsave(p_corr_pos, filename = "Correlation substitutions between AB42 and IAPP per position.pdf", width = 14, height = 6)

#Panel E: Plot correlation of nucleation scores of AB42 and IAPP per aligned position - one plot per position
p_corr_pos <- ggplot(insertions_AB_IAPP, aes(y = nscore_c_AB, x = nscore_c_IAPP)) + geom_point() + facet_wrap(~as.numeric(ins_pos_IAPP), nrow = 4)
p_corr_pos <- p_corr_pos + 
  geom_hline(yintercept = 0, linetype="dotted", color="black")+
  geom_vline(xintercept = 0,linetype="dotted", color="black")+
  stat_cor(size = 3) 
p_corr_pos <- p_corr_pos + ylab("Nucleation score AB42") + xlab("Nucleation score IAPP") 
p_corr_pos <- p_corr_pos + geom_text(aes(label=labelIAPP, x=-1, y=-4), color = "#F69415")
p_corr_pos <- p_corr_pos + geom_text(aes(label=labelAB, x=-1, y=-5.5), color = "#8127BB")
p_corr_pos <- p_corr_pos + theme_bw() + theme(panel.grid = element_blank(), legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 8))
p_corr_pos <- p_corr_pos + ylim(-7, 3.5) + xlim(c(-7, 3.5))
ggsave(p_corr_pos, filename = "Correlation insertions between AB42 and IAPP per position.pdf", width = 14, height = 6)

