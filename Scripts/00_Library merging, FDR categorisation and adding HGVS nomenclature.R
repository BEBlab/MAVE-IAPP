
require(ggplot2)
require(dplyr)
require(deming)
library(Biostrings)
library(ggpubr)
library(readxl)

load("Mut_types_IAPP.RData")
load("INDEL_datasets.RData")
load("MB_IAPP_NNK_fitness_replicates.RData")
NNK_dataset <- all_variants
NNK_dataset$DS <- "NNK"
NNK_synonymous <- synonymous
NNK_synonymous$DS <- "NNK"
NNK_synonymous$Mutation_type <- "Synonymous"
load("MB_IAPP_2024_combined_fitness_replicates.RData")
Twist_dataset <- all_variants[all_variants$mean_count > 100, ]
Twist_dataset <- Twist_dataset[Twist_dataset$aa_seq %in% AllVariantsDesigned.df$aa_seq, ]
Twist_dataset$DS <- "Twist"
Twist_synonymous <- synonymous
Twist_synonymous$DS <- "Twist"
Twist_synonymous$Mutation_type <- "Synonymous"
doubles_synonymous <- synonymous


#I want to normalize each dataset to its own synonymous distribution and plot how they look before and after
NNK_dataset$Mutation_type <- "Singles"
NNK_dataset[NNK_dataset$indel == FALSE & NNK_dataset$STOP == TRUE, ]$Mutation_type <- "STOP"
NNK_dataset <- rbind(NNK_synonymous, NNK_dataset)
p_NNK <- ggplot(NNK_dataset, aes(x = fitness, group = Mutation_type, fill = Mutation_type)) + geom_histogram() + theme_bw() + xlab("Nucleation score before synonymous centering")
p_NNK <- p_NNK + scale_fill_manual(values = c("purple", "orange", "darkgreen")) + ylab("Count")
p_NNK
ggsave(p_NNK, filename = "Distribution_pre_merged_NNK.pdf", width = 3, height = 3)


Twist_dataset$Mutation_type <- "Singles, deletions and insertions"
Twist_dataset[Twist_dataset$aa_seq %in% Truncations.df$aa_seq & startsWith(x = Twist_dataset$aa_seq, prefix = "KC"), ]$Mutation_type <- "STOP"
Twist_dataset <- rbind(Twist_synonymous, Twist_dataset)
p_Twist <- ggplot(Twist_dataset, aes(x = fitness, group = Mutation_type, fill = Mutation_type)) + geom_histogram() + theme_bw() + xlab("Nucleation score before synonymous centering")
p_Twist <- p_Twist + scale_fill_manual(values = c("purple", "orange", "darkgreen")) + ylab("Count")
p_Twist
ggsave(p_Twist, filename = "Distribution_pre_merged_Twist.pdf", width = 3, height = 3)


mean_syn_NNK <- weighted.mean(NNK_synonymous$fitness, NNK_synonymous$sigma^-2, na.rm = T)
mean_syn_Twist <- weighted.mean(Twist_synonymous$fitness, Twist_synonymous$sigma^-2, na.rm = T)

NNK_dataset$fitness_merged <- NNK_dataset$fitness +(-mean_syn_NNK)
Twist_dataset$fitness_merged <- Twist_dataset$fitness +(-mean_syn_Twist)

NNK_dataset <- NNK_dataset[NNK_dataset$Mutation_type != "Synonymous", ]
Twist_dataset <- Twist_dataset[Twist_dataset$Mutation_type != "Synonymous", ]

joined_NNKTwist <- inner_join(NNK_dataset, Twist_dataset, by = "aa_seq")



p_NNKTwist <- ggplot(joined_NNKTwist, aes(x = fitness.x, y = fitness.y)) + geom_point() + stat_cor() + theme_bw() + xlab("Nucleation score of the IAPP NNK library") + ylab("Nucleation score of the IAPP Twist library")
p_NNKTwist <- p_NNKTwist + annotate("text", 
                                      x = Inf, y = -Inf,  # bottom-right corner
                                      hjust = 1.1, vjust = -0.5,
                                      label = paste0("n = ", nrow(joined_NNKTwist)))
p_NNKTwist
ggsave(p_NNKTwist, filename = "NNKTwist_correlation.pdf", width = 3, height = 3)


demmodelIAPPTwistNNK <- deming(fitness.y ~ fitness.x, data = joined_NNKTwist, xstd = sigma.x, ystd = sigma.y, model = TRUE)
Twist_dataset$fitness_merged <- Twist_dataset$fitness_merged * demmodelIAPPTwistNNK$coefficients[2] + demmodelIAPPTwistNNK$coefficients[1]
Twist_dataset$sigma_scaled <- Twist_dataset$sigma * demmodelIAPPTwistNNK$coefficients[2] 

NNK_dataset$sigma_scaled <- NNK_dataset$sigma
outer_Twist <- Twist_dataset[!(Twist_dataset$aa_seq %in% NNK_dataset$aa_seq), ]

all_ds <- rbind(NNK_dataset, outer_Twist)
insertions_out <- all_ds[nchar(all_ds$aa_seq) == 38 & grepl(pattern = "\\*", x = all_ds$aa_seq), ]
all_ds <- all_ds[!(all_ds$aa_seq %in% insertions_out$aa_seq), ]

all_ds$zscore <- all_ds$fitness_merged/all_ds$sigma_scaled
all_ds$p.adjust<-p.adjust(2*pnorm(-abs(all_ds$zscore)), method = "BH")
all_ds$sig_fdr<-FALSE
all_ds[all_ds$p.adjust<0.1,]$sig_fdr<-TRUE

all_ds$category_fdr<-"WT-like"
all_ds$aa_seq <- lapply(all_ds$aa_seq, function(i) strsplit(x = i, split = "\\*")[[1]][1])
all_ds$aa_seq <- unlist(all_ds$aa_seq)
all_ds[all_ds$sig_fdr==T & all_ds$fitness_merged < 0,]$category_fdr<-"NS_dec"
all_ds[all_ds$sig_fdr==T & all_ds$fitness_merged > 0,]$category_fdr<-"NS_inc"
all_ds[(all_ds$count_e3_s1 == 0 | all_ds$count_e2_s1 == 0 | all_ds$count_e1_s1 == 0) & all_ds$fitness_merged < 0, ]$category_fdr<-"NS_dec"
all_ds$Mutation_type <- "Multi aa deletions"
all_ds[!(startsWith(all_ds$aa_seq, prefix = "K")), ]$Mutation_type <- "N-terminal truncations"
all_ds[!(endsWith(all_ds$aa_seq, suffix = "Y")), ]$Mutation_type <- "C-terminal truncations"
all_ds[nchar(all_ds$aa_seq) == 36, ]$Mutation_type <- "Single aa deletions"
all_ds[nchar(all_ds$aa_seq) == 37 & all_ds$Nham_aa == 1, ]$Mutation_type <- "Singles"
all_ds[nchar(all_ds$aa_seq) == 37 & all_ds$Nham_aa > 1, ]$Mutation_type <- "Multiple aa substitutions"
all_ds[nchar(all_ds$aa_seq) == 38, ]$Mutation_type <- "Insertions"
all_ds[nchar(all_ds$aa_seq) == 39, ]$Mutation_type <- "Polymerase slip"
dupl_trunc <- all_ds[duplicated(all_ds$aa_seq), ]$aa_seq
all_ds_stop <- all_ds[all_ds$aa_seq %in% dupl_trunc, ]
all_ds <- all_ds[!(all_ds$aa_seq %in% dupl_trunc & all_ds$DS == "NNK"), ]
all_ds[all_ds$aa_seq == "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY", ]$Mutation_type <- "WT"
all_ds <- all_ds[all_ds$aa_seq != "", ]
# Install dependencies if needed
# install.packages("BiocManager")
# BiocManager::install(c("Biostrings", "pwalign"))


ref <- "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"

names_hgvs <- lapply(all_ds$aa_seq, function(i) generate_HGVS_aligned(ref_seq = ref, var_seq = i, mutation_type = all_ds$Mutation_type[i]))
all_ds$HGVS_nomenclature <- unlist(names_hgvs)

shift_positions <- function(hgvs_string, offset = 33) {
  # Extract all numbers
  matches <- gregexpr("\\d+", hgvs_string, perl = TRUE)
  numbers <- regmatches(hgvs_string, matches)[[1]]
  
  if (length(numbers) == 0) return(hgvs_string)
  
  # Shift numbers
  shifted <- as.character(as.integer(numbers) + offset)
  
  # Insert back into the string
  regmatches(hgvs_string, matches) <- list(shifted)
  hgvs_string
}

all_ds$HGVS_non_relative <- vapply(
  all_ds$HGVS_nomenclature,
  shift_positions,
  FUN.VALUE = character(1)
)


save(all_ds, file = "All_mutants.RData")
all_ds$HGVS_preproIAPPcontext <- all_ds$HGVS_non_relative
clean_df <- all_ds[, c("aa_seq", "HGVS_nomenclature", "HGVS_preproIAPPcontext", "Mutation_type", "mean_count", "fitness_merged", "sigma_scaled", "category_fdr")]

require(xlsx)

write.xlsx(clean_df, file = "Supplementary Table S2. IAPP nucleation scores and associated error estimates.xlsx")
