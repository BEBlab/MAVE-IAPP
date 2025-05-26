
#Script to filter out all variants with not enough input reads, center nucleation scores and set the nucleation score of non-nucleating variants - variants with no reads in the output - as the mode of the stops.

#Load required packages
require(dplyr)
require(seqinr)

#Load files from DimSum pipeline
load("MAVE_IAPP_fitness_replicates.RData")
load("MAVE_IAPP_variant_data_merge.RData")
load("Designed_variants.RData")

#Filter variants with low number of reads - according to DiMSum thresholds
all_variants_df <- all_variants
variant_data_merge_df <- variant_data_merge
variant_data_merge_df_1 <- variant_data_merge_df[!(is.na(variant_data_merge_df$Nham_nt)), ]
variant_data_merge_df_1 <- variant_data_merge_df_1[variant_data_merge_df_1$Nham_nt == 1, ]
variant_data_merge_df_1 <- variant_data_merge_df_1[variant_data_merge_df_1$input1_e1_s0_bNA_count > 69 | variant_data_merge_df_1$input2_e2_s0_bNA_count > 69 | variant_data_merge_df_1$input3_e3_s0_bNA_count > 69, ]

variant_data_merge_df_2 <- variant_data_merge_df[!(is.na(variant_data_merge_df$Nham_nt)), ]
variant_data_merge_df_2 <- variant_data_merge_df_2[variant_data_merge_df_2$Nham_nt == 2, ]
variant_data_merge_df_2 <- variant_data_merge_df_2[variant_data_merge_df_2$input1_e1_s0_bNA_count > 4 | variant_data_merge_df_2$input2_e2_s0_bNA_count > 4 | variant_data_merge_df_2$input3_e3_s0_bNA_count > 4, ]

variant_data_merge_df_3 <- variant_data_merge_df[!(is.na(variant_data_merge_df$Nham_nt)), ]
variant_data_merge_df_3 <- variant_data_merge_df_3[variant_data_merge_df_3$Nham_nt == 3, ]
variant_data_merge_df_3 <- variant_data_merge_df_3[variant_data_merge_df_3$input1_e1_s0_bNA_count > 1 | variant_data_merge_df_3$input2_e2_s0_bNA_count > 1 | variant_data_merge_df_3$input3_e3_s0_bNA_count > 1, ]

variant_data_merge_df_indels <- variant_data_merge_df[variant_data_merge_df$indel == TRUE, ]
variant_data_merge_df_indels <- variant_data_merge_df_indels[variant_data_merge_df_indels$input1_e1_s0_bNA_count > 1 | variant_data_merge_df_indels$input2_e2_s0_bNA_count > 1 | variant_data_merge_df_indels$input3_e3_s0_bNA_count > 1, ]

#Remerge variant data merge including only variants with enough reads. 
Variants_filtered <- rbind(variant_data_merge_df_1, variant_data_merge_df_2, variant_data_merge_df_3, variant_data_merge_df_indels)
Variants_filtered <- Variants_filtered[Variants_filtered$aa_seq %in% AllVariantsDesigned.df$aa_seq, ]
all_variants_df <- all_variants_df[all_variants_df$aa_seq %in% Variants_filtered$aa_seq, ]

#Identify dead (non-nucleating) variants

Non_nucleating.df <- Variants_filtered[Variants_filtered$output1_e1_s1_b1_count==0  &
                               Variants_filtered$output2_e2_s1_b1_count==0  &
                               Variants_filtered$output3_e3_s1_b1_count==0 ,]

Non_nucleating.df <- as.data.frame(Non_nucleating.df[ ,c("aa_seq")])
colnames(Non_nucleating.df) <- "aa_seq"

#Centering all nucleation scores with the fitness of synonymous mutations
mean_syn_1codon <- weighted.mean(synonymous$fitness, synonymous$sigma^-2, na.rm = T)

indels_library <- rbind(Insertions.df, Deletions.df)
indels_library$nt_seq <- tolower(indels_library$nt_seq)
nt_seq_indel <- inner_join(variant_data_merge_df_indels, indels_library, by = "nt_seq")

all_variants_df %>% select(-nt_seq, -WT, -indel, -Nham_nt, -Nmut_codons)

SinglesStops.df$aa_seq <- lapply(SinglesStops.df$nt_seq, function(i) paste(seqinr::translate(strsplit(as.character(i), split = "")[[1]]), collapse = ""))
Deletions.df$aa_seq <- lapply(Deletions.df$nt_seq, function(i) paste(seqinr::translate(strsplit(as.character(i), split = "")[[1]]), collapse = ""))
Insertions.df$aa_seq <- lapply(Insertions.df$nt_seq, function(i) paste(seqinr::translate(strsplit(as.character(i), split = "")[[1]]), collapse = ""))


non_nuc_df_singles <- SinglesStops.df[SinglesStops.df$aa_seq %in% Non_nucleating.df$aa_seq, ]
non_nuc_df_deletion <- Deletions.df[Deletions.df$aa_seq %in% Non_nucleating.df$aa_seq, ]
non_nuc_df_insertion <- Insertions.df[Insertions.df$aa_seq %in% Non_nucleating.df$aa_seq, ]
non_nuc_df <- rbind(non_nuc_df_singles, non_nuc_df_deletion, non_nuc_df_insertion)
non_nuc_df$aa_seq <- as.character(non_nuc_df$aa_seq)

#Set the stop mode as the nucleation score for the non-nucleating variants

stops <- all_variants_df[all_variants_df$STOP == TRUE, ]
dens <- density(stops$fitness)
estimated_mode <- dens$x[which.max(dens$y)]
non_nuc_df$fitness <- estimated_mode
non_nuc_df <- non_nuc_df[!(non_nuc_df$aa_seq %in% all_variants_df$aa_seq), ]

all_var_df <- full_join(all_variants_df, non_nuc_df)
all_var_df$nscore_c <- all_var_df$fitness+(-mean_syn_1codon)

write.table(all_var_df, file="MBG_indels_processed_data.tsv", sep="\t", quote = F, row.names = F)


