##############################################
## Load required packages and define colours for plotting
##############################################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(stringr)
library(ggrepel)
colors_aa <- c("darkgrey", "#9A703EFF", "#EE0011FF", "#FEC10BFF", "#15983DFF", "#0C5BB0FF")

##############################################
## Define alignment dictionaries (IAPP ↔ Aβ42)
## Alignment performed with T-COFFEE
##############################################

# Substitutions
equiv_dict <- data.frame(
  Pos_IAPP = c(1, 2, rep("-", 5), 3:36, "-", 37),
  Pos_AB   = c(1:34, "-", 35:42),
  plotting_number = 1:43
)

# Full positional vectors
complete_pos             <- as.character(c(1:34, "-", 35:42))
complete_pos_for_plotting <- 1:43

# Insertions
equiv_dict_ins <- data.frame(
  Pos_IAPP = c(1, 2, rep("-", 5), 3:35, "-", 36),
  Pos_AB   = c(1:34, "-", 35:41),
  plotting_number = 1:42
)

##############################################
## Load + parse IAPP dataset
##############################################

load("All_mutants.RData")

# Single substitutions
singles_IAPP <- as.data.frame(all_ds[all_ds$Mutation_type == "Singles", ])
singles_IAPP$Pos <- as.numeric(str_extract(singles_IAPP$HGVS_nomenclature, "\\d+"))
singles_IAPP$WT_AA <- str_extract(singles_IAPP$HGVS_nomenclature, "(?<=p\\.)[A-Z][a-z]{2}")
singles_IAPP$Mut   <- str_extract(singles_IAPP$HGVS_nomenclature, "[A-Z][a-z]{2}$")

# Insertions
insertions_IAPP <- as.data.frame(all_ds[all_ds$Mutation_type == "Insertions", ])
insertions_IAPP$ins_pos <- as.numeric(str_extract(insertions_IAPP$HGVS_nomenclature, "\\d+"))
insertions_IAPP$ins_aa  <- sapply(insertions_IAPP$HGVS_nomenclature, function(x) strsplit(x, "ins")[[1]][2])

##############################################
## Map 3-letter → 1-letter amino-acid codes
##############################################

aa3_to_1 <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", 
  Gln="Q", Glu="E", Gly="G", His="H", Ile="I",
  Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
  Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V"
)

singles_IAPP$Mut     <- aa3_to_1[singles_IAPP$Mut]
singles_IAPP$WT_AA   <- aa3_to_1[singles_IAPP$WT_AA]
insertions_IAPP$ins_aa <- aa3_to_1[insertions_IAPP$ins_aa]

##############################################
## Format parsed IAPP data
##############################################

singles_IAPP <- singles_IAPP[, c("fitness_merged", "sigma_scaled", "WT_AA", "Mut", "Pos", "category_fdr")]
colnames(singles_IAPP) <- c("nscore_c_IAPP", "sigma_IAPP", "WT_AA_IAPP", "Mut", "Pos_IAPP", "category_fdr_IAPP")
singles_IAPP$SEQ <- "IAPP"
singles_IAPP$pos_plotting <- equiv_dict$plotting_number[
  match(singles_IAPP$Pos_IAPP, equiv_dict$Pos_IAPP)
]

insertions_IAPP <- insertions_IAPP[, c("fitness_merged", "sigma_scaled", "ins_aa", "ins_pos", "category_fdr")]
colnames(insertions_IAPP) <- c("nscore_c_IAPP", "sigma_IAPP", "ins_aa", "ins_pos_IAPP", "category_fdr_IAPP")
insertions_IAPP$SEQ <- "IAPP"
insertions_IAPP$pos_plotting <- equiv_dict_ins$plotting_number[
  match(insertions_IAPP$ins_pos_IAPP, equiv_dict_ins$Pos_IAPP)
]

##############################################
## Load + parse Aβ42 dataset (Seuma et al.)
##############################################

load("Required files/INDEL_datasets_AB.Rdata")

# Substitutions
singles_AB <- singles.df[, c("nscore_c", "sigma", "WT_AA", "Mut", "Pos", "category_fdr")]
colnames(singles_AB) <- c("nscore_c_AB", "sigma_AB", "WT_AA_AB", "Mut", "Pos_AB", "category_fdr_AB")
singles_AB$SEQ <- "AB"
singles_AB$Pos_WT_AB <- singles_AB$Pos_AB
singles_AB$pos_plotting <- equiv_dict$plotting_number[
  match(singles_AB$Pos_AB, equiv_dict$Pos_AB)
]
singles_AB$ID_Equiv <- as.character(singles_AB$Pos_AB)

# Insertions
insertions_reps_AB <- insertions_reps[, c("nscore_c", "sigma", "ins_aa", "ins_pos", "category_fdr")]
colnames(insertions_reps_AB) <- c("nscore_c_AB", "sigma_AB", "ins_aa", "ins_pos_AB", "category_fdr_AB")
insertions_reps_AB$SEQ <- "AB"
insertions_reps_AB$pos_plotting <- equiv_dict_ins$plotting_number[
  match(insertions_reps_AB$ins_pos_AB, equiv_dict_ins$Pos_AB)
]

##############################################
## Merge IAPP ↔ Aβ42 substitutions
##############################################

singles_IAPP$ID_Equiv <- equiv_dict$Pos_AB[
  match(singles_IAPP$Pos_IAPP, equiv_dict$Pos_IAPP)
]

singles_AB_IAPP <- merge(
  singles_AB, singles_IAPP,
  by = c("pos_plotting", "Mut", "ID_Equiv")
)

singles_AB_IAPP$labelAB   <- paste("AB42 =", singles_AB_IAPP$WT_AA_AB, singles_AB_IAPP$Pos_WT_AB)
singles_AB_IAPP$labelIAPP <- paste("IAPP =", singles_AB_IAPP$WT_AA_IAPP, singles_AB_IAPP$Pos_IAPP)

##############################################
## Merge IAPP ↔ Aβ42 insertions
##############################################

insertions_AB_IAPP <- merge(
  insertions_reps_AB, insertions_IAPP,
  by = c("pos_plotting", "ins_aa")
)

insertions_AB_IAPP$labelAB   <- paste("AB42 =", insertions_AB_IAPP$ins_pos_AB)
insertions_AB_IAPP$labelIAPP <- paste("IAPP =", insertions_AB_IAPP$ins_pos_IAPP)

##############################################
## Summary of the effect of the insertion of the same for IAPP and AB
##############################################

sum_ins <- insertions_AB_IAPP %>%
  group_by(ins_aa, category_fdr_IAPP) %>%
  summarise(
    mean_value_IAPP = mean(nscore_c_IAPP),
    mean_value_AB   = mean(nscore_c_AB),
    .groups = "drop"
  )

# AA classification
sum_ins$type <- dplyr::case_when(
  sum_ins$ins_aa %in% c("A","L","M","V","I") ~ "Aliphatic",
  sum_ins$ins_aa %in% c("K","R")            ~ "Positive",
  sum_ins$ins_aa %in% c("D","E")            ~ "Negative",
  sum_ins$ins_aa %in% c("Q","N","S","T","H","C") ~ "Polar",
  sum_ins$ins_aa %in% c("G","P")            ~ "P/G",
  sum_ins$ins_aa %in% c("F","Y","W")        ~ "Aromatic",
  TRUE ~ "Other"
)


##############################################
## Plot: mean nucleation scores (insertions)
##############################################

p_ins <- ggplot(sum_ins, aes(mean_value_IAPP, mean_value_AB, colour = type)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_text_repel(aes(label = ins_aa), max.overlaps = 30,
                  arrow = arrow(length = unit(0.015, "npc"))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12)) +
  labs(x = "IAPP mean nucleation score",
       y = expression(A*beta*"42 mean nucleation score"),
       colour = "AA insertion type") + scale_color_manual(values = colors_aa) +
  xlim(-4, 2) + ylim(-4, 2) +
  facet_wrap(~category_fdr_IAPP) + stat_cor(aes(x = mean_value_IAPP, y = mean_value_AB), inherit.aes = FALSE)

ggsave("p_corr_facetted_ns_insertions.pdf", p_ins, width = 10, height = 4.5)

##############################################
## Summary of the effect of the same substitutions for IAPP and AB
##############################################

sum_sub <- singles_AB_IAPP %>%
  group_by(Mut, category_fdr_IAPP) %>%
  summarise(
    mean_value_IAPP = mean(nscore_c_IAPP),
    mean_value_AB   = mean(nscore_c_AB),
    .groups = "drop"
  )

sum_sub$type <- dplyr::case_when(
  sum_sub$Mut %in% c("A","L","M","V","I") ~ "Aliphatic",
  sum_sub$Mut %in% c("K","R")            ~ "Positive",
  sum_sub$Mut %in% c("D","E")            ~ "Negative",
  sum_sub$Mut %in% c("Q","N","S","T","H","C") ~ "Polar",
  sum_sub$Mut %in% c("G","P")            ~ "P/G",
  sum_sub$Mut %in% c("F","Y","W")        ~ "Aromatic",
  TRUE ~ "Other"
)

##############################################
## Plot: mean nucleation scores (substitutions)
##############################################

p_sub <- ggplot(sum_sub, aes(mean_value_IAPP, mean_value_AB, colour = type)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_text_repel(aes(label = Mut), max.overlaps = 30,
                  arrow = arrow(length = unit(0.015, "npc"))) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12)) +
  labs(x = "IAPP mean nucleation score",
       y = expression(A*beta*"42 mean nucleation score"),
       colour = "AA substitution type") +
  scale_color_manual(values = colors_aa) + xlim(-4.2, 3) + ylim(-4.2, 3) + stat_cor(aes(x = mean_value_IAPP, y = mean_value_AB), inherit.aes = FALSE) + 
  facet_wrap(~category_fdr_IAPP)

ggsave("p_corr_facetted_ns_substitutions.pdf", p_sub, width = 10, height = 4.5)

##############################################
## Summary of the effect of the mutations in the same equivalent position for IAPP and AB
##############################################

sum_pos <- singles_AB_IAPP %>%
  group_by(ID_Equiv, category_fdr_IAPP) %>%
  summarise(
    mean_value_IAPP = mean(nscore_c_IAPP),
    mean_value_AB   = mean(nscore_c_AB),
    .groups = "drop"
  )

##############################################
## Plot
##############################################

p_pos <- ggplot(sum_pos, aes(mean_value_IAPP, mean_value_AB)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12)) +
  labs(x = "IAPP mean nucleation score",
       y = expression(A*beta*"42 mean nucleation score")) +
  xlim(-4.2, 3) + ylim(-4.2, 3) + stat_cor(aes(x = mean_value_IAPP, y = mean_value_AB), inherit.aes = FALSE) + 
  facet_wrap(~category_fdr_IAPP)

ggsave("p_facetted_pos_mutfrom.pdf", p_pos, width = 10, height = 4)

