# ---- Load libraries ----
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(purrr)

# ---- Load datasets ----
load("All_mutants.RData")
load("aggregation_tools_results_Camsol_IAPP.RData")
load("TANGO_coef.RData")
load("zgg_scores_IAPP.RData")
load("s4pred_IAPP.RData")
load("popEVE_seq.RData")
load("PCs_IAPP_20.RData")
load("s4pred_helical_prop.RData")
load("Alphamissense.RData")

# ---- Normalize predictors ----
result_list_final_zyggregator$Zyggregator_norm <- result_list_final_zyggregator$Zyggregator / nchar(result_list_final_zyggregator$aa_seq)
helix_prop$Numb_helix <- helix_prop$Numb_helix / nchar(helix_prop$aa_seq)
helix_prop$Sum_helix <- helix_prop$Sum_helix
df_camsol$Camsol_score7 <- df_camsol$Camsol_score7 / nchar(df_camsol$aa_seq)
tango_pred$Tango_helix <- tango_pred$Tango_helix / nchar(tango_pred$aa_seq)
tango_pred$Tango_beta <- tango_pred$Tango_beta / nchar(tango_pred$aa_seq)
s4pred_df$Mean_Coil_Prob <- s4pred_df$Mean_Coil_Prob / nchar(s4pred_df$aa_seq)
s4pred_df$Mean_Strand_Prob <- s4pred_df$Mean_Strand_Prob / nchar(s4pred_df$aa_seq)
s4pred_df$Mean_Helix_Prob <- s4pred_df$Mean_Helix_Prob / nchar(s4pred_df$aa_seq)

# ---- Process PCsIAPP ----
PCsIAPP <- PCsIAPP %>%
  mutate(across(where(is.list), ~map_dbl(.x, ~ if(length(.x) == 0) NA_real_ else mean(as.numeric(.x), na.rm=TRUE)))) %>%
  mutate(aa_seq = as.character(aa_seq))

# ---- Merge all datasets ----
merged_df <- all_ds %>%
  mutate(aa_seq = as.character(aa_seq)) %>%
  left_join(PCsIAPP, by = "aa_seq") %>%
  left_join(df_camsol, by = "aa_seq") %>%
  left_join(tango_pred, by = "aa_seq") %>%
  left_join(result_list_final_zyggregator, by = "aa_seq") %>%
  left_join(s4pred_df %>% distinct(aa_seq, .keep_all = TRUE), by = "aa_seq") %>%
  left_join(popeve_df, by = "aa_seq") %>%
  left_join(helix_prop %>% distinct(aa_seq, .keep_all = TRUE), by = "aa_seq") %>%
  left_join(alphamissense, by = "aa_seq")

# ---- Parse mutations for category and core ----
core_positions <- 15:32
parse_mutation <- function(mut_str, core_positions) {
  mut_str <- str_trim(mut_str)
  if(mut_str=="WT") return(data.frame(Mutation_type_full=mut_str, Mut_category="WT", Mut_positions=NA, Core=NA, stringsAsFactors=FALSE))
  if(!str_detect(mut_str,"ins|del")) {
    pos <- as.numeric(str_extract(mut_str,"\\d+"))
    return(data.frame(Mutation_type_full=mut_str, Mut_category="substitutions", Mut_positions=paste(pos,collapse=","), Core=ifelse(any(pos %in% core_positions),"Core","Non-core"), stringsAsFactors=FALSE))
  }
  if(str_detect(mut_str,"del")) {
    pos <- as.numeric(str_extract_all(mut_str,"\\d+")[[1]])
    return(data.frame(Mutation_type_full=mut_str, Mut_category="deletion", Mut_positions=paste(pos,collapse=","), Core=ifelse(any(pos %in% core_positions),"Core","Non-core"), stringsAsFactors=FALSE))
  }
  if(str_detect(mut_str,"ins")) {
    pos <- as.numeric(str_extract_all(mut_str,"\\d+")[[1]])
    inserted_res <- str_extract(mut_str,"ins([A-Za-z]+)$") %>% str_remove("ins")
    type <- ifelse(nchar(inserted_res)>3,"polymerase_slip","insertion")
    return(data.frame(Mutation_type_full=mut_str, Mut_category=type, Mut_positions=paste(pos,collapse=","), Core=ifelse(any(pos %in% core_positions),"Core","Non-core"), stringsAsFactors=FALSE))
  }
}

mut_info <- merged_df %>%
  rowwise() %>%
  mutate(parsed=list(parse_mutation(HGVS_nomenclature, core_positions))) %>%
  unnest(cols=c(parsed))

# ---- Select numeric predictors ----
predictors_order <- c("Camsol_score7","Tango_helix","Tango_beta","Zyggregator_norm",
                      "Mean_Coil_Prob","Mean_Helix_Prob","Mean_Strand_Prob",
                      "EVE","pop.adjusted_EVE","ESM1v","Numb_helix",
                      paste0("PC",1:20),"pathogenicity.score")
numeric_cols <- mut_info %>% select(all_of(predictors_order)) %>% select(where(is.numeric)) %>% colnames()

# ---- Compute correlations ----
safe_cor <- function(x,y){if(all(is.na(x))|all(is.na(y))) return(NA_real_); cor(x,y,use="complete.obs")}

cor_results <- mut_info %>%
  filter(Mut_category!="WT") %>%
  group_by(Mut_category, Core) %>%
  summarise(across(all_of(numeric_cols), ~safe_cor(.x, fitness_merged)), .groups="drop") %>%
  pivot_longer(cols=all_of(numeric_cols), names_to="Predictor", values_to="Correlation") %>%
  mutate(Mut_category_Core=paste(Mut_category, Core, sep="_"))

# ---- Add "All mutants" column ----
all_mut_col <- mut_info %>%
  summarise(across(all_of(numeric_cols), ~safe_cor(.x, fitness_merged))) %>%
  pivot_longer(cols=everything(), names_to="Predictor", values_to="Correlation") %>%
  mutate(Mut_category_Core="All")
cor_results_all <- bind_rows(cor_results, all_mut_col)

# ---- Create figure-friendly names ----
predictor_labels <- c(
  "Camsol_score7"="Cansol",
  "Tango_helix"="Helix propensity (AGADIR)",
  "Tango_beta"="Beta strand propensity (TANGO)",
  "Zyggregator_norm"="Zyggregator",
  "Mean_Coil_Prob"="Coil Probability (s4pred)",
  "Mean_Helix_Prob"="Helix Probability (s4pred)",
  "Mean_Strand_Prob"="Strand Probability (s4pred)",
  "EVE"="EVE",
  "pop.adjusted_EVE"="PopEVE",
  "ESM1v"="ESM1v",
  "Numb_helix"="Helix Count (s4pred)",
  paste0("PC",1:20)=paste0("PC",1:20),
  "pathogenicity.score"="AlphaMissense"
)
x_axis_labels <- c(
  "All"="All mutants",
  "deletion_Core"="Single and multiple aa deletions in core",
  "deletion_Non-core"="Single and multiple aa deletions outside core",
  "insertion_Core"="Single aa insertions in core",
  "insertion_Non-core"="Single aa insertions outside core",
  "polymerase_slip_Core"="Polymerase slippage variants in core",
  "polymerase_slip_Non-core"="Polymerase slippage variants outside core",
  "substitutions_Core"="Single aa substitutions in core",
  "substitutions_Non-core"="Single aa substitutions outside core"
)

# ---- Recode and set factor levels ----
cor_results_all <- cor_results_all %>%
  mutate(Predictor=recode(Predictor, !!!predictor_labels),
         Mut_category_Core=recode(Mut_category_Core, !!!x_axis_labels))
cor_results_all$Predictor <- factor(cor_results_all$Predictor, levels=predictor_labels)
cor_results_all$Mut_category_Core <- factor(cor_results_all$Mut_category_Core, levels=x_axis_labels)

# ---- Identify strong correlations for marking ----
more_than_hm <- cor_results_all %>%
  filter(Correlation < -0.5 | Correlation > 0.5) %>%
  # Make sure x and y are factors matching the plot
  mutate(Mut_category_Core = factor(Mut_category_Core, levels = levels(cor_results_all$Mut_category_Core)),
         Predictor = factor(Predictor, levels = levels(cor_results_all$Predictor)))

# ---- Plot heatmap with asterisks ----
p_heatmap <- ggplot(cor_results_all, aes(x=Mut_category_Core, y=Predictor, fill=Correlation)) +
  geom_tile() +
  geom_text(data = more_than_hm, aes(label="*"), color="black", size=5) +
  scale_fill_gradient2(low="maroon", mid="gray80", high="darkgreen", midpoint=0, limits=c(-1,1), name="R", na.value=NA) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=13),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line()) +
  xlab("") +
  ylab("")

p_heatmap
ggsave(p_heatmap, filename = "Correlation_matrix_all.pdf", width = 7, height = 9)


