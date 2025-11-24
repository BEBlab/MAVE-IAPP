
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to generate all possible deletions and plot heatmap
plot_deletion_heatmap <- function(WT, mut_data, fitness_col = "fitness_merged",
                                  fdr_col = "category_fdr", show_legend_inside = TRUE) {
  # Split WT sequence into characters
  WT_chars <- unlist(strsplit(WT, ""))
  WT_length <- length(WT_chars)
  
  # ---------------------------
  # Generate internal deletions
  # ---------------------------
  internal_deletions <- data.frame(
    del_start = integer(),
    del_end = integer(),
    del_length = integer(),
    mutant_seq = character()
  )
  
  for(start_pos in 2:(WT_length-1)) {
    for(end_pos in start_pos:(WT_length-1)) {
      mutant_seq <- paste0(c(WT_chars[1:(start_pos-1)], WT_chars[(end_pos+1):WT_length]), collapse = "")
      internal_deletions <- rbind(internal_deletions, data.frame(
        del_start = start_pos,
        del_end = end_pos,
        del_length = end_pos - start_pos + 1,
        mutant_seq = mutant_seq
      ))
    }
  }
  
  # ---------------------------
  # Generate truncations
  # ---------------------------
  truncations <- data.frame(
    del_start = integer(),
    del_end = integer(),
    del_length = integer(),
    mutant_seq = character()
  )
  
  # N-terminal truncations
  for(start_pos in 2:WT_length) {
    mutant_seq <- paste0(WT_chars[start_pos:WT_length], collapse = "")
    truncations <- rbind(truncations, data.frame(
      del_start = 1,
      del_end = start_pos - 1,
      del_length = start_pos - 1,
      mutant_seq = mutant_seq
    ))
  }
  
  # C-terminal truncations
  for(end_pos in 1:(WT_length-1)) {
    mutant_seq <- paste0(WT_chars[1:end_pos], collapse = "")
    truncations <- rbind(truncations, data.frame(
      del_start = end_pos + 1,
      del_end = WT_length,
      del_length = WT_length - end_pos,
      mutant_seq = mutant_seq
    ))
  }
  
  truncations <- unique(truncations)
  
  # Combine all deletions
  all_deletions <- rbind(internal_deletions, truncations)
  
  # Merge with actual experimental data
  merged_del <- left_join(all_deletions, mut_data, by = c("mutant_seq" = "aa_seq"))
  
  # Create upper-triangle grid for heatmap
  full_grid <- expand.grid(del_start = 1:WT_length, del_end = 1:WT_length) %>%
    filter(del_start <= del_end)
  
  # Ensure correct start and end
  merged_del$del_start <- merged_del$del_start.x %||% merged_del$del_start
  merged_del$del_end <- merged_del$del_end.x %||% merged_del$del_end
  
  heatmap_df <- full_grid %>%
    left_join(merged_del, by = c("del_start", "del_end"))
  
  # ---------------------------
  # Define colors
  # ---------------------------
  min_val <- abs(min(mut_data[[fitness_col]], na.rm = TRUE))
  max_val <- abs(max(mut_data[[fitness_col]], na.rm = TRUE))
  cols <- c(
    colorRampPalette(c("brown3", "grey95"))((min_val/(min_val+max_val)*100)-0.5),
    colorRampPalette("grey95")(1),
    colorRampPalette(c("grey95", "darkblue"), bias=1)((max_val/(min_val+max_val)*100)-0.5)
  )
  
  # X and Y labels
  IAPPseq_pos <- paste0(WT_chars, "\n", 1:WT_length)
  IAPPseq_pos_y <- paste0(WT_chars, 1:WT_length)
  
  # ---------------------------
  # Plot fitness heatmap
  # ---------------------------
  p_fit <- ggplot(heatmap_df, aes(x = factor(del_start, levels = 1:WT_length),
                                  y = factor(del_end, levels = 1:WT_length),
                                  fill = .data[[fitness_col]])) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = cols, limits = c(-min_val, max_val), na.value = "grey60") +
    labs(x = "First deleted position", y = "Last deleted position", fill = "Nucleation\nscore") +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_discrete(labels = IAPPseq_pos) +
    scale_y_discrete(labels = IAPPseq_pos_y, position = "left")
  
  if(show_legend_inside) {
    p_fit <- p_fit + theme(legend.position = c(0.8, 0.2),
                           legend.background = element_rect(fill = "white", color = "black"))
  }
  
  # ---------------------------
  # Plot FDR category heatmap
  # ---------------------------
  colors_fdr <- c("darkred", "darkblue", "white")
  p_fdr <- ggplot(heatmap_df, aes(x = factor(del_start, levels = 1:WT_length),
                                  y = factor(del_end, levels = 1:WT_length),
                                  fill = .data[[fdr_col]])) +
    geom_tile(color = "white") +
    scale_fill_manual(values = colors_fdr, na.value = "grey60", name = "Nucleation\nscore (FDR=0.1)") +
    labs(x = "First deleted position", y = "Last deleted position") +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_discrete(labels = IAPPseq_pos) +
    scale_y_discrete(labels = IAPPseq_pos_y, position = "left")
  
  if(show_legend_inside) {
    p_fdr <- p_fdr + theme(legend.position = c(0.8, 0.2),
                           legend.background = element_rect(fill = "white", color = "black"))
  }
  
  return(list(fitness_heatmap = p_fit, fdr_heatmap = p_fdr, heatmap_data = heatmap_df))
}

# ---------------------------
# Example usage
# ---------------------------
WT = "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"

load("All_mutants.RData")

multiaa <- all_ds[all_ds$Mutation_type %in% c("C-terminal truncations", "N-terminal truncations", "Single aa deletions", "Multi aa deletions")]

result <- plot_deletion_heatmap(WT, multiaa)
result$fitness_heatmap
result$fdr_heatmap

# Save plots
ggsave(result$fitness_heatmap, filename = "Multi_amino_acid_deletions_heatmap.pdf", width = 9, height = 9)
ggsave(result$fdr_heatmap, filename = "Multi_amino_acid_deletions_heatmap_by_FDR.pdf", width = 9, height = 9)
