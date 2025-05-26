## Load required packages
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(tidyr)
library(RColorBrewer)

## Load data
load("INDEL.df.RData")
load("INDEL_datasets.RData")

## Define WT sequence and position labels
IAPP_wt <- "KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY"
IAPPseq <- strsplit(IAPP_wt, "")[[1]]
IAPPseq_pos <- paste0(IAPPseq, "\n", 1:37)

## Add position and insertions to PolSlip.df
PolSlip.df$Pos <- sapply(strsplit(PolSlip.df$name, "_"), function(x) as.numeric(x[3]))
PolSlip.df$ins_2aa <- apply(PolSlip.df, 1, function(i) {
  substring(i["aa_seq"], as.numeric(i["Pos"]) + 1, as.numeric(i["Pos"]) + 2)
})

## Subset and complete 2i and 2c data
polslip2i <- PolSlip.df[grep("PolSlip_2i", PolSlip.df$name), ] %>%
  complete(Pos = 1:37) %>% mutate(y = 1)

polslip2c <- PolSlip.df[grep("PolSlip_2c", PolSlip.df$name), ] %>%
  mutate(y = ifelse(Pos %% 2 == 0, 3, 1)) %>%
  complete(Pos = 1:37)

## Helper to categorize data
assign_categories <- function(df) {
  df <- df %>%
    mutate(
      p.adjust = as.numeric(p.adjust),
      category = "WT-like",
      category = case_when(
        !is.na(nscore_c) & p.adjust < 0.01 & nscore_c < 0 ~ "NS- 1%",
        !is.na(nscore_c) & p.adjust < 0.05 & nscore_c < 0 ~ "NS- 5%",
        !is.na(nscore_c) & p.adjust < 0.1  & nscore_c < 0 ~ "NS- 10%",
        !is.na(nscore_c) & p.adjust < 0.25 & nscore_c < 0 ~ "NS- 25%",
        !is.na(nscore_c) & p.adjust < 0.01 & nscore_c > 0 ~ "NS+ 1%",
        !is.na(nscore_c) & p.adjust < 0.05 & nscore_c > 0 ~ "NS+ 5%",
        !is.na(nscore_c) & p.adjust < 0.1  & nscore_c > 0 ~ "NS+ 10%",
        !is.na(nscore_c) & p.adjust < 0.25 & nscore_c > 0 ~ "NS+ 25%",
        is.na(nscore_c) ~ NA_character_,
        TRUE ~ "WT-like"
      ),
      category_fdr = case_when(
        is.na(nscore_c) ~ NA_character_,
        p.adjust < 0.1 & nscore_c < 0 ~ "NS_dec",
        p.adjust < 0.1 & nscore_c > 0 ~ "NS_inc",
        TRUE ~ "WT-like"
      )
    )
  return(df)
}

## Apply to both datasets
fdr_polslip_i <- assign_categories(polslip2i[, c("aa_seq", "Pos", "ins_2aa", "nscore_c", "sigma", "p.adjust", "sig_fdr", "y")])
fdr_polslip_c <- assign_categories(polslip2c[, c("aa_seq", "Pos", "ins_2aa", "nscore_c", "sigma", "p.adjust", "sig_fdr", "y")])

## Define plotting colors
levels <- c("NS- 1%", "NS- 5%", "NS- 10%", "NS- 25%", "WT-like", "NS+ 25%", "NS+ 10%", "NS+ 5%", "NS+ 1%")
colors <- c("darkred", "#D45B5B", "#DC8484", "#E8C2C2", "grey90", "#C3C3DE", "#A1A1CF", "#7E7EC0", "darkblue")
myColors <- setNames(colors, levels)

gray_colors <- c("white", rep("black", 7), "white")
myColors_gray <- setNames(gray_colors, levels)

## Plot function
plot_polslip <- function(df, filename, rect_plot = FALSE, height = 1) {
  if (rect_plot) {
    p <- ggplot(df, aes(xmin = Pos - 0.45, xmax = Pos + 1.45, ymin = y - 1, ymax = y + 1, fill = category)) +
      geom_rect(size = 1) +
      geom_rect(data = df[df$sig_fdr == TRUE, ], colour = "black") +
      geom_text(aes(x = Pos + 0.5, y = y, label = ins_2aa, colour = category_fdr), size = 5, show.legend = FALSE)
  } else {
    p <- ggplot(df, aes(x = Pos, fill = category, y = y)) +
      geom_tile() +
      geom_tile(data = df[df$sig_fdr == TRUE, ], colour = "black", linejoin = "round", size = 0.4) +
      geom_text(aes(label = ins_2aa, colour = category), size = 5)
  }
  
  p <- p +
    theme_bw() +
    scale_fill_manual(values = myColors, name = "Nucleation\nscore", na.value = "white") +
    scale_colour_manual(values = myColors_gray, na.value = "white") +
    scale_x_continuous(breaks = 1:37, labels = 1:37, expand = c(0, 0)) +
    theme(
      axis.ticks.y = element_blank(),
      plot.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 13),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      panel.border = element_blank(),
      legend.position = if (rect_plot) "bottom" else "none"
    )
  
  ggsave(filename, plot = p, height = height, width = 13)
}

## Generate plots
plot_polslip(fdr_polslip_i, "PolSlip_2i_fdr.pdf", rect_plot = FALSE, height = 1)
plot_polslip(fdr_polslip_c, "PolSlip_double_fdr.pdf", rect_plot = TRUE, height = 2)
