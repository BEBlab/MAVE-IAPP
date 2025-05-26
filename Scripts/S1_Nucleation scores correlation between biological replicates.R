
#Script required to generate Supplementary Figure 1. 

#Load required packages
require(ggplot2)
require(ggpubr)
require(ggpointdensity)

#Load required data 
load("INDEL.df.RData")

#Tech replicate correlations 
p_corr_replicate_12 <- ggplot(INDEL.df, aes(x = fitness1_uncorr, y = fitness2_uncorr)) + stat_binhex(geom = "hex", bins = 50) + xlab("NS of replicate 1") + ylab("NS of replicate 2")
p_corr_replicate_12 <- p_corr_replicate_12 + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                                                     axis.ticks.x.top = element_line(),       
                                                                     plot.title = element_blank(),
                                                                     panel.grid.major = element_blank(), 
                                                                     panel.grid.minor = element_blank(),
                                                                     axis.text = element_text(size=8), 
                                                                     axis.line = element_line(size = 0.2),
                                                                     axis.title = element_text(size = 8), 
                                                                     legend.position = "none")
p_corr_replicate_12 <- p_corr_replicate_12 + scale_x_continuous(limits = c(-8, 5)) + scale_y_continuous(limits = c(-8, 5)) + stat_cor(size = 2) + scale_fill_gradient(high="black", low="grey90")
p_corr_replicate_12 <- p_corr_replicate_12 + annotate("text", x = -6.3, y = 3.2, label = paste("n = ", nrow(INDEL.df[!(is.na(INDEL.df$fitness2_uncorr)) & !(is.na(INDEL.df$fitness1_uncorr)), ]), sep = ""), size = 2)
p_corr_replicate_12
ggsave(p_corr_replicate_12, filename = "Rep12.pdf", width = 2, height = 2)


p_corr_replicate_13 <- ggplot(INDEL.df, aes(x = fitness1_uncorr, y = fitness3_uncorr)) + stat_binhex(geom = "hex", bins = 50) + xlab("NS of replicate 1") + ylab("NS of replicate 3")
p_corr_replicate_13 <- p_corr_replicate_13 + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                                                     axis.ticks.x.top = element_line(),       
                                                                     plot.title = element_blank(),
                                                                     panel.grid.major = element_blank(), 
                                                                     panel.grid.minor = element_blank(),
                                                                     axis.text = element_text(size=8), 
                                                                     axis.line = element_line(size = 0.2),
                                                                     axis.title = element_text(size = 8), 
                                                                     legend.position = "none")
p_corr_replicate_13 <- p_corr_replicate_13 + scale_x_continuous(limits = c(-8, 5)) + scale_y_continuous(limits = c(-8, 5)) + stat_cor(size = 2) + scale_fill_gradient(high="black", low="grey90")
p_corr_replicate_13 <- p_corr_replicate_13 + annotate("text", x = -6.3, y = 3.2, label = paste("n = ", nrow(INDEL.df[!(is.na(INDEL.df$fitness3_uncorr)) & !(is.na(INDEL.df$fitness1_uncorr)), ]), sep = ""), size = 2)
ggsave(p_corr_replicate_13, filename = "Rep13.pdf", width = 2, height = 2)

p_corr_replicate_23 <- ggplot(INDEL.df, aes(x = fitness2_uncorr, y = fitness3_uncorr)) + stat_binhex(geom = "hex", bins = 50) + xlab("NS of replicate 2") + ylab("NS of replicate 3")
p_corr_replicate_23 <- p_corr_replicate_23 + theme_minimal() + theme(axis.ticks.y=element_blank(),
                                                                     axis.ticks.x.top = element_line(),       
                                                                     plot.title = element_blank(),
                                                                     panel.grid.major = element_blank(), 
                                                                     panel.grid.minor = element_blank(),
                                                                     axis.text = element_text(size=8), 
                                                                     axis.line = element_line(size = 0.2),
                                                                     axis.title = element_text(size = 8))
p_corr_replicate_23 <- p_corr_replicate_23 + scale_x_continuous(limits = c(-8, 5)) + scale_y_continuous(limits = c(-8, 5)) + stat_cor(size = 2) + scale_fill_gradient(high="black", low="grey90", name = "Count")
p_corr_replicate_23 <- p_corr_replicate_23 + annotate("text", x = -6.3, y = 3.2, label = paste("n = ", nrow(INDEL.df[!(is.na(INDEL.df$fitness3_uncorr)) & !(is.na(INDEL.df$fitness3_uncorr)), ]), sep = ""), size = 2)
p_corr_replicate_23
ggsave(p_corr_replicate_23, filename = "Rep23.pdf", width = 2.7, height = 2)


