
#Script to generate the text file required for Chimera to paint the structures with the mean nucleation score for each position. 
#This script outputs the file that will be used as an input in the Chimera software. 

require(dplyr)
load("INDEL_datasets.RData")
df <- as.data.frame(Singles.df %>% group_by(Pos) %>% dplyr::summarize(Mean = mean(nscore_c, na.rm = TRUE), na.rm=TRUE))
res_number <- paste("\t", ":", seq(1, 37, 1), "\t", df$Mean)
write_text_mean <-c("attribute: median\t", "recipient: residues\t", res_number)
write.table(write_text_mean, file="attr_file_IAPP_singles_mean", quote = F, row.names = F, col.names = F)

