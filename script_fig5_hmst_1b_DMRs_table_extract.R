######
# extract_5hmC_table
# 2021_12_03
#####

library(tidyverse)
library(stringr)

#!/bin/sh
# extract HMST DMR table, including info chr_position, transcripts, geneName, median, p and rratio
# filter geneList in CG by difference >= 20 #test_2021_01_16

# module load R/4.0.0-foss-2020a, 
# run for multiple files, 2021_12_03 on Mac: Given input file in out/data/. work directory: inside data fold
# for file in [ND]*_vs_WT_*5*C_[Tg]*_imputedWith_zeros_DMRs_hyp*[ro].csv; do Rscript ../../../script_fig5_hmst_1b_DMRs_table_extract.R $file table_$file filter_$file; done

# set arguments, input_file= args[1], out_file= args[2], out_filter = args[3]
args=commandArgs(trailingOnly = T)

# read file
df <- read.csv(file = args[1])

# test
#df <- read.csv(file = "/Volumes/WD_Y4/GWBS_katja2_backup/hmst_katja2_hmC/out/data/Neil1_vs_WT_5mC_TSS_imputedWith_zeros_DMRs_hypo.csv")
#df <- read.csv(file = "./dataIn_fig5/out/data/Neil1_vs_WT_5mC_TSS_imputedWith_zeros_DMRs_hypo.csv")

# remove string
df2 <- df %>% 
  separate(name_comb, c("location", "gene"), sep = "(:[-+]:)", extra = "merge") %>% 
  separate(gene, c("transcript", "gene"), sep = "(:0)", extra = "merge") %>%    
  separate(gene, c("gene1", NA), sep = "(:[-+]:)", extra = "merge") %>%                                   
  separate(gene1, c("feature", "gene"), sep = "(:1000:1000)", extra = "merge") 
  
# remove all characters such as: || (only keep letter and numbers)
df2$gene <- str_replace_all(df2$gene, "[^[:alnum:]]", "")
df2$feature <- str_replace_all(df2$feature, "[^[:alnum:]]", "")
df2 <- df2[, -2]

# gene filtered by mean difference:
diff_cutoff = 0.2
df2_diff <- df2 %>% 
  filter(abs(median_KO - median_WT) >= diff_cutoff) %>% 
  dplyr::select(gene) %>%
  distinct(gene) %>% 
  arrange(gene)

write.table(df2, file= args[2], sep = "\t", quote = F, row.names = F)
write.table(df2_diff, file= args[3], sep = "\t", quote = F, row.names = F, col.names = F)

print("Extraction table and filter done!")
# ----------------------------------------


