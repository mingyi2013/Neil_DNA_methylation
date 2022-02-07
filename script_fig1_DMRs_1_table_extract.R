###############
# HMST-Seq-Analyzer data process part I: DMRs table extraction
###############
# extract HMST DMR table, including info chr_position, transcripts, geneName, median, p and rratio
# filter geneList in CG by difference >= 20 
# original script file: extract_hmst_DMRs_table_CG.R
# original run location: SAGA /cluster/projects/nn9383k/mingyiy/hmst_katja2/*/*/data/

# Example here: 
# setwd("~/project_katja2/manuscript_Neil/supplementary/dataIn_fig1/data_DMRs_1_table_extract")
  # df <- read.csv(file = "./KO1.me_vs_WT1.me_5mC_geneBody_imputedWith_zeros_DMRs_hyper.csv") # test inside Rstudio
# run in terminal: for file in KO*5mC_[Tg]*_imputedWith_zeros_DMRs_hyp*[ro].csv; do Rscript ../../script_fig1_DMRs_1_table_extract.R  $file table_$file filter_$file; done

# module load R/4.0.0-foss-2020a
#!/bin/sh

# set arguments, input_file= args[1], out_file= args[2], out_filter = args[3]
args=commandArgs(trailingOnly = T)

df <- read.csv(file = args[1])

# remove string
library(tidyverse)
library(stringr)

df2 <- df %>% 
  separate(name_comb, c("location", "gene"), sep = "(:[-+]:)", extra = "merge") %>% 
  separate(gene, c("transcript", "gene"), sep = "(:0)", extra = "merge") %>%    
  separate(gene, c("gene1", NA), sep = "(:[-+]:)", extra = "merge") %>%                                    
  separate(gene1, c("feature", "gene"), sep = "(:1000:1000)", extra = "merge")

# remove all characters such as: || (only keep letter and numbers)
df2$gene <- str_replace_all(df2$gene, "[^[:alnum:]]", " ")
df2$feature <- str_replace_all(df2$feature, "[^[:alnum:]]", " ")
df2 <- df2[, -2]

# gene filtered by mean difference:
df2_diff10 <- df2 %>% 
  filter(abs(median_KO - median_WT) >= 0.20) %>% 
  select(gene) %>%
  distinct(gene) %>% 
  arrange(gene)

write.table(df2, file= args[2], sep = "\t", quote = F, row.names = F)
write.table(df2_diff10, file= args[3], sep = "\t", quote = F, row.names = F, col.names = F)

print("End!")
# -------------
