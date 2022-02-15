####################
# HMST-Seq-Analyzer data process part III: format merged table
####################
# including format DMRs table, filter, count and plot 
  # extract methylation_levle (me = original median-1) 
  # mean methylation (mean_me), mean_rratio, mean_pvalue(pvalue)
# filter by abs(rratio) >= 0.15 for CG and >= 0.06 for CHX

# original path: SAGA /cluster/projects/nn9383k/mingyiy/hmst_katja2/hmst_Neil1_vs_WT_6M
# original file name: script_DMR_table_format.R 
# run Rscript ../script_DMR_table_format.R


# module load R/4.0.0-foss-2020a
#!/bin/sh
library(tidyverse)

# example here: Neil1_vs_WT
# setwd("~/project_katja2/manuscript_Neil/supplementary/dataIn_fig1/data_DMRs_3_table_format/")
# work path: dataIn_fig1/data_DMRs_3_table_format/
out = "out_DMRs_table"
dir.create(out)

# input data
#context = "CHG" #test 
# for_loop for context CG, CHG and CHH:
for (context in c("CG", "CHG", "CHH")) {print(context)

File = paste("out_table/allChr_table_group_KO_vs_WT_", context,".txt", sep = "")
total_df <- read.delim(file= File)
#total_df <- read.delim("out_table/allChr_table_group_KO_vs_WT_CG.txt")
tibble(total_df)
colnames(total_df)

##############
# Step 1. table format
##############
total_df2 <- total_df %>% 
      mutate(mean_KO = (median_KO_1.x + median_KO_1.y + median_KO_2.x + median_KO_2.y)/4 -1,
         mean_WT = (median_WT_1.x + median_WT_1.y + median_WT_2.x + median_WT_2.y)/4 -1,
         mean_dif = (mean_KO - mean_WT),
         mean_pval = (pvals_1 + pvals_2 + pvals_3 + pvals_4)/4, 
         mean_rratio = (rratio_1 + rratio_2 +rratio_3 +rratio_4)/4) %>% 
      mutate(KO_1a = median_KO_1.x -1,       
          KO_1b = median_KO_1.y -1,
          KO_2a = median_KO_2.x -1,       
          KO_2b = median_KO_2.y -1,
          
          WT_1a = median_WT_1.x -1,       
          WT_1b = median_WT_1.y -1,
          WT_2a = median_WT_2.x -1,       
          WT_2b = median_WT_2.y -1) %>% 
    select(-c(4:5, 8:9, 12:13, 16:17)) %>% 
    select(c(1:3,12:24,4,6,8,10))
summary(total_df2)         
head(total_df2, n=3)
colnames(total_df2)
# save table_format: 
outFile = paste(out, "/table_format_DMRs_group_KO_vs_WT_", context, ".txt", sep = "")
write.table(total_df2, file = outFile, row.names = F, quote = F, sep = "\t")
print("Done step 1. table_format")

##############
# Step 2. table_format_filter
##############
# filter by mean diff
dif1 = 0.15 # for CG
dif2 = 0.06 # for CHX
if (context == "CG") {d = dif1
  print(paste("Filter CG by abs diff >= ", d, sep = ""))
} else {
  d = dif2
  print(paste("Filter CHX by abs diff >= ", d, sep = ""))
}
df3 <- total_df2 %>% 
  filter(abs(KO_1a - WT_1a) >= d &
           abs(KO_1b - WT_2a) >= d &
           abs(KO_2a - WT_1b) >= d &
           abs(KO_2b - WT_2b) >= d )
#summary(df3$mean_dif[df3$mean_dif > 0])
#summary(df3$mean_dif[df3$mean_dif < 0])
#length(total_df2$mean_dif[total_df2$mean_dif < - 0.06])

# save table_format_filter: 
outFile = paste(out, "/table_format_DMRs_group_KO_vs_WT_", context, "_filter_", d, "_dif.txt", sep = "")
write.table(df3, file = outFile, row.names = F, quote = F, sep = "\t")

# count appearance times (frequency) for each gene and export genes with >= 3 times:
df3b <- table(df3$gene) %>% as.data.frame() %>% arrange(desc(Freq)) %>% filter(Freq >=3) %>% arrange(Var1) %>% rename(gene_name="Var1")
outFile = paste(out, "/table_format_DMRs_group_KO_vs_WT_", context, "_filter_", d, "_dif_gene_count_over_three_times_geneList.txt", sep = "")
write.table(df3b, file = outFile, row.names = F, quote = F, sep = "\t")

# extract gene table with count over 3 times 
df_total_count <- data.frame()
for (i in 1:nrow(df3b)) { print(df3b$gene_name[i])
  df3c <- df3 %>% filter(gene== as.character(df3b$gene_name[i]))
  df_total_count <- rbind(df_total_count, df3c)
}

df3d <- df_total_count %>% select(1:7)
outFile = paste(out, "/table_format_DMRs_group_KO_vs_WT_", context, "_filter_", d, "_dif_gene_count_over_three_times_matrix.txt", sep = "")
write.table(df3d, file = outFile, row.names = F, quote = F, sep = "\t")

# filter unique genes and export geneList for enrichment analysis:
df4 <- df3 %>%   
      distinct(gene, .keep_all = T) %>% 
      select(3:6,7,2,1) %>% 
      arrange(gene)
head(df4)
outFile = paste(out, "/geneList_DMRs_", context, "_filter_", d, "_dif.txt", sep = "")
write.table(df4, file = outFile, row.names = F, quote = F, sep = "\t")
print("Done step 2. table_format_filter")


##############
# Step 3. count filtered DMRs genes by feature
##############
# create a empty dataframe for gene count
df_total = data.frame()
for (Feature in c("gene", "TSS", "TES")) { print(Feature)
  df_up <- df4 %>% 
    filter(feature == Feature) %>% 
    filter(mean_dif > 0) %>% 
    arrange(gene) %>% 
    distinct(gene, .keep_all=T)
  df_down <- df4 %>% 
    filter(feature == Feature) %>% 
    filter(mean_dif < 0) %>% 
    arrange(gene) %>% 
    distinct(gene, .keep_all=T)
  head(df_down,3)
  
  genelist_up <- df_up[1]
  genelist_down <- df_down[1]
  
  outF1= paste(out, "/DMR_table_", context, "_", Feature, "_hyper.txt", sep="")
  outF2= paste(out, "/DMR_table_", context, "_", Feature, "_hypo.txt", sep="")
  outF3= paste(out, "/geneList_", context, "_", Feature, "_hyper.txt", sep="")
  outF4= paste(out, "/geneList_", context, "_", Feature, "_hypo.txt", sep="")
  
  write.table(df_up, file=outF1, col.names = T, row.names = F, quote = F)
  write.table(df_down, file=outF2, col.names = T, row.names = F, quote = F)
  write.table(genelist_up, file=outF3, col.names = F, row.names = F, quote = F)
  write.table(genelist_down, file=outF4, col.names = F, row.names = F, quote = F)
  
  # count genes and append data frame
  df <- data.frame(feature = Feature, 
                   type = c("hyper", "hypo"), 
                   count = c(nrow(df_up),nrow(df_down)))
  df_total <- rbind(df_total,df)
  print(df_total)
  outF5= paste(out, "/sum_gene_count_", context, ".txt", sep="")
  write.table(df_total, file=outF5, col.names = T, row.names = F, quote = F)
}  # end of for_loop in Feature
print("Done step 3. count DMRs genes")

############################
# Step 4. plot for each context (CG, CHG, CHH)
############################
# plot for all
library("RColorBrewer")
df1 = df_total 
df1[df1$feature == "gene", "feature"] <- "gene_Body"
print(df1)


# set new ylim, 2021_11_05 for manuscript
if (context == "CG") {y = 300
  print(paste("set ylim for CG at ", y, sep = ""))
} else if (context == "CHG") { y = 2000
  print(paste("set ylim for CHG at ", y, sep = ""))
} else {
  y = 4000
  print(paste("set ylim for CHH at ", y, sep = ""))
}

Ylim = y
print(Ylim)

p <- ggplot(data=df1, aes(x=feature, y=count, fill=type)) +
  geom_bar(stat="identity", color="black",position=position_dodge())+
  geom_text(aes(label=count), vjust=1.6, color="white",
            position = position_dodge(0.9), size=6)+
  scale_fill_manual(values=c('#E69F00', '#999999'))+
  theme_minimal(base_size = 18) +
  labs(y = "DMRs gene count") +
  ggtitle(context) +
  ylim(0, Ylim) 

outF6= paste(out, "/plot_sum_gene_count_", context, ".jpg", sep="")
print(outF6)
jpeg(file = outF6, 
     units="in",
     width=6,
     height=6,
     res=300, type="cairo") # in saga with: type="cairo"
print(p)
dev.off()
print("Done step 4. plot count DMRs genes")

}  # end of for_loop in context CG, CHG and CHH.


# sum total count for CG, CHG and CHH
df_sum <- data.frame()
for ( context in c("CG", "CHG", "CHH") ) { print(context)  
  File = paste(out, "/sum_gene_count_", context, ".txt", sep = "")
  df_in <- read.table(File, header = T)
  df_temp <- data.frame(context = context, DMRs_genes = sum(df_in$count) )
  df_sum <- rbind(df_sum, df_temp)
  }
df_sum_total <- data.frame(context = "Total", DMRs_genes = sum(df_sum$DMRs_genes))
df_sum <- rbind(df_sum, df_sum_total)

outFile = paste(out, "/sum_gene_count_total.txt", sep = "")
write.table(df_sum, file = outFile, row.names = F, quote = F, sep = "\t")

print("Done DMR table_format, filter, gene count and plot !")
# -------