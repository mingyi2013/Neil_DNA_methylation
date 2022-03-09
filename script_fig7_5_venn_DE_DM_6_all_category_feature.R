###
# overlap DE genes and DM, total
# ggplot and extract count table by catalog, correlation, feature, mC_type, DE_type.
###
###
# overlap DE genes among comparisions
###
library(tidyverse)

#setwd
setwd("/Users/mingyi/project_katja2_RNAseq")
out = "venn_DE_DM_all_test"
dir.create(out)

for (group in c("Neil1", "Neil2", "DK")) {print(group)
  
# read file
Catalog = c("DMC_CG", "DMP_CG","DMR_CG_tile", 
            "DMR_CCGG_5hmC_gene", "DMR_CCGG_5hmC_TSS","DMR_CCGG_5hmC_TES",
            "DMR_CCGG_gene","DMR_CCGG_TSS","DMR_CCGG_TES",
            "DMRs_CG","DMRs_CHG","DMRs_CHH")
#Catalog= Catalog[1:3]   # test    
df_total= data.frame()
df_count_total = data.frame()
for (catalog in Catalog) { print(catalog)  
inFile = paste("dataIn_fig7/venn_DE_",catalog, "/commonGenes_matrix_", group, "_vs_WT.txt", sep = "")
df <- read.delim(file = inFile, header = T) 
print(nrow(df))
df_total <- rbind(df_total, df)

df_count <- data.frame(catalog = catalog, 
                       positive = nrow(df[df$correlation == "positive", ]),
                       negative = nrow(df[df$correlation == "negative", ]))
df_count_total <- rbind(df_count_total, df_count)
}

head(df_total)
df_total$catalog <- as.factor(df_total$catalog)
df_total$DE_dif <- round(df_total$DE_dif, 2) # keep 2 decimal places
df_total$mC_dif <- round(df_total$mC_dif, 2)

df2 <- df_total %>% arrange(gene_name, catalog, correlation)
outFile = paste(out, "/DE_DM_overlap_all_table_", group, "_vs_WT.txt", sep = "")
write.table(df2, file=outFile, sep = "\t", quote = F, row.names = F)

gene_count <- df2 %>% group_by(gene_name) %>% summarise(count= n()) %>% arrange(desc(count))
outFile = paste(out, "/gene_count_", group, "_vs_WT.txt", sep = "")
write.table(gene_count, file=outFile, sep = "\t", quote = F, row.names = F)

# gene count large than 2 times: 21 genes
gene_top <- gene_count %>% filter(count > 2) %>% select(gene_name)
df2_top <- df2 %>% filter(gene_name %in% gene_top$gene_name )

outFile = paste(out, "/DE_DM_overlap_all_table_top_count_over_3times_", group, "_vs_WT.txt", sep = "")
write.table(df2_top, file=outFile, sep = "\t", quote = F, row.names = F)

# count by coorelation: 
df2_coorelation <- df2 %>%  
        group_by(correlation) %>% summarise(count= n())
outFile = paste(out, "/DE_DM_overlap_count_correlation_", group, "_vs_WT.txt", sep = "")
write.table(df2_coorelation, file=outFile, sep = "\t", quote = F, row.names = F)

# count by catalog: 
df2_catalog <- df2 %>%  
  group_by(catalog) %>% summarise(count= n())
outFile = paste(out, "/DE_DM_overlap_count_catalog_", group, "_vs_WT.txt", sep = "")
write.table(df2_catalog, file=outFile, sep = "\t", quote = F, row.names = F)

# count by others
#X = c("DE_type", "mC_type", "feature")
# count by DE_type: 
df2_x <- df2 %>%  
  group_by(DE_type) %>% summarise(count= n())
outFile = paste(out, "/DE_DM_overlap_count_DE_type_", group, "_vs_WT.txt", sep = "")
write.table(df2_x, file=outFile, sep = "\t", quote = F, row.names = F)

# count by mC_type: 
df2_x <- df2 %>%  
  group_by(mC_type) %>% summarise(count= n())
outFile = paste(out, "/DE_DM_overlap_count_mC_type_", group, "_vs_WT.txt", sep = "")
write.table(df2_x, file=outFile, sep = "\t", quote = F, row.names = F)

# count by feature: 
df2_x <- df2 %>%  
  group_by(feature) %>% summarise(count= n())
outFile = paste(out, "/DE_DM_overlap_count_feature_", group, "_vs_WT.txt", sep = "")
write.table(df2_x, file=outFile, sep = "\t", quote = F, row.names = F)


### plot by catalog
# reshape data by gather()
#df3 <- df_count_total %>% gather(key = "coorelation", value= "count", - catalog)
g <- ggplot(df2, aes(catalog)) +
    geom_bar(aes(fill = correlation)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual("legend", values = c("negative" = "blue", "positive" = "orange"))

 
jpeg(file=paste(out, "/plot_DE_DM_overlap_count_by_catalog_", group, "_vs_WT.jpg", sep = ""),
     units="in", 
     width=7, 
     height=3, 
     res=300)
print(g)
dev.off()


### plot by feature
g <- ggplot(df2, aes(feature)) +
  geom_bar(aes(fill = correlation)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual("legend", values = c("negative" = "blue", "positive" = "orange"))
  

jpeg(file=paste(out, "/plot_DE_DM_overlap_count_by_feature_", group, "_vs_WT.jpg", sep = ""),
     units="in", 
     width=3, 
     height=3, 
     res=300)
print(g)
dev.off()

}  # end of for_loop in group.


