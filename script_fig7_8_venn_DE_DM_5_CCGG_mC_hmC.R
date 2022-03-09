###
# overlap DE and 5mC_hypo/5hmC_hyper pattern:
###
library(tidyverse)
library("ggvenn")

#setwd
#setwd("/Users/mingyi/project_katja2_RNAseq")

Catalog = "DMR_CCGG_TSS_5mC_5hmC"
out = paste("venn_DE_", Catalog, sep = "")
dir.create(out)

# read file
Dir <- "./DE_groups/"

for (group in c("Neil1", "Neil2", "DK")) {print(group)
  #for (group in c("Neil1")) {print(group)} # test
  
  Dir <- "./DE_groups/"
  # DE all data
  inFile = paste(Dir, group, "_vs_WT/DE_data_sigDiff_condition_", group, "_vs_WT.txt", sep = "")
  DE_data <- read.delim(file = inFile, header = T) %>% arrange(gene_name) %>% distinct(gene_name, .keep_all = T)
  
  # DE_up
  inFile = paste(Dir, group, "_vs_WT/DE_geneList_up_condition_", group, "_vs_WT.txt", sep = "")
  DE_up <- read.delim(file = inFile, header = F) %>% 
    rename(gene_name = V1) %>% 
    arrange(gene_name) %>% 
    distinct(gene_name)
  df_1 <-  DE_up %>% unlist()
  
  inFile = paste(Dir, group, "_vs_WT/DE_geneList_down_condition_", group, "_vs_WT.txt", sep = "")
  DE_down <- read.delim(file = inFile, header = F) %>% 
    rename(gene_name = V1) %>% 
    arrange(gene_name) %>% 
    distinct(gene_name)
  df_2 <-  DE_down %>% unlist()
  
  
  
  # pattern in 5mC_hypo TSS
  #driver <- "/Volumes/WD_Y4/GWBS_katja2_backup/hmst_katja2_hmC/"
  driver <- "dataIn_fig7/"
  pattern = "TSS"
  inFile = paste(driver, "out/data/table_", group, "_vs_WT_5mC_", pattern, "_imputedWith_zeros_DMRs_hypo.csv", sep = "")
  df_5mC_hypo <- read.delim(file = inFile, header = T, strip.white = T, sep = "\t")  # strip white space in gene column !  
  df_5mC_hypo$mC_type <- "mC_hypo_hmC_hyper"
  
  df_5mC_hypo$catalog <- Catalog
  df_5mC_hypo$mC_dif <- df_5mC_hypo$median_KO - df_5mC_hypo$median_WT
  Location <- substring(df_5mC_hypo$location, 4, 300)
  df_5mC_hypo$location <- Location
  df_5mC_hypo$distance_to_gene <- 0
  names(df_5mC_hypo)[names(df_5mC_hypo) == "gene"] <- "gene_name"
  
  # pattern in 5hmC_hyper TSS
  inFile = paste(driver, "out/data/table_", group, "_vs_WT_5hmC_", pattern, "_imputedWith_zeros_DMRs_hyper.csv", sep = "")
  df_5hmC_hyper <- read.delim(file = inFile, header = T, strip.white = T, sep = "\t")  # strip white space in gene column ! 
  df_5hmC_hyper$mC_type <- "mC_hypo_hmC_hyper"
  
  
  df_5hmC_hyper$catalog <- Catalog
  df_5hmC_hyper$mC_dif <- df_5hmC_hyper$median_KO - df_5hmC_hyper$median_WT
  Location <- substring(df_5hmC_hyper$location, 4, 300)
  df_5hmC_hyper$location <- Location
  df_5hmC_hyper$distance_to_gene <- 0
  names(df_5hmC_hyper)[names(df_5hmC_hyper) == "mC_dif"] <- "hmC_dif"  #rename
  names(df_5hmC_hyper)[names(df_5hmC_hyper) == "gene"] <- "gene_name"
  
  
  # common genes in 5mC_hypo/5hmC_hyper
  common_mC_hmC <- inner_join(df_5mC_hypo, df_5hmC_hyper, 
                            by= c("gene_name", "location", "feature", "catalog", "distance_to_gene", "mC_type")) %>% 
                  arrange(gene_name) %>% distinct(gene_name, .keep_all = T) %>% 
                  select(c(gene_name, mC_type, catalog, feature, mC_dif, hmC_dif,location))
  
  # save common gene table for 5mC_hypo/5hmC hyper
  outFile = paste(out, "/commonGenes_table_mC_hypo_hmC_hyper_", group, "_vs_WT.txt", sep = "")
  write.table(common_mC_hmC, file=outFile, sep = "\t", quote = F, row.names = F)        
  
  
  # common genes in DE and 5mC_hypo/5hmC hyper
  df_commonGene_1 <- inner_join(DE_down, common_mC_hmC, by= "gene_name") %>%  
    mutate(correlation = "positive") %>% 
    mutate(DE_type = "down")
  
  df_commonGene_2 <- inner_join(DE_up, common_mC_hmC, by= "gene_name") %>%  
    mutate(correlation = "negative") %>% 
    mutate(DE_type = "up")
  
  df_merge <- rbind(df_commonGene_1, df_commonGene_2) %>% 
    arrange(correlation, gene_name) 
  
  # append DE_df 
  df_merge$DE_df = DE_data$log2FoldChange[DE_data$gene_name %in% df_merge$gene_name]  
  
  df_merge <- df_merge %>% 
            select(c(gene_name, correlation, DE_type, mC_type, 
           catalog, feature, DE_df, mC_dif, hmC_dif, location))    
    
  # save common gene matrix for DE and 5mC_hypo/5hmC hyper:
  outFile = paste(out, "/commonGenes_matrix_", group, "_vs_WT.txt", sep = "")
  write.table(df_merge, file=outFile, sep = "\t", quote = F, row.names = F)     
    
    
  
  # plot list for common 5mC_hypo and 5hmC_hyper:  
  df_3 <- common_mC_hmC %>% select(gene_name) %>% unlist()
  
  
  # plot
  D <- list("DE_up"= df_1, "DE_down"= df_2, "5mC_hypo/5hmC_hyper"= df_3)
  jpeg(file=paste(out, "/venn_", group, "_DE_down_and_mC_hypo.jpg", sep = ""),
       units="in", 
       width=3, 
       height=3, 
       res=300)
  p <- ggvenn(D,show_percentage=F,
              fill_color=c("pink","orange", "blue"),
              stroke_color = "white",
              stroke_size = 0,
              set_name_size = 4,
              text_size = 6)
  print(p)
  dev.off()
}

print("End!")
# -------------