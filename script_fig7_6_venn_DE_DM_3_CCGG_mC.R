###
# overlap DE genes and DM in CCGG context
###

library(tidyverse)
library("ggvenn")

#setwd
#setwd("/Users/mingyi/project_katja2_RNAseq")

# read file
#driver <- "/Volumes/WD_Y4/GWBS_katja2_backup/hmst_katja2_hmC/"
#driver <- "/Volumes/WD_Y5/GWBS_Neil_HC_6M_HMSTseq/hmst_katja2_hmC/"
driver <- "dataIn_fig7/"

### input data
## loop_1 in Feature
Feature = c("TSS", "gene", "TES")
#Feature = c("gene")
for (feature in Feature){print(feature)

Catalog = paste("DMR_CCGG_", feature, sep = "")
out = paste("venn_DE_", Catalog, sep = "")
dir.create(out)

# hyper
if (feature == "gene") {
  pattern = paste(feature, "Body", sep = "")
} else {
  pattern = feature
}

# for_loop_2 in groups: 
for (group in c("Neil1", "Neil2", "DK")) {print(group)
 
inFile = paste(driver, "out/data/table_", group, "_vs_WT_5mC_", pattern, "_imputedWith_zeros_DMRs_hyper.csv", sep = "")
df_5mC_hyper <- read.delim(file = inFile, header = T, strip.white = T, sep = "\t")  # strip white space in gene column ! 
df_5mC_hyper$mC_type <- "hyper"

# hypo
inFile = paste(driver, "out/data/table_", group, "_vs_WT_5mC_", pattern, "_imputedWith_zeros_DMRs_hypo.csv", sep = "")
df_5mC_hypo <- read.delim(file = inFile, header = T, strip.white = T, sep = "\t")  # strip white space in gene column !  
df_5mC_hypo$mC_type <- "hypo"

# rbind
df_5mC <- rbind(df_5mC_hyper, df_5mC_hypo)
df_5mC$catalog <- Catalog
df_5mC$mC_dif <- df_5mC$median_KO - df_5mC$median_WT
Location <- substring(df_5mC$location, 4, 300)
df_5mC$location <- Location
df_5mC$distance_to_gene <- 0


## overlap between DE and DMRs:
  # read DE file
  Dir <- "./DE_groups/"
  df_commonGene_list <- list()
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
    
    
    # pattern in 5mC_hyper
    mC_data <- df_5mC %>% arrange(gene) %>% distinct(gene, .keep_all =T) # OBS, keep 2D data frame by distinct(df, .keep_all =T).
    
    mC_hyper <- mC_data %>% filter(mC_dif >= 0) %>% 
      select(gene) %>% rename(gene_name=gene) %>% 
      arrange(gene_name) %>% 
      distinct(gene_name) 
    df_3 <- mC_hyper %>% unlist()

    # pattern in 5mC_hypo
    mC_hypo <- mC_data %>% filter(mC_dif <= 0) %>% 
      select(gene) %>% rename(gene_name=gene) %>% 
      arrange(gene_name) %>% 
      distinct(gene_name) 
    df_4 <- mC_hypo %>% unlist()
    
    
    # common genes in postive and negative correlation
    row_L1 <- nrow(inner_join(DE_down, mC_hypo, by= "gene_name")) # check length in data
    
    df_commonGene_1 <- inner_join(DE_down, mC_hypo, by= "gene_name") %>%  
      mutate(correlation = "positive") %>% 
      mutate(DE_type = "down", mC_type = "hypo", catalog = Catalog)
    
    df_commonGene_2 <- inner_join(DE_up, mC_hyper, by= "gene_name") %>%  
      mutate(correlation = "positive") %>% 
      mutate(DE_type = "up", mC_type = "hyper", catalog = Catalog)
    
    df_commonGene_3 <- inner_join(DE_up, mC_hypo, by= "gene_name") %>%  
      mutate(correlation = "negative") %>% 
      mutate(DE_type = "up", mC_type = "hypo", catalog = Catalog)
    
    df_commonGene_4 <- inner_join(DE_down, mC_hyper, by= "gene_name") %>%  
      mutate(correlation = "negative") %>% 
      mutate(DE_type = "down", mC_type = "hyper", catalog = Catalog)
    
    df_merge <- rbind(df_commonGene_1, df_commonGene_2, df_commonGene_3, df_commonGene_4) %>% 
      arrange(gene_name)
    
    # append feature and location
    inFile = paste(Dir, group, "_vs_WT/DE_data_sigDiff_condition_", group, "_vs_WT.txt", sep = "")
    DE_data <- read.delim(file = inFile, header = T) %>% arrange(gene_name)
    
    df_merge$feature = mC_data$feature[mC_data$gene %in% df_merge$gene_name]
    df_merge$DE_dif = DE_data$log2FoldChange[DE_data$gene_name %in% df_merge$gene_name]
    
    df_merge$mC_dif = mC_data$mC_dif[mC_data$gene %in% df_merge$gene_name]
    df_merge$location = mC_data$location[mC_data$gene %in% df_merge$gene_name]
    
    # create a new column of distance to gene
    df_merge$distance_to_gene = mC_data$distance_to_gene[mC_data$gene %in% df_merge$gene_name]
    # sort 
    df_merge <- df_merge %>% 
      arrange(correlation, DE_type, gene_name, feature) 
    
    # check if the df_merge data is in the right order in column mC_dif and location as original data.
    # it is important to sort by arrange( ) before merging data by using %in% func.
    #mC_data[mC_data$gene %in% "Rabep1", ]
    
    # save common gene matrix
    outFile = paste(out, "/commonGenes_matrix_", group, "_vs_WT.txt", sep = "")
    write.table(df_merge, file=outFile, sep = "\t", quote = F, row.names = F)
    
    
    # plot
    D <- list("DE_up"= df_1, "DE_down"= df_2, "hyper"= df_3, "hypo"= df_4)
    jpeg(file=paste(out, "/venn_", group, "_vs_WT.jpg", sep = ""),
         units="in", 
         width=3, 
         height=3, 
         res=300)
    p <- ggvenn(D,show_percentage=F,
                fill_color=c("pink","orange", "blue", "green"),
                stroke_color = "white",
                stroke_size = 0,
                set_name_size = 4,
                text_size = 6)
    print(p)
    dev.off()
  } # end of for_lop in group
} # end of for_loop in Feature

print("End!")
# -------------

