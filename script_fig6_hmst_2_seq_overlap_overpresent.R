#####
# GO over-representation analysis
# input: commonGenes in Neil1_vs_WT and Neil2_vs_WT with pattern 5mC_hypo and 5hmC_hyper
#####

# set work directory
# setwd("/Users/mingyi/project_katja2/hmst_5hmC")

library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(purrr)
library(org.Mm.eg.db)  # annotation using Entrez Gene identifiers
library(ggnewscale)
library(ggridges)
library(pathview) # bioconductor
set.seed(188)

#library("org.Dm.eg.db", character.only = TRUE)
organism = org.Mm.eg.db  # OBS, no quote here, if you want use quote, you can do by:organism = get("org.Mm.eg.db")
out = "overpresent_commonN1N2"
dir.create(out)

# input commonGenes file: dataIn_fig6/commonGenes_Neil1_Neil2_TSS* (generated by script_fig6_hmst_1_*)
    File_in <- "dataIn_fig6/commonGenes_Neil1_Neil2_TSS_in_pattern_5mC_hypo_5hmC_hyper.txt"
    data <- read.delim(file = File_in)
    head(data,3)
    
    # set and input mean_dif = 1, and rename colname of "gene_name" to "gene"
    
    data$mean_dif <- 1
    names(data)[names(data)=="gene_name"] <- "gene"
    df <- data
   
     # keep length of df <= 2000 in over-represent analysis
    if (nrow(df) > 2000) {
      df <- df %>%  mutate(abs_mean_dif = abs(mean_dif)) %>%  arrange(desc(abs_mean_dif)) %>% head(2000)
    }
    if (nrow(df) == 0) {
      print("common gene is null, quit for_loop")
      break
    } 
    
    print(paste("The df length is:", nrow(df), sep = ""))
    
    
    # format gene_list as vector with names:
    original_gene_list <- df$mean_dif
    
    # name the vector
    names(original_gene_list) <- df$gene
    
    # Convert gene IDs
    # We will lose some genes here because not all IDs will be converted
    keytypes(organism)
    ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
    
    # remove duplicate IDS 
    dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
    head(dedup_ids)
    head(df,3)
    
    # Create a new dataframe, which has only the genes successfully mapped using the bitr function above
    df2 = df[df$gene %in% dedup_ids$SYMBOL,] %>% distinct(gene, .keep_all = T)
    head(df2,3)
    
    # Create a new column in df2 with the corresponding ENTREZ IDs
    df2$Y = dedup_ids$ENTREZID
    
    # Create a vector of the gene unuiverse
    kegg_gene_list <- df2$mean_dif
    
    # Name vector with ENTREZ ids
    names(kegg_gene_list) <- df2$Y
    head(kegg_gene_list,3)
    
    # omit any NA values 
    kegg_gene_list<-na.omit(kegg_gene_list)
    
    # sort the list in decreasing order 
    kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
    
    # GO Over-representation test (Boyle et al. 2004):
    head(kegg_gene_list,3)
   
    # filter gene and create geneList df 
      gene <- names(kegg_gene_list)[abs(kegg_gene_list) > 0]
      gene_df <- bitr(gene, fromType = "ENTREZID",
                      toType = c("ENSEMBL", "SYMBOL"),
                      OrgDb = org.Mm.eg.db)
      head(gene_df)
      
## GO classfication
      type = c("BP", "MF", "CC")
      for (j in type) {print(j)
      go_class <- groupGO(gene     = gene,
                          OrgDb    = organism,
                          ont      = j,
                          level    = 3,
                          readable = TRUE)
      head(go_class)
      length(go_class$ID)
      
## GO over-representation test:
      Cutoff = 0.05
      ego <- enrichGO(gene          = gene,
                      OrgDb         = organism,
                      ont           = j,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = Cutoff,
                      qvalueCutoff  = Cutoff,
                      readable      = TRUE,
                      keyType       = 'ENTREZID')
      head(ego,3)
      length(ego$ID)
      
      # check length(ego$ID)
      if (length(ego$ID) == 0) { print(" null line in result, quit for_loop")
        next
      } 
      
      print("not null in result, continue ...")
      str = paste("The total length of enrichment ID is : ", length(ego$ID), sep = "")
      print(str)
      
      # filter GO to restrict the result at specific GO level
      ego_filter <- gofilter(ego, level = 6)
      head(ego_filter)
      length(ego_filter$ID)
      
      # check length(ego_filter$ID):
      if (length(ego_filter$ID) == 0) { print(" null line in ego_filter, set it equal to ego")
        ego_filter = ego
      }
      
      # simplify, to remove redundant term:
      ego_simplify <- clusterProfiler::simplify(ego_filter, measure='Wang', semData=NULL) 
      head(ego_simplify)
      length(ego_simplify$ID)
      colnames(as.data.frame(ego_simplify))
      cat("The total filtered and simplified enrichment ID is: ", length(ego_simplify$ID), "\n")
      
      
      # export over-representation result matrix
      x <- as.data.frame(ego_simplify)
      df_ego <- x[, -c(5)]
      File = paste(out, "/GO_matrix_over_representation_", j, ".txt", sep = "")
      write.table(df_ego, file = File, sep = "\t", quote = F, row.names = F, col.names = T)
      
      # input GO_matrix_over_representation file:
      #File = paste(out, "/GO_matrix_over_representation_", j, ".txt", sep = "")
      
      # alternative plot by wrap x lab text:
      #https://stackoverflow.com/questions/21878974/wrap-long-axis-labels-via-labeller-label-wrap-in-ggplot2
      # set output number of GO terms by arg m, default m = 10:
      x <- as.data.frame(ego_simplify) %>% arrange(p.adjust)
      if (nrow(x) == 0) {cat("enrichment matrix is:", nrow(x), "quit! \n")
         next
      } else if (nrow(x) < 10){ m = nrow(x)
      } else { m = 10 }
      
      data_x = x$Description[1:m]
      data_y = x$Count[1:m]
      data_p <- x$p.adjust[1:m]
      df3 <- data.frame(x= data_x, y = data_y, p_adjust = data_p)
      df3$newx = str_wrap(df3$x, width = 50) 
      df3$newx = factor(df3$newx, levels = df3$newx)
      
      p = ggplot(df3, aes(newx, y, fill = p_adjust)) + 
        xlab("") + ylab("Number of genes") +
        geom_bar(stat = "identity") + theme_bw() +
        coord_flip() +
        labs(title=paste("GO terms in ", j, sep = "" )) +
        theme(axis.text.x = element_text(color = "black", size = 16, angle = 0),
              axis.text.y = element_text(color = "black", size = 16, angle = 0)) 
      
      
      File = paste(out, "/boxplot_GO_over_represent_format_", j, ".jpg", sep = "")
      jpeg(file=File,
           units="in", 
           width=8, 
           height=8, 
           res=300)
      print(p)
      dev.off()
      
      
      ## plot 2. net-plot
      ## convert gene ID to Symbol
      ego_net <- setReadable(ego_simplify, 'org.Mm.eg.db', 'ENTREZID')
      p1 <- cnetplot(ego_net, foldChange=kegg_gene_list,
                     showCategory = 5,
                     node_label = "all", # or "gene", "category"
                     node_label_size = 1.4,
                     cex_category = 1.4,
                     cex_gene = 1,
                     cex_label_category = 1.2,
                     cex_label_gene = 1.2)
      
      File = paste(out, "/netplot_GO_over_represent_format_", j, ".jpg", sep = "")
      jpeg(file=File,
           units="in", 
           width=8, 
           height=8, 
           res=300)
      print(p1)
      dev.off()
      
      ## plot 3.heatmap
      p <- heatplot(ego_simplify, foldChange=kegg_gene_list)
      File = paste(out, "/heatplot_GO_over_represent_", j, ".jpg", sep = "")
      jpeg(file=File,
           units="in", 
           width=12, 
           height=4, 
           res=300)
      print(p)
      dev.off()
      
} # end of for_loop for j in BP, MF and CC
  
print("End!")
######################
