#####################
# GO enrichment analysis: over_representation
# update_2022_03_3
#####################
#setwd("~/project_katja2_RNAseq/")

library(clusterProfiler)
library(enrichplot)
library(tidyverse)

library(org.Mm.eg.db)  # annotation using Entrez Gene identifiers
library(stringr)  # for wrap text in barplot
library(ggnewscale)
library(ggridges)

organism = org.Mm.eg.db

# create out fold
Contrast = c("Neil1_vs_WT", 
             "Neil2_vs_WT", 
             "DK_vs_WT")

### for_loop for Contrast:
for (i in 1:length(Contrast)) { print(i)
  cat("Create over_represent fold for Contrast:", Contrast[i])
  out =  paste("DE_over_representation_", Contrast[i], sep = "")
  dir.create(out)
  
# input data
  #Dir <- file.path("~/project_katja2_RNAseq/DE_groups/", Contrast[i])
  Dir <- file.path("dataIn_fig7/DE_groups/", Contrast[i])
  File = paste(Dir, "/DE_data_sigDiff_condition_", Contrast[i], ".txt", sep = "")
  df = read.delim(File)
  data <- df
  head(data,3)
  
# transfer gene name to entrez id:
  keytypes(org.Mm.eg.db)
  ids<- bitr(data$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)
  head(ids,3)
  
  data2 <- data %>% left_join(ids, by= c("gene_name" = "SYMBOL")) %>% 
    distinct(ENTREZID, .keep_all = T) %>% 
    filter(ENTREZID != "NA")
# format gene_list as vector + names:
  gene_list <- data2$log2FoldChange 
  names(gene_list) <- data2$ENTREZID 
  
# sort gene_list
  gene_list <- sort(gene_list, decreasing = TRUE)
  head(gene_list,3)
  
  
# enrichGO analysis, enrichment by only gene_list info:
# set qvalueCutoff= 0.05    
# option in for_loop of j (ont = j): one of “BP”, “MF”, “CC” or “ALL”

for (j in c("BP", "MF", "CC")) {print(j) 
#j = "BP"
# gene: a vector of entrez gene ID.

## GO classfication, set levle = 3
go_class <- groupGO(gene     = names(gene_list),
                    OrgDb    = organism,
                    ont      = j,
                    level    = 3,
                    readable = TRUE) 
head(go_class,2)
length(go_class$ID)

# top GO term table
go_class_top <- data.frame(go_class) %>% dplyr::select(c(1:4)) %>% arrange(desc(Count)) %>% head(10)
print(" GO classification top 10 terms:")
print(go_class_top) 


## GO over-representation test(Boyle et al. 2004):
ego <- enrichGO(gene          = names(gene_list),
                OrgDb         = org.Mm.eg.db,  # OR: OrgDb = organism
                keyType = "ENTREZID",          # required type, default
                ont           = j,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.08,
                qvalueCutoff  = 0.08,
                readable      = TRUE)

head(ego)  
length(ego$ID)

# check length(ego$ID)
if (length(ego$ID) == 0) { print("The enrichGO is null, quit for_loop")
} else {
  print("The enrichGO is  not null, continue ...")
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
  
  # export over-representation result matrix
  x <- as.data.frame(ego_simplify) %>% arrange(p.adjust)
  df_ego <- x[, c(1:7,9,8)]
  File = paste(out, "/GO_matrix_over_representation_", j, ".txt", sep = "")
  write.table(df_ego, file = File, sep = "\t", quote = F, row.names = F, col.names = T)
  
  ## plot 1. barplot by package enrichplot
   p <- barplot(ego_simplify, 
               drop = TRUE, 
               showCategory = 10, 
               title = paste("GO enrichment in ", j, sep = "" ),
               font.size = 12)

  # alternative plot by wrap x lab text:
  # set output number of GO terms by arg m, default m = 10:
   if (nrow(x) == 0) {cat("enrichment matrix is:", nrow(x), "quit! \n")
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
    labs(title=paste("GO enrichment in ", j, sep = "" )) +
    theme(axis.text.x = element_text(color = "black", size = 16, angle = 0),
          axis.text.y = element_text(color = "black", size = 16, angle = 0)) 
  
  File = paste(out, "/boxplot_GO_over_represent_", j, ".jpg", sep = "")
  jpeg(file=File,
       units="in", 
       width=8, 
       height=4, 
       res=300)
  print(p)
  dev.off()
  
  ## plot 2. net-plot
  ## convert gene ID to Symbol,  # default set: showCategory = 4
  ego_net <- setReadable(ego_simplify, 'org.Mm.eg.db', 'ENTREZID')
  p1 <- cnetplot(ego_net, 
                 cex_label_category = 0,
                 cex_label_gene = 2,
                 showCategory = 4,
                 foldChange=gene_list)
  
  File = paste(out, "/netplot_GO_over_represent_", j, ".jpg", sep = "")
  jpeg(file=File,
       units="in", 
       width=10, 
       height=8, 
       res=300)
  print(p1)
  dev.off()
  
  ## plot 3.heatmap
  p <- heatplot(ego_simplify, 
                showCategory = 30,
                foldChange=gene_list)
  File = paste(out, "/heatplot_GO_over_represent_", j, ".jpg", sep = "")
  jpeg(file=File,
       units="in", 
       width=25, 
       height=4, 
       res=300)
  print(p)
  dev.off()
   
   } # end of if_loop: checking length(ego$ID).
  }  # end of for_loop: j in BP, MF and CC
}    # end of for_loop:  i in Contrast of groups

print("End!")
# -----------
