######
# script figure 4a: venn diagarm
######

setwd("~/project_katja2/manuscript_Neil/supplementary/")
library(VennDiagram)
library(RColorBrewer)
library(tidyverse)

out = "venn_commonGene_N1_N2_DK"
dir.create(out)

# set input file Path: given input data in: ./dataOut_fig3/geneList_DMRs_*
Path = "./dataOut_fig3/"
#read.delim("dataOut_fig3/geneList_DMRs_CG_filter_0.15_dif_Neil1_vs_WT.txt") # test

Context = c("CG", "CHG", "CHH")
#Context = "CHH" # test

for (context in Context) {print(context)
  
  if (context == "CG") {
    parameter1 = 0.15} else{parameter1= 0.06}
  inFile <- paste(Path, "geneList_DMRs_", context, "_filter_", parameter1, "_dif_Neil1_vs_WT.txt" , sep = "")
  N1 <- read.delim(inFile)
  
  inFile <- paste(Path, "geneList_DMRs_", context, "_filter_", parameter1, "_dif_Neil2_vs_WT.txt" , sep = "")
  N2 <- read.delim(inFile)
  
  inFile <- paste(Path, "geneList_DMRs_", context, "_filter_", parameter1, "_dif_DK_vs_WT.txt" , sep = "")
  DK <- read.delim(inFile)
  
  # Generate 3 sets:
  set1 <- N1 %>% dplyr::select(location) %>% arrange(location) %>% distinct(location) %>% unlist()
  set2 <- N2 %>% dplyr::select(location) %>% arrange(location) %>% distinct(location) %>% unlist()
  set3 <- DK %>% dplyr::select(location) %>% arrange(location) %>% distinct(location) %>% unlist()
  
  # Prepare a palette of 3 colors with R colorbrewer:
  myCol <- brewer.pal(3, "Pastel2")
  
  ## Chart for overlap among Neil1_vs_WT, Neil2_vs_WT and DK_vs_WT:
  venn.diagram(
    x = list(set1 , 
             set2 , 
             set3),
    category.names = c("" , "" , ""),
    filename = paste(out, "/venn_", context, ".png", sep = ""),
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )

# extract matrix of common genes 
common1 <-  N1 %>% 
            inner_join(N2, by = c("location", "gene", "feature")) %>% 
            inner_join(DK, by = c("location", "gene", "feature") ) 
# re-order columns and rename colname as "*.z" for data from DK_vs_WT.  
        # "*.x": data from Neil1_vs_WT1;  "*.y" data from Neil2_vs_WT.
common1b <- common1 %>% 
            dplyr::select(c(location, gene, feature),everything()) %>% 
            dplyr::rename(mean_KO.z = mean_KO) %>% 
            dplyr::rename(mean_WT.z = mean_WT) %>%
            dplyr::rename(mean_dif.z = mean_dif) %>%
            dplyr::rename(mean_pval.z = mean_pval) 
          
common1b %>% arrange(gene) %>% distinct(gene) %>% nrow()

print( paste("total common genes in contex of", context, " : ", nrow(common1),sep= "") )
         
File_out <- paste(out, "/commonGenes_", context, ".txt", sep = "")
write.table(common1b, file = File_out, row.names = F, quote = F, sep = "\t" )

}  # end for-loop in Context

print("End!")