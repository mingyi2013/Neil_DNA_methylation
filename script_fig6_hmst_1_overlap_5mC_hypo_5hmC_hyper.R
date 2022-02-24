#####
# overlap genes between 5mC_hypo and 5hmC_hyper, venn plot for figure 6A-B
#####

library(VennDiagram)
library(RColorBrewer)
library(tidyverse)
library("ggvenn")

# ref: https://www.r-graph-gallery.com/14-venn-diagramm.html
# set work directory:
#setwd("/Users/mingyi/project_katja2/hmst_5hmC")

### read file
#driver <- "/Volumes/WD_Y4/GWBS_katja2_backup/hmst_katja2_hmC/out/data/" #original

driver="dataIn_fig6/" # data in github

# Neil1
inFile = paste(driver, "Neil1_vs_WT_5mC_TSS_imputedWith_zeros_DMRs_hypo_geneNames.csv", sep = "")
Neil1_mC_TSS_hypo <- read.csv(file = inFile, header = T)
head(Neil1_mC_TSS_hypo)

inFile = paste(driver, "Neil1_vs_WT_5hmC_TSS_imputedWith_zeros_DMRs_hyper_geneNames.csv", sep = "")
Neil1_hmC_TSS_hyper <- read.csv(file = inFile, header = T)

# Neil2
inFile = paste(driver, "Neil2_vs_WT_5mC_TSS_imputedWith_zeros_DMRs_hypo_geneNames.csv", sep = "")
Neil2_mC_TSS_hypo <- read.csv(file = inFile, header = T)

inFile = paste(driver, "Neil2_vs_WT_5hmC_TSS_imputedWith_zeros_DMRs_hyper_geneNames.csv", sep = "")
Neil2_hmC_TSS_hyper <- read.csv(file = inFile, header = T)

# DK
inFile = paste(driver, "DK_vs_WT_5mC_TSS_imputedWith_zeros_DMRs_hypo_geneNames.csv", sep = "")
DK_mC_TSS_hypo <- read.csv(file = inFile, header = T)

inFile = paste(driver, "DK_vs_WT_5hmC_TSS_imputedWith_zeros_DMRs_hyper_geneNames.csv", sep = "")
DK_hmC_TSS_hyper <- read.csv(file = inFile, header = T)

# Generate 3 sets:
set1 <- Neil1_mC_TSS_hypo %>% arrange(gene_name) %>% distinct(gene_name)
set2 <- Neil2_mC_TSS_hypo %>% arrange(gene_name) %>% distinct(gene_name)
set3 <- DK_mC_TSS_hypo %>% arrange(gene_name) %>% distinct(gene_name)

# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(3, "Pastel2")


### Chart for overlap among Neil1,2 and DK in 5mC_hypo TSS
venn.diagram(
  x = list(set1 %>% select(gene_name) %>% unlist(), set2 %>% unlist(), set3 %>% unlist()),
  category.names = c("Neil1" , "Neil2" , "DK"),
  filename = 'venn_diagramm_5mC_TSS_hypo.png',
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


### repeat for overlap among Neil1,2 and DK in 5hmC_hyper TSS
# Generate 3 sets:
set1 <- Neil1_hmC_TSS_hyper
set2 <- Neil2_hmC_TSS_hyper
set3 <- DK_hmC_TSS_hyper

# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(set1 %>% select(gene_name) %>% unlist(), set2 %>% unlist(), set3 %>% unlist()),
  category.names = c("Neil1" , "Neil2" , "DK"),
  filename = 'venn_diagramm_5hmC_TSS_hyper.png',
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


### overlap between 5mC_hypo and 5hmC_hyper in TSS, Figure 6B
df_commonGene_list <- list()
for (group in c("Neil1", "Neil2", "DK")) {print(group)
# 5mC hypo
  inFile = paste(driver, group, "_vs_WT_5mC_TSS_imputedWith_zeros_DMRs_hypo_geneNames.csv", sep = "")
  mC_TSS_hypo <- read.csv(file = inFile, header = T)
  head(mC_TSS_hypo)
  df_1 <-  mC_TSS_hypo %>% 
          arrange(gene_name) %>% 
          distinct(gene_name) %>% unlist()
# 5hmC hyper
  inFile = paste(driver, group, "_vs_WT_5hmC_TSS_imputedWith_zeros_DMRs_hyper_geneNames.csv", sep = "")
  hmC_TSS_hyper <- read.csv(file = inFile, header = T)
  head(hmC_TSS_hyper)
 
  df_2 <- hmC_TSS_hyper %>% 
  arrange(gene_name) %>% 
  distinct(gene_name) %>% unlist()

# common TSS genes between 5mC_hypo and 5hmC_hyper
df_commonGene <- inner_join(mC_TSS_hypo, hmC_TSS_hyper, by= "gene_name") %>% arrange(gene_name) %>% distinct(gene_name)
colnames(df_commonGene) <- group
head(df_commonGene)

# merged df of common TSS genes in each group
newList <- list(df_commonGene[[1]])
names(newList) <- group
df_commonGene_list <- append(df_commonGene_list, newList)


# plot
D <- list("5mC_hypo"= df_1, "5hmC_hyper"= df_2)
jpeg(file=paste("venn_", group, "_TSS_5mC_hypo_vs_5hmC_hyper.jpg", sep = ""),
       units="in", 
       width=4, 
       height=4, 
       res=300)
p <- ggvenn(D,show_percentage=T,
       fill_color=c("pink","orange"))
print(p)
dev.off()
}


### common genes between Neil1_vs_WT and Neil2_vs_WT in TSS with pattern of 5mC_hyper/5hmC_hypo 
N1 <- df_commonGene_list[names(df_commonGene_list)== "Neil1"] %>% as.data.frame() %>% rename(gene_name = Neil1)
N2 <- df_commonGene_list[names(df_commonGene_list)== "Neil2"] %>% as.data.frame() %>% rename(gene_name = Neil2)
DK <- df_commonGene_list[names(df_commonGene_list)== "DK"] %>% as.data.frame() %>% rename(gene_name = DK)

common_N1N2 <- inner_join(N1, N2)

write.table(common_N1N2, file="commonGenes_Neil1_Neil2_TSS_in_pattern_5mC_hypo_5hmC_hyper.txt", sep = "\t", quote = F, row.names = F)

print("End!")
# ------------







