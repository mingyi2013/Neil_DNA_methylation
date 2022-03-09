##########
# venn plot for overlap DE genes among comparisions
##########
library(VennDiagram)
library(RColorBrewer)
library(tidyverse)
 
# setwd
setwd("/Users/mingyi/project_katja2_RNAseq")
out = "venn_DE"
dir.create(out)

# read file
Dir <- "dataIn_fig7/DE_groups/"
# Neil1
inFile = paste(Dir, "Neil1_vs_WT/DE_geneList_down_condition_Neil1_vs_WT.txt", sep = "")
Neil1_down <- read.delim(file = inFile, header = F)
head(Neil1_down)

inFile = paste(Dir, "Neil1_vs_WT/DE_geneList_up_condition_Neil1_vs_WT.txt", sep = "")
Neil1_up <- read.delim(file = inFile, header = F)
head(Neil1_up)

# Neil2
inFile = paste(Dir, "Neil2_vs_WT/DE_geneList_down_condition_Neil2_vs_WT.txt", sep = "")
Neil2_down <- read.delim(file = inFile, header = F)
head(Neil1_down)

inFile = paste(Dir, "Neil2_vs_WT/DE_geneList_up_condition_Neil2_vs_WT.txt", sep = "")
Neil2_up <- read.delim(file = inFile, header = F)
head(Neil2_up)

# DK
inFile = paste(Dir, "DK_vs_WT/DE_geneList_down_condition_DK_vs_WT.txt", sep = "")
DK_down <- read.delim(file = inFile, header = F)
head(DK_down)

inFile = paste(Dir, "DK_vs_WT/DE_geneList_up_condition_DK_vs_WT.txt", sep = "")
DK_up <- read.delim(file = inFile, header = F)
head(DK_up)


### venn for DE_up:
# Generate 3 sets:
library(magrittr)
set1 <- Neil1_down %>% rename(gene_name = V1) %>% arrange(gene_name) %>% distinct(gene_name)
set2 <- Neil2_down %>% rename(gene_name = V1) %>% arrange(gene_name) %>% distinct(gene_name)
set3 <- DK_down %>% rename(gene_name = V1) %>% arrange(gene_name) %>% distinct(gene_name)

# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(set1 %>% unlist(), set2 %>% unlist(), set3 %>% unlist()),
  category.names = c("Neil1" , "Neil2" , "DK"),
  filename = paste(out, '/venn_diagramm_DE_down.png', sep = ""),
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
  cex = .8,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.8,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


### venn for DE_up, repeat above as for DE_down:
# Generate 3 sets:
set1 <- Neil1_up %>% rename(gene_name = V1) %>% arrange(gene_name) %>% distinct(gene_name)
set2 <- Neil2_up %>% rename(gene_name = V1) %>% arrange(gene_name) %>% distinct(gene_name)
set3 <- DK_up %>% rename(gene_name = V1) %>% arrange(gene_name) %>% distinct(gene_name)

# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(3, "Pastel2")

#
venn.diagram(
  x = list(set1 %>% unlist(), set2 %>% unlist(), set3 %>% unlist()),
  category.names = c("Neil1" , "Neil2" , "DK"),
  filename = paste(out, '/venn_diagramm_DE_up.png', sep = ""),
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
  cex = .8,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.8,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

## remove *.log file
# list file
file_clean <- list.files(pattern = ".log", recursive = TRUE)
# remove files
for (j in 1:length(file_clean)) {file.remove(file_clean[j])}

print("End!")
# -------------