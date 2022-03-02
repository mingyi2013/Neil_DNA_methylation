#!/bin/sh/Rscript
##########################
# DE analysis with DESeq2
# out: normalized count, PCA, DEGs and heatmap
##########################

library(tximport)
library(tximportData)
library(DESeq2)
library(tidyverse)
library(apeglm) # Approximate posterior estimation for GLM coefficients
library(pheatmap)
library(dendextend)

setwd("~/project_katja2_RNAseq/")
out= "DE_count" # read count for each transcript or gene, before and after normalization by DESeq.
dir.create(out)

# set DE adj_p threhold
p.threshold <- 0.05

##########
# 1. import data from stringTie
##########
#df <- read.delim("~/project_katja2_RNAseq/ballgown/WT_1/t_data.ctab") # test
dir <- "~/project_katja2_RNAseq/ballgown"
list.files(dir)

## create a vector of multiple files (path and name)
#sample_list <- read.table("input_alignment_count/sum_align_sample_list.txt")
sample_list <- list.files(dir) %>% sort()  # sort in ascend order

files <- file.path(dir, sample_list, "t_data.ctab") # create a vector of file list with path.
names(files) <- sample_list   # give the name for correlated files

# create a reference df of tx2gene containing transcript name (t_name) and gene_name, from one of the StringTie file:
tmp <- read.delim(files[1]) # use info from the first file.
tx2gene <- tmp[, c("t_name", "gene_name")]


## -- create and save count_byGene file:
# The coverage value is reversed into original count by formula: count = cov * (gene_length/readLength)
txi <- tximport(files, type = "stringtie", tx2gene = tx2gene, readLength = 100) # default readlength = 75
names(txi)
# show counts:
head(txi$counts, 3) # estimate gene count from the transcript abundance.
head(txi$abundance, 3) # calculated from coverage value.
head(txi$length, 3)

# remove zero and dot in all samples:
x <- as.data.frame(txi$counts)
txi2 <- x[!rowSums(x)== 0 & ! row.names(x) == ".", ] 

# convert rownames to first column:
txi2_b <- tibble::rownames_to_column(txi2, "gene_name") 

# save count_byGene file:
File = paste(out, "/count_byGene_all_samples.txt", sep = "")
write.table(txi2_b, file = File, sep = "\t", quote = F, row.names = F)
# --


## -- create and save count_byTranscript file:
txi_tr <- tximport(files, type = "stringtie",  
                   txOut = TRUE, readLength = 100) # default readlength = 75
# remove zero and dot in all samples:
x <- as.data.frame(txi_tr$counts)
txi3 <- x[!rowSums(x)== 0 & ! row.names(x) == ".", ] 
# convert rownames to first column:
txi3_b <- tibble::rownames_to_column(txi3, "trans_name") 
# save count_byTranscript file:
File = paste(out, "/count_byTranscript_all_samples.txt", sep = "")
write.table(txi3_b, file = File, sep = "\t", quote = F, row.names = F)
# --

## input group info: condition
pheno_data <- read.csv("data_in/list_phenodata.csv", sep = ",") %>% arrange(ids)
pheno_data$ids

# check colnames of input matrix txi, and rename them coordinatedly to make sure the same order as pheno_data$ids.
colnames(txi$counts)
all.equal(colnames(txi$counts), pheno_data$ids)  # True


##########
# 2. DESeq2, extract normalized count
##########
# create a colData file (here: pData_1) from phenotype info:  
  samples <- pheno_data$ids
  condition <- pheno_data$phenotype
  pData_1 <-data.frame(samples= samples, condition= condition)
  pData_1$condition <- factor(pData_1$condition, 
                              levels = c("WT", "Neil1", "Neil2","DK"))
  pData_1$condition
  # OBS, the default ref group is the first item in "Levels: ...". Set the order by levels()
  
## input data into DESeq, method 1: from Tximport:  
  dds <- DESeqDataSetFromTximport(txi,
                                  colData = pData_1,
                                  design = ~ condition)
  
# remove low count of less than 10 in sum:
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
# export count after filtering and normalization by DESeq
dds <- estimateSizeFactors(dds) # correcting for library size
df_normCount <- counts(dds, normalized = T) # return counts after normalizing by size factor and dipersion factor.
head(df_normCount, 3)
df_normCount[rownames(df_normCount)== "Gapdh", ] # check example gene

df_normCount2 <- data.frame(df_normCount) %>% 
  rownames_to_column("gene_name") %>% 
  filter(!gene_name==".") %>% 
  arrange(gene_name) %>% 
  distinct(gene_name, .keep_all=T)

head(df_normCount2, 3)

File=  paste(out,"/data_normalization_count.txt", sep='')
write.table(df_normCount2, file = File, sep = "\t", col.names = T, row.names = F, quote = F )


##########
# 3. DESeq2, PCA for normalized count in all groups
##########
df1 <- df_normCount2[, -1]

# make sure colnames is correlated with condition in order:
colnames(df1)
condition <- pheno_data$phenotype

# convert to log10 scale
df1_log <- log10(df1+1)

# PCA
pca_df1=prcomp(t(df1_log))  # Prinical Comp-onents analysis.
head(pca_df1$x, 3) # check out

pca_df1_perc=round(100*pca_df1$sdev^2/sum(pca_df1$sdev^2),1)
pca_df1_perc

df_pca_df1=data.frame(PC1 = pca_df1$x[,1], PC2 = pca_df1$x[,2], sample = colnames(df1), condition=condition)
df_pca_df1

# ggplot labeled by sample number:
jpeg(file=  paste(out,"/","plot_PCA_allSamples_labeled.jpg", sep=''),
     units="in", 
     width=5, 
     height=4, 
     res=300)
ggplot(df_pca_df1, aes(PC1,PC2, color = condition))+
  geom_point(size=0.5)+
  labs(x=paste0("PC1 (",pca_df1_perc[1],"%)"), y=paste0("PC2 (",pca_df1_perc[2],"%)")) +
  geom_text(
    label=df_pca_df1$sample, 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T, cex=4)
dev.off()

# ggplot with color by conditions:
jpeg(file=  paste(out,"/","plot_PCA_allSamples.jpg", sep=''),
     units="in", 
     width=3, 
     height=2, 
     res=300)
ggplot(df_pca_df1, aes(PC1,PC2, color = condition))+
  geom_point(size=3)+
  labs(x=paste0("PC1 (",pca_df1_perc[1],"%)"), y=paste0("PC2 (",pca_df1_perc[2],"%)")) 
dev.off()


##########
# 4. DESeq2, DE betwen groups
##########
# Create DESeq2 datasets, transfer to interger
data_count <- as.data.frame(lapply(txi2, as.integer))
rownames(data_count) <- rownames(txi2)
  #summary(data_count)

pData_1$condition <- relevel(pData_1$condition, ref = "WT") # set WT as control by manual.

# input data into DESeq, method 2:from matrix
dds <- DESeqDataSetFromMatrix(countData = data_count, 
                              colData = pData_1, 
                              design = ~ condition)
dds2 <-DESeq(dds)

## -- design contrast and run all group comparisons  
# design pair comparision:
control = c("WT", "WT", "WT")
test= c("Neil1", "Neil2", "DK")

for (m in 1:length(control)) {
  info = paste("compare: ", test[m], "_vs_",control[m], sep = "" )
  print(info)
  
  # creat out fold
  out_group = paste(test[m], "_vs_", control[m], sep = "")
  dir.create(out_group)
  #dir.create(file.path("out_DE",out_group), recursive = TRUE)

# set group and ref (for merging data in sigDiff step, and in pheatmap)
group_test <- pheno_data$ids[pheno_data$phenotype == test[m]]
group_ref <- pheno_data$ids[pheno_data$phenotype == control[m]]

# extract DE by DESeq2
resultsNames(dds2) #check resultsNames, and set the contrast with the same name in above list.
Condition = paste("condition_", test[m], "_vs_", control[m], sep = "")  # the same format as resultsNames!

contrast.deseq2 <- list(Condition) 
deseq2_results <- results(dds2, contrast=contrast.deseq2) # "condition_treated_vs_untreated"
summary(deseq2_results, alpha=0.05) # quickly check how many genes with adjusted P value < 0.05.

# shrink and ranking:
resLFC <- lfcShrink(dds2, coef=Condition, type="apeglm")
resLFC
DESeq2::plotMA(resLFC, ylim=c(-2, 2),xlim = c(1, 1000) )

# MA plot
File = paste(out_group,"/plot_MA_", Condition, ".jpg", sep='')
jpeg(file= File ,
     units="in", 
     width=6, 
     height=5, 
     res=300)
DESeq2::plotMA(deseq2_results, ylim=c(-2, 2), main="DESeq2", xlim = c(1, 1000))
dev.off()

# Filter DE genes ID with padj<0.05 
res <- deseq2_results
res$threshold <- as.logical(res$padj < p.threshold)
genes.deseq <- row.names(res)[which(res$threshold)]
length(genes.deseq)  # 169 genes in DK_vs_WT
head(res, 3)

# extract DE data and geneList: padj<0.05 and abs(log2FoldChange)>= 1
de <- as.data.frame(res)
de2 <- de %>% 
  filter(padj <= 0.05 & abs(log2FoldChange)>= 1) %>% 
  tibble::rownames_to_column("gene_name") %>% 
  select(-8) 

# info on matrix of de2 (DE_data)
de2_readme_DE_data <- mcols(res)$description

de2_geneList <- de2 %>% 
  filter(!is.na(gene_name)) %>% 
  distinct(gene_name, .keep_all=T) %>% 
  select(1)

# sigDiff: merge de2, de2_geneList with df_temp2:
df_temp2 <- df_normCount2
df_test <- df_temp2[, group_test]
df_ref <- df_temp2[, group_ref]
df_test_ref <-cbind(df_test, df_ref)
df_test_ref$gene_name <- df_temp2$gene_name

df_merged <- left_join(de2, df_test_ref, by= "gene_name") %>% 
  select(-c("lfcSE", "pvalue", "stat")) %>% 
  distinct(gene_name, .keep_all=T)

# up and down genes
de2_geneList_up <- de2 %>% 
  filter(!is.na(gene_name)) %>%
  distinct(gene_name, .keep_all=T) %>%
  filter(log2FoldChange > 0) %>% 
  select(1)
de2_geneList_down <- de2 %>% 
    filter(!is.na(gene_name)) %>%
    distinct(gene_name, .keep_all=T) %>%
    filter(log2FoldChange < 0) %>% 
    select(1)
File = paste(out_group, "/DE_data_", Condition, ".txt", sep = "")
write.table(de2, file = File, sep = "\t", quote = F, row.names = F)

File = paste(out_group, "/DE_data_", Condition, "_readme.txt", sep = "")
write.table(de2_readme_DE_data, file = File, sep = "\t", quote = F, row.names = F, col.names = F)

File = paste(out_group, "/DE_data_sigDiff_", Condition, ".txt", sep = "")
write.table(df_merged, file = File, sep = "\t", quote = F, row.names = F)


File = paste(out_group, "/DE_geneList_", Condition, ".txt", sep = "")
write.table(de2_geneList, file = File, sep = "\t", quote = F, row.names = F, col.names = F)

File = paste(out_group, "/DE_geneList_up_", Condition, ".txt", sep = "")
write.table(de2_geneList_up, file = File, sep = "\t", quote = F, row.names = F, col.names = F)

File = paste(out_group, "/DE_geneList_down_", Condition, ".txt", sep = "")
write.table(de2_geneList_down, file = File, sep = "\t", quote = F, row.names = F, col.names = F)

## plot counts for genes
# plot for top 5 up genes:
de2_geneList_up <- de2 %>% 
  filter(!is.na(gene_name)) %>%
  filter(log2FoldChange > 0) %>% 
  distinct(gene_name, .keep_all=T) %>% 
  arrange(padj) %>% 
  select(1,3,7)
head(de2_geneList_up,3)

if (nrow(de2_geneList_up) < 6) {
  cat("DE_up gene length:",nrow(de2_geneList_up))
     n =  nrow(de2_geneList_up)
} else {
    n = 5
  }

for (i in 1:n) {print(i)
  padj <- de2_geneList_up$padj[i]
  geneID <- de2_geneList_up$gene_name[i]
  File = paste(out_group, "/plot_count_gene_topUp_", i, "_", geneID, ".jpg", sep = "")
  jpeg(file=  File,
       units="in", 
       width=8, 
       height=6, 
       res=300)
  plotCounts(dds, 
             gene=geneID, 
             intgroup="condition",
             returnData = F, 
             #ylim = c(0.1, 1000), 
             col = c("red"), cex.main=1.3, cex.axis=1.3, cex.lab=1.4, cex.sub=1.3,
             main = paste(test[m], "_vs_", control[m], "; Gene:", geneID, "; p_adj:", padj, sep = ""))
 dev.off()
}

# plot for top 5 down:
de2_geneList_down <- de2 %>% 
  filter(!is.na(gene_name)) %>%
  filter(log2FoldChange < 0) %>% 
  distinct(gene_name, .keep_all=T) %>% 
  arrange(padj) %>% 
  select(1,3,7)
head(de2_geneList_down)

if (nrow(de2_geneList_down) < 6) {
  cat("DE_up gene length:",nrow(de2_geneList_down))
  n =  nrow(de2_geneList_down)
} else {
  n = 5
}

for (i in 1:n) {print(i)
  padj <- de2_geneList_up$padj[i]
  geneID <- de2_geneList_down$gene_name[i]
  File = paste(out_group, "/plot_count_gene_topDown_", i, "_", geneID, ".jpg", sep = "")
  jpeg(file=  File,
       units="in", 
       width=8, 
       height=6, 
       res=300)
  plotCounts(dds, gene=geneID, 
             intgroup="condition",
             returnData = F, 
             #ylim = c(min_y, 1000), 
             col = c("red"), cex.main=1.3, cex.axis=1.3, cex.lab=1.4, cex.sub=1.3,
             main = paste(test[m], "_vs_", control[m], "; Gene:", geneID, "; p_adj:", padj, sep = ""))
  dev.off()
}

#######################
# heatmap for df_merged data
########################
# format data
data_subset1 <- df_merged %>% 
      tibble::column_to_rownames(var="gene_name") %>%
      select(-c("log2FoldChange", "padj", "baseMean"))
head(data_subset1)
# pheatmap with absolute value
pheatmap(data_subset1, show_rownames = F, cluster_cols = F, cluster_rows = F)

# Convert to z score:
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)}
data_subset_norm <- t(apply(data_subset1, 1, cal_z_score)) %>% 
                    na.omit()
  
# simple plot 
File = paste(out_group, "/heatmap_without_label.jpg", sep="")
jpeg(File)
pheatmap(data_subset_norm, show_rownames = F, cluster_cols = T, cluster_rows = T,
         units="in", 
         width= 2, 
         height=3,
         cex= 1.6,
         res=300)
dev.off()

# hierarchical clustering for annaotation
my_hclust_gene <- hclust(dist(data_subset1), method = "complete")

as.dendrogram(my_hclust_gene) %>%
  plot(horiz = T, cex=0.3)

# set cluster K=2 (up and down)
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)
my_gene_col

#
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, 
                                           yes = "cluster 1", no = "cluster 2"))
head(my_gene_col)

# add multiple row annotations to a heatmap
set.seed(1984)
my_random <- as.factor(sample(x = 1:2, size = nrow(my_gene_col), replace = TRUE))
my_gene_col$random <- my_random

# crate col-data:
#my_sample_col <- data.frame(sample = rep(c("Neil1", "WT"), c(4,4))) # set col data by manual
my_sample_col <- data.frame(sample = rep(c(test[m], control[m]),c(2,2))) # set replicates = 2 by manual !
row.names(my_sample_col) <- colnames(data_subset1)

# plot with col and row annoation:
my_heatmap <- pheatmap(data_subset_norm, annotation_row = my_gene_col, 
                       annotation_col = my_sample_col, fontsize_row = 7)  
save_pheatmap_png <- function(x, filename, width=800, height=3000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
 }

File = paste(out_group, "/heatmap_label_col_samples.png", sep="")
save_pheatmap_png(my_heatmap, File)

# without row_annotation
data_subset_norm2 <- data_subset_norm
row.names(data_subset_norm2) <- NULL
my_heatmap2 <- pheatmap(data_subset_norm2,  
                       annotation_col = my_sample_col, fontsize_row = 7)  
save_pheatmap_png <- function(x, filename, width=500, height=700, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
  }
File = paste(out_group, "/heatmap_label_col.png", sep="")
save_pheatmap_png(my_heatmap2, File)

} # end for_loop in DEseq group comparisons.

print("End!")
#-----------
