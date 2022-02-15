#!/usr/bin/sh
####################################
# DMC_genes_extraction  
# v2_2021_11_08
####################################

# set work location:
#setwd("~/project_katja2/methylkit_2th/6M_Neil1_vs_WT")

# create output fold:
dir.create("./out3_methylkit_DM_genes")
library("biomaRt")

# database for mouse
host = c('http://jul2019.archive.ensembl.org') 
ensembl_mouse = useMart(biomart="ensembl", dataset = "mmusculus_gene_ensembl",
                        host = host)

# input data
DMC <- read.table(file = "./out2_methylkit_DMR/DMC_diffAnn_TSS_plus_Diff25_KO_vs_WT_chrAll.txt", 
                  sep = "\t", header = T)

transcriptID_list <- DMC$feature.name
transcriptID_list <- as.character(transcriptID_list)

# unique transcriptID in file:
uniq_transID <- unique(transcriptID_list)
geneName_list <- getBM(attributes=c("ensembl_transcript_id_version", "ensembl_gene_id",
                                    "mgi_symbol", "entrezgene_description"),
                       filters = "ensembl_transcript_id_version",
                       values = transcriptID_list,
                       mart = ensembl_mouse)
# rename col "ensembl_transcript_id_version" as "feature.name":
colnames(geneName_list)[1] <- "feature.name"
head(geneName_list,3)

DMC_merge <- merge(DMC,geneName_list, by = "feature.name") #before/after merge: 10174/10183

DMC_merge_format <- DMC_merge[, c(1:4, 6:8,10,12:14)]

DMC_uniq_transcripts <- DMC_merge_format[!duplicated(DMC_merge_format$feature.name), ]
DMC_uniq_transcripts_hyper <- DMC_uniq_transcripts[DMC_uniq_transcripts$meth.diff > 0, ]
DMC_uniq_transcripts_hypo <- DMC_uniq_transcripts[DMC_uniq_transcripts$meth.diff < 0, ]

DMC_uniq_genes <- DMC_merge_format[!duplicated(DMC_merge_format$mgi_symbol), ]
DMC_uniq_genes_hyper <- DMC_uniq_genes[DMC_uniq_genes$meth.diff > 0, ]
DMC_uniq_genes_hypo <- DMC_uniq_genes[DMC_uniq_genes$meth.diff < 0, ]

write.table(DMC_uniq_transcripts, file="./out3_methylkit_DM_genes/DMC_transcripts_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMC_uniq_transcripts_hyper, file="./out3_methylkit_DM_genes/DMC_transcripts_hyper_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMC_uniq_transcripts_hypo, file="./out3_methylkit_DM_genes/DMC_transcripts_hypo_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

write.table(DMC_uniq_genes, file="./out3_methylkit_DM_genes/DMC_genes_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMC_uniq_genes_hyper, file="./out3_methylkit_DM_genes/DMC_genes_hyper_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMC_uniq_genes_hypo, file="./out3_methylkit_DM_genes/DMC_genes_hypo_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

# sum merge:
t1 <- nrow(DMC)  # total DMC transcirptID
x1<- nrow(DMC_uniq_transcripts)
x2 <- nrow(DMC_uniq_transcripts_hyper)
x3 <- nrow(DMC_uniq_transcripts_hypo)

x4 <- nrow(DMC_uniq_genes)
x5 <- nrow(DMC_uniq_genes_hyper)
x6 <- nrow(DMC_uniq_genes_hypo)

# sum_DMC <- data.frame(transcripts_total = x1, transcripts_hyper = x2, transcripts_hypo = x3, genes_total = x4, gene_hyper = x5, genes_hypo = x6)
names <- c("transcripts_total", "transcripts_hyper", "transcripts_hypo", "genes_total", "genes_hyper", "genes_hypo")
counts <- c(x1, x2, x3, x4, x5, x6)
sum_DMC <- data.frame(names, counts)
write.table(sum_DMC, file="./out3_methylkit_DM_genes/sum_DMC_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

# plot
sum_DMC
barplot(counts, names.arg = names)

## stack barplot
# Create the input vectors.
library("RColorBrewer")
colors = brewer.pal(n = 3, name = "RdBu")
y_scale <- as.integer((counts[1])*1.3)

# y_scale for manuscript, 2021_11_08
y_scale = 4000

subtype <- c("hyper", "hypo")
feature <- c("transcripts", "genes")

# Create the input vectors, the data in the order of row by row (each row represent one subtype!)
values <- matrix(c(x2, x5, x3, x6), nrow = 2, ncol = 2, byrow = TRUE)
values
class(values)

# plot
jpeg(file = "./out3_methylkit_DM_genes/barplot_stacked_DMC__KO_vs_WT_chrAll.jpg", 
	units="in",
        width=4,
        height=5,
        res=300,
	type = "cairo")

barplot(values, main = "DMC", names.arg = feature, cex.axis = 1.4, cex.names = 1.4, 
        ylab = "Counts", cex.lab =1.2, col = colors,
        ylim = c(0, y_scale))
legend("topright", subtype, cex = 1.2, fill = colors, bty = "n") # no border
dev.off()

# format above plot for manuscript_2021_11_08
jpeg(file = "./out3_methylkit_DM_genes/barplot_stacked_DMC__KO_vs_WT_chrAll_format.jpg", 
     units="in",
     width=4,
     height=5,
     res=300,
     type = "cairo")
barplot(values, main = "", names.arg = feature, cex.axis = 1.4, cex.names = 1.4, 
        ylab = "", cex.lab =1.2, col = colors,
        ylim = c(0, y_scale))
dev.off()

print("End!")
# ---------------
