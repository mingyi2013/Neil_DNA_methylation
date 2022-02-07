#!/usr/bin/sh
####################################
# methylKit gene annotation for DMP
####################################

#dir.create("./out3_methylkit_DM_genes")  #create output fold
library("biomaRt")

# database for mouse
ensembl_mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# input data
DMP <- read.table(file = "./out2_methylkit_DMR/promoter_diffAnn_TSS_plus_Diff25_KO_vs_WT_chrAll.txt", 
                  sep = "\t", header = T)

transcriptID_list <- DMP$feature.name
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


DMP_merge <- merge(DMP,geneName_list, by = "feature.name") 
colnames(DMP_merge)
DMP_merge_format <- DMP_merge[, c(1:4, 6:8,10,12:14)]
colnames(DMP_merge_format)

DMP_uniq_transcripts <- DMP_merge_format[!duplicated(DMP_merge_format$feature.name), ]
DMP_uniq_transcripts_hyper <- DMP_uniq_transcripts[DMP_uniq_transcripts$meth.diff > 0, ]
DMP_uniq_transcripts_hypo <- DMP_uniq_transcripts[DMP_uniq_transcripts$meth.diff < 0, ]

DMP_uniq_genes <- DMP_merge_format[!duplicated(DMP_merge_format$mgi_symbol), ]
DMP_uniq_genes_hyper <- DMP_uniq_genes[DMP_uniq_genes$meth.diff > 0, ]
DMP_uniq_genes_hypo <- DMP_uniq_genes[DMP_uniq_genes$meth.diff < 0, ]

write.table(DMP_uniq_transcripts, file="./out3_methylkit_DM_genes/DMP_transcripts_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMP_uniq_transcripts_hyper, file="./out3_methylkit_DM_genes/DMP_transcripts_hyper_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMP_uniq_transcripts_hypo, file="./out3_methylkit_DM_genes/DMP_transcripts_hypo_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

write.table(DMP_uniq_genes, file="./out3_methylkit_DM_genes/DMP_genes_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMP_uniq_genes_hyper, file="./out3_methylkit_DM_genes/DMP_genes_hyper_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMP_uniq_genes_hypo, file="./out3_methylkit_DM_genes/DMP_genes_hypo_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

# sum merge:
t1 <- nrow(DMP)  # total DMP transcirptID

x1<- nrow(DMP_uniq_transcripts)
x2 <- nrow(DMP_uniq_transcripts_hyper)
x3 <- nrow(DMP_uniq_transcripts_hypo)

x4 <- nrow(DMP_uniq_genes)
x5 <- nrow(DMP_uniq_genes_hyper)
x6 <- nrow(DMP_uniq_genes_hypo)

# sum_DMP <- data.frame(transcripts_total = x1, transcripts_hyper = x2, transcripts_hypo = x3, genes_total = x4, gene_hyper = x5, genes_hypo = x6)
names <- c("transcripts_total", "transcripts_hyper", "transcripts_hypo", "genes_total", "genes_hyper", "genes_hypo")
counts <- c(x1, x2, x3, x4, x5, x6)
sum_DMP <- data.frame(names, counts)
write.table(sum_DMP, file="./out3_methylkit_DM_genes/sum_DMP_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

# stack barplot
# Create the input vectors.
library("RColorBrewer")
colors = brewer.pal(n = 2, name = "RdBu")
y_scale <- as.integer((counts[1])*1.3)

subtype <- c("hyper", "hypo")
feature <- c("transcripts", "genes")

# Create the input vectors, the data in the order of row by row (each row represent one subtype!)
values <- matrix(c(x2, x5, x3, x6), nrow = 2, ncol = 2, byrow = TRUE)

# Give the chart file a name
jpeg(file = "./out3_methylkit_DM_genes/barplot_stacked_DMP__KO_vs_WT_chrAll.jpg", 
     units="in",
     width=4,
     height=6,
     res=300,
     type="cairo") 

# Create the barplot
barplot(values, main = "DMP", names.arg = feature, cex.axis = 1.4, cex.names = 1.4, 
        ylab = "Counts", cex.lab =1.4, col = colors,
        ylim = c(0, y_scale))

# Add the legend to the chart
legend("topright", subtype, cex = 1.4, fill = colors, bty = "n") # no border
dev.off()

print("End !")
# ----------
