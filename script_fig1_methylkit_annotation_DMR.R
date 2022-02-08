#!/usr/bin/sh
####################################
# methylKit gene annotation for DMR
####################################

#dir.create("./out3_methylkit_DM_genes")  # create output fold
library("biomaRt")

# database for mouse
ensembl_mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# input data
DMR <- read.table(file = "./out2_methylkit_DMR/tiles_diffAnn_TSS_plus_Diff25_KO_vs_WT_chrAll.txt", 
                  sep = "\t", header = T)

transcriptID_list <- DMR$feature.name
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

DMR_merge <- merge(DMR,geneName_list, by = "feature.name") #before/after merge: 242/237

DMR_merge_format <- DMR_merge[, c(1:4, 6:8,10,12:14)]


DMR_uniq_transcripts <- DMR_merge_format[!duplicated(DMR_merge_format$feature.name), ]
DMR_uniq_transcripts_hyper <- DMR_uniq_transcripts[DMR_uniq_transcripts$meth.diff > 0, ]
DMR_uniq_transcripts_hypo <- DMR_uniq_transcripts[DMR_uniq_transcripts$meth.diff < 0, ]

DMR_uniq_genes <- DMR_merge_format[!duplicated(DMR_merge_format$mgi_symbol), ]
DMR_uniq_genes_hyper <- DMR_uniq_genes[DMR_uniq_genes$meth.diff > 0, ]
DMR_uniq_genes_hypo <- DMR_uniq_genes[DMR_uniq_genes$meth.diff < 0, ]

write.table(DMR_uniq_transcripts, file="./out3_methylkit_DM_genes/DMR_transcripts_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMR_uniq_transcripts_hyper, file="./out3_methylkit_DM_genes/DMR_transcripts_hyper_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMR_uniq_transcripts_hypo, file="./out3_methylkit_DM_genes/DMR_transcripts_hypo_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

write.table(DMR_uniq_genes, file="./out3_methylkit_DM_genes/DMR_genes_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMR_uniq_genes_hyper, file="./out3_methylkit_DM_genes/DMR_genes_hyper_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)
write.table(DMR_uniq_genes_hypo, file="./out3_methylkit_DM_genes/DMR_genes_hypo_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

# sum merge:
t1 <- nrow(DMR)  # total DMR transcirptID

x1<- nrow(DMR_uniq_transcripts)
x2 <- nrow(DMR_uniq_transcripts_hyper)
x3 <- nrow(DMR_uniq_transcripts_hypo)

x4 <- nrow(DMR_uniq_genes)
x5 <- nrow(DMR_uniq_genes_hyper)
x6 <- nrow(DMR_uniq_genes_hypo)


# sum_DMR <- data.frame(transcripts_total = x1, transcripts_hyper = x2, transcripts_hypo = x3, genes_total = x4, gene_hyper = x5, genes_hypo = x6)
names <- c("transcripts_total", "transcripts_hyper", "transcripts_hypo", "genes_total", "genes_hyper", "genes_hypo")
counts <- c(x1, x2, x3, x4, x5, x6)
sum_DMR <- data.frame(names, counts)
write.table(sum_DMR, file="./out3_methylkit_DM_genes/sum_DMR_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

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
jpeg(file = "./out3_methylkit_DM_genes/barplot_stacked_DMR__KO_vs_WT_chrAll.jpg", 
     units="in", 
     width=4, 
     height=6, 
     res=300,
     type="cairo")

# Create the barplot
barplot(values, main = "DMR", names.arg = feature, cex.axis = 1.4, cex.names = 1.4, 
        ylab = "Counts", cex.lab =1.4, col = colors,
        ylim = c(0, y_scale))
    legend("topright", subtype, cex = 1.4, fill = colors, bty = "n") 
dev.off()

print("End !")
# ----------