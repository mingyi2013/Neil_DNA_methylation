#!/usr/bin/sh

###########
# script2_DMR, modified_2021_02_13.
##########
# Annotating DM based on promoters and tiling window
# before running, copy required ref file in current path: refseq.mm10.bed.txt
# set cov.bases = 10 for promters
# set cov.bases = 5 for tiling window

library("methylKit")
#setwd("~/project_katja2/methylkit_2th/6M_Neil1_vs_WT") 
# original sample name for 6M: WT_3,4 and KO_3,4.

# create out fold:
dir.create("./out2_methylkit_DMR")

# load session RData file from out1, including meth, myDiff and myobj_filter:
load(file="./out1_methylkit_DMC/session1_DMC.RData")

#check load files
ls()

# get myDiff25
myobj <- myobj_filter
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

# load library and check columns:
library("genomation")
gene.obj=readTranscriptFeatures("refseq.mm10.bed.txt")  #Ensemble_seq
# gene.obj$introns
# gene.obj$promoters
# gene.obj$exons
# gene.obj$TSSes

# check gene.obj format, column chr is string, chr1-19, chrX, chrY, ChrM
# class(gene.obj)
# seqnames(gene.obj)
# seqlevels(gene.obj)

# check myDiff25p, myDiff25p_GRanges format, column in chr should be in string, ok 
# if it is numeric! convert to chr* by seqlevels.
head(myDiff25p,3)

myDiff25p_GRanges <- as(myDiff25p,"GRanges")
seqnames(myDiff25p_GRanges)

# method for converting numeric chr to string chr*:
#seqlevels(myDiff25p_GRanges) <- c("1"="chr1", "2"="chr2", "3"="chr3","4"="chr4", "5"="chr5",
# "6"="chr6", "7"="chr7", "8"="chr8", "9"="chr9", "10"="chr10",
#"11"="chr11", "12"="chr12", "13"="chr13", "14"="chr14", "15"="chr15", 
#"16"="chr16", "17"="chr17", "18"="chr18", "19"="chr19",
#Y="chrY")  

# percentage of target features (DMC) overlapping with annotation
annotateWithGeneParts(myDiff25p_GRanges, gene.obj)


### Step 1. Regional analysis in promoter (TSS +/- 1000bp):
promoters=regionCounts(myobj, gene.obj$promoters, cov.bases = 10)
head(promoters[[1]])
class(promoters)  #methylRawlist as myobj

# promoter (p) merge:
meth_p=unite(promoters, destrand=T) #destrand=T only for CpG of merging reads in both strands, not for CH.
head(meth_p,3)

## calculate MDRs_p as the difference of treatment (group: 1) - control (group: 0)
myDiff_p=calculateDiffMeth(meth_p)

head(myDiff_p)
write.table(myDiff_p, file="./out2_methylkit_DMR/promoter_Diff_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

bedgraph(myDiff_p, file.name="./out2_methylkit_DMR/promoter_bg_Diff_KO_vs_WT_chrAll.bed", 
         col.name="meth.diff", unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

# Get hyper DM promoter:
myDiff25p.hyper_p=getMethylDiff(myDiff_p,difference=25,qvalue=0.01,type="hyper")
write.table(myDiff25p.hyper_p, file="./out2_methylkit_DMR/promoter_Diff25p_hyper_KO_vs_WT_chrAll.txt", 
            sep = "\t", quote = F, row.names = F)

# get hypo DM promoter
myDiff25p.hypo_p=getMethylDiff(myDiff_p,difference=25,qvalue=0.01,type="hypo")
write.table(myDiff25p.hypo_p, file="./out2_methylkit_DMR/promoter_Diff25p_hypo_KO_vs_WT_chrAll.txt",
            sep = "\t", quote = F, row.names = F)

# get all DM promter
myDiff25p_p=getMethylDiff(myDiff_p,difference=25,qvalue=0.01)
write.table(myDiff25p_p, file="./out2_methylkit_DMR/promoter_Diff25p_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

# get myDiff25p_p bedgarph file:
bedgraph(myDiff25p_p, file.name="./out2_methylkit_DMR/promoter_bg_Diff25p_KO_vs_WT_chrAll.bed", col.name="meth.diff",unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

## count hyper hypo DM promoter:
diff25p_count_p <- diffMethPerChr(myDiff_p, plot=F,qvalue.cutoff=0.01, meth.cutoff=25)
write.table(diff25p_count_p, file="./out2_methylkit_DMR/promoter_count_Diff25p_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)


### Step 2. regional analysis for random tiling window (DMR):
tiles=tileMethylCounts(myobj,win.size=1000,step.size=1000, cov.bases = 5) 
head(tiles[[1]],3)
summary(tiles[[1]])

# the regions covered by all samples are merged to one object
meth_tiles=unite(tiles, destrand=T) #destrand=T only for CpG of merging reads in both strands, not for CpH.
class(meth_tiles)
meth_tiles

# calculate MDRs as the difference of treatment (group: 1) - control (group: 0)
# running about 4 min per chr.
myDiff_tiles=calculateDiffMeth(meth_tiles)
head(myDiff_tiles) 
write.table(myDiff_tiles, file="./out2_methylkit_DMR/tiles_Diff_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

bedgraph(myDiff_tiles, file.name="./out2_methylkit_DMR/tiles_bg_Diff_KO_vs_WT_chrAll.bed", col.name="meth.diff", 
         unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

# Get hyper DMRs:
myDiff25p.hyper_tiles=getMethylDiff(myDiff_tiles,difference=25,qvalue=0.01,type="hyper")
write.table(myDiff25p.hyper_tiles, file="./out2_methylkit_DMR/tiles_Diff25p_hyper_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

# get hypo DMRs
myDiff25p.hypo_tiles=getMethylDiff(myDiff_tiles,difference=25,qvalue=0.01,type="hypo")
write.table(myDiff25p.hypo_tiles, file="./out2_methylkit_DMR/tiles_Diff25p_hypo_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

# get all DMRs
myDiff25p_tiles=getMethylDiff(myDiff_tiles,difference=25,qvalue=0.01)
write.table(myDiff25p_tiles, file="./out2_methylkit_DMR/tiles_Diff25p_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

## count hyper hypo DMRs:
diff25p_count_tiles <- diffMethPerChr(myDiff_tiles, plot=F,qvalue.cutoff=0.01, meth.cutoff=25)
write.table(diff25p_count_tiles, file="./out2_methylkit_DMR/tiles_count_Diff25p_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

# Get bedgrapg for myDiff25p_titles
bedgraph(myDiff25p_tiles, file.name="./out2_methylkit_DMR/tiles_bg_Diff25p_KO_vs_WT_chrAll.bed", col.name="meth.diff", 
         unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")


### Step 3. Convenience functions for annotation objects,  methylation associated with transcripts_ID,promoter/exon/intron/intergenic
## 4a. annotation for DMC myDiff25p;  output file: diffAnn_TSS
diffAnn=annotateWithGeneParts(myDiff25p_GRanges, gene.obj)

# target.row is the row number in myDiff25p
TSS <- getAssociationWithTSS(diffAnn)
write.table(TSS, file="./out2_methylkit_DMR/DMC_diffAnn_TSS_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

TSS_plus_diff25 <- cbind(myDiff25p, TSS)
write.table(TSS_plus_diff25, file="./out2_methylkit_DMR/DMC_diffAnn_TSS_plus_Diff25_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

# It is also desirable to get percentage/number of differentially methylated regions 
diff_stats <- getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

diff2 <- as.data.frame(diff_stats)
colnames(diff2) <- c("Percentage")
diff2$count <- nrow(myDiff25p)*diff2$Percentage/100 +0.5
diff2$count <- as.integer(diff2$count)

write.table(diff2, file="./out2_methylkit_DMR/DMC_diffAnn_stats_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = T)

## 4b Repeat 4a by annotation for DMR myDiff25p_tiles;  output file: tiles_diffAnn_TSS
myDiff25p_tiles_GRanges <- as(myDiff25p_tiles,"GRanges")

diffAnn_tiles=annotateWithGeneParts(myDiff25p_tiles_GRanges, gene.obj)

TSS_tiles <- getAssociationWithTSS(diffAnn_tiles)
write.table(TSS_tiles, file="./out2_methylkit_DMR/tiles_diffAnn_TSS_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

TSS_plus_diff25_tiles <- cbind(myDiff25p_tiles, TSS_tiles)
write.table(TSS_plus_diff25_tiles, file="./out2_methylkit_DMR/tiles_diffAnn_TSS_plus_Diff25_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

diff_stats_tiles <- getTargetAnnotationStats(diffAnn_tiles,percentage=TRUE,precedence=TRUE)

diff2_tiles <- as.data.frame(diff_stats_tiles)
colnames(diff2_tiles) <- c("Percentage")
diff2_tiles$count <- nrow(myDiff25p_tiles)*diff2_tiles$Percentage/100 +0.5
diff2_tiles$count <- as.integer(diff2_tiles$count)

write.table(diff2_tiles, file="./out2_methylkit_DMR/tiles_diffAnn_stats_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = T)

## 4c Repeat 4a by annotation for DMR myDiff25p_p;  output file: promoter_diffAnn_TSS
# generate GRange data, chr is numeric:
myDiff25p_p_GRanges <- as(myDiff25p_p,"GRanges")

# use genmonation package func:
diffAnn_p=annotateWithGeneParts(myDiff25p_p_GRanges, gene.obj)

TSS_p <- getAssociationWithTSS(diffAnn_p)
write.table(TSS_p, file="./out2_methylkit_DMR/promoter_diffAnn_TSS_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

TSS_plus_diff25_p <- cbind(myDiff25p_p, TSS_p)
write.table(TSS_plus_diff25_p, file="./out2_methylkit_DMR/promoter_diffAnn_TSS_plus_Diff25_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = F)

diff_stats_p <- getTargetAnnotationStats(diffAnn_p,percentage=TRUE,precedence=TRUE)

diff2_p <- as.data.frame(diff_stats_p)
colnames(diff2_p) <- c("Percentage")
diff2_p$count <- nrow(myDiff25p_p)*diff2_p$Percentage/100 +0.5
diff2_p$count <- as.integer(diff2_p$count)
write.table(diff2_p, file="./out2_methylkit_DMR/promoter_diffAnn_stats_KO_vs_WT_chrAll.txt", sep = "\t", quote = F, row.names = T)


### Step 4. Plot the percentage of differentially methylated bases overlapping with exon/intron/promoters
## 4a plot for DMC diffAnn, base resolution
jpeg(file="./out2_methylkit_DMR/DMC_plot_DM_annotation_KO_vs_WT_chrAll.jpg",
     units="in", 
     width=4, 
     height=4, 
     res=300, 
     type="cairo")
plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="DMC annotation by genes",
                     cex.legend = 0.5,
                     radius = 0.5,
                     border = F)
dev.off()

## 4b plot for DMR diffAnn_tiles, per 1000bp window
jpeg(file="./out2_methylkit_DMR/tiles_plot_DM_annotation_KO_vs_WT_chrAll.jpg",
     units="in", 
     width=4, 
     height=4, 
     res=300, 
     type="cairo")
plotTargetAnnotation(diffAnn_tiles,precedence=TRUE,
                     main="DMR annotation by genes",
                     cex.legend = 0.5,
                     radius = 0.5,
                     border = F)
dev.off()

### save session
save(myobj,myDiff,
     myDiff_tiles, myDiff_p,
     diffAnn, diffAnn_tiles, diffAnn_p,
     file="./out2_methylkit_DMR/session2_DMR.RData")

print("End !")
# -----------------
