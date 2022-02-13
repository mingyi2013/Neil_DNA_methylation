#!/usr/bin/sh
### script1_DMC, 4samples, KO_1/2 vs WT_1/2. 
library("methylKit")
#setwd("/cluster/projects/nn9383k/mingyiy/methylkit_6M_WTvsNeil1")

# create out fold
dir.create("./out1_methylkit_DMC")

### step 1. Input sorted bam and index file, filter and plot mC profile:

## Input data:
# choose one of "CpG", "CHH" or "CHG" string in read.context
file.list1 = list("/cluster/work/users/mingyiy/index_WT_1/data2_input/WT_1_sorted.bam",
                  "/cluster/work/users/mingyiy/index_WT_2/data2_input/WT_2_sorted.bam",
                  "./data2_input/Neil1_1_sorted.bam",
                  "./data2_input/Neil1_2_sorted.bam")

myobj=processBismarkAln(file.list1, 
                             sample.id= list("WT_1","WT_2","Neil1_1","Neil1_2"),
                             assembly="mm10",
                             treatment=c(0,0,1,1),
                             read.context="CpG",
                             save.folder= "./out1_methylkit_DMC/data_myobj",
                             mincov = 10,
                             minqual = 20)

## Filtering sample based on coverage:
myobj_filter <- filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                   hi.count=NULL,hi.perc=99.9)
myobj_filter 
summary(myobj_filter)

## Plot for mC profile and coverage:
# plot for mC stats:
jpeg(file="./out1_methylkit_DMC/plot_methylation_stats_WT_1.jpg",
        units="in", 
        width=6, 
        height=8, 
        res=300, type='cairo')
getMethylationStats(myobj_filter[[1]],plot=T,both.strands=FALSE)
dev.off()

jpeg(file="./out1_methylkit_DMC/plot_methylation_stats_WT_2.jpg",
        units="in",
        width=6,
        height=8, 
        res=300, type='cairo')
getMethylationStats(myobj_filter[[2]],plot=T,both.strands=FALSE)
dev.off()

jpeg(file="./out1_methylkit_DMC/plot_methylation_stats_KO_1.jpg",
        units="in", 
        width=6, 
        height=8, 
        res=300, type='cairo')
getMethylationStats(myobj_filter[[3]],plot=T,both.strands=FALSE)
dev.off()

jpeg(file="./out1_methylkit_DMC/plot_methylation_stats_KO_2.jpg",
        units="in",
        width=6,
        height=8, 
        res=300, type='cairo')
getMethylationStats(myobj_filter[[4]],plot=T,both.strands=FALSE)
dev.off()

# plot fo coverage:
jpeg(file="./out1_methylkit_DMC/plot_CpG_coverage_WT_1.jpg",
         units="in", 
         width=6, 
         height=8, 
         res=300, type='cairo')
getCoverageStats(myobj_filter[[1]],plot=TRUE,both.strands=FALSE)
dev.off()

jpeg(file="./out1_methylkit_DMC/plot_CpG_coverage_WT_2.jpg",
         units="in",
         width=6,
         height=8, 
         res=300, type='cairo')
getCoverageStats(myobj_filter[[2]],plot=TRUE,both.strands=FALSE)
dev.off()

jpeg(file="./out1_methylkit_DMC/plot_CpG_coverage_KO_1.jpg",
         units="in", 
         width=6, 
         height=8, 
         res=300,type='cairo' )
getCoverageStats(myobj_filter[[3]],plot=TRUE,both.strands=FALSE)
dev.off()

jpeg(file="./out1_methylkit_DMC/plot_CpG_coverage_KO_2.jpg",
         units="in",
         width=6,
         height=8, 
         res=300,type='cairo' )
getCoverageStats(myobj_filter[[4]],plot=TRUE,both.strands=FALSE)
dev.off()

### step2. Make bedgrapg for myObj with mC percentages:
bedgraph(myobj_filter[[1]], file.name="./out1_methylkit_DMC/bg_WT_1.bed", col.name="perc.meth", 
         unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

bedgraph(myobj_filter[[2]], file.name="./out1_methylkit_DMC/bg_WT_2.bed", col.name="perc.meth",
         unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

bedgraph(myobj_filter[[3]], file.name="./out1_methylkit_DMC/bg_KO_1.bed", col.name="perc.meth", 
         unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

bedgraph(myobj_filter[[4]], file.name="./out1_methylkit_DMC/bg_KO_2.bed", col.name="perc.meth",
         unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")


### step 3. Merging samples
# the bases covered by all samples are merged to one object, called methylBase
meth=unite(myobj_filter, destrand=T) #destrand=T only for CpG of merging reads in both strands, not for CH.

# PCA plot for meth:
jpeg(file="./out1_methylkit_DMC/plot_PCA_KO_vs_WT.jpg",
         units="in", 
         width=8, 
         height=8, 
         res=300, type='cairo')
PCASamples(meth, adj.lim = c(0.4, 0.4) )
dev.off()

# plot for clustering:
jpeg(file="./out1_methylkit_DMC/plot_clustering_KO_vs_WT.jpg",
         units="in", 
         width=6, 
         height=4, 
         res=300, type='cairo')
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off()

### step 4. Get myDiff, need bout 1.5 hours for 4 samples. calculate methylation difference as the difference of treatment(group:1) and control(group:0)
myDiff=calculateDiffMeth(meth, mc.cores=2)
write.table(myDiff, file="./out1_methylkit_DMC/DMC_KO_vs_WT.txt", sep = "\t", quote = F, row.names = F)

# myDiff_input <- read.table("./out1_methylkit_DMC/DMC_KO_vs_WT.txt", header = T, sep = "\t")

bedgraph(myDiff, file.name="./out1_methylkit_DMC/bg_DMC_KO_vs_WT.bed", col.name="meth.diff", 
         unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

### step 5. Get DMC, myDiff_hyper and hypo:
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
write.table(myDiff25p, file="./out1_methylkit_DMC/DMC_diff25p_KO_vs_WT.txt", sep = "\t", quote = F, row.names = F)

# get hyper methylated bases
myDiff25p_hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
write.table(myDiff25p_hyper, file="./out1_methylkit_DMC/DMC_diff25p_hyper_KO_vs_WT.txt", sep = "\t", quote = F, row.names = F)

# get hypo methylated bases
myDiff25p_hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
write.table(myDiff25p_hypo, file="./out1_methylkit_DMC/DMC_diff25p_hypo_KO_vs_WT.txt", sep = "\t", quote = F, row.names = F)

### step 6. Get bedgrapg for myDiff25p
bedgraph(myDiff25p, file.name="./out1_methylkit_DMC/bg_DMC_diff25p_KO_vs_WT.bed", col.name="meth.diff", 
         unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

##save session
save(myobj_filter, meth, myDiff, file="./out1_methylkit_DMC/session1_DMC.RData")

print("End !")
# -------------