######
# Format figure 2 pie_plot
######
library("methylKit")

# set work location: in group comparisions: eg.
#setwd("~/project_katja2/methylkit_2th/6M_Neil1_vs_WT") 

# load session RData file from out2, including diffAnn and diffAnn_tiles:
load(file="./out2_methylkit_DMR/session2_DMR.RData")
ls()

 # modified step 4 for pie_plot
# plot for DMC diffAnn, base resolution
jpeg(file="./out2_methylkit_DMR/DMC_plot_DM_annotation_KO_vs_WT_chrAll.jpg",
     units="in", 
     width=4, 
     height=4, 
     res=300, 
     type="cairo")

plotTargetAnnotation(diffAnn,precedence=T,
                     main="DMC annotation by genes",
                     cex.legend = 0.5,
                     radius = 0.5,
                     border = F)
dev.off()

# plot for DMR diffAnn_tiles, per 1000bp window
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

# plot for DMP
jpeg(file="./out2_methylkit_DMR/promoter_plot_DM_annotation_KO_vs_WT_chrAll.jpg",
     units="in", 
     width=4, 
     height=4, 
     res=300, 
     type="cairo")
plotTargetAnnotation(diffAnn_p, precedence=TRUE,
                     main="DMP annotation by genes",
                     cex.legend = 0.5,
                     radius = 0.5,
                     border = F)
dev.off()

print("End!")
# ------------
