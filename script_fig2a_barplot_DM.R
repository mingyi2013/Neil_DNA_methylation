#########
# Figure 2a, barplot for DM count in CG context
#########
library("RColorBrewer")
colors = brewer.pal(n = 2, name = "RdBu")


# data in
# fig3a/DMP_genes_KO_vs_WT_chrAll.txt
work_path = "/Users/mingyi/project_katja2/manuscript_Neil/supplementary/fig2a"
setwd(work_path)
Out = "out"
dir.create(Out)

# DMP
Context = c("DMC", "DMR", "DMP")
for (context in Context) {print(context)

Group = c("Neil1", "Neil2", "DK")
for (group in Group ) {print(group)
  inFile = paste0(group, "_", context,".txt")
  data <- read.table(file = inFile,  sep = "\t", header = T)
  
  count_hyper <- nrow(data[data$meth.diff > 0, ])
  count_hypo <-  nrow(data[data$meth.diff < 0, ])
  
  jpeg(file = paste0(Out, "/barplot_Count_", context, "_", group, "_vs_WT.jpg"), 
       units="in", 
       width=4, 
       height=5, 
       res=300)
  
  # set y-axis scale
  if (context == "DMC") {
    y_scale = 20000
  } else if (context == "DMR") {
    y_scale = 1100
  } else {
    y_scale = 50
  }
  
  barplot(c(count_hyper, count_hypo), main = "", names.arg = c("hyper", "hypo"), cex.axis = 1.4, cex.names = 1.4, 
          ylab = "", cex.lab =1.2, col = colors,
          ylim = c(0, y_scale))
    dev.off() 
  }
}

print("End")
# ----------
