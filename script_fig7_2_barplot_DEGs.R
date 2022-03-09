#!/usr/local/bin/Rscript
###############
# count of DE genes and plot
##############

# work path
#setwd("~/project_katja2_RNAseq/")

library(tidyverse)
library(qdap) # remove strings

# create out fold
out = "DE_gene_count"
dir.create(out)

# generate count file in terminal and move into fold: dataIn_fig7
  # $ wc -l DE_groups/*/DE_geneList* > sum_transcript_diffExp.txt
  
## input above generated data:
de <- read.table("dataIn_fig7/sum_transcript_diffExp.txt", 
                 quote="\"", comment.char="")

# remove string out_expression_ and /gene_results_diffExp.csv:
de1 <- de[c(1,4,7), ] # Extract DE gene count:
de1[[2]] <- unlist(genXtract(de1[[2]], "DE_groups/", "/DE_geneList_condition")) 
  #pick up str in boundary between "xx" and "yy".

# extract DE up/down gene:
de2 <- de[-c(1,4,7,10), ]
de2[[3]] <- unlist(genXtract(de2[[2]], "DE_groups/", "/DE_geneList"))
de2[[4]] <- unlist(genXtract(de2[[2]], "DE_geneList_", "_condition"))
de3 <- de2[, c(3,1,4)] 
colnames(de3) <- c("group", "count", "type")

# save file:
write.csv(de1, paste(out, "DE_gene_count.csv", sep = '/'), 
          row.names=FALSE)
write.csv(de3, paste(out, "DE_gene_count_by_up_down.csv", sep = '/'), 
          row.names=FALSE)


## plot for total DE genes:
data = de1
y_max= ceiling(max(data$V1)/10)*12

fig <- ggplot(data, aes(V2, V1)) + 
  geom_col(colour = "black")+
  xlab("") +
  ylab("DE gene count") +
  theme_gray(base_size = 10) +
  ylim(c(0,y_max)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_text(aes(label = V1), vjust = -0.5) +
  theme_bw()
fig   
# save plot 
File = paste(out, "/plot_bar_DE_gene_count.jpg", sep="")
ggsave(filename = File, plot = fig_enrichment, 
       width = 4, 
       height = 6, 
       dpi = 300, 
       units = "cm")
dev.off()


## plot for DE genes by up/down:
data <- de3
  y_max= ceiling(max(data$count)/10)*12
  dodge <- position_dodge(width = 0.5)

  fig <- ggplot(data, aes(x= group, y=count, fill = type)) + 
    geom_bar(stat = "identity", position = position_dodge()) +
    xlab(NULL) +
    ylab("DE genes count") +
    ylim(c(0,y_max)) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    theme(text = element_text(size = 12)) +
    theme_bw( ) +
    scale_fill_manual("legend", values = c("down" = "blue", "up" = "orange"))
# save plot 
  print(fig)
  File = paste(out, "/plot_bar_DE_genes_count_by_up_down.jpg", sep="")
  ggsave(filename = File, plot = fig_enrichment, 
       width = 6, 
       height = 6, 
       dpi = 300, 
       units = "cm")

print("End !")
#-------------