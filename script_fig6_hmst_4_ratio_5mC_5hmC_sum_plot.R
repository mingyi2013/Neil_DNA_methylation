####
# ratio of 5mC/5hmC in TSS, boxplot for sum Neil1_vs_WT, Neil2_vs_WT for common 453 TSS, plus DK_vs_WT.
###

library(tidyverse)

#setwd("/Users/mingyi/project_katja2/hmst_5hmC")
setwd("./dataIn_fig6/")  # in github

out = "ratio_5mC_5hmC"
# read file
df_N1 <- read.delim("ratio_5mC_5hmC/matrix_ratio_5mC_5hmC_Neil1_vs_WT.txt")
df_N2 <- read.delim("ratio_5mC_5hmC/matrix_ratio_5mC_5hmC_Neil2_vs_WT.txt")
df_DK <- read.delim("ratio_5mC_5hmC/matrix_ratio_5mC_5hmC_DK_vs_WT.txt")
head(df_N2,3)

# extract common index of location
Index <- inner_join(df_N1, df_N2, by = "location") %>% distinct(location) %>% select(location)

# extract df inside index
df_WT <- df_N1 %>% filter(group=="mC_vs_hmC_WT") 
df_N1b <- df_N1 %>% filter(group=="mC_vs_hmC_KO")
df_N2b <- df_N2 %>% filter(group=="mC_vs_hmC_KO")
df_DKb <- df_DK %>% filter(group=="mC_vs_hmC_KO")

WT <- df_WT[df_WT$location %in% Index$location, ]
Neil1 <- df_N1b[df_N1b$location %in% Index$location, ]
Neil2 <- df_N2b[df_N2b$location %in% Index$location, ]
DK <- df_DKb[df_DKb$location %in% Index$location, ]


WT$group <- "WT"
Neil1$group <- "Neil1"
Neil2$group <- "Neil2"
DK$group <- "DK"
  
df_m <- rbind(WT, Neil1, Neil2, DK)
# Basic box plot
  p <- ggplot(df_m, aes(x=group, y=ratio)) + 
    geom_boxplot(aes(colour = group), outlier.shape = NA)
  
# compute lower and upper whiskers, to remove outliers:
# ref:https://stackoverflow.com/questions/5677885/ignore-outliers-in-ggplot2-boxplot
  ylim1 = boxplot.stats(df_m$ratio)$stats[c(1, 5)]
# scale y limits based on ylim1
  p2 = p + coord_cartesian(ylim = c(0, ylim1[2]*1.05) )
  
# Box plot with jittered points
  p3 <- p2 + geom_jitter(shape=16, position=position_jitter(0.25), size = 0.8, aes(colour = group)) +
    labs(title="",x="", y = "Ratio 5mC/5hmC")+
    theme_classic() +
    theme(legend.position = "None")
  jpeg(file=paste(out,"/boxplot_ratio_5mC_by_5hmC_TSS_sum_Neil1_and_Neil2.jpg", sep = ""),
       units="in", 
       width=3, 
       height=3, 
       res=300)
  print(p3)
  dev.off()

# t test for Neil1
  t_test1 <- t.test(Neil1$ratio, WT$ratio)
  chars <- capture.output(print(t_test1))
  File = paste(out, "/t_test_common_5mC_hypo_and_5hmC_hyper_TSS_sum_Neil1_vs_WT.txt", sep = "")
  writeLines(chars, con = file(File))
  
# t test for Neil2
  t_test1 <- t.test(Neil2$ratio, WT$ratio)
  chars <- capture.output(print(t_test1))
  File = paste(out, "/t_test_common_5mC_hypo_and_5hmC_hyper_TSS_sum_Neil2_vs_WT.txt", sep = "")
  writeLines(chars, con = file(File)) 

# t test for DK
  t_test1 <- t.test(DK$ratio, WT$ratio)
  chars <- capture.output(print(t_test1))
  File = paste(out, "/t_test_common_5mC_hypo_and_5hmC_hyper_TSS_sum_DK_vs_WT.txt", sep = "")
  writeLines(chars, con = file(File))   

print("End !")
# --------------

