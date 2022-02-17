##########
# ggplot for mC profile
# 2021_01_03 
##########

### step1: ggplot for each sample
#setwd("~/project_katja2/hmst_5hmC/")
#tss1 <- read.csv("./out_chrAll/plots/plotData/plotData_TSSgeneTES_smoothed_5mC_allMRs.csv", header=T) # original

tss1 <- read.csv("./dataIn_fig5/out_chrAll/plots/plotData/plotData_TSSgeneTES_smoothed_5mC_allMRs.csv", header=T)

library(tidyverse)

head(tss1)
as.tibble(tss1)
tss1 <- tss1[, -1]

# samples:
dim(tss1)
tss1 <- tss1[c(1:8001), ]

x <- c(1:8001)
y1 <- tss1$WT
y2 <- tss1$Neil1
y3 <- tss1$Neil2
y4 <- tss1$DK

# plot:
ggplot(tss1) +
  #geom_point() + 
  geom_smooth(aes(x, y1, colour = "WT"), method = "loess", span = 0.3, se = F) + 
  geom_smooth(aes(x, y2, colour = "Neil1"), method = "loess", span = 0.3, se = F) +
  geom_smooth(aes(x, y3, colour = "Neil2"), method = "loess", span = 0.3, se = F) +
  geom_smooth(aes(x, y4, colour = "DK"), method = "loess", span = 0.3, se = F) +
  labs( y = "Average methylation Level", x = " ", title = "5mC profile") +
  # Add a vertical line at x = 3
  geom_vline(xintercept = c(0, 1000, 2000, 6000, 7000, 8000), 
             linetype="dotted", 
             color = "black", size=0.5) +
  scale_color_manual("", breaks = c("WT", "Neil1", "Neil2", "DK"),
                     values = c("WT" = "gold2", "Neil1" = "red3", "Neil2" = "springgreen4","DK"="royalblue3" )) +
  #theme(axis.text.x = element_blank()) +
  scale_x_continuous(breaks = c(1, 1000, 2000, 4000, 6000, 7000, 8000), 
                     labels = c("-1Kb", "TSS", "+1Kb", "50%", "-1Kb", "TES", "+1Kb")) +
  theme(text = element_text(size=14))+
  theme_bw()

dev.copy(jpeg, file="ggplot_5mC_profile.jpg",
         units="in", 
         width=6, 
         height=4, 
         res=300)
dev.off()


### step2: T-test
WT <- y1
Neil1 <- y2
Neil2 <- y3
DK <- y4

data1_tgt <- tibble(WT= WT, Neil1 = Neil1, Neil2 = Neil2, DK=DK)
data2_tss <-data1_tgt[c(1:2000), ]
data3_gene <-data1_tgt[c(2001:6000), ]
data4_tes <-data1_tgt[c(6001:8001), ]

# boxplot for Figure 5C
boxplot(data1_tgt,
        main = "5mC in TSS_gene_TES",
        xlab = "",
        ylab = "Average Methylation Level",
        ylim = c(0.6, 2.3) ) 

boxplot(data2_tss,
        main = "5mC in TSS",
        xlab = "",
        ylab = "Average Methylation Level",
        ylim = c(0.5, 1.6) ) 
dev.copy(jpeg, file="boxplot_6M_5mC_TSS.jpg",
         units="in", 
         width=4, 
         height=4, 
         res=300)
dev.off()


boxplot(data3_gene,
        main = "5mC in gene_body",
        xlab = "",
        ylab = "Average Methylation Level",
        ylim = c(1.0, 2.2) )  
dev.copy(jpeg, file="boxplot_6M_5mC_geneBody.jpg",
         units="in", 
         width=4, 
         height=4, 
         res=300)
dev.off()


boxplot(data4_tes,
        main = "5mC in TES",
        xlab = "",
        ylab = "Average Methylation Level",
        ylim = c(1.0, 2.1) )  
dev.copy(jpeg, file="boxplot_6M_5mC_TES.jpg",
         units="in", 
         width=4, 
         height=4, 
         res=300)
dev.off()

## two sided test, all are sig.
t.test(data1_tgt$WT, data1_tgt$Neil1, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)
t.test(data1_tgt$WT, data1_tgt$Neil2, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)
t.test(data1_tgt$WT, data1_tgt$DK, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)


t.test(data2_tss$WT, data2_tss$Neil1, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)
t.test(data2_tss$WT, data2_tss$Neil2, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)
t.test(data2_tss$WT, data2_tss$DK, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99) 

t.test(data3_gene$WT, data3_gene$Neil1, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)
t.test(data3_gene$WT, data3_gene$Neil2, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)
t.test(data3_gene$WT, data3_gene$DK, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)

t.test(data4_tes$WT, data4_tes$Neil1, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)
t.test(data4_tes$WT, data4_tes$Neil2, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)
t.test(data4_tes$WT, data4_tes$DK, mu= 0, alt= "two.sided", paired = T, conf.level = 0.99)



