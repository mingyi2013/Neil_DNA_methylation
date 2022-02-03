SAGA /cluster/projects/nn9383k/mingyiy/hmst_katja2]$ cat script_merge_table.R 
#!/bin/sh/Rscript

###############################
# New Part II: merge tables of output in hmst, by R
# mingyi yang,  2021_07_06, correct for CHX input file pattern in step1
# script_merge_table.R
# /cluster/projects/nn9383k/mingyiy/hmst_katja1_CHX/
# sbatch under: /cluster/projects/nn9383k/mingyiy/hmst_katja2/Neil1_vs_WT_10d/

# test: module load R/4.0.0-foss-2020a
# Rscript ../script_merge_table.R

# batch job: sbatch ../job_sbatch_merge_table.sm
# skip data in chrY for CHX. skip chrY for KO2_vs_WT2 in CG
###############################

###########
# step 1. For KO_vs_WT: rbind of TSS, TES and genes in each chr:
###########
dir.create("out_table") 
library(tidyverse)

r_threshold = 0.1 # set value for filter rratio

### merged chr for CG:
context = "CG"
print(context)
## KO_1_vs_WT_1
for (i in c(1:19,"X"))
  { print(i)
  Path= paste ("./out_", context, "_chr",i,"/data", sep = "")
  merged_df <- dir(path=Path,
                   pattern = "table_KO1.me_vs_WT1.me_5mC_", full.names = T) %>%
                map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE )
  outFile = paste("out_table/merged_table_KO1_vs_WT1_", context, "_chr_", i, ".txt", sep = "")
  print(outFile)
  write.table(merged_df, file = outFile, row.names = F, quote = F, sep = "\t")
}

## KO_1_vs_WT_2
for (i in c(1:19,"X"))
{print(i)
  Path= paste ("./out_", context, "_chr",i,"/data", sep = "")
  merged_df <- dir(path=Path,
                   pattern = "table_KO1.me_vs_WT2.me_5mC_", full.names = T) %>%
    map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE )
  outFile = paste("out_table/merged_table_KO1_vs_WT2_", context, "_chr_", i, ".txt", sep = "")
  print(outFile)
  write.table(merged_df, file = outFile, row.names = F, quote = F, sep = "\t")
}

## KO_2_vs_WT_1
for (i in c(1:19,"X"))
{print(i)
  Path= paste ("./out_", context, "_chr",i,"/data", sep = "")
  merged_df <- dir(path=Path,
                   pattern = "table_KO2.me_vs_WT1.me_5mC_", full.names = T) %>%
    map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE )
  outFile = paste("out_table/merged_table_KO2_vs_WT1_", context, "_chr_", i, ".txt", sep = "")
  print(outFile)
  write.table(merged_df, file = outFile, row.names = F, quote = F, sep = "\t")
}

## KO_2_vs_WT_2
for (i in c(1:19,"X"))
{print(i)
  Path= paste ("./out_", context, "_chr",i,"/data", sep = "")
  merged_df <- dir(path=Path,
                   pattern = "table_KO2.me_vs_WT2.me_5mC_", full.names = T) %>%
    map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE )
  outFile = paste("out_table/merged_table_KO2_vs_WT2_", context, "_chr_", i, ".txt", sep = "")
  print(outFile)
  write.table(merged_df, file = outFile, row.names = F, quote = F, sep = "\t")
 }


### merged chr for CHX: data in chrX chrY are empty, skip.
# check file name pattern, it is diff between CG and CHX: in CHX: out_CHG_chr1/data/table_KO1.me_vs_WT1*.csv
# for loop for context: until end of step1.
for (context in c("CHG", "CHH")) {
print(context) 
## KO_1_vs_WT_1
for (i in c(1:19)) 
  { print(i)
  Path= paste ("./out_", context, "_chr",i,"/data", sep = "")
  merged_df <- dir(path=Path, 
                   pattern = "table_KO1.me_vs_WT1", full.names = T) %>% 
                map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE )
  outFile = paste("out_table/merged_table_KO1_vs_WT1_", context, "_chr_", i, ".txt", sep = "")
  print(outFile)
  write.table(merged_df, file = outFile, row.names = F, quote = F, sep = "\t")
}

## KO_1_vs_WT_2
for (i in c(1:19)) 
{print(i)
  Path= paste ("./out_", context, "_chr",i,"/data", sep = "")
  merged_df <- dir(path=Path, 
                   pattern = "table_KO1.me_vs_WT2", full.names = T) %>% 
    map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE )
  outFile = paste("out_table/merged_table_KO1_vs_WT2_", context, "_chr_", i, ".txt", sep = "")
  print(outFile)
  write.table(merged_df, file = outFile, row.names = F, quote = F, sep = "\t")
}

## KO_2_vs_WT_1
for (i in c(1:19)) 
{print(i)
  Path= paste ("./out_", context, "_chr",i,"/data", sep = "")
  merged_df <- dir(path=Path, 
                   pattern = "table_KO2.me_vs_WT1", full.names = T) %>% 
    map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE )
  outFile = paste("out_table/merged_table_KO2_vs_WT1_", context, "_chr_", i, ".txt", sep = "")
  print(outFile)
  write.table(merged_df, file = outFile, row.names = F, quote = F, sep = "\t")
}

## KO_2_vs_WT_2
for (i in c(1:19)) 
{print(i)
  Path= paste ("./out_", context, "_chr",i,"/data", sep = "")
  merged_df <- dir(path=Path, 
                   pattern = "table_KO2.me_vs_WT2", full.names = T) %>% 
    map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE )
  outFile = paste("out_table/merged_table_KO2_vs_WT2_", context, "_chr_", i, ".txt", sep = "")
  print(outFile)
  write.table(merged_df, file = outFile, row.names = F, quote = F, sep = "\t")
 }
} # end of for_loop in step1 for CHX merge.

############# 
# step 2. Merged all chr* for KO_vs_WT
############
# for_loop for context: until end in script.
for (context in c("CG","CHG", "CHH")) {
print(context)

## for KO1_vs_WT1 all_chr*
Path= "out_table"
Pattern = paste("merged_table_KO1_vs_WT1_", context, sep = "")
merged_df1 <- dir(path=Path, 
                 pattern = Pattern, full.names = T) %>% 
  map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE ) %>% 
  select(c(1:7)) %>% 
  rename(median_WT_1 = median_WT, median_KO_1 = median_KO,rratio_1 = rratio, pvals_1= pvals ) %>% 
  arrange(location) %>% 
  distinct(location, .keep_all = TRUE)
head(merged_df1, n=3)
outFile = paste("out_table/allChr_table_KO1_vs_WT1_", context, ".txt", sep = "")
print(outFile)
write.table(merged_df1, file = outFile, row.names = F, quote = F, sep = "\t")

## for KO1_vs_WT2 all_chr*
Path= "out_table"
Pattern = paste("merged_table_KO1_vs_WT2_", context, sep = "")
merged_df2 <- dir(path=Path, 
                 pattern = Pattern, full.names = T) %>% 
  map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE ) %>% 
  select(c(1:7)) %>% 
  rename(median_WT_2 = median_WT, median_KO_1 = median_KO,rratio_2 = rratio, pvals_2= pvals ) %>% 
  arrange(location) %>% 
  distinct(location, .keep_all = TRUE)
outFile = paste("out_table/allChr_table_KO1_vs_WT2_", context, ".txt", sep = "")
print(outFile)
write.table(merged_df2, file = outFile, row.names = F, quote = F, sep = "\t")

## for KO2_vs_WT1 all_chr*
Path= "out_table"
Pattern = paste("merged_table_KO2_vs_WT1_", context, sep = "")
merged_df3 <- dir(path=Path, 
      pattern = Pattern, full.names = T) %>% 
      map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE ) %>% 
      select(c(1:7)) %>% 
      rename(median_WT_1 = median_WT, median_KO_2 = median_KO,rratio_3 = rratio, pvals_3= pvals ) %>% 
      arrange(location) %>% 
      distinct(location, .keep_all = TRUE)
outFile = paste("out_table/allChr_table_KO2_vs_WT1_", context, ".txt", sep = "")
print(outFile)
write.table(merged_df3, file = outFile, row.names = F, quote = F, sep = "\t")

## for KO2_vs_WT2 all_chr*
Path= "out_table"
Pattern = paste("merged_table_KO2_vs_WT2_", context, sep = "")
merged_df4 <- dir(path=Path, 
        pattern = Pattern, full.names = T) %>% 
        map_df(read_delim, "\t", escape_double = FALSE, trim_ws = TRUE ) %>% 
        select(c(1:7)) %>% 
        rename(median_WT_2 = median_WT, median_KO_2 = median_KO,rratio_4 = rratio, pvals_4= pvals ) %>% 
        arrange(location) %>% 
        distinct(location, .keep_all = TRUE)
outFile = paste("out_table/allChr_table_KO2_vs_WT2_", context, ".txt", sep = "")
print(outFile)
write.table(merged_df4, file = outFile, row.names = F, quote = F, sep = "\t")


############# 
# step 3.  merged 4 pair-wised comparisons
############
total_df <- merged_df1 %>% 
        inner_join(merged_df2, by=c("location", "feature", "gene")) %>% 
        inner_join(merged_df3, by=c("location", "feature", "gene")) %>%
        inner_join(merged_df4, by=c("location", "feature", "gene"))
# limited the number of dicimals in defined columns:
col_list= c(4:6, 8:10, 12:14, 16:18)
total_df[, col_list] <- round(total_df[, col_list], digits = 2) 
# filter for keeping only co-up or co-down rows:
total_df2 <- total_df %>% 
      filter((rratio_1>0 & rratio_2>0 & rratio_3>0 & rratio_4 >0) | (rratio_1<0 & rratio_2<0 & rratio_3<0 & rratio_4<0))

head(total_df2, n=10)
# write and save df: 
outFile = paste("out_table/allChr_table_group_KO_vs_WT_", context, ".txt", sep = "")
write.table(total_df2, file = outFile, row.names = F, quote = F, sep = "\t")
}   #For loop for context, from start in step2 to end in script.

print("Done_merge_table!")
##############################################