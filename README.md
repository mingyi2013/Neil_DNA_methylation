# Neil_DNA_methylation


This sesion contains the code information and data analysis processs for manuscript: Epigenetic landscape of 5mC and 5hmC in hippocampus of Neil1- and Neil2 deficient mice

# 1. Workflow for GWBS

This part contains the workflow for genome wide bisulfite sequence (GWBS) data analysis, codes for raw data quality control, sequence alignment and mC extraction with Bismark package.

## 1.1 Workflow for GWBS
Outline (**Figure 1**):
+ Quality control by FasteQC
+ Alignment by Bismark package
+ Extraction of 5mC in CG, CHG and CHH context by Bismark
+ Identification of Differential Methylation (**DM**):
    - CG in genome-wide by methyKit 
        - Differentially methylated cytosine (DMC) in base resolution
        - Differentially methylated region (DMR) by tiled window (1000 bp per window)
        - Differentially methylate promoter region (DMP) between 0 and - 1000 bp of TSS. 
    - CG, CHG and CHH in gene feature regions, by HMST-Seq-Analyzer
        - TSS +/- 1000 bp
        - TES +/- 1000 bp
        - gene body (between TSS and TES)
+ data visualization and Gene Ontology (GO) enrichment analysis
    - plot by ggplot2
    - GO enrichment analysis by clusterProfile

## 1.2 Code for quanlity control
```console
#!/usr/bin/sh
mkdir ./out_fastqc
fastqc -o ./out_fastqc *_1.fq.gz
fastqc -o ./out_fastqc *_2.fq.gz
```
## 1.3 Code for sequence alignment and 5mC extraction
Step 1. Sequence alignment
```console
#!/usr/bin/sh
mkdir ./out_align
bismark -p 4 -o ./out_align --gzip --genome /cluster/projects/nn9383k/mingyiy/bismark_saga/genomes/mouse/GRCm38_ncbi/ \
-1 *_1.fq.gz -2 *_2.fq.gz

Step 2. Deduplication
```console
#!/usr/bin/sh
deduplicate_bismark -p --bam ./out_align/*_1_bismark_bt2_pe.bam

mkdir out_deduplicate 
mv *deduplicated.bam ./out_deduplicate
mv *deduplication_report.txt ./out_deduplicate
```

Step 4. Extraction of 5mC in CG, CHG and CHH 
```console
#!/usr/bin/sh
bismark_methylation_extractor -o out_mExtractor_non_CpG --multicore 3 -p --ignore_r2 2 --comprehensive --cytosine_report --bedGraph --CX_context --ample_memory --split_by_chromosome --gzip --genome /cluster/projects/nn9383k/mingyiy/bismark_saga/genomes/mouse/GRCm38_ncbi/ ./out_deduplicate/*1_bismark_bt2_pe.deduplicated.bam
```

Step 5. Bismark report (**Table 1**)
```console
bismark2report --dir ./out_report --alignment_report ./out_align/*report.txt \
--dedup_report ./out_deduplicate/*report.txt --splitting_report ./out_mExtractor_non_CpG/*_splitting_report.txt \
--mbias_report ./out_mExtractor_non_CpG/*.M-bias.txt 
```

# 2. DM in genome-wide CG by methylKit tool
(**Figure 2A-B**)
## 2.1 Differentially methylated cytosine (DMC)
```console
Rscript script_fig1_methylkit_DMC.R
```
## 2.2 Differentially methylated region (DMR)
```console
Rscript script_fig1_methylkit_DMR.R
```
## 2.3 Differential methylation gene annotation

- Differentially methylated cytosine (DMC)
```console
Rscript script_fig1_methylkit_annotation_DMC.R
```
- Differentially methylated region (DMR)
```console
Rscript script_fig1_methylkit_annotation_DMR.R
```
- Differentially methylated promoter (DMP)
```console
Rscript script_fig1_methylkit_annotation_DMP.R
```

# 3. DM in gene feature region by HMST-Seq Analyzer
- DMRs process data generated by HMST-Seq-Analyzer
```console
bash script_fig1_hmst_CG.sh
bash script_fig1_hmst_CHG.sh
bash script_fig1_hmst_CHH.sh
```
- Extract, merge process data by R packages
```console
Rscript script_fig1_DMRs_1_table_extract.R 
Rscript script_fig1_DMRs_2_table_merge.R 
```
- DMRs gene count and plot (**Figure 3**)
```console
Rscript script_fig1_DMRs_3_table_format_filter.R
```
out data for next step: geneList_DMRs_*

# 4. GO enrichment in DM associated genes
## 4.1 Extract common DM genes among Nei1,2 and double knockouts
input data: dataOut_fig3/geneList_DMRs_*  
out data: commonGenes_C*.txt  
out plot: venn diagram (**figure 4A**)

```console
Rscript script_fig4a_venn.R
```

## 4.2 GO enrichment for GWBS data
-  GO enrichment in CG
-  GO enrichment in CHG
-  GO enrichment in CHH  

input data: dataIn_fig4/commonGenes_C*.txt (generated from step 4.1)  
out plot: **Figure 4B-E** in fold: GO_commonGenes_C* 

```console
Rscript script_fig4b_GO_enrich.R
```

# 5. HMST-sequencing data analysis
## 5.1 Data processing by HMST-Seq-Analyzer
Extraction and gene annotation of DMRs in 5mC and 5hmC (**Figure 5E-H**)
```console
bash script_fig5_hmst_1_DMR.sh
```

## 5.2 Down-stream analysis by R package
- Extract DMRs table
```console
Rscript script_fig5_hmst_1b_DMRs_table_extract.R  <file_in> <table_out> <filtered_talbe_out>
```
- Methylation profile and boxplot of 5mC and 5hmC (**Figure 5A-D**)  
input data for 5mC:  
        dataIn_fig5/out_chrAll/plots/plotData/plotData_TSSgeneTES_smoothed_5mC_allMRs.csv  
input data for 5hmC:  
    dataIn_fig5/out_chrAll/plots/plotData/plotData_TSSgeneTES_smoothed_5hmC_allMRs.csv
```console
Rscript script_fig5_hmst_2_plot_mC_profile.R
Rscript script_fig5_hmst_3_plot_5hmC_profile.R
```
- venn diagram plot for overlap DMRs between groups (**Figure 6A**) and common DMRs genes in 5mC_hypo and 5hmC_hyper (**Figure 6B**)  
input data for Figure 6: dataIn_fig6/
```console
Rscript script_fig6_hmst_1_overlap_5mC_hypo_5hmC_hyper.R
```

- GO enrichment (**Figure 6C-F**)
```console
Rscript script_fig6_hmst_2_seq_overlap_overpresent.R
```

- Boxplot of 5mC/5hmC ratio (**Figure 6G**)
```console
Rscript script_fig6_hmst_4_ratio_5mC_5hmC_sum_plot.R
```

# 6. RNA-seq analysis
## 6.1 QC and Sequence alignment
- Quality control (QC)
```console
mkdir ./out_fastqc
for i in N1_1 N1_2 N2_1 N2_2 N12_1 N12_2 WT_1 WT_2
do
  echo "QC for sample:"
  echo ${i}
  fastqc -o ./out_fastqc data_in/${i}_1.fq.gz    
  fastqc -o ./out_fastqc data_in/${i}_2.fq.gz 
done
````

- Sequence alignment
```console
#!/usr/bin/sh
# alignment by Hisat2
mkdir ./out_align
SAMPLES="N1_1 N1_2 N2_1 N2_2 N12_1 N12_2 WT_1 WT_2"
for file in $SAMPLES
do
	echo 'start alignment for sample:'
	echo $file
 	hisat2 -p 11 --dta -x /cluster/projects/nn9383k/mingyiy/RNAseq_Neil/index_grcm38/grcm38_tran/genome_tran \
	-1 data_in/${file}_1.fq.gz \
	-2 data_in/${file}_2.fq.gz \
	-S out_align/${file}.sam
done
echo 'Done!'
```
## 6.2 sam file sorting and bam file indexing
- sam file sorting
```console
#!/usr/bin/sh
# sort by Samtools
mkdir ./out_sort
SAMPLES="N1_1 N1_2 N2_1 N2_2 N12_1 N12_2 WT_1 WT_2"
for file in $SAMPLES
do
	echo $file
	samtools sort -@ 11 -o out_sort/${file}.bam out_align/${file}.sam
done
```
- bam file indexing
```console
#!/usr/bin/sh
# index Bam file by Samtools
SAMPLES="N1_1 N1_2 N2_1 N2_2 N12_1 N12_2 WT_1 WT_2"
for file in $SAMPLES
do
	echo run_bam: ${file}
	samtools index out_sort/${file}.bam out_sort/${file}.bam.bai
done
```

## 6.3 Assemble, merge and count transcripts by stringTie
- Assemble transcripts
Reference mm10.ncbiRefSeq.gtf download from: http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/  
It contains gene_name and transcript id as NCBI NM*.  The file was formated by removing chr in column 1 (unix: sed -E 's/^chr//g').  
```console
#!/usr/bin/sh
# assemble transcripts by stringTie
mkdir ./out_stringTie
SAMPLES="N1_1 N1_2 N2_1 N2_2 N12_1 N12_2 WT_1 WT_2"

for file in $SAMPLES; do
echo run_stringTie: $file
stringtie -p 11 -G /cluster/projects/nn9383k/mingyiy/RNAseq_Neil/index_grcm38/grcm38_tran/mm10.ncbiRefSeq_format.gtf \
-o out_stringTie/${file}.gtf \
-l ${file} out_sort/${file}.bam
done
```

- Merge transcripts
```console
#!/usr/bin/sh
# merge transcripts by stringTie
# out directory: ./out_stringTie

# create a mergelist.txt file
echo 'create a mergelist file:'
ls out_stringTie/*.gtf | sort > mergelist.txt
cat mergelist.txt

echo 'merge transcripts:'
stringtie --merge -p 11 -G /cluster/projects/nn9383k/mingyiy/RNAseq_Neil/index_grcm38/grcm38_tran/mm10.ncbiRefSeq_format.gtf \
-o out_stringTie/stringtie_merged.gtf mergelist.txt

echo 'DONE!'
```

- Count transcripts
```console
#!/usr/bin/sh
# coverage of transcripts by stringTie -eB
# set output fold: ballgown

SAMPLES="N1_1 N1_2 N2_1 N2_2 N12_1 N12_2 WT_1 WT_2"

for SAMPLE in $SAMPLES; do
    echo $SAMPLE
    stringtie -e -B -p 11 -G out_stringTie/stringtie_merged.gtf -o ballgown/${SAMPLE}/${SAMPLE}.gtf out_sort/${SAMPLE}.bam
done

echo 'Done!'
```

## 6.4 DEGs identification by DEseq2
Data input: the count table file of "t_data.ctab" in fold ballgown/*/, generated by stringTie  
- PCA plot (**Figure 7A**)
- heatmap of DE genes (**Figure S7**)
```console
Rscript script_fig7_1_DESeq.R
```

- barplot of DE gene cout (**Figure 7B**)
Data input: dataIn_fig7/sum_transcript_diffExp.txt  

```console
Rscript script_fig7_2_barplot_DEGs.R
```

- venn diagarm of overlap DE genes among group comparisons (**Figure 7C**)  
Data input: dataIn_fig7/DE_groups/*/DE_geneList_*.txt

```console
Rscript script_fig7_3_venn.R
```

## 6.5 GO enrichment by clusterProfile
- GO enrichment in DEGs in Neil1 vs WT (**Figure 7D**)  
Data input: dataIn_fig7/DE_groups/*/DE_data_sigDiff_condition*.txt
```console
Rscript script_fig7_4_GO_enrich.R
```
# 7. Correlation between differential methylaiton (DM) and differential expression (DE)
- Overlap genes sorted by catalog

```console
Rscript 
```
- Overlap genes sorted by features
```console
Rscript 
```

- venn diagarm for overlap genes
```console
Rscript 
```

