# Neil_DNA_methylation


This sesion contains the code information and data analysis processs for manuscript: Epigenetic landscape of 5mC and 5hmC in hippocampus of Neil1- and Neil2 deficient mice

# 1. Workflow for GWBS

This part contains the workflow for genome wide bisulfite sequence (GWBS) data analysis, codes for raw data quality control, sequence alignment and mC extraction with Bismark package.

## 1.1 Workflow for GWBS
**Outline:**
+ Quality control by FasteQC
+ Alignment by Bismark package
+ Extraction of 5mC in CG, CHG and CHH context by Bismark
+ Identification of differential, by methylation (DM):
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

```


# 2. DM in genome-wide CG by methylKit tool

## 2.1 DMC
## 2.2 DMR
## 2.3 DMP

# 3. DM in gene feature region by HMST-Seq Analyzer
## 3.1 DMRs in CG

## 3.2 DMRs in CHG
## 3.3 DMRs in CHH

# 4. GO enrichment in DM associated genes
## 4.1 GO enrichment in CG
## 4.1 GO enrichment in CHG
## 4.1 GO enrichment in CHH

# 5. Workflow for HMST-sequencing

