# run location:
# SAGA /cluster/projects/nn9383k/mingyiy/hmst_katja2/hmst_Neil1_vs_WT_6M 
#-----------------------
#!/bin/bash
# HMST GWBS_mouse v4b

# original name:script2a_CG_hmst_m_v4.sh, update 2020_12_28
# New setting: hmst_seq_analyzer find_MRs -W yes (try to use equal bin to remove low density mC region) -a 200.


# define string
string="run HMST GWBS_mouse..."
# print variable on screen
echo $string

set -x
which python

Type='CG'
#CHR='chr9'
CHR="$1"

echo $Type
echo $CHR


hmst_seq_analyzer gene_annotation -F out_${Type}_${CHR} -hu no \
-r ../in_data/refFlat.txt \
-n yes -X 1000 -Y 1000 -M 10000 -N 50000 -l 2000 -xL 5000 \
-g ../in_data/mm10_chrom_length_sorted
echo gene_annotation-DONE

# if inpute data in zip: -z no
hmst_seq_analyzer data_preprocessing -F out_${Type}_${CHR} -z no -m yes -hu no -n yes -c 0.4999 \
-fko ../mC_count_Neil1_3/me_${CHR}_${Type}.bed \
-fko ../mC_count_Neil1_4/me_${CHR}_${Type}.bed \
-fwt ../mC_count_WT_3/me_${CHR}_${Type}.bed \
-fwt ../mC_count_WT_4/me_${CHR}_${Type}.bed \
-g ../in_data/mm10_chrom_length_sorted
echo data_preprocessing-DONE

# GWBS: -W yes, adjust -a as new default.
hmst_seq_analyzer find_MRs \
-W yes -a 200 -mc1 3 -mc2 5 -mc3 3 \
-F out_${Type}_${CHR} -p 5 \
-fko out_${Type}_${CHR}/list_mC_hmC_files_KO.txt \
-fwt out_${Type}_${CHR}/list_mC_hmC_files_WT.txt \
-ref out_${Type}_${CHR}/data/refFlat_clean_sorted.bed \
-reg out_${Type}_${CHR}/list_region_files.txt \
-e ../in_data/all_tissueSpecific_enhancers4mm10_sorted2.bed
echo find_MRs-DONE

# set if it is same trend yes: -isST 1; -mc1 10 -mc2 10 -mc3 10; -a 2000(default); -i zeros(default)
hmst_seq_analyzer prepare_for_DMR_finding \
-isST 1 -a 2000 -mc1 10 -mc2 10 -mc3 10 -i zeros \
-F out_${Type}_${CHR} -p 4 -i zeros \
-ko out_${Type}_${CHR}/list_all_filtered_formatted_MRs_KO.txt \
-wt out_${Type}_${CHR}/list_all_filtered_formatted_MRs_WT.txt
echo prepare_for_DMR_finding-DONE

# set FDR=0.05, method = Ttest
hmst_seq_analyzer DMR_search \
-T Ttest \
-F out_${Type}_${CHR} -p 5 \
-pl no -fdr 0.05 \
-f out_${Type}_${CHR}/list_prepared_for_DMR_finding_imputed.txt
echo DMR_search-DONE

hmst_seq_analyzer prep4plot \
-F out_${Type}_${CHR} -c 0.4999 \
-ko out_${Type}_${CHR}/list_all_filtered_formatted_MRs_KO.txt \
-wt out_${Type}_${CHR}/list_all_filtered_formatted_MRs_WT.txt
echo prep4plot-DONE

hmst_seq_analyzer plot_all \
-F out_${Type}_${CHR} \
-reg out_${Type}_${CHR}/list_count_allMRs_regions_files.txt \
-sit out_${Type}_${CHR}/list_count_sites_files.txt \
-cmc out_${Type}_${CHR}/counts_DMR_hypo_hyper_imputed_*_5mC.csv \
-chmc out_${Type}_${CHR}/counts_DMR_hypo_hyper_imputed_*_5hmC.csv \
-aMR out_${Type}_${CHR}/list_TSS_genebody_TES_enhancer_allMRs.txt
#-oMR out_${Type}_${CHR}/list_overlapping_MRs.txt #used only for 2 samples, one WT and one KO.
echo plot_all-DONE

#---------------------------------------------------
