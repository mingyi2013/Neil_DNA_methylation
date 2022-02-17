#!/bin/bash
echo starting HMST-Seq-Analyzer

TEST_METHOD='Ttest'

echo ${TEST_METHOD}
# Do the work
hmst_seq_analyzer gene_annotation -l 2000 -xL 1000000 -M 10000 -N 100000 \
-F  out \
-hu no -n no \
-r in_data/refFlat.txt \
-g in_data/mm10.chrom.sizes.clear.sorted
echo gene_annotation-DONE

hmst_seq_analyzer data_preprocessing \
-F out \
-z no -m no -hu no -n no \
-fko data_input/Neil1.txt \
-fko data_input/Neil2.txt \
-fko data_input/DK.txt \
-fwt data_input/WT.txt \
-g in_data/mm10.chrom.sizes.clear.sorted
echo data_preprocessing-DONE

#
hmst_seq_analyzer find_MRs -W no -a 2000 -mc1 3 -mc2 5 -mc3 3 \
-F out -p 5 \
-fko out/list_mC_hmC_files_KO.txt \
-fwt out/list_mC_hmC_files_WT.txt \
-ref out/data/refFlat_clean_sorted.bed \
-reg out/list_region_files.txt
echo find_MRs-DONE

#
hmst_seq_analyzer prepare_for_DMR_finding \
-isST 1 -a 2000 -mc1 3 -mc2 5 -mc3 3 \
-F out \
-p 8 \
-ko out/list_all_filtered_formatted_MRs_KO.txt \
-wt out/list_all_filtered_formatted_MRs_WT.txt
echo prepare_for_DMR_finding-DONE

hmst_seq_analyzer DMR_search \
-T ${TEST_METHOD} \
-F out \
-p 8 \
-f out/list_prepared_for_DMR_finding_imputed.txt
echo DMR_search-DONE

hmst_seq_analyzer prep4plot \
-F out \
-ko out/list_all_filtered_formatted_MRs_KO.txt \
-wt out/list_all_filtered_formatted_MRs_WT.txt
echo prep4plot-DONE

hmst_seq_analyzer plot_all \
-F out \
-reg out/list_count_allMRs_regions_files.txt \
-sit out/list_count_sites_files.txt \
-cmc out/counts_DMR_hypo_hyper_imputed_${TEST_METHOD}_5mC.csv \
-chmc out/counts_DMR_hypo_hyper_imputed_${TEST_METHOD}_5hmC.csv \
-aMR out/list_TSS_genebody_TES_enhancer_allMRs.txt 
echo plot_all-DONE