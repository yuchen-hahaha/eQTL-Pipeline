

1. Update PLEXUS_SPARC_ALLRNASEQ_METAFILE_n1322_TIN_20260520.txt 
$Rscript add_TIN_values.R PLEXUS_SPARC_ALLRNASEQ_METAFILE_n1322.txt PLEXUS_SPARC_ALLRNASEQ_METAFILE_n1322_TIN_20260520.txt 

2. Copy raw ileum meta data and save to a txt file RNAseq_ileum_metadata_raw.txt (https://adminliveunc-my.sharepoint.com/:x:/r/personal/away_ad_unc_edu/_layouts/15/doc2.aspx?sourcedoc=%7B54506B86-AC7E-4F3F-A715-B3F87C0834C7%7D&file=ALL_ILEUM_COLDATA.xlsx&action=default&mobileredirect=true&DefaultItemOpen=1)

3. Generate/Execute R script to update RNAseq_ileum_metadata_raw.txt into UNC_ileum_n156_meta.txt 

$Rscript Table.add_TIN_to_UNC_ileum_meta.R /work/users/y/u/yuchenh/sQTL_PLEXUS_UNC/00.input/metafiles/RNAseq_ileum_metadata_raw.txt /work/users/y/u/yuchenh/sQTL_PLEXUS_UNC/00.input/metafiles/RNAseq_ileum_metadata_20260521.txt

4. Generate temporary UNC ileum Cutoff 0.8 file 
$ awk -F'/' '{print $NF}' UNC_ileum_n156_bampaths.txt | awk -F'.' '{print $1}' > /proj/fureylab/data/RNA-seq/human/verifyBamID/202605_PLEXUS_R2_0.3_UNC_R2_0.3_Ileum/UNC_ILEUM_hg38v49_R2_0.3_RNA_seq_Cutoff_0.8_Pass_SampleID_tmp.txt  

5. Generate temporary PLEXUS ileum Cutoff 0.8 file 
Download the UNC_ILEUM_hg38v49r20.3_best_results_uniq.txt extract values in the BAM_ID column where MATCH_RESULTS column said yes.

Rscript /work/users/y/u/yuchenh/sQTL_PLEXUS_UNC/scripts/00.input/Table.generate_PLEXUS_UNC_merged_ileum_meta.R \
  /proj/fureylab/projects/PLEXUS_metafiles_preprocessing/PLEXUS_SPARC_ALLRNASEQ_METAFILE_n1322_TIN_20260520.txt \
  /work/users/y/u/yuchenh/sQTL_PLEXUS_UNC/00.input/metafiles/RNAseq_ileum_metadata_20260521.txt \
  /proj/fureylab/data/Genotypes/human/PLEXUS/vcf/ancestry_PCA_PLEXUS_UNC_R2_0.3_shared/plink_pca.eigenvec \
  /proj/fureylab/data/RNA-seq/human/verifyBamID/202605_PLEXUS_R2_0.3_UNC_R2_0.3_Ileum/PLEXUS_ILEUM_hg38v49_R2_0.3_RNA_seq_Cutoff_0.8_Pass_SampleID.txt \
  /proj/fureylab/data/RNA-seq/human/verifyBamID/202605_PLEXUS_R2_0.3_UNC_R2_0.3_Ileum/UNC_ILEUM_hg38v49_R2_0.3_RNA_seq_Cutoff_0.8_Pass_SampleID_tmp.txt \
  output_test3 \
  output_test3 \
  827 \
  uninflamed 
  
