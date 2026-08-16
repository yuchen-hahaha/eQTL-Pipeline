#!/bin/bash
#SBATCH --job-name=gene
#SBATCH --time=12:00:00
#SBATCH --mem=5G
#SBATCH -n 1
#SBATCH -N 1


while read bamfile; do
    sample=$(basename "$bamfile" | cut -d. -f1)
    QTLtools quan \
        --bam "$bamfile" \
        --gtf /proj/fureylab/genomes/human/hg38_gencodeV47_reference/RNA_annotation_v47/gencode.v47.annotation.gtf \
        --sample "$sample" \
        --out-prefix UNC_colon_n326/quant/"$sample" \
        --filter-mapping-quality 150 \
        --filter-mismatch 5 \
        --filter-mismatch-total 5 \
        --rpkm --no-merge
#done < input/UNC_n255_colon_bamfiles.txt
done < input/UNC_colon_nonIBD_bampaths4.txt

 