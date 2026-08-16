#!/bin/bash
#SBATCH --job-name=eqtl_nom
#SBATCH --output=logs/eqtl_nom_%A_%a.out
#SBATCH --error=logs/eqtl_nom_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1


module load qtltools
PROJECT="PLEXUS_colon_n1149" #Change this
K="28" #Change this 


mkdir -p logs

  
QTLtools cis \
  --vcf ${PROJECT}/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz \
  --bed ${PROJECT}/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz \
  --cov ${PROJECT}/input/cov/cov_${BEST_K}.txt \
  --seed 12345 \
  --normal \
  --nominal 1 \
  --std-err \
  --out ${PROJECT}/RUVr/nominals/nominals_cov${BEST_K}.txt
  
  