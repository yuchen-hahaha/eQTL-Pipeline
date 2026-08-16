#!/bin/bash
#SBATCH --job-name=eqtl_permutation
#SBATCH --output=logs/eqtl_perm_%A_%a.out
#SBATCH --error=logs/eqtl_perm_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-30

module load qtltools
PROJECT="PLEXUS_colon_n1149" # Change this
BEST_K="28" # Change this 

CHUNK_ID=${SLURM_ARRAY_TASK_ID}
TOTAL_CHUNKS=23



QTLtools cis \
  --vcf ${PROJECT}/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz \
  --bed ${PROJECT}/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz \
  --cov ${PROJECT}/input/cov/cov_${BEST_K}.txt \
  --seed 12345 \
  --permute 1000 \
  --normal \
  --std-err \
  --chunk ${CHUNK_ID} ${TOTAL_CHUNKS} \
  --out ${PROJECT}/RUVr/permutations/permutations_cov${BEST_K}_chunk${CHUNK_ID}.txt