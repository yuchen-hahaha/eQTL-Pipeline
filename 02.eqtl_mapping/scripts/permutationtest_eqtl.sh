#!/bin/bash
#SBATCH --job-name=eqtl_permutationtest
#SBATCH --output=logs/eqtl_permtest_%A_%a.out
#SBATCH --error=logs/eqtl_permtest_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-30

module load qtltools

COV=${SLURM_ARRAY_TASK_ID}
PROJECT="PLEXUS_colon_n1149"


mkdir -p logs

QTLtools cis \
  --vcf ${PROJECT}/input/vcf/${PROJECT}_R2_0.3_rsidsnp_maf001.vcf.gz \
  --bed ${PROJECT}/input/bed/${PROJECT}_eqtl_qtltools.bed.gz \
  --cov ${PROJECT}/input/cov/cov_${COV}.txt \
  --seed 12345 \
  --permute 100 \
  --normal \
  --std-err \
  --out  ${PROJECT}/RUVr/permutationtest/perm100/permute_cov${COV}.txt