#!/bin/bash
#SBATCH --job-name=eqtl_cond
#SBATCH --output=logs/eqtl_cond_%A_%a.out
#SBATCH --error=logs/eqtl_cond_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-24

module load qtltools
PROJECT="PLEXUS_colon_n1149" # Change this
K="28" # Change this 

CHUNK=${SLURM_ARRAY_TASK_ID}
TOTAL_CHUNKS=24


OUTDIR="${PROJECT}/RUVr/cond"
mkdir -p ${OUTDIR}
mkdir -p logs

QTLtools cis \
  --vcf ${PROJECT}/input/vcf/${PROJECT}_R2_0.3_rsidsnp_maf001.vcf.gz \
  --bed ${PROJECT}/input/bed/${PROJECT}_eqtl_qtltools.bed.gz \
  --cov ${PROJECT}/input/cov/qtltools/cov_${K}.txt \
  --mapping ${PROJECT}/RUVr/permutations/${PROJECT}.eqtl.permutations_cov${K}.fdr0.05.thresholds.txt \
  --normal \
  --std-err \
  --chunk ${CHUNK} ${TOTAL_CHUNKS} \
  --out ${OUTDIR}/cond_cov${K}_chunk${CHUNK}.txt