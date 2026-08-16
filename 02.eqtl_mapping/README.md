# eQTL QTLtools Pipeline
Update on 2026/06/03 
1. Change the n_chunks for conditional analysis to 24 
2. Add header suggestion for the output files. 

## Purpose

This document summarizes the **eQTL-only** workflow.

The workflow has four main stages:

1.  **Covariate optimization**: test different numbers of covariates and choose the best `K`.
2.  **Final eQTL mapping**: run nominal associations and permutation testing using the selected `K`.
3.  **Post-processing and statistics**: merge chunks, run FDR, generate significant pairs, and run conditional analysis.
4.  **Organization**: copy the final outputs into shared summary-statistics folders.

------------------------------------------------------------------------

## Required software

Load or install the following tools before running the workflow:

``` bash
module load qtltools
module load R
```

Required command-line tools:

-   `QTLtools`
-   `Rscript`
-   `bgzip`
-   `tabix`
-   `awk`
-   `cat`
-   `mkdir`
-   `cp`

------------------------------------------------------------------------

## Required scripts

The following scripts are used by the workflow.

| Script                          | Purpose                                                                                                      |
|------------------|------------------------------------------------------|
| `qtltools_runFDR_cis.R`         | Takes a QTLtools permutation output file and produces significant genes plus FDR thresholds.                 |
| `plot_RUVK_SigGene_eqtl.R`      | Plots number of significant eGenes across covariate values and writes the selected best `K`.                 |
| `rename_nominals_singlefile.sh` | Splits or renames the large nominal output into chromosome-level files.                                      |
| `generate_sigpairs.sh`          | Uses chromosome-level nominal results and the FDR threshold file to generate significant variant-gene pairs. |
| `cond_eqtl.sh`                  | Example to run conditional analysis                                                                          |

Recommended shared script folder:

``` bash
scripts
```

------------------------------------------------------------------------

## Required input files

### 1. Genotype VCF

Input:

``` bash
PLEXUS_colon_n1149/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz
```

Expected format:

-   bgzipped VCF
-   indexed with `.tbi`
-   I set SNP ID as rs151190501|chr1:727233:G:A using command

``` bash
input_shared_vcf_gz="/proj/fureylab/data/Genotypes/human/PLEXUS/vcf/PLEXUS_UNC_colon_R2_0.3_shared/PLEXUS_R2_0.3_MAF_0.01_UNC_R2_0.3_shared.MAF0.01.vcf.gz"
final_vcf_gz="PLEXUS_colon_n1149/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz"
THREADs="16"

bcftools annotate \
    --set-id '%ID|%CHROM:%POS:%REF:%ALT' \
    --threads ${THREADS} \
    -Oz \
    -o "${final_vcf_gz}" \
    "${input_shared_vcf_gz}"
```
-   samples should match the expression/covariate files
-   Example: /proj/fureylab/data/Genotypes/human/PLEXUS/vcf/PLEXUS_UNC_colon_R2_0.3_shared/PLEXUS_UNC_R2_0.3_colon_n1498_rsidsnp_maf001.vcf.gz

### 2. eQTL BED file

Input:

``` bash
PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed
```

Expected format:

-   QTLtools-compatible BED
-   one row per gene
-   gene expression values across samples
-   sample order should match the covariate file and VCF sample order

This file will be compressed and indexed before QTLtools is run.

### 3. Covariate files

Input directory:

``` bash
PLEXUS_colon_n1149/input/cov
```

Expected files:

``` bash
cov_0.txt
cov_1.txt
...
cov_30.txt
```

The covariate optimization step tests these files and chooses the best covariate number `K`.

------------------------------------------------------------------------

## Key configuration values

These are the main values used in the current eQTL analysis.

``` yaml
project_name: "PLEXUS_colon_n1149"
qtl: "eqtl"
RUV: "RUVr"

n_perm_test: 100
n_permutations: 1000
fdr_threshold: 0.05

n_chunks: 24
n_chunks_test: 24

cov_start: 0
cov_indices: 30
```

------------------------------------------------------------------------

# Stage 0. Prepare the BED file

## Goal

Compress and index the eQTL BED file for QTLtools.

## Input

``` bash
PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed
```

## Command

``` bash
bgzip -f PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed
tabix -p bed PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz
```

## Output

``` bash
PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz
PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz.tbi
```

## Make directories

``` bash
mkdir -p PLEXUS_colon_n1149/RUVr/permutationtest/perm100
mkdir -p PLEXUS_colon_n1149/RUVr/permutations
mkdir -p PLEXUS_colon_n1149/RUVr/nominals
mkdir -p PLEXUS_colon_n1149/RUVr/plots
mkdir -p PLEXUS_colon_n1149/RUVr/sgenes_sigpairs
```

------------------------------------------------------------------------

# Stage 1. Covariate optimization

## Goal

Test a range of covariate files and select the best covariate number `K`.

In this example, the tested covariates are:

``` bash
cov_0.txt
cov_1.txt
...
cov_30.txt
```

Each covariate value is tested using a small number of permutations, here `100`.

------------------------------------------------------------------------

## Step 1.1 Run QTLtools permutation test for each covariate value

## Inputs

For each covariate value `K`:

``` bash
PLEXUS_colon_n1149/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz
PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz
PLEXUS_colon_n1149/input/cov/cov_K.txt
```

## Command template

``` bash
module load qtltools
 
QTLtools cis \
  --vcf PLEXUS_colon_n1149/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz \
  --bed PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz \
  --cov PLEXUS_colon_n1149/input/cov/qtltools/cov_K.txt \
  --seed 12345 \
  --permute 100 \
  --normal \
  --std-err \
  --out PLEXUS_colon_n1149/RUVr/permutationtest/perm100/permute_covK.txt
```

#### Send sbatch directly

``` bash
sbatch scripts/permutationtest_eqtl.sh
```

## Output

``` bash
PLEXUS_colon_n1149/RUVr/permutationtest/perm100/permute_covK.txt
```

------------------------------------------------------------------------

## Step 1.2 Run FDR for each covariate value

## Input

``` bash
PLEXUS_colon_n1149/RUVr/permutationtest/perm100/permute_covK.txt
```

## Command

``` bash
Rscript scripts/qtltools_runFDR_cis.R \
  PLEXUS_colon_n1149/RUVr/permutationtest/perm100/permute_covK.txt \
  0.05 \
  PLEXUS_colon_n1149/RUVr/permutationtest/perm100/PLEXUS_colon_n1149.eqtl.covK.fdr0.05
```

## Outputs

``` bash
PLEXUS_colon_n1149/RUVr/permutationtest/perm100/PLEXUS_colon_n1149.eqtl.covK.fdr0.05.significant.txt
PLEXUS_colon_n1149/RUVr/permutationtest/perm100/PLEXUS_colon_n1149.eqtl.covK.fdr0.05.thresholds.txt
```

------------------------------------------------------------------------

## Step 1.3 Plot number of significant eGenes and choose best K

## Inputs

All files matching:

``` bash
PLEXUS_colon_n1149/RUVr/permutationtest/perm100/PLEXUS_colon_n1149.eqtl.cov*.fdr0.05.significant.txt
```

## Command

``` bash
Rscript scripts/plot_RUVK_SigGene_eqtl.R \
  PLEXUS_colon_n1149/RUVr/permutationtest/perm100/PLEXUS_colon_n1149.eqtl.cov \
  PLEXUS_colon_n1149/RUVr/plots \
  PLEXUS_colon_n1149_eqtl_sigGene_perm100 \
  30 \
  0.05
```

## Outputs

``` bash
PLEXUS_colon_n1149/RUVr/plots/PLEXUS_colon_n1149_eqtl_sigGene_perm100.png
PLEXUS_colon_n1149/RUVr/plots/PLEXUS_colon_n1149_eqtl_sigGene_perm100_best_k.txt
```

## Important check

``` bash
cat PLEXUS_colon_n1149/RUVr/plots/PLEXUS_colon_n1149_eqtl_sigGene_perm100_best_k.txt
```

This value is the selected covariate number used in later stages.

------------------------------------------------------------------------

# Stage 2. Final eQTL mapping using best K

## Goal

Using the selected `best_k`, run:

1.  nominal eQTL association testing
2.  final permutation testing with more permutations

For this analysis:

``` bash
n_permutations = 1000
n_chunks = 24
```

------------------------------------------------------------------------

## Step 2.1 Save best K as a shell variable

``` bash
BEST_K=$(cat PLEXUS_colon_n1149/RUVr/plots/PLEXUS_colon_n1149_eqtl_sigGene_perm100_best_k.txt)
echo $BEST_K
```

------------------------------------------------------------------------

## Step 2.2 Run nominal association testing

## Inputs

``` bash
PLEXUS_colon_n1149/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz
PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz
PLEXUS_colon_n1149/input/cov/cov_${BEST_K}.txt
```

## Command

``` bash
module load qtltools

QTLtools cis \
  --vcf PLEXUS_colon_n1149/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz \
  --bed PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz \
  --cov PLEXUS_colon_n1149/input/cov/cov_${BEST_K}.txt \
  --seed 12345 \
  --normal \
  --nominal 1 \
  --std-err \
  --out PLEXUS_colon_n1149/RUVr/nominals/nominals_cov${BEST_K}.txt
```
#### Send sbatch directly

``` bash 
sbatch scripts/nominals_eqtl.sh
```

## Output

``` bash
PLEXUS_colon_n1149/RUVr/nominals/nominals_cov${BEST_K}.txt
```

------------------------------------------------------------------------

## Step 2.3 Run final permutation testing in chunks

## Inputs

``` bash
PLEXUS_colon_n1149/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz
PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz
PLEXUS_colon_n1149/input/cov/cov_${BEST_K}.txt
```

## Command template

Submit one job per chunk, where `CHUNK_ID` ranges from `1` to `24`.

``` bash
module load qtltools

mkdir -p PLEXUS_colon_n1149/RUVr/permutations
BEST_K= ${BEST_K}

QTLtools cis \
  --vcf PLEXUS_colon_n1149/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz \
  --bed PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz \
  --cov PLEXUS_colon_n1149/input/cov/cov_${BEST_K}.txt \
  --seed 12345 \
  --permute 1000 \
  --normal \
  --std-err \
  --chunk CHUNK_ID 24 \
  --out PLEXUS_colon_n1149/RUVr/permutations/permutations_cov${BEST_K}_chunk${CHUNK_ID}.txt
```

#### Send sbatch directly

``` bash
sbatch scripts/permutation_eqtl.sh
```

## Output for each chunk

``` bash
PLEXUS_colon_n1149/RUVr/permutations/permutations_cov${BEST_K}_chunk${CHUNK_ID}.txt
```

------------------------------------------------------------------------

# Stage 3. Statistics and post-processing

## Goal

Merge the final permutation chunks, run FDR, split nominal results by chromosome, generate significant pairs, and run conditional analysis.

------------------------------------------------------------------------

## Step 3.1 Merge final permutation chunks

## Input

``` bash
PLEXUS_colon_n1149/RUVr/permutations/permutations_cov${BEST_K}_chunk*.txt
```

## Command

``` bash
cat PLEXUS_colon_n1149/RUVr/permutations/permutations_cov${BEST_K}_chunk*.txt \
  > PLEXUS_colon_n1149/RUVr/permutations/permutations_cov${BEST_K}.txt
```

## Output

``` bash
PLEXUS_colon_n1149/RUVr/permutations/permutations_cov${BEST_K}.txt
```

------------------------------------------------------------------------

## Step 3.2 Run final FDR

## Input

``` bash
PLEXUS_colon_n1149/RUVr/permutations/permutations_cov${BEST_K}.txt
```

## Command

``` bash
Rscript scripts/qtltools_runFDR_cis.R \
  PLEXUS_colon_n1149/RUVr/permutations/permutations_cov${BEST_K}.txt \
  0.05 \
  PLEXUS_colon_n1149/RUVr/permutations/PLEXUS_colon_n1149.eqtl.permutations_cov${BEST_K}.fdr0.05
```

## Outputs

``` bash
PLEXUS_colon_n1149/RUVr/permutations/PLEXUS_colon_n1149.eqtl.permutations_cov${BEST_K}.fdr0.05.significant.txt
PLEXUS_colon_n1149/RUVr/permutations/PLEXUS_colon_n1149.eqtl.permutations_cov${BEST_K}.fdr0.05.thresholds.txt
```

The `.significant.txt` file is the final list of significant eGenes.

The `.thresholds.txt` file is used to define significant variant-gene pairs from the nominal results and is also used for conditional analysis.

------------------------------------------------------------------------

## Step 3.3 Split nominal results by chromosome

## Input

``` bash
PLEXUS_colon_n1149/RUVr/nominals/nominals_cov${BEST_K}.txt
```

## Command

``` bash
mkdir -p PLEXUS_colon_n1149/RUVr/nominals/splitchr

chmod +x scripts/rename_nominals_singlefile.sh

scripts/rename_nominals_singlefile.sh \
  PLEXUS_colon_n1149/RUVr/nominals/nominals_cov${BEST_K}.txt \
  PLEXUS_colon_n1149/RUVr/nominals/splitchr
```

## Outputs

``` bash
PLEXUS_colon_n1149/eqtl/RUVr/nominals/splitchr/chr1.txt
PLEXUS_colon_n1149/eqtl/RUVr/nominals/splitchr/chr2.txt
...
PLEXUS_colon_n1149/eqtl/RUVr/nominals/splitchr/chr22.txt
PLEXUS_colon_n1149/eqtl/RUVr/nominals/splitchr/chrX.txt
```

------------------------------------------------------------------------

## Step 3.4 Generate significant variant-gene pairs

## Inputs

``` bash
PLEXUS_colon_n1149/RUVr/nominals/splitchr/chr*.txt
PLEXUS_colon_n1149/RUVr/permutations/PLEXUS_colon_n1149.eqtl.permutations_cov${BEST_K}.fdr0.05.thresholds.txt
```

## Command

``` bash
chmod +x scripts/generate_sigpairs.sh

scripts/generate_sigpairs.sh \
  PLEXUS_colon_n1149/RUVr/nominals/splitchr \
  PLEXUS_colon_n1149/RUVr/permutations/PLEXUS_colon_n1149.eqtl.permutations_cov${BEST_K}.fdr0.05.thresholds.txt \
  PLEXUS_colon_n1149/RUVr/nominals
```

## Output

``` bash
PLEXUS_colon_n1149/RUVr/nominals/sigpairs.txt
```

This file contains significant eQTL variant-gene pairs.

------------------------------------------------------------------------

## Step 3.5 Run conditional analysis in chunks

## Inputs

``` bash
PLEXUS_colon_n1149/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz
PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz
PLEXUS_colon_n1149/input/cov/cov_${BEST_K}.txt
PLEXUS_colon_n1149/RUVr/permutations/PLEXUS_colon_n1149.eqtl.permutations_cov${BEST_K}.fdr0.05.thresholds.txt
```

## Command template

Submit one job per chunk, where `CHUNK_ID` ranges from `1` to `24`.

``` bash
module load qtltools

mkdir -p PLEXUS_colon_n1149/RUVr/cond

QTLtools cis \
  --vcf PLEXUS_colon_n1149/input/vcf/PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz \
  --bed PLEXUS_colon_n1149/input/bed/PLEXUS_colon_n1149_eqtl_qtltools.bed.gz \
  --cov PLEXUS_colon_n1149/input/cov/cov_${BEST_K}.txt \
  --mapping PLEXUS_colon_n1149/RUVr/permutations/PLEXUS_colon_n1149.eqtl.permutations_cov${BEST_K}.fdr0.05.thresholds.txt \
  --normal \
  --std-err \
  --chunk CHUNK_ID 24 \
  --out PLEXUS_colon_n1149/RUVr/cond/cond_cov${BEST_K}_chunk${CHUNK_ID}.txt
```
#### Send sbatch directly

``` bash
sbatch scripts/cond_eqtl.sh
``` 

## Output for each chunk

``` bash
PLEXUS_colon_n1149/RUVr/cond/cond_cov${BEST_K}_chunk${CHUNK_ID}.txt
```

------------------------------------------------------------------------

## Step 3.6 Merge and filter conditional results

## Inputs

``` bash
PLEXUS_colon_n1149/RUVr/cond/cond_cov${BEST_K}_chunk*.txt
```

## Commands

``` bash
cat PLEXUS_colon_n1149/RUVr/cond/cond_cov${BEST_K}_chunk*.txt \
  > PLEXUS_colon_n1149/RUVr/cond/combined_cond_cov${BEST_K}.txt

awk -F' ' '$23 == 1' \
  PLEXUS_colon_n1149/RUVr/cond/combined_cond_cov${BEST_K}.txt \
  > PLEXUS_colon_n1149/RUVr/cond/cond_cov${BEST_K}_sgenes.txt
```

## Outputs

``` bash
PLEXUS_colon_n1149/RUVr/cond/combined_cond_cov${BEST_K}.txt
PLEXUS_colon_n1149/RUVr/cond/cond_cov${BEST_K}_sgenes.txt
```

The filtered file keeps rows where column 23 equals `1`.

------------------------------------------------------------------------

# Stage 4. Organize final outputs

## Goal

Copy final files into shared summary-statistics folders.

## Final output directories

``` bash
PLEXUS_colon_n1149/RUVr/nominals/splitchr
PLEXUS_colon_n1149/RUVr/sgenes_sigpairs/sgenes_sigpairs
```

## Step 4.1 Create final significant pair file

## Input

``` bash
PLEXUS_colon_n1149/RUVr/nominals/sigpairs.txt
```

## Command

``` bash
mkdir -p PLEXUS_colon_n1149/RUVr/sgenes_sigpairs
cp PLEXUS_colon_n1149/RUVr/nominals/sigpairs.txt PLEXUS_colon_n1149/RUVr/sgenes_sigpairs/sigpairs.txt
```

## Output

``` bash
PLEXUS_colon_n1149/RUVr/sgenes_sigpairs/sigpairs.txt
```

------------------------------------------------------------------------

## Step 4.2 Copy final significant eGenes

## Input

``` bash
PLEXUS_colon_n1149/RUVr/permutations/PLEXUS_colon_n1149.eqtl.permutations_cov${BEST_K}.fdr0.05.significant.txt
```

## Command

``` bash
cp PLEXUS_colon_n1149/RUVr/permutations/PLEXUS_colon_n1149.eqtl.permutations_cov${BEST_K}.fdr0.05.significant.txt \
   PLEXUS_colon_n1149/RUVr/sgenes_sigpairs/IBD.eqtl.sgenes.txt
```

## Output

``` bash
PLEXUS_colon_n1149/RUVr/sgenes_sigpairs/sgenes_sigpairs/IBD.eqtl.sgenes.txt
```

------------------------------------------------------------------------

## Step 4.3 Copy conditional eGenes

## Input

``` bash
PLEXUS_colon_n1149/RUVr/cond/cond_cov${BEST_K}_sgenes.txt
```

## Command

``` bash
cp PLEXUS_colon_n1149/RUVr/cond/cond_cov${BEST_K}_sgenes.txt \
   PLEXUS_colon_n1149/RUVr/sgenes_sigpairs//sgenes_sigpairs/IBD.eqtl.sgenes.cond.txt
```

## Output

``` bash
PLEXUS_colon_n1149/RUVr/sgenes_sigpairs/IBD.eqtl.sgenes.cond.txt
```

------------------------------------------------------------------------

# Recommended final files to share with lab members

The final organized files are:

| File                                       | Description                                               |
|-------------------------------|-----------------------------------------|
| `nominals/IBD/eqtl/chr*.txt`               | Chromosome-level nominal eQTL results.                    |
| `sgenes_sigpairs/IBD.eqtl.sigpairs.txt`    | Significant eQTL variant-gene pairs.                      |
| `sgenes_sigpairs/IBD.eqtl.sgenes.txt`      | Significant eGenes from permutation FDR.                  |
| `sgenes_sigpairs/IBD.eqtl.sgenes.cond.txt` | Conditional eQTL results filtered to independent signals. |

------------------------------------------------------------------------

# Minimal folder structure

A clean project folder should look like this:

``` text
scripts/
├── permutationtest_eqtl.sh
├── nominals_eqtl.sh
├── permutation_eqtl.sh
├── cond_eqtl.sh 
├── generate_sigpairs.sh 
├── rename_nominals_singlefile.sh 
├── plot_RUVK_SigGene_eqtl.R
├── qtltools_runFDR_cis.R
└── README_eqtl_qtltools_pipeline.md


PLEXUS_colon_n1149/
├── input/
    ├── bed/
    │   ├── PLEXUS_colon_n1149_eqtl_qtltools.bed.gz
    │   └── PLEXUS_colon_n1149_eqtl_qtltools.bed.gz.tbi
    ├── cov/
    │   ├── cov_0.txt
    │   ├── ...
    │   └── cov_30.txt
    ├── vcf/
    │   └── PLEXUS_colon_n1149_R2_0.3_rsidsnp_maf001.vcf.gz
└── RUVr/
    ├── permutationtest/
    │   └── perm100/
    ├── plots/
    │   ├── PLEXUS_colon_n1149_eqtl_sigGene_perm100.png
    │   └── PLEXUS_colon_n1149_eqtl_sigGene_perm100_best_k.txt
    ├── nominals/
    │   ├── nominals_covK.txt
    │   ├── sigpairs.txt
    │   └── splitchr/
    │       ├── chr1.txt
    │       ├── ...
    │       └── chrX.txt
    ├── permutations/
    │   ├── permutations_covK_chunk1.txt
    │   ├── ...
    │   ├── permutations_covK_chunk24.txt
    │   ├── permutations_covK.txt
    │   ├── PLEXUS_colon_n1149.eqtl.permutations_covK.fdr0.05.significant.txt
    │   └── PLEXUS_colon_n1149.eqtl.permutations_covK.fdr0.05.thresholds.txt
    ├── cond/
    │   ├── cond_covK_chunk1.txt
    │   ├── ...
    │   ├── cond_covK_chunk24.txt
    │   ├── combined_cond_covK.txt
    │   └── cond_covK_sgenes.txt
    └── sgenes_sigpairs/
    │   ├── IBD.eqtl.sigpairs.txt
    │   ├── IBD.eqtl.sgenes.txt
    │   └── IBD.eqtl.sgenes.cond.txt
```

Replace `K` with the selected best covariate number.

------------------------------------------------------------------------

#Headers to output files (Optional):

To create header-added copies without changing the original files:

```bash
# Nominal/Sigpairs output header 
printf '%s\n' 'phe_id phe_chr phe_from phe_to phe_strd n_var_in_cis dist_phe_var var_id var_chr var_from var_to nom_pval r_squared slope slope_se best_hit' \
  | cat - PLEXUS_colon_n1149/RUVr/nominals/nominals_cov${BEST_K}.txt \
  > PLEXUS_colon_n1149/RUVr/nominals/nominals_cov${BEST_K}.with_header.txt

# Final permutation output header
printf '%s\n' 'phe_id phe_chr phe_from phe_to phe_strd n_var_in_cis dist_phe_var var_id var_chr var_from var_to dof1 dof2 bml1 bml2 nom_pval r_squared slope slope_se adj_emp_pval adj_beta_pval' \
  | cat - PLEXUS_colon_n1149/RUVr/permutations/permutations_cov${BEST_K}.txt \
  > PLEXUS_colon_n1149/RUVr/permutations/permutations_cov${BEST_K}.with_header.txt


# Final FDR significant eGenes header
printf '%s\n' 'phe_id phe_chr phe_from phe_to phe_strd n_var_in_cis dist_phe_var var_id var_chr var_from var_to dof1 dof2 bml1 bml2 nom_pval r_squared slope slope_se adj_emp_pval adj_beta_pval qval sig.thresholds' \
  | cat - PLEXUS_colon_n1149/RUVr/sgenes_sigpairs/IBD.eqtl.sgenes.txt \
  > PLEXUS_colon_n1149/RUVr/sgenes_sigpairs/IBD.eqtl.sgenes.with_header.txt

# Conditional output header
printf '%s\n' 'phe_id phe_chr phe_from phe_to phe_strd n_var_in_cis dist_phe_var var_id var_chr var_from var_to rank fwd_pval fwd_r_squared fwd_slope fwd_slope_se fwd_best_hit fwd_sig bwd_pval bwd_r_squared bwd_slope bwd_slope_se bwd_best_hit bwd_sig' \
  | cat - PLEXUS_colon_n1149/RUVr/sgenes_sigpairs/IBD.eqtl.sgenes.cond.txt \
  > PLEXUS_colon_n1149/RUVr/sgenes_sigpairs/IBD.eqtl.sgenes.cond.with_header.txt
```


# Workflow summary

This eQTL pipeline first chooses the best number of covariates by running a small permutation test across several covariate files. After choosing the best `K`, it runs the final eQTL analysis in two parts: nominal testing to get all variant-gene association statistics, and permutation testing to identify significant eGenes and empirical FDR thresholds. The nominal results are then split by chromosome, and the threshold file is used to extract significant variant-gene pairs. Finally, QTLtools conditional analysis is run to identify independent eQTL signals, and all final files are copied into shared summary-statistics folders.
