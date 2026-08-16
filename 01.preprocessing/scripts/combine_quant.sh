#!/bin/bash
#SBATCH --job-name=gene
#SBATCH --time=12:00:00
#SBATCH --mem=5G
#SBATCH -n 1
#SBATCH -N 1

project=$1
sample_list=$2  # text file with one sample name per line

quant_dir="${project}/quant"
tmp_file="${quant_dir}/${project}_tmp.txt"
output_file="${quant_dir}/${project}_gene_quant.txt"

# Get the first sample in the list and use its file for the first 6 columns
first_sample=$(head -n 1 "$sample_list")
first_file="${quant_dir}/${first_sample}.gene.count.bed"   # <--- changed

cut -f1-6 "$first_file" > "$tmp_file"

while read sample; do
    file="${quant_dir}/${sample}.gene.count.bed"           # <--- changed
    if [[ ! -f "$file" ]]; then
        echo "Warning: file not found: $file"
        continue
    fi

    cut -f7 "$file" > "${quant_dir}/quant_${sample}.txt"
    paste "$tmp_file" "${quant_dir}/quant_${sample}.txt" > tmp && mv tmp "$tmp_file"
    rm "${quant_dir}/quant_${sample}.txt"
done < "$sample_list"

mv "$tmp_file" "$output_file"
