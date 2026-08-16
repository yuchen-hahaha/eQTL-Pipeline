#!/bin/bash
#SBATCH --job-name=sigpairs
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --output="logs/eqtlsigpairs.%A.out"
#SBATCH --error="logs/eqtlsigpairs.%A.err"

 

# Usage:./scripts/generate_sigpairs.sh <input_dir> <threshold_file> <output_dir>
# Ex: ./scripts/generate_sigpairs.sh UNC_colon_n255/sqtl/nominals/RUVg UNC_colon_n255/sqtl/permutations/RUVg/UNC.sqtl.cov2.fdr005.thresholds.txt UNC_colon_n255/sqtl/nominals/RUVg 
 
input_dir="$1"
threshold_file="$2"
output_dir="$3"



if [[ -z "$input_dir" || -z "$threshold_file" || -z "$output_dir" ]]; then
  echo "Usage: $0 <input_dir> <threshold_file> <output_dir><project><data>"
  exit 1
fi

#mkdir -p "$output_dir"

# Step 1: Filter each chr*.txt by $12 < 0.05
for file in "$input_dir"/chr*.txt; do
  base=$(basename "$file" .txt)
  echo "Processing $base"
  awk '$12 < 0.05' "$file" > "${output_dir}/${base}.p005.txt"
done

# Step 2: Concatenate all .p005.txt files
cat "$output_dir"/*.p005.txt > "${output_dir}/nominals.p005.txt"

# Step 3: Filter by thresholds
output_sigpairs="${output_dir}/sigpairs.txt"

awk 'BEGIN {FS=OFS=" "}
FNR==NR {threshold[$1] = $2; next}
$1 in threshold && $12 <= threshold[$1]' \
"$threshold_file" \
"${output_dir}/nominals.p005.txt" \
> "$output_sigpairs"

# Optional cleanup
rm "$output_dir"/chr*.p005.txt
rm "$output_dir"/nominals.p005.txt
