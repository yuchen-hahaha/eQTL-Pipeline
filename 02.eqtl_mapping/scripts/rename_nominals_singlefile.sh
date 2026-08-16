#!/bin/bash
# for eqtl
# Usage: ./scripts/rename_nominals_singlefile.sh <input_file> <output_dir>
# Example: ./scripts/rename_nominals_singlefile.sh UNC_colon_n255/eqtl/nominals/qtltools/nominals_cov18.txt UNC_colon_n255/eqtl/nominals/qtltools

input_file="$1"
output_dir="$2"

if [[ -z "$input_file" || -z "$output_dir" ]]; then
  echo "Usage: $0 input_file output_directory"
  exit 1
fi

mkdir -p "$output_dir"

# Split by column 2 (chromosome), assuming tab-delimited file
awk -v outdir="$output_dir" '
  {
    chr = $2
    print >> outdir "/" chr ".txt"
  }
' "$input_file"

echo "Done: Split completed with only present chromosomes in ${output_dir}."