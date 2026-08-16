import os

# Base directory
base_dir = '/proj/fureylab/data/RNA-seq/human/colon_tissue/nonIBD'

# Output file
output_file = 'input/UNC_colon_nonIBD_bampaths.txt'

# Collect BAM file paths
bam_paths = []

# Go through each item in the base directory
for sample in os.listdir(base_dir):
    sample_dir = os.path.join(base_dir, sample)
    if not os.path.isdir(sample_dir):
        continue  # Skip non-directories

    star_dir = os.path.join(sample_dir, 'snakemakeRNA_WASP_hg38_v47', 'star')
    if os.path.isdir(star_dir):
        # Look for BAM file named like: <sample>.Aligned.sortedByCoord.out.bam
        bam_filename = f"{sample}.Aligned.sortedByCoord.out.bam"
        bam_path = os.path.join(star_dir, bam_filename)
        if os.path.isfile(bam_path):
            bam_paths.append(bam_path)

# Write to output file
with open(output_file, 'w') as f:
    for path in sorted(bam_paths):
        f.write(path + '\n')

print(f"Collected {len(bam_paths)} BAM paths into {output_file}")
