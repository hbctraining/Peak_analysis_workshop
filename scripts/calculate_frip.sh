#!/bin/bash

# Arguments check, requires 3 arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bam_dir> <peak_dir> <output_file>"
    exit 1
fi

# Assiging variables
bam_dir=$1
peak_dir=$2
output_file=$3

# create output file
echo -e "Sample\tFripscore" > "$output_file"


for bam_file in "$bam_dir"/*.bam; do
    # Extract the sample name (remove .mLb.clN.sorted.bam extension)
    sample_name=$(basename "$bam_file" .mLb.clN.sorted.bam)

    # peak file
    peak_file="$peak_dir/${sample_name}_peaks.narrowPeak"

    # Check if peak file exists
    if [ ! -f "$peak_file" ]; then
        echo "Peak file not found for $sample_name"
        continue
    fi

    # bedtools intersect between the BAM and peak file and save it (temporarily)
    bedtools intersect -abam "$bam_file" -b "$peak_file" > "${sample_name}.bedtools.out"

    # Count reads in the intersect file produced above, this will be the numerator for frip calculation
    num_reads_in_peaks=$(samtools view -c "${sample_name}.bedtools.out")

    # Count total reads in the BAM file this will be the denominator for frip calculation
    total_reads=$(samtools view -c "$bam_file")

    # Calculate FRiP 
    if [ "$total_reads" -gt 0 ]; then
        frip_score=$(echo "scale=5; $num_reads_in_peaks / $total_reads" | bc)
    else
        frip_score="NA"  # If total reads is zero, mark FRiP score as not available
    fi

    #  Save the result 
    echo -e "${sample_name}\t${frip_score}" >> "$output_file"

    # remove the temporary intersect file
    rm "${sample_name}.bedtools.out"
done

echo "FRiP scores saved to $output_file"
