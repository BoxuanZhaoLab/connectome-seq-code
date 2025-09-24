#!/bin/bash

input_dir="/home/zhouwan2/expansion_link/RNA_seq_processed_data/cseq1_virus_lib_results_080325"
output_dir="$input_dir"

for hist_file in "$input_dir"/hist_*.txt; do
    base_name=$(basename "$hist_file" .txt)
    output_file="$output_dir/complexity_${base_name#hist_}.txt"

    preseq lc_extrap \
        -H \
        -o "$output_file" \
        "$hist_file"
done
