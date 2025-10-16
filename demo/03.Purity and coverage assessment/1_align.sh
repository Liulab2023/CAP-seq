#!/bin/bash

WORK_DIR="/data/disk3-2/code-test"
cd $WORK_DIR

dir="./bbmap_output_many_seqs" 
mkdir minimap2_output_many_seqs

for file in "$dir"/*.fastq
do 
    base=$(basename $file)
    echo "Processing $base" 

    minimap2 -ax map-ont --secondary=no mock_community.fasta "$file" > "./minimap2_output_many_seqs/${base%.fastq}_minimap2_output.txt"
done

#After align, you can swicth to purity.ipynb to calculate Purity

