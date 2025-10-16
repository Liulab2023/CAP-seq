#!/bin/bash
#SBATCH -J align
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00

WORK_DIR="dir_to_set_working_dir"
cd $WORK_DIR

dir="./bbmap_output_many_seqs" 
mkdir minimap2_output_many_seqs

for file in "$dir"/*.fastq
do 
    base=$(basename $file)
    echo "Processing $base" 

    minimap2 -ax map-ont --secondary=no {dir_to_reference_genomes} "$file" > "./minimap2_output_many_seqs/${base%.fastq}_minimap2_output.txt"
done

#After align, you can swicth to purity.ipynb to calculate Purity

