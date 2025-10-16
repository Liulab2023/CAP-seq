#!/bin/bash
#SBATCH -J quast
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00


WORK_DIR="dir_to_set_working_dir"
cd $WORK_DIR

mkdir -p quast_output

for fasta_file in ./miniasm_fasta_output/*.fasta; do
    base_name=$(basename "$fasta_file" .fasta)

    output_dir=./quast_output/"$base_name"

    nanopore_file=./bbmap_output_many_seqs/"$base_name".fastq

    if [ -f "$nanopore_file" ]; then
        echo "Processing $fasta_file with Nanopore file $nanopore_file..."
        
        quast "$fasta_file" -o "$output_dir" -t 32 -f --nanopore "$nanopore_file"
    else
        echo "Warning: Nanopore file $nanopore_file not found. Skipping $fasta_file."
    fi
done
