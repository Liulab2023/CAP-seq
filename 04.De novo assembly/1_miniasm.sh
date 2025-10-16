#!/bin/bash
#SBATCH -J miniasm
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00


WORK_DIR="dir_to_set_working_dir"
cd $WORK_DIR

mkdir -p miniasm_paf_output
mkdir -p miniasm_fasta_output
dir="./bbmap_output_many_seqs" 
for fastq_file in "$dir"/*.fastq
do 
	base_name=$(basename "$fastq_file" .fastq)
    echo "Processing $base" 
	
    minimap2 -x ava-ont -t 32 "$fastq_file" "$fastq_file" | gzip -1 > ./miniasm_paf_output/"$base_name".paf.gz

    miniasm -f "$fastq_file" ./miniasm_paf_output/"$base_name".paf.gz -i 0.03 -m 50 -s 50 -e 1 -g 2000 -F 0.5 -c 1 -n 0 -1 -2 > ./miniasm_paf_output/"$base_name".gfa
    awk '/^S/{print ">"$2"\n"$3}' ./miniasm_paf_output/"$base_name".gfa | fold > ./miniasm_fasta_output/"$base_name".fasta
    rm ./miniasm_paf_output/"$base_name".paf.gz
    rm ./miniasm_paf_output/"$base_name".gfa
done

