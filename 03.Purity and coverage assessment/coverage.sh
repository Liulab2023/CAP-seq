#!/bin/bash
#SBATCH -J coverage
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00

WORK_DIR="dir_to_set_working_dir"
cd $WORK_DIR
mkdir coverage
mkdir coverage/species

cd ./coverage/species
FILE_PATH="path_to_fastq_list_species.txt" 

REF="species.fasta"
samtools faidx species.fasta
awk '{sum+=$2} END{print sum}' species.fasta.fai > total_bases.txt


while IFS= read -r line
do
    echo "Processing file: $line"
    mkdir ${line%_minimap2_output.txt}
    minimap2 -ax map-ont --secondary=no -t 32 $REF $WORK_DIR/bbmap_output_many_seqs/${line%_minimap2_output.txt}.fastq > ${line%_minimap2_output.txt}/aln.txt
    samtools view -bS ${line%_minimap2_output.txt}/aln.txt > ${line%_minimap2_output.txt}/aln.bam
    samtools sort ${line%_minimap2_output.txt}/aln.bam -o ${line%_minimap2_output.txt}/aln_sort.bam
    samtools depth ${line%_minimap2_output.txt}/aln_sort.bam > ${line%_minimap2_output.txt}/aln_sort_depth
    awk '{if($3>0) count++} END{print count}' ${line%_minimap2_output.txt}/aln_sort_depth > ${line%_minimap2_output.txt}/covered_bases.txt
    paste ${line%_minimap2_output.txt}/covered_bases.txt total_bases.txt | awk '{print $1/$2*100}' > ${line%_minimap2_output.txt}/coverage_percentage.txt

done < "$FILE_PATH"
