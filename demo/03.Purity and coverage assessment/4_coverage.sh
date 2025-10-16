#!/bin/bash

WORK_DIR="/data/disk3-2/code-test"
cd $WORK_DIR
mkdir coverage
mkdir coverage/Bacillus

cd ./coverage/Bacillus
FILE_PATH="../../fastq_list_Bacillus.txt" 

REF="Bacillus.fasta"
samtools faidx Bacillus.fasta
awk '{sum+=$2} END{print sum}' Bacillus.fasta.fai > total_bases.txt


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
