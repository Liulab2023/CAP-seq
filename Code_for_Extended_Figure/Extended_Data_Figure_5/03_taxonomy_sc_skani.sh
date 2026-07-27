cd /path/to/sc
mkdir -p miniasm_paf_output miniasm_fasta_output
dir="./sc_split_bc"

for fastq_file in "$dir"/*.fq.gz; do
  base_name=$(basename "$fastq_file" .fq.gz)
  echo "Processing $base_name"

  minimap2 -x ava-ont -t 32 "$fastq_file" "$fastq_file" | gzip -1 > ./miniasm_paf_output/"$base_name".paf.gz

  miniasm -f "$fastq_file" ./miniasm_paf_output/"$base_name".paf.gz \
          -i 0.03 -m 50 -s 50 -e 1 -g 2000 -F 0.5 -c 1 -n 0 -1 -2 \
          > ./miniasm_paf_output/"$base_name".gfa

  awk '/^S/{print ">"$2"\n"$3}' ./miniasm_paf_output/"$base_name".gfa | fold \
      > ./miniasm_fasta_output/"$base_name".fasta

  rm ./miniasm_paf_output/"$base_name".paf.gz ./miniasm_paf_output/"$base_name".gfa
done

skani search ./miniasm_fasta_output/*.fasta \