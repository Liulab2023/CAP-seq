dir="sc_split_bc" 

for fastq_file in "$dir"/*.fq.gz; do
  base_name=$(basename "$fastq_file" .fq.gz)
  echo "Processing $base_name"

  metaphlan "$fastq_file" --nproc 20 --input_type fastq \
    -o ./metaout/"$base_name".txt --mapout ./metaout/"$base_name"_map.txt \
    --db_dir /path/to/metaphlan --ignore_eukaryotes --split_reads \
    --minimap2_exe /path/to/minimap2 --stat_q 0.01 --avoid_disqm
done