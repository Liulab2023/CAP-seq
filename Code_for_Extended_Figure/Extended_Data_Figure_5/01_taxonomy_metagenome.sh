# Quality trimming
fastp -i metagenome_raw_1.fq.gz -I metagenome_raw_2.fq.gz \
      -o metagenome_1.trim.fastq.gz -O metagenome_2.trim.fastq.gz \
      --length_required 50 --thread 8 \
      -j metagenome_fastp.json -h metagenome_fastp.html

# Remove human reads (hg38)
bowtie2 -x /path/to/ref/hg38_index \
        -1 metagenome_1.trim.fastq.gz -2 metagenome_2.trim.fastq.gz \
        -p 16 --very-sensitive-local \
        --un-conc-gz metagenome_clean_%.fastq.gz \
        > /dev/null

# Run MetaPhlAn4
metaphlan metagenome_clean_1.fastq.gz,metagenome_clean_2.fastq.gz \
         -s metagenome.sam.bz2 --mapout metagenome.bowtie2.bz2 \
         --nproc 16 --input_type fastq -o metagenome_tax.txt \
         --db_dir /path/to/metaphlan

# Convert SGB to GTDB nomenclature
sgb_to_gtdb_profile.py -i metagenome_tax.txt -o metagenome_tax_gtdb.txt