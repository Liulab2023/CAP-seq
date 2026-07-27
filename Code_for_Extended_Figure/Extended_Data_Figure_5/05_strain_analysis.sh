#!/usr/bin/env bash

##StrainPhlAn4 (Phylogenomic markers)
mkdir -p consensus_markers db_markers output

sample2markers.py -i meta.sam.bz2 -o consensus_markers -n 16 \
                  -d /path/to/metaphlan/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl

extract_markers.py -c t__SGB1855 -o db_markers/ \
                   -d /path/to/metaphlan/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl

strainphlan -s consensus_markers/*.json \
            -m db_markers/t__SGB1855.fna \
            -r reference_genome/GCF_000009925.1_ASM992v1_genomic.fna \
            -o output -n 8 -c t__SGB1855 --mutation_rates \
            -d /path/to/metaphlan/mpa_vJan25_CHOCOPhlAnSGB_202503.pkl

add_metadata_tree.py -t output/RAxML_bestTree.t__SGB1855.StrainPhlAn4.tre \
                     -f metadata.txt -m subjectID --string_to_remove .fastq.gz

plot_tree_graphlan.py -t output/RAxML_bestTree.t__SGB1855.StrainPhlAn4.tre.metadata -m subjectID

##inStrain (Read mapping & variant profiling)
# Metagenome
bowtie2 -p 10 -x ./ref/ref -1 metagenome_clean_1.fastq.gz -2 metagenome_clean_2.fastq.gz | \
  samtools sort -O bam -@ 10 -o - > metagenome.bam

inStrain profile metagenome.bam ./ref_genomes/GCF_000009925.1_ASM992v1_genomic.fna \
              -o metagenome.IS -p 16

# Pooled single-cell
minimap2 -ax map-ont ../ref_genomes/GCF_000009925.1_ASM992v1_genomic.fna sc.fq.gz | \
  samtools sort -O bam -@ 10 -o - > sc.bam

inStrain profile sc.bam ./ref_genomes/GCF_000009925.1_ASM992v1_genomic.fna \
              -o sc.IS -p 16


##Snippy (SNP calling in contigs)
INPUT="./query_list.tsv"
OUTDIR="./contig_snp_total"
mkdir -p $OUTDIR

run_snippy() {
    query="$1"; ref="$2"
    sample_name=$(basename "$query" .fasta)
    sample_dir="$OUTDIR/${sample_name}"
    snippy --cpus 2 --outdir "$sample_dir" --ref "$ref" --ctgs "$query" --quiet
    rm -rf "$OUTDIR/${sample_name}/reference"
}
export -f run_snippy

tail -n +2 "$INPUT" | parallel --colsep '\t' -j $(nproc) run_snippy {1} {2}