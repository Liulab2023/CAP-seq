#!/bin/bash
#SBATCH -J demux
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00
#SBATCH --exclude=dm-c02n03

# code path
CODE_DIR="dir_to_set_code_dir"
WORK_DIR="dir_to_set_working_dir"

# search for linker between barcode and genomic DNA
python "$CODE_DIR"/search_fastq.py -i "$WORK_DIR" -o "$WORK_DIR"/search_output -s GTCTCGTGGGCTCGG -t 26

# search for valid barcode
python "$CODE_DIR"/search_barcode.py -i "$WORK_DIR"/search_output -e "$CODE_DIR"/barcode_seq.xlsx -o "$WORK_DIR"/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode.py -i "$WORK_DIR"/search_output -e "$CODE_DIR"/barcode_seq.xlsx -o "$WORK_DIR"/barcode_match2_output -s1 ACAG -s2 GTCA -r

cd $WORK_DIR
find ./barcode_match*_output/ -type f -name "*.fastq" -print0 | xargs -0 cat >> ./bbmap_input_nocut.fastq

cutadapt -g AGATGTGTATAAGAGACAG -O 9 -a CTGTCTCTTATACACATCT -O 10 -m 50 --cores 32 -e 0.15 -o ./bbmap_input.fastq ./bbmap_input_nocut.fastq 

OUTPUT_DIR="$WORK_DIR/bbmap_output"
demuxbyname.sh in=bbmap_input.fastq out=$OUTPUT_DIR/out_%.fastq prefixmode=f delimiter=_ -Xmx30g ow=f

python $CODE_DIR/statistical_bbmap.py -i $WORK_DIR/bbmap_output -o $WORK_DIR/ -s 1000
