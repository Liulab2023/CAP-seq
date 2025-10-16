#!/bin/bash
#SBATCH -J basecalling_qc
#SBATCH -p dm_pub_gpu
#SBATCH -N 1 
#SBATCH -n 17
#SBATCH --time=3-00:00:00
#SBATCH --gres=gpu:1 

# dorado path
DORADO_PATH="dir_to_dorado"
KIT_PATH="dir_to_{dna_r10.4.1_e8.2_400bps_sup@v5.0.0/}"
POD5_DIR="dir_to_pod5"

# work path
WORK_DIR="dir_to_set_working_dir"

# code path
CODE_DIR="dir_to_set_code_dir"

$DORADO_PATH basecaller $KIT_PATH $POD5_DIR --emit-fastq > "$WORK_DIR/calls.fastq"
NanoFilt -q 15 -l 20 --maxlength 60000 calls.fastq > calls_filt.fastq