![GitHub watchers](https://img.shields.io/github/watchers/Liulab2023/CAP-seq?style=flat&logo=Github&color=blue) [![Stars](https://img.shields.io/github/stars/Liulab2023/CAP-seq?style=flat&logo=GitHub&color=blue)](https://github.com/Liulab2023/CAP-seq/stargazers)
# Analysis Workflow for **Large-scale single-cell long-read genomics enables high-resolution microbiome profiling using cap-seq**

## Description

This is a comprehensive bioinformatics pipeline for single-cell genomics data analysis, including basecalling, barcode demultiplexing, quality assessment, and de novo assembly of single amplified genomes (SAGs).

## Contents

- [System Requirements](#System-Requirements)
- [Installation Guide](#Installation-Guide)
- [Demo](#Demo)
- [Instructions for Use](#Instructions-for-Use)
- [Troubleshooting](#Troubleshooting)
- [Citation](#Citation)
- [Support](#Support)


## System Requirements

### Software Dependencies

- **Operating Systems**: Linux (tested on ![Static Badge](https://img.shields.io/badge/Ubuntu-24.04-blue))
- **Python**: ![Static Badge](https://img.shields.io/badge/Version-3.12.3-blue)
- **Key Tools and Versions**:
  - Dorado ![Static Badge](https://img.shields.io/badge/Version-0.8.3-blue) (basecalling)
  - NanoFilt  ![Static Badge](https://img.shields.io/badge/Version-2.8.0-blue)(quality filtering)
  - Cutadapt ![Static Badge](https://img.shields.io/badge/Version-5.1-blue) (adapter trimming)
  - BBMap ![Static Badge](https://img.shields.io/badge/Version-39.37-blue) (demultiplexing)
  - Minimap2 ![Static Badge](https://img.shields.io/badge/Version-2.30-blue) (alignment)
  - Samtools ![Static Badge](https://img.shields.io/badge/Version-1.22.1-blue) (BAM processing)
  - Miniasm ![Static Badge](https://img.shields.io/badge/Version-r179-blue) (de novo assembly)
  - QUAST ![Static Badge](https://img.shields.io/badge/Version-5.3.0-blue) (assembly quality assessment)
  - CheckM ![Static Badge](https://img.shields.io/badge/Version-1.0.13-blue) (genome completeness assessment)

### GPU Requirements

- **Dorado basecaller** requires a GPU with proper VRAM (recommended 16GB+ for large datasets)
- **CUDA**: ![Static Badge](https://img.shields.io/badge/Version-11.0-blue) or higher (compatible with Dorado ![Static Badge](https://img.shields.io/badge/Version-0.8.3-blue))
- **NVIDIA drivers**: Latest stable version recommended

### Python Dependencies

The pipeline requires the following Python packages (automatically installed via environment.yaml):

```
argparse
pandas
biopython
multiprocessing
glob
shutil
collections
matplotlib
pysam
numpy
```

*Complete dependency list available in `environment.yaml`*

### Tested Environments

- Ubuntu 24.04 LTS with Conda environment
- Python 3.12.3 environment
- The pipeline has been validated with the provided Conda environment specification

### Hardware Requirements

- **Standard hardware**: No specialized hardware required
- **Minimum**: 32GB RAM, 4-core CPU, GPU
- **Recommended**: 64GB+ RAM, 16+ core CPU for efficient processing
- **Storage**: Sufficient space for raw POD5 files and intermediate files (typically 500-1500GB depending on dataset size)



## Installation Guide

### Instructions

1. **Clone the repository**:

   ```
   git clone https://github.com/Liulab2023/CAP-seq
   cd CAP-seq
   ```

2. **Create and activate Conda environment from environment.yaml**:

   ```
   conda env create -f environment.yaml
   conda activate captain
   ```

3. **Verify installation**:

   ```
   python -c "import pandas, Bio, pysam; print('All packages imported successfully')"
   ```

4. **Download and install external tools**:

   - Dorado: Download from Oxford Nanopore Technologies
   - Ensure all tools in the dependency list are accessible in PATH

### Typical Install Time

- **Conda environment setup**: 2-3h minutes (due to comprehensive bioinformatics packages)
- **Complete installation**: 3-4h minutes on a standard desktop computer



## Demo

### Instructions to Run on Demo Data
You can get complete demo data and results from https://zenodo.org/records/17385053.

1. **Prepare demo data**:

   - Place POD5 files in `dir_to_pod5/` directory
   - Update paths in configuration scripts

2. **Execute the pipeline sequentially**

   - 01-basecalling_and_qc

   - 02-Barcode identification and demultiplexing
   - 03-Purity_and_coverage_assessment

### Expected Output

1. **Basecalling**:

   - `calls.fastq`: Raw basecalled reads

   - `calls_filt.fastq`: Quality-filtered reads

     Since the raw `pod5 files` is too big, we just upload a demo `fastq file calls.fastq` and `calls_filt.fastq` can be get through command `NanoFilt -q 15 -l 200 --maxlength 60000 calls.fastq > calls_filt.fastq`.

2. **Demultiplexing**:

   - Individual FASTQ files for each SAG in `bbmap_output_many_seqs/` directory
   - Barcode matching statistics (see `statistical.txt` and `summary.txt`)
   - Quality control plots (see `histogram_and_curve_gt_1_log_scale`.png)

3. **Quality Assessment**:

   - Alignment files in SAM/BAM format (processed with pysam)

   - Coverage statistics and purity metrics (see `2_purity.ipynb` and `purity_minimap.txt`)

   - Coverage percentage files for each SAG (see `fastq_list_Bacillus.txt` and `3_get_aligned_result.ipynb`)

     note: we just use Bacillus subtilis as an example

   - Visualization plots (see `coverage_Bacillus.txt` and `5_coverage_distribution.ipynb`)

### Expected Run Time

- **Small dataset (1-5GB)**: 2-4 hours on a standard desktop computer
- **Medium dataset (5-20GB)**: 6-12 hours
- **Large dataset (20GB+)**: 12+ hours (recommend cluster execution)



## Instructions for Use

### Running on Your Own Data

1. **Environment Setup**:

   ```
   conda activate captain
   ```

2. **Configuration**:

   - Update all directory paths in the shell scripts:
     - `DORADO_PATH`: Path to Dorado executable
     - `KIT_PATH`: Path to Dorado kit
     - `POD5_DIR`: Directory containing POD5 files
     - `WORK_DIR`: Working directory for output files
     - `CODE_DIR`: Directory containing Python scripts

3. **Barcode customization**:

   - Modify `barcode_seq.xlsx` with your barcode sequences
   - Update adapter sequences in `search_fastq.py` if different from default (GTCTCGTGGGCTCGG)

4. **Reference genomes**:

   - Provide reference genomes for purity assessment in step 3
   - Update the reference path in the coverage assessment script

5. **Parameter tuning**:

   - Adjust quality thresholds in NanoFilt command as needed
   - Modify alignment parameters in Minimap2 commands for different applications
   - Adjust Miniasm assembly parameters based on read characteristics

### Key Python Scripts and Functions

- **`search_fastq.py`**: Uses Bio.Align.PairwiseAligner for adapter sequence identification
- **`search_barcode.py`**: Processes barcode matching with pandas for data handling
- **`statistical_bbmap.py`**: Generates quality metrics and visualizations using matplotlib and numpy



## Troubleshooting

### Common Issues

1. **Environment activation fails**:

   ```
   conda deactivate
   conda activate captain
   ```

2. **Missing Python packages**:

   ```
   conda env update -f environment.yaml
   ```

3. **Permission denied for scripts**:

   ```
   chmod +x *.sh
   ```

## Citation

If you use this pipeline in your research, please cite the corresponding publication.

Large-scale single-cell long-read genomics enables high-resolution microbiome profiling

Mengzhen Li, Shiyan Li, Hengxin Liu, Xuanpei Zhai, Jie Li, Yanan Du, Xiaolu Li, Jiawei Tong, Jian Zhang, Rong Zhang, Zhangyue Song, Quanjiang Ji, Yuan Luo, Ting Zhang, Wu Wei, Yifan Liu

bioRxiv 2024.09.10.612220; doi: https://doi.org/10.1101/2024.09.10.612220



## Support


For technical support or to report issues, please open an issue on the GitHub repository.










