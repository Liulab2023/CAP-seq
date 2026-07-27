

# Comparison of CAP-seq with Shotgun Metagenomics

[![Pipeline](https://img.shields.io/badge/Workflow-MetaPhlAn4--StrainPhlAn--inStrain-blue)](LICENSE)
[![Language](https://img.shields.io/badge/Language-Shell--R-yellow)](README.md)

This repository contains a complete bioinformatics pipeline for comparing **CAP-seq (single-cell/amplicon)** data against **shotgun metagenomics** data. It covers:

- **Taxonomic profiling** using MetaPhlAn4 and skani.
- **Strain-level variant analysis** using StrainPhlAn4, inStrain, and Snippy.
- **Comprehensive visualization** in R, including Alluvial/Sankey diagrams, UpSet plots, Venn diagrams, and hierarchical clustering of SNP profiles.

The primary goal is to systematically evaluate the concordance and discrepancies between single-cell-derived taxonomic/SNP calls and those obtained from bulk shotgun metagenomes.

---

## Table of Contents

- [Dependencies & Installation](#dependencies--installation)
- [Pipeline Overview](#pipeline-overview)
- [Usage](#usage)
  - [1. Taxonomy Assignment (Metagenome)](#1-taxonomy-assignment-metagenome)
  - [2. Taxonomy Assignment (Single-cell)](#2-taxonomy-assignment-single-cell)
  - [3. Taxonomy Assignment using skani (Single-cell)](#3-taxonomy-assignment-using-skani-single-cell)
  - [4. Taxonomic Profile Comparison (R)](#4-taxonomic-profile-comparison-r)
  - [5. Strain-level Analysis](#5-strain-level-analysis)
  - [6. SNP Comparison & Clustering (R)](#6-snp-comparison--clustering-r)
- [Outputs](#outputs)
- [Citation](#citation)

---

## Dependencies & Installation

### Required Software Tools
| Tool                                                    | Usage                                |
| :------------------------------------------------------ | :----------------------------------- |
| [fastp](https://github.com/OpenGene/fastp)              | Quality control and trimming         |
| [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/)   | Host read removal (hg38)             |
| [MetaPhlAn 4](https://github.com/biobakery/MetaPhlAn)   | Taxonomic profiling                  |
| [Minimap2](https://github.com/lh3/minimap2)             | Long-read/assembly alignment         |
| [Miniasm](https://github.com/lh3/miniasm)               | Long-read assembly for single cells  |
| [skani](https://github.com/bluenote-1577/skani)         | Fast ANI-based taxonomic search      |
| [StrainPhlAn 4](https://github.com/biobakery/MetaPhlAn) | Strain-level phylogenetic analysis   |
| [inStrain](https://github.com/MrOlm/inStrain)           | Microbial strain profiling from BAMs |
| [Snippy](https://github.com/tseemann/snippy)            | Rapid SNP detection                  |
| [samtools](http://www.htslib.org/)                      | BAM processing                       |

### R Environment
Run the following in R to install required packages:

```r
install.packages(c("tidyverse", "ggvenn", "eulerr", "ggplot2", "reshape2",
                   "pheatmap", "vegan", "ggpubr", "ggsci", "ggalluvial",
                   "ggbreak", "UpSetR", "ape", "ggtree", "data.table", 
                   "tidydr", "readxl"))
```

> **Note:** Please ensure `ggtree` and `tidydr` are installed via Bioconductor if standard installation fails:
> ```r
> BiocManager::install(c("ggtree", "tidydr"))
> ```

---

## Pipeline Overview

1. **Preprocessing**: Adapter trimming and host (human) contamination removal.
2. **Profiling**: 
   - Metagenomes → MetaPhlAn4.
   - Single-cells → MetaPhlAn4  or skani (via miniasm assemblies).
3. **Comparison**: R scripts for taxon overlap, abundance flows, and set intersections.
4. **Strain tracking**:
   - Mapping reads to specific reference genomes.
   - Variant calling via inStrain/Snippy.
   - Phylogenetic clustering based on SNP presence/absence.

---

## Usage

### 1. Taxonomy Assignment (Metagenome)

```
01_taxonomy_metagenome.sh
```

### 2. Taxonomy Assignment (Single-cell)

```shell
02_taxonomy_sc_metaphlan.sh
```

### 3. Taxonomy Assignment using skani (Single-cell)

*Uses miniasm for de novo assembly of long-read single-cell data.*

```shell
03_taxonomy_sc_skani.sh
```

### 4. Taxonomic Profile Comparison (R)

*Compare profiles from Metagenome, MetaPhlAn4 (Single-cell), and skani (Single-cell). Generates:*
- *Alluvial/Sankey diagrams per taxonomic level*
- *UpSet plots & abundance bar charts per sample*

```r
04_compare_taxonomy.R
```

> **Note:** The full R script integrates `extract_level()`, automatic color palette selection (NPG, scico, or HCL), and batch export of PDF plots to `venn/` and `upset/` directories.

### 5. Strain-level Analysis

#### a) StrainPhlAn4 (Phylogenomic markers)

#### b) inStrain (Read mapping & variant profiling)

#### c) Snippy (SNP calling in contigs)

```shell
05_strain_analysis.sh
```

### 6. SNP Comparison & Clustering (R)

*Compares SNP positions/alleles across Metagenome, Single-cell, and SAG datasets.*
- *Generates Venn diagrams and UpSet plots for shared/unique variants.*
- *Performs hierarchical clustering (binary Jaccard distance) and visualizes sample relationships via `ggtree`.*

```r
06_compare_snp_cluster.R
```

---

## Outputs

All analysis outputs are saved in organized directories:

| Directory            | Contents                                                     |
| :------------------- | :----------------------------------------------------------- |
| `venn/`              | Venn diagrams for taxonomic overlaps.                        |
| `bar/`               | Abundance bar charts.                                        |
| `upset/`             | UpSet plots and accompanying abundance bar plots.            |
| `consensus_markers/` | Marker gene JSON files for StrainPhlAn.                      |
| `contig_snp_total/`  | Snippy per-sample SNP call directories.                      |
| Root directory       | Final comparison PDFs: `venn.pdf`, `upset.pdf`, `clustering_tree_*.pdf`. |

