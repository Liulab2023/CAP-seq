

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

```shell
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
```

### 2. Taxonomy Assignment (Single-cell)

```shell
dir="sc_split_bc" 

for fastq_file in "$dir"/*.fq.gz; do
  base_name=$(basename "$fastq_file" .fq.gz)
  echo "Processing $base_name"

  metaphlan "$fastq_file" --nproc 20 --input_type fastq \
    -o ./metaout/"$base_name".txt --mapout ./metaout/"$base_name"_map.txt \
    --db_dir /path/to/metaphlan --ignore_eukaryotes --split_reads \
    --minimap2_exe /path/to/minimap2 --stat_q 0.01 --avoid_disqm
done
```

### 3. Taxonomy Assignment using skani (Single-cell)

*Uses miniasm for de novo assembly of long-read single-cell data.*

```shell
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
     -d /path/to/SGB-database/SGB_skani -o sc_SGB.txt -t 32
```

### 4. Taxonomic Profile Comparison (R)

*Compare profiles from Metagenome, MetaPhlAn4 (Single-cell), and skani (Single-cell). Generates:*
- *Alluvial/Sankey diagrams per taxonomic level*
- *UpSet plots & abundance bar charts per sample*

```r
required_packages <- c("tidyverse", "ggvenn", "eulerr", 
                       "ggplot2", "reshape2", "pheatmap", 
                       "vegan", "ggpubr", "ggsci", "ggalluvial")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}
dir.create("venn")
dir.create("bar")
# Read raw data (assuming identical sample names)
metaphlan <- read.table("metaphlan_merged.tsv", header = TRUE, sep = "\t", 
                        row.names = 1, check.names = FALSE)
skani <- read.table("skani_profile.tsv", header = TRUE, sep = "\t", 
                    row.names = 1, check.names = FALSE)
sag  <- read.table("sc_profile.tsv", header = TRUE, sep = "\t", 
                   row.names = 1, check.names = FALSE)
# ------ Separate Unclassified rows ------
# Identify common variants: UNCLASSIFIED / unclassified / Unclassified etc.
metaphlan_uncl <- metaphlan[grepl("unclassified", rownames(metaphlan), ignore.case = TRUE), , drop = FALSE]
skani_uncl <- skani[grepl("unclassified", rownames(skani), ignore.case = TRUE), , drop = FALSE]
sag_uncl = sag[grepl("unclassified", rownames(sag), ignore.case = TRUE), , drop = FALSE]

# Remove these rows from the main table (subsequent taxonomic level extraction will not include them)
metaphlan <- metaphlan[!grepl("unclassified", rownames(metaphlan), ignore.case = TRUE), , drop = FALSE]
skani <- skani[!grepl("unclassified", rownames(skani), ignore.case = TRUE), , drop = FALSE]
sag <- sag[!grepl("unclassified", rownames(sag), ignore.case = TRUE), , drop = FALSE]

# =============================================
# 2. Extract abundance at a specific taxonomic level
# =============================================
# Format: "k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus|s__Lactobacillus_acidophilus"
# Function: extract abundances for a given level from a taxonomic abundance matrix
# df: matrix/data.frame with full taxonomic paths (pipe-separated) as rows, samples as columns
# level: level number (1~7) or shorthand name
#        1 = "k" (Kingdom), 2 = "p" (Phylum), 3 = "c" (Class), 4 = "o" (Order), 5 = "f" (Family),
#        6 = "g" (Genus), 7 = "s" (Species)
# Incomplete paths (less than 7 levels) will be filled with NA and filtered out.

extract_level <- function(df, level) {
  # Support input like 'k', 'p' etc. by converting to number
  if (is.character(level)) {
    level <- match(tolower(level), c("k","p","c","o","f","g","s"))
    if (is.na(level)) stop("level must be one of k,p,c,o,f,g,s")
  }
  
  df <- as.data.frame(df)
  df$full_taxonomy <- rownames(df)
  
  df_long <- df %>%
    pivot_longer(-full_taxonomy, names_to = "sample", values_to = "abundance") %>%
    filter(abundance > 0) %>%
    # Split taxonomic path by pipe
    mutate(levels = str_split(full_taxonomy, "\\;")) %>%
    # Extract target level name (remove prefix like p__)
    mutate(taxon = map_chr(levels, function(x) {
      if (length(x) >= level) {
        str_replace(x[level], "^[a-z]__", "")
      } else {
        NA_character_
      }
    })) %>%
    filter(!is.na(taxon)) %>%
    group_by(sample, taxon) %>%
    summarise(abundance = sum(abundance), .groups = "drop")
  
  # Convert back to wide format
  df_wide <- df_long %>%
    pivot_wider(names_from = sample, values_from = abundance, values_fill = 0) %>%
    column_to_rownames("taxon")
  return(df_wide)
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)   # install.packages("ggalluvial")
library(grDevices)
if (!require("scico")) install.packages("scico")
library(scico)
install.packages("randtoolbox")
library(qualpalr)
get_nature_palette <- function(n) {
  if (n <= 10) {
    # For small counts, use NPG directly
    require(ggsci)
    return(pal_npg("nrc")(n))
  } else if (n <= 20) {
    # Between 10-20, use scico "roma" or "batlow" (colourblind-friendly, high contrast)
    return(scico::scico(n, palette = "batlow"))
  } else {
    # Over 20, use hcl.colors to generate evenly spaced palette (more scientific)
    # Options: "Dynamic", "Harmonic", "Tropic" – all suitable for journals
    return(hcl.colors(n, palette = "Dynamic", rev = FALSE))
  }
}

levels <- c("k", "p", "c", "o", "f", "g", "s")
level_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

for(i in 1:length(levels)){
  metaphlan_level <- extract_level(metaphlan, levels[i])
  index = which(rownames(metaphlan_level) == "")
  rownames(metaphlan_level)[index] = "Others"
  
  skani_level <- extract_level(skani, levels[i])
  
  sag_level <- extract_level(sag, levels[i])
  index = which(rownames(sag_level) == "")
  rownames(sag_level)[index] = "Others"
  
  # Combine data from two sources and label method
  metaphlan_level_long <- metaphlan_level %>% 
    rownames_to_column("level") %>%
    pivot_longer(-level, names_to = "sample", values_to = "abundance") %>%
    mutate(method = "Metagenome")
  
  skani_level_long <- skani_level %>% 
    rownames_to_column("level") %>%
    pivot_longer(-level, names_to = "sample", values_to = "abundance") %>%
    mutate(method = "Single cell")
  
  sag_level_long <- sag_level %>% 
    rownames_to_column("level") %>%
    pivot_longer(-level, names_to = "sample", values_to = "abundance") %>%
    mutate(method = "SAG")
  
  uncl_meta_long <- metaphlan_uncl %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "sample", values_to = "abundance") %>%
    mutate(level = "Unclassified", method = "Metagenome") %>%
    select(-taxon)
  
  uncl_skani_long <- skani_uncl %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "sample", values_to = "abundance") %>%
    mutate(level = "Unclassified", method = "Single cell") %>%
    select(-taxon)
  
  uncl_sag_long <- sag_uncl %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "sample", values_to = "abundance") %>%
    mutate(level = "Unclassified", method = "SAG") %>%
    select(-taxon)
  # Merge with level-specific data
  combined_level <- bind_rows(metaphlan_level_long, skani_level_long,sag_level_long,
                              uncl_meta_long, uncl_skani_long,uncl_sag_long)
  
  # Filter low abundance (<1% total) and label as "Others"
  level_threshold <- combined_level %>% 
    group_by(level) %>% 
    summarise(total = sum(abundance), .groups = "drop") %>%
    mutate(keep = ifelse(total >= 1, level, "Others"))
  
  combined_level <- combined_level %>% 
    left_join(level_threshold, by = "level") %>%
    mutate(level = keep) %>%
    group_by(sample, method, level) %>%
    summarise(abundance = sum(abundance), .groups = "drop")
  
  # Key modification: summarise by method and compute relative abundance (equivalent to position="fill")
  combined_summary <- combined_level %>%
    group_by(method, level) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    group_by(method) %>%
    mutate(prop = abundance / sum(abundance)) %>%
    ungroup()
  combined_summary <- combined_summary %>%
    complete(method, level, fill = list(prop = 0, abundance = 0))
  
  combined_summary$method = factor(combined_summary$method,levels = c("Metagenome","Single cell","SAG"))
  npg_base <- pal_npg()(10)           # NPG base 10 colours
  n_levels <- length(unique(combined_summary$level))
  
  if (n_levels <= 10) {
    # Within 10, use NPG directly
    fill_scale <- scale_fill_npg()
  } else if (n_levels <= 20) {
    # 10-20 use D3 category20 (similar style)
    fill_scale <- scale_fill_manual(values = scico::scico(n_levels, palette = "batlow"))
  } else {
    # 20-50 use IGV or NPG interpolated extension
    extended <- colorRampPalette(npg_base)(n_levels)
    fill_scale <- scale_fill_manual(values = extended)
  } 
  
  # Draw alluvial / Sankey-style plot
  p <- ggplot(combined_summary,
              aes(x = method, y = prop, 
                  alluvium = level, stratum = level,
                  fill = level)) +
    # Draw flow lines connecting the same level across two bars
    geom_flow(stat = "alluvium", 
              lode.guidance = "frontback",
              color = "white", 
              alpha = 0.7, 
              width = 0.25) +
    # Draw strata (bars)
    geom_stratum(stat = "alluvium", 
                 width = 0.4, 
                 color = "white", 
                 size = 0.5) +
    # Optional: add level labels on strata (uncomment if few taxa)
    # geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
    fill_scale +
    theme_minimal() +
    labs(y = "Relative Abundance", 
         x = "", 
         fill = level_names[i]) +
    ggtitle(paste0(level_names[i], " Level")) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # x-axis tick label size
      axis.text.y = element_text(size = 16),                        # y-axis tick label size
      axis.title.x = element_text(size = 14),                       # x-axis title size
      axis.title.y = element_text(size = 14),                      # y-axis title size
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(paste0(level_names[i], "_comparison_sankey.pdf"), 
         plot = p, width = 10, height = 6)
}

library(UpSetR)
if(!dir.exists("./upset")) dir.create("./upset")

levels <- c("k", "p", "c", "o", "f", "g", "s")
level_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

for(i in seq_along(levels)) {
  lev <- levels[i]
  meta_l <- extract_level(metaphlan, lev)
  skani_l <- extract_level(skani, lev)
  sag_l   <- extract_level(sag, lev)
  
  # Remove empty row names
  meta_l <- meta_l[rownames(meta_l) != "", , drop = FALSE]
  skani_l <- skani_l[rownames(skani_l) != "", , drop = FALSE]
  sag_l   <- sag_l[rownames(sag_l) != "", , drop = FALSE]
  
  samples <- intersect(colnames(meta_l), colnames(skani_l))
  samples <- intersect(samples, colnames(sag_l))
  
  for (samp in samples) {
    # Get taxa present (>0) in each set for this sample
    meta_set <- rownames(meta_l)[meta_l[[samp]] > 0]
    skani_set <- rownames(skani_l)[skani_l[[samp]] > 0]
    sag_set   <- rownames(sag_l)[sag_l[[samp]] > 0]
    
    # Union of all taxa that appear
    all_taxa <- unique(c(meta_set, skani_set, sag_set))
    if(length(all_taxa) == 0) next
    
    # Build boolean presence matrix (rows: taxa, columns: three sets)
    presence <- data.frame(
      Metagenome = as.integer(all_taxa %in% meta_set),   # TRUE → 1, FALSE → 0
      Single_cell = as.integer(all_taxa %in% skani_set),
      SAG_metagenome = as.integer(all_taxa %in% sag_set),
      row.names = all_taxa,
      stringsAsFactors = FALSE
    )
    
    # Draw UpSet plot
    # Adjust parameters: bar colours, set name sizes, intersection size fonts, point sizes, etc.
    set_order <- rev(c("Metagenome", "Single_cell", "SAG_metagenome"))
    p <- upset(presence[, set_order, drop = FALSE],
               sets = set_order,
               keep.order = TRUE,
               sets.bar.color = c("#009E73", "#D55E00", "#0072B2"),
               main.bar.color = "black",
               mainbar.y.label = "Intersection Size",
               sets.x.label = "Set Size",
               text.scale = c(2.5, 2.5, 2.2, 2.2, 2.5, 2.2),
               point.size = 2.5,
               line.size = 0.8,
               mb.ratio = c(0.6, 0.4),
               number.angles = 0,
               show.numbers = TRUE,
               nsets = 3,
               nintersects = NA,
               order.by = "freq",
               decreasing = TRUE)
    
    # Save as PDF (width can be adjusted; for 3 sets 6~7 inches is enough)
    outfile <- paste0("./upset/upset_", level_names[i], "_", samp, ".pdf")
    pdf(outfile, width = 12, height = 6)
    print(p)
    dev.off()
    message("Saved: ", outfile)
  }
}

library(UpSetR)
library(ggplot2)
library(dplyr)
library(tidyr)

if(!dir.exists("./upset")) dir.create("./upset")

levels <- c("k", "p", "c", "o", "f", "g", "s")
level_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

sets <- c("Metagenome", "Single_cell", "SAG_metagenome")

set_cols <- c(
  Metagenome = "#0072B2",
  Single_cell = "#D55E00",
  SAG_metagenome = "#009E73"
)

# Safely convert to numeric
to_num <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}

# Extract abundance for given taxa and sample from a matrix
get_abund <- function(mat, taxa, samp) {
  v <- setNames(rep(0, length(taxa)), taxa)
  common <- intersect(taxa, rownames(mat))
  
  if(length(common) > 0) {
    v[common] <- to_num(mat[common, samp, drop = TRUE])
  }
  
  unname(v)
}

for(i in seq_along(levels)) {
  
  lev <- levels[i]
  
  meta_l <- extract_level(metaphlan, lev)
  skani_l <- extract_level(skani, lev)
  sag_l   <- extract_level(sag, lev)
  
  # Remove empty row names
  meta_l  <- meta_l[rownames(meta_l) != "", , drop = FALSE]
  skani_l <- skani_l[rownames(skani_l) != "", , drop = FALSE]
  sag_l   <- sag_l[rownames(sag_l) != "", , drop = FALSE]
  
  samples <- intersect(colnames(meta_l), colnames(skani_l))
  samples <- intersect(samples, colnames(sag_l))
  
  for (samp in samples) {
    
    meta_val  <- to_num(meta_l[, samp, drop = TRUE])
    skani_val <- to_num(skani_l[, samp, drop = TRUE])
    sag_val   <- to_num(sag_l[, samp, drop = TRUE])
    
    # Get taxa present (>0) in each set
    meta_set  <- rownames(meta_l)[meta_val > 0]
    skani_set <- rownames(skani_l)[skani_val > 0]
    sag_set   <- rownames(sag_l)[sag_val > 0]
    
    # Union of all taxa
    all_taxa <- unique(c(meta_set, skani_set, sag_set))
    if(length(all_taxa) == 0) next
    
    # Build 0/1 presence matrix for UpSetR
    presence <- data.frame(
      Metagenome = as.integer(all_taxa %in% meta_set),
      Single_cell = as.integer(all_taxa %in% skani_set),
      SAG_metagenome = as.integer(all_taxa %in% sag_set),
      row.names = all_taxa,
      stringsAsFactors = FALSE
    )
    
    # Build abundance information
    abundance_df <- data.frame(
      Taxon = all_taxa,
      Abund_Metagenome = get_abund(meta_l, all_taxa, samp),
      Abund_Single_cell = get_abund(skani_l, all_taxa, samp),
      Abund_SAG_metagenome = get_abund(sag_l, all_taxa, samp),
      stringsAsFactors = FALSE
    )
    
    # Determine exact intersection for each taxon
    presence_tmp <- presence
    presence_tmp$Taxon <- rownames(presence_tmp)
    
    presence_tmp$Intersection <- apply(
      presence_tmp[, sets, drop = FALSE],
      1,
      function(x) {
        paste(sets[as.logical(as.integer(x))], collapse = "&")
      }
    )
    
    plot_df <- presence_tmp %>%
      left_join(abundance_df, by = "Taxon")
    
    # Summarise counts and total abundance per intersection
    inter_sum <- plot_df %>%
      group_by(Intersection) %>%
      summarise(
        Count = n(),
        Metagenome = sum(Abund_Metagenome, na.rm = TRUE),
        Single_cell = sum(Abund_Single_cell, na.rm = TRUE),
        SAG_metagenome = sum(Abund_SAG_metagenome, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(Count))
    
    if(nrow(inter_sum) == 0) next
    
    # Order for both UpSet and bar plot – keeping consistency
    inter_order <- inter_sum$Intersection
    
    # For UpSetR intersections parameter (list format)
    intersections_list <- lapply(inter_order, function(x) {
      strsplit(x, "&", fixed = TRUE)[[1]]
    })
    
    safe_samp <- gsub("[^A-Za-z0-9_.-]+", "_", samp)
    
    # ============================================================
    # 1. Draw UpSetR plot (optional, shown earlier)
    # ============================================================
    
    # ============================================================
    # 2. Draw separate abundance bar plot
    #    x-axis order follows the same intersection order as UpSet
    # ============================================================
    
    abund_long <- inter_sum %>%
      mutate(
        Intersection = factor(Intersection, levels = inter_order)
      ) %>%
      select(Intersection, all_of(sets)) %>%
      pivot_longer(
        cols = all_of(sets),
        names_to = "Group",
        values_to = "Abundance"
      )
    group_order <- c("Metagenome", "Single_cell", "SAG_metagenome")
    abund_long <- inter_sum %>%
      mutate(
        Intersection = factor(Intersection, levels = inter_order)
      ) %>%
      select(Intersection, all_of(group_order)) %>%
      pivot_longer(
        cols = all_of(group_order),
        names_to = "Group",
        values_to = "Abundance"
      ) %>%
      mutate(
        Group = factor(Group, levels = group_order)
      )
    p_abund <- ggplot(
      abund_long,
      aes(x = Intersection, y = -Abundance, fill = Group)
      ) +
      scale_y_continuous(
        limits = c(-100, 0),
        breaks = seq(0, -100, -25),
        labels = seq(0, 100, 25)
      ) +
      geom_col(
        position = position_dodge(width = 0.8),
        width = 0.7
      ) +
      scale_fill_manual(values = set_cols) +
      theme_bw() +
      theme(
        panel.border = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    outfile_bar <- paste0(
      "./upset/abundance_bar_",
      level_names[i],
      "_",
      safe_samp,
      ".pdf"
    )
    
    ggsave(
      filename = outfile_bar,
      plot = p_abund,
      width = 8,
      height = 4
    )
    
    message("Saved abundance bar plot: ", outfile_bar)
  }
}
```

> **Note:** The full R script integrates `extract_level()`, automatic color palette selection (NPG, scico, or HCL), and batch export of PDF plots to `venn/` and `upset/` directories.

### 5. Strain-level Analysis

#### a) StrainPhlAn4 (Phylogenomic markers)

```shell
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
```

#### b) inStrain (Read mapping & variant profiling)

```shell
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
```

#### c) Snippy (SNP calling in contigs)

```shell
#!/usr/bin/env bash
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
```

### 6. SNP Comparison & Clustering (R)

*Compares SNP positions/alleles across Metagenome, Single-cell, and SAG datasets.*
- *Generates Venn diagrams and UpSet plots for shared/unique variants.*
- *Performs hierarchical clustering (binary Jaccard distance) and visualizes sample relationships via `ggtree`.*

```r
library(tidydr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggvenn)
inStrain = readxl::read_xlsx("merge_snv.xlsx")
inStrain = as.data.frame(inStrain)
Single_cell = fread("snp_all.tsv")
Single_cell = as.data.frame(Single_cell)
Single_cell = Single_cell[which(Single_cell$count >= 5),]
sc_instrain = readxl::read_xlsx("sc_SNVs.xlsx")
sc_instrain = as.data.frame(sc_instrain)
inStrain = inStrain[which(inStrain$scaffold == "NC_006347.1"),]
Single_cell = Single_cell[which(Single_cell$ref == "NC_006347.1"),]
sc_instrain = sc_instrain[which(sc_instrain$scaffold == "NC_006347.1"),]
#--------------------------------------------------------------------------
Meta = paste(inStrain$POS,inStrain$REF,inStrain$ALT,sep = "-")
SAG = paste(sc_instrain$POS,sc_instrain$REF,sc_instrain$ALT,sep = "-")
sc = paste(Single_cell$POS,Single_cell$REF,Single_cell$ALT,sep = "-")

venn_list <- list(Metagenome = Meta, Single_cell = sc,
                  SAG_metagonme = SAG)

# ggvenn plot
p <- ggvenn(venn_list, fill_color = c("#0072B2","#D55E00","#009E73"),
            stroke_size = 0.5, set_name_size = 4)
ggsave(paste0("venn.pdf"), p, width = 6, height = 6)

library(UpSetR)
all_taxa <- unique(c(Meta, SAG, sc))
presence <- data.frame(
  Metagenome = as.integer(all_taxa %in% Meta),   # TRUE → 1, FALSE → 0
  Single_cell = as.integer(all_taxa %in% sc),
  SAG_metagenome = as.integer(all_taxa %in% SAG),
  row.names = all_taxa,
  stringsAsFactors = FALSE
)

# Draw UpSet plot
# Adjust parameters: bar colours, set name sizes, intersection size fonts, point sizes, etc.
set_order <- rev(c("Metagenome", "Single_cell", "SAG_metagenome"))
p <- upset(presence[, set_order, drop = FALSE],
           sets = set_order,
           keep.order = TRUE,
           sets.bar.color = c("#009E73", "#D55E00", "#0072B2"),
           main.bar.color = "black",
           mainbar.y.label = "Intersection Size",
           sets.x.label = "Set Size",
           text.scale = c(2.2, 2.2, 1.9, 1.9, 2.2, 1.9),
           point.size = 2.5,
           line.size = 0.8,
           mb.ratio = c(0.6, 0.4),
           number.angles = 0,
           show.numbers = TRUE,
           nsets = 3,
           order.by = "degree",
           decreasing = TRUE,
           nintersects = NA)

# Save as PDF (width can be adjusted; for 3 sets 6~7 inches is enough)
outfile <- "upset.pdf"
pdf(outfile, width = 12, height = 6)
print(p)
dev.off()


library(dplyr)
library(ggplot2)
library(ggbreak)
library(UpSetR)
set_order <- c("Metagenome", "Single_cell", "SAG_metagenome")
presence_use <- presence[, set_order, drop = FALSE]
presence_use[] <- lapply(presence_use, as.integer)

inter_df <- presence_use %>%
  count(across(all_of(set_order)), name = "n")

# Compute how many sets each intersection involves
inter_df$degree <- rowSums(inter_df[, set_order, drop = FALSE])

# Remove the case where none of the three sets are present
inter_df <- inter_df[inter_df$degree > 0, ]

# Sort by degree then by count (can adjust as needed)
inter_df <- inter_df %>%
  arrange(desc(degree), desc(n))

# Add x-axis indices

# Generate label for each intersection
inter_df$intersection_label <- apply(
  inter_df[, set_order, drop = FALSE],
  1,
  function(x) paste(set_order[as.logical(x)], collapse = "&")
)

inter_df = inter_df[c(1,4,2,3,7,5,6),]
inter_df$xid <- seq_len(nrow(inter_df))

inter_list <- lapply(seq_len(nrow(inter_df)), function(i) {
  set_order[as.logical(inter_df[i, set_order])]
})

inter_list

break_low  <- 30000
break_high <- 300000
gap        <- 1000       # visual gap height for the broken axis
compress   <- 0.08      # compression factor for the high-value part

inter_df$y_plot <- ifelse(
  inter_df$n <= break_low,
  inter_df$n,
  break_low + gap + (inter_df$n - break_high) * compress
)

bottom_breaks <- pretty(c(0, break_low), n = 4)
bottom_breaks <- bottom_breaks[bottom_breaks <= break_low]

top_raw_breaks <- pretty(c(break_high, max(inter_df$n)), n = 4)
top_raw_breaks <- top_raw_breaks[top_raw_breaks >= break_high]

top_plot_breaks <- break_low + gap + 
  (top_raw_breaks - break_high) * compress

y_breaks <- c(bottom_breaks, top_plot_breaks)
y_labels <- c(bottom_breaks, top_raw_breaks)

library(ggplot2)

p_bar_manual <- ggplot(inter_df, aes(x = factor(xid), y = n)) +
  geom_col(width = 0.7, fill = "black") +
  
  # Add actual value labels
  geom_text(
    aes(label = n),
    vjust = -0.3,
    size = 4
  ) +
  scale_y_continuous(limits = c(0, 325000),
                     breaks = c(0,10000,20000,30000))+
  scale_y_break(c(30000, 300000),ticklabels=c(300000,325000),scales = "free",expand=expansion(add = c(0, 0)))+
  scale_x_discrete(labels = inter_df$intersection_label) +
  
  labs(
    x = NULL,
    y = "Intersection Size"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )

ggsave("upset_top_bar.pdf", p_bar_manual, width = 8, height = 6)

##hierarchical clustering
library(data.table)
data1 <- fread("snp.tsv",sep = "\n",header = F) 
data = as.data.frame(data1)
index = which(grepl("==",data1$V1) == TRUE)
index2 = c((index[2:length(index)]-1),nrow(data1))
pos = fread("snp_all.tsv") 
pos = pos[which(pos$count >= 5),]

result = pos[,2,drop = F]
pos_tmp = result
for(i in 1:length(index)){
  tmp = as.data.frame(data[((index[i]+1):index2[i]),])
  tmp = tmp[which(grepl("snp",tmp[,1]) == TRUE),,drop = F]
  tmp = tmp[which(grepl("NC_006347.1",tmp[,1]) == TRUE),,drop = F]
  sag = strsplit(data[index[i],1],split = "/")[[1]][3]
  p = c()
  for(j in 1:nrow(tmp)){
    tt = strsplit(tmp[j,1],split = "\t")[[1]]
    p = c(p,as.numeric(tt[2]))

  }
  tmp_tmp = pos_tmp
  tmp_tmp[[sag]] = 0
  p_ind = which(tmp_tmp$POS %in% p)
  tmp_tmp[p_ind,2] = 1
  result = cbind(result,tmp_tmp[,2,drop = F])
}

result = as.data.frame(result)
rownames(result) = result[,1]
result = result[,-1]
save(result,file = "snp_dist.RData")
#-------------------------------------------------------------------
library(ape)      # for as.phylo
library(ggtree)   # for tree visualisation (optional)
library(tidyverse)
load("snp_dist.RData")
result = as.data.frame(result)
rownames(result) = result[,1]
result = result[,-1]
mat_t <- t(result)
rm(result)
# 2. Compute Jaccard distance (binary distance)
cat("Computing Jaccard distance matrix...\n")
col_dist <- dist(mat_t, method = "binary")   # outputs dist object


# 3. Hierarchical clustering
cat("Performing hierarchical clustering...\n")
hc <- hclust(col_dist, method = "ward.D2")

# 4. Cut tree into 2 clusters
cat("Cutting tree into 2 clusters...\n")
cluster_assignment <- cutree(hc, k = 2)   # returns named vector, values 1 or 2

# 5. Format output
result_df <- data.frame(
  Column_Name = names(cluster_assignment),
  Cluster = cluster_assignment
)

head(result_df)

# Save results to CSV
#write.csv(result_df, "clusters_k2.csv", row.names = FALSE)
cat("Results saved to 'column_clusters_k2.csv'\n")
result_df = read.csv("./clusters_k2.csv")
# 6. (Optional) View number of columns per cluster
table(result_df$Cluster)

# 7. (Optional) Visualise tree with cluster colours
# Convert hclust object to phylo
tree <- as.phylo(hc)

# Build annotation data frame (must have 'label' column)
annot_df <- data.frame(
  label = result_df$Column_Name,   # original column names
  cluster = factor(cluster_assignment)
)

# Draw tree
p <- ggtree(tree, aes(color = cluster), size = 0.3) %<+% annot_df +
  geom_tippoint(aes(color = cluster), size = 0.1) +
  #geom_tiplab(size = 1.5, offset = 0.01, hjust = 0) +
  scale_color_manual(values = c("1" = "#E41A1C", "2" = "#2E78B4"), name = "Cluster") +
  theme(legend.position = "right") +
  labs(title = "Hierarchical Clustering of Columns (k=2)")

print(p)

# Save plot (width/height may need adjustment due to many columns)
ggsave("clustering_tree_k2.pdf", p, width = 4, height = 10, limitsize = TRUE)

annot_df <- data.frame(
  label = result_df$Column_Name,   # original column names
  cluster = factor(result_df$group)
)

my_cols <- c(
  "FMT-2_Before_treatment" = "#B10026",
  "FMT-2_Week2"            = "#E31A1C",
  "FMT-2_Week4"            = "#FC4E2A",
  "VAN-2_Before_treatment" = "#08519C",
  "VAN-2_Week2"            = "#3182BD",
  "VAN-2_Week4"            = "#9ECAE1"
)

# Draw tree
p <- ggtree(tree, aes(color = cluster), size = 0.3) %<+% annot_df +
  geom_tippoint(aes(color = cluster), size = 0.1) +
  #geom_tiplab(size = 1.5, offset = 0.01, hjust = 0) +
  scale_color_manual(values = my_cols, na.value = "grey70")+
  theme(legend.position = "right") +
  labs(title = "Hierarchical Clustering of Columns (k=2)")

print(p)

# Save plot
ggsave("clustering_tree_group.pdf", p, width = 4, height = 10, limitsize = TRUE)
#-------------------------------------------------------------------------------------------
cluster_vec <- setNames(result_df$group, result_df$Column_Name)

# Average profile per cluster
cluster_mat <- sapply(split(names(cluster_vec), cluster_vec), function(cols) {
  rowMeans(result[, cols, drop = FALSE], na.rm = TRUE)
})

# Cluster the 6 clusters again
d6 <- dist(t(cluster_mat))
hc6 <- hclust(d6, method = "complete")
tree6 <- as.phylo(hc6)
my_cols <- c(
  "FMT-2_Before_treatment" = "#E41A1C",
  "FMT-2_Week2" = "#E41A1C",
  "FMT-2_Week4" = "#E41A1C",
  "VAN-2_Before_treatment" = "#2E78B4",
  "VAN-2_Week2" = "#2E78B4",
  "VAN-2_Week4" = "#2E78B4"
)
#tree6$edge.length <- tree6$edge.length * 0.3
p <- ggtree(tree6, branch.length = "none", aes(color = label), size = 1) +
  geom_tippoint(aes(color = label), size = 3, show.legend = FALSE) +
  geom_tiplab(
    aes(color = label),
    size = 4,
    offset = 0.2,
    hjust = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(values = my_cols, na.value = "grey70") +
  labs(title = "Cluster-level tree, k = 6") +
  theme_tree2() +
  theme(
    legend.position = "none",
    plot.margin = margin(5.5, 100, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")
print(p)

# Save plot
ggsave("clustering_tree_6sample.pdf", p, width = 11.34, height = 3.85, limitsize = TRUE)
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

