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
    mutate(method = "CAP-seq in house pipeline")
  
  sag_level_long <- sag_level %>% 
    rownames_to_column("level") %>%
    pivot_longer(-level, names_to = "sample", values_to = "abundance") %>%
    mutate(method = "CAP-seq SAG MetaPhlan4")
  
  uncl_meta_long <- metaphlan_uncl %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "sample", values_to = "abundance") %>%
    mutate(level = "Unclassified", method = "Metagenome") %>%
    select(-taxon)
  
  uncl_skani_long <- skani_uncl %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "sample", values_to = "abundance") %>%
    mutate(level = "Unclassified", method = "CAP-seq in house pipeline") %>%
    select(-taxon)
  
  uncl_sag_long <- sag_uncl %>%
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon, names_to = "sample", values_to = "abundance") %>%
    mutate(level = "Unclassified", method = "CAP-seq SAG MetaPhlan4") %>%
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
  
  combined_summary$method = factor(combined_summary$method,levels = c("Metagenome","CAP-seq in house pipeline","CAP-seq SAG MetaPhlan4"))
  stat <- combined_level %>%
    group_by( level)%>%
    summarise(abundance = sum(abundance), .groups = "drop")
  stat = stat[order(-stat$abundance),]
  index = which(stat$level %in% c("Others","Unclassified"))
  valid = stat$level[c(setdiff(1:nrow(stat),index),index[1],index[2])]
  
  npg_base <- pal_npg()(10)           # NPG base 10 colours
  n_levels <- length(unique(combined_summary$level))
  
  extended <- colorRampPalette(npg_base)(n_levels)
  set.seed(1)
  color_mapping <- setNames(extended[sample(x = n_levels)], valid)
  color_mapping[n_levels] = "darkgrey"
  fill_scale <- scale_fill_manual(values = color_mapping)
  
  combined_summary$level <- factor(combined_summary$level, levels = valid)
  
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
#-------------------------------------------------------------------------------
##diversity
library(vegan)

df <- readxl::read_xlsx("./merge_dist_Species.xlsx")  
df = as.data.frame(df);rownames(df) = df[,1]
df = df[,-1]
df_t <- t(df)

shannon <- diversity(df_t, index = "shannon")  

simpson <- diversity(df_t, index = "simpson")  

inv_simpson <- diversity(df_t, index = "invsimpson")

result <- data.frame(
  Sample = names(shannon),
  Shannon = shannon,
  Simpson_1_minus_D = simpson,
  InvSimpson = inv_simpson
)
print(result)
