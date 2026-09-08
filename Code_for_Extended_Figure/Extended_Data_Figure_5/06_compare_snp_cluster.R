#############################################################################
##--------Process snp calling data for SAG using Snippy----------------------
library(data.table)
data1 <- fread("snp_sag_snippy.tsv",sep = "\n",header = F) 
data = as.data.frame(data1)
index = which(grepl("==",data1$V1) == TRUE)
index2 = c((index[2:length(index)]-1),nrow(data1))
snp = data.frame(ref = "NC_006347.1",pos = 1:5277274,count = 0,REF = NA,A = 0,T = 0,C = 0,G = 0,TYPE = "snp",ALT = NA)
for(i in 1:length(index)){
  sample = strsplit(data[index[i],1],split = "/")[[1]][3]
  tmp = as.data.frame(data[((index[i]+1):index2[i]),])
  tmp = tmp[which(grepl("snp",tmp[,1]) == TRUE),,drop = F]
  ind = c()
  ref = c()
  ind_a = c()
  ind_t = c()
  ind_c = c()
  ind_g = c()
  for(j in 1:nrow(tmp)){
    tt = strsplit(tmp[j,1],split = "\t")[[1]]
    if(tt[1] == "NC_006347.1"){
      ind = c(ind,as.numeric(tt[2]))
      ref = c(ref,tt[4])
      if(tt[5] == "A"){ind_a = c(ind_a,as.numeric(tt[2]))}
      if(tt[5] == "T"){ind_t = c(ind_t,as.numeric(tt[2]))}
      if(tt[5] == "C"){ind_c = c(ind_c,as.numeric(tt[2]))}
      if(tt[5] == "G"){ind_g = c(ind_g,as.numeric(tt[2]))}
    }
  }
  snp$count[ind] = snp$count[ind] + 1
  snp$REF[ind] = ref
  snp$A[ind_a] = snp$A[ind_a] + 1
  snp$T[ind_t] = snp$T[ind_t] + 1
  snp$C[ind_c] = snp$C[ind_c] + 1
  snp$G[ind_g] = snp$G[ind_g] + 1
}
snp = snp[which(snp$REF != ""),]
for(i in 1:nrow(snp)){
  tmp = snp[i,5:8]
  snp$ALT[i] = names(which.max(tmp))
}
write.table(snp,file = "SNP_CAP-seq_snippy",sep = "\t",quote = F,row.names = F)

#############################################################################
##----------------Comparison for different methods---------------------------
#----------------------------------------------------------------------------
library(tidydr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggvenn)
inStrain = readxl::read_xlsx("SNP_metagenomics_inStrain.xlsx")
inStrain = as.data.frame(inStrain)
Single_cell = fread("SNP_CAP-seq_snippy.tsv")
Single_cell = as.data.frame(Single_cell)
Single_cell = Single_cell[which(Single_cell$count >= 5),]
sc_instrain = readxl::read_xlsx("SNP_CAP-seq-pooled_inStrain.xlsx")
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

#############################################################################
##-------------------------hierarchical clustering---------------------------
library(data.table)
data <- fread("snp_sag_snippy.tsv",sep = "\n",header = F) 
data = as.data.frame(data)
index = which(grepl("==",data$V1) == TRUE)
index2 = c((index[2:length(index)]-1),nrow(data))
pos = fread("SNP_CAP-seq_snippy.tsv") 
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
save(result,file = "snp.RData")
#-------------------------------------------------------------------
library(ape)      # for as.phylo
library(ggtree)   # for tree visualisation (optional)
library(tidyverse)
load("snp.RData")
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
