# from QIIME2 to r - see qiime2_r_tutorial_script (latest date...)

library(tidyverse)
library(vegan)
library(BiocManager)
library(phyloseq)
library(ANCOMBC)
library(DESeq2)
library(ComplexHeatmap)
library(viridis)
library(pheatmap)
otu <- read.table(file = "/home/oded/Documents/post/articles/springs/NGS/heatmap/Na/filter2000/filter_combined/qiime2_convert_tsv/reordered/biofilm/last4/feature-table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 1, comment.char = "")
# clean the taxonomy, Greengenes format (I need to see it is the same with SILVA format). 
#if the format is different, I think it will just be in the taxonomy.tsv - so k__ will be different...)



metadata <- read.table(file = "/home/oded/Documents/post/articles/springs/NGS/heatmap/Na/filter2000/filter_combined/qiime2_convert_tsv/reordered/biofilm/sample-metadata.tsv", sep = "\t", header = T, row.names = 1)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata)
# merge the data (this is the phyloseq object)
ps1 <- phyloseq(OTU, SAMPLE)

# for now I omit all the rerafaction, rarefying, alpha, beta...

# I will continue from: "Differential abundance analysis"

sample_data(ps1)$origin_new <- as.factor(sample_data(ps1)$origin_new) # factorize for DESeq2

# Let’s see what features are enriched differentially between tongue and gut samples.
# now DESeq2, make pseudocount + 1, to omit 0.


ps1.pse <- ps1
otu_table(ps1.pse) <- otu_table(ps1) + 1

# heatmap using complexheatmap

matrix <- as.matrix(data.frame(otu_table(ps1.pse)))
metadata <- data.frame(sample_data(ps1.pse))



# Define the annotation color for columns and rows
annotation_var = data.frame(
  `Spring cluster` = as.factor(metadata$origin_new),
  Na = as.numeric(metadata$Na),
  check.names = FALSE
)

# ann_color should be named vectors
ann_colors = list(
  `Spring cluster` = c(Tiberias_Hot_Springs = "purple", Tabgha = "cyan", Fuliya = "red"),
  Na = c("white", "blue4", "green", "deeppink2"))

# in the code below, with cluster_cols = FALSE - the order is that specified in the feature table (.tsv) supplied. 
# So this may need to be edited

ComplexHeatmap::pheatmap(matrix, scale= "row", cluster_rows = TRUE, clustering_method = "average",
                           treeheight_row = 100,
                           col=viridis(365), cluster_cols = FALSE, cellheight = 2.4,
                           annotation_col = annotation_var,
                           border_color = NA, cutree_rows = 2,
                           annotation_colors = ann_colors, fontsize_row = 5)

# t() - transposes the matrix. So we have to change the lol to wors... 
  p.heatmap.biofilm <- ComplexHeatmap::pheatmap(mat = t(matrix), scale = "column", cluster_col = TRUE, clustering_method = "average",
                         treeheight_col = 100, cluster_row = FALSE, col=viridis(365), border_color = NA,
                         annotation_colors = ann_colors, annotation_row = annotation_var, fontsize_col = 5, cutree_cols = 8)