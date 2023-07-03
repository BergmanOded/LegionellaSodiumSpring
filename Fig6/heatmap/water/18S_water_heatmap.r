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
otu <- read.table(file = "/home/oded/Documents/post/articles/springs/18S/NGS/heatmap/Na/filter500/filter-separate/qiime2_convert_tsv/reordered/water/relative_abundance/last5/feature-table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 1, comment.char = "")
# clean the taxonomy, Greengenes format (I need to see it is the same with SILVA format). 
#if the format is different, I think it will just be in the taxonomy.tsv - so k__ will be different...)



metadata <- read.table(file = "/home/oded/Documents/post/articles/springs/18S/NGS/heatmap/Na/filter500/filter-separate/qiime2_convert_tsv/reordered/water/relative_abundance/sample-metadata.tsv", sep = "\t", header = T, row.names = 1)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata)

# merge the data (this is the phyloseq object)
ps1 <- phyloseq(OTU, SAMPLE)

# for now I omit all the rerafaction, rarefying, alpha, beta...


# heatmap using complexheatmap

matrix <- as.matrix(data.frame(otu_table(ps1)))
metadata <- data.frame(sample_data(ps1))



# Define the annotation color for columns and rows
annotation_var = data.frame(
  `Spring cluster` = as.factor(metadata$origin_new),
  Na = as.numeric(metadata$Na),
  check.names = FALSE
)

# ann_color should be named vectors
ann_colors = list(
  `Spring cluster` = c(Tiberias_Hot_Springs = "purple", Haon = "dark green", Tabgha = "cyan", Fuliya = "red"),
  Na = c("white", "blue4", "green", "deeppink2"))

# t() - transposes the matrix. So we have to change the lol to wors... 
  p.heatmap.water <- ComplexHeatmap::pheatmap(mat = t(matrix), scale = "column", cluster_col = TRUE, clustering_method = "average",
                         treeheight_col = 100, cluster_row = FALSE, col=viridis(365), border_color = NA,
                         annotation_colors = ann_colors, annotation_row = annotation_var, fontsize_col = 5, cutree_cols = 5)

  pdf("/home/oded/Documents/post/articles/springs/18S/NGS/heatmap/Na/filter500/filter-separate/qiime2_convert_tsv/reordered/water/relative_abundance/last5/Legionella_host_water500_Na.pdf",  width = 20, height = 10, family="ArialMT")
  p.heatmap.water
  dev.off()
  
  
  png("/home/oded/Documents/post/articles/springs/18S/NGS/heatmap/Na/filter500/filter-separate/qiime2_convert_tsv/reordered/water/relative_abundance/last4/Legionella_host_water500_Na.png",  width = 14800, height = 9000, res = 1200)
  