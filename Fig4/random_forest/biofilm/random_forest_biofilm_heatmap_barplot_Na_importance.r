  # from QIIME2 to r - see qiime2_r_tutorial_script (latest date...)

library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(phyloseq)
#make heatmap of relative abundance for each of the 20 features in the spring clusters
# the files are 20 top features importance (spring-cluster or Na) and feature relative abundance 

Na_feature_table <- read.table(file = "/home/oded/Documents/post/articles/springs/random_forest/biofilm20/last5/feature-table.tsv", sep = "\t", row.names = 1, header = T, 
                               check.names = FALSE,  skip = 1, comment.char = "")
imp.50 <- read.table(file = "/home/oded/Documents/post/articles/springs/random_forest/biofilm20/last5/heatmap/imp.50.csv", sep = ",", header = T, 
                     check.names = FALSE, comment.char = "",  row.names = 1)
imp.40 <- read.table(file = "/home/oded/Documents/post/articles/springs/random_forest/biofilm20/last5/heatmap/imp.40.csv", sep = ",", header = T, 
                     check.names = FALSE, comment.char = "",  row.names = 1)
imp.30 <- read.table(file = "/home/oded/Documents/post/articles/springs/random_forest/biofilm20/last5/heatmap/imp.30.csv", sep = ",", header = T, 
                     check.names = FALSE, comment.char = "",  row.names = 1)
imp.20 <- read.table(file = "/home/oded/Documents/post/articles/springs/random_forest/biofilm20/last5/heatmap/imp.20.csv", sep = ",", header = T, 
                     check.names = FALSE, comment.char = "",  row.names = 1)

Na_feature_table50 <- subset(Na_feature_table, rownames(Na_feature_table) %in% rownames(imp.50))
Na_feature_table40 <- subset(Na_feature_table, rownames(Na_feature_table) %in% rownames(imp.40))
Na_feature_table30 <- subset(Na_feature_table, rownames(Na_feature_table) %in% rownames(imp.30))
Na_feature_table20 <- subset(Na_feature_table, rownames(Na_feature_table) %in% rownames(imp.20))

# order Na_feature_table by Na_feature_importance
Na_feature_table50 <- Na_feature_table50[match(rownames(imp.50), rownames(Na_feature_table50)), ]
Na_feature_table40 <- Na_feature_table40[match(rownames(imp.40), rownames(Na_feature_table40)), ]
Na_feature_table30 <- Na_feature_table30[match(rownames(imp.30), rownames(Na_feature_table30)), ]
Na_feature_table20 <- Na_feature_table20[match(rownames(imp.20), rownames(Na_feature_table20)), ]


#new df with mean feature relative abundance per spring-clustures (choose hohw many features)
feature_table_Fuliya <- Na_feature_table30 %>%
  select(contains('Fuliya')) %>% # 'select' function to choose the columns you want 
  mutate(Fuliya = rowMeans(.))

feature_table_Tiberias <- Na_feature_table30 %>%
  select(contains('Tiberias')) %>% # 'select' function to choose the columns you want 
  mutate(Tiberias_Hot_Springs = rowMeans(.))

feature_table_Tabgha <- Na_feature_table30 %>%
  select(contains(c("Barbutim", "Maayan", "Gesher_hana", "Ein7", "Nur"))) %>% # 'select' function to choose the columns you want 
  mutate(Tabgha = rowMeans(.))

#feature_table_reference_becteria <- feature_table %>%
# select(contains(c("reference"))) %>% # 'select' function to choose the columns you want 
#mutate(Reference_Bacteria = rowSums(.))


feature_table_cbind <- cbind(feature_table_Tiberias, feature_table_Fuliya, feature_table_Tabgha)
str(feature_table_cbind)

feature_table_cluster <- feature_table_cbind %>%
  select("Tiberias_Hot_Springs", "Tabgha", "Fuliya")
str(feature_table_cluster)


OTU = otu_table(as.matrix(feature_table_cluster), taxa_are_rows = TRUE)

# merge the data (this is the phyloseq object)
ps1 <- phyloseq(OTU)

# heatmap using complexheatmap

matrix <- as.matrix(data.frame(ps1))

#order 

importance30 = imp.30$IncMSE

column_anno = function(index) {
  n = length(index)
  # since middle of columns are in 1, 2, ..., n and each column has width 1
  # then the most left should be 1 - 0.5 and the most right should be n + 0.5
  pushViewport(viewport(xscale = c(0.5, n + 0.5), yscale = range(importance)))
  # since order of columns will be adjusted by clustering, here we also 
  # need to change the order by `[index]`
  grid.points(index, importance[index], pch = 16, default.unit = "native")
  # this is very important in order not to mess up the layout
  upViewport() 
}

barplot.imp = HeatmapAnnotation("%IncMSE"= row_anno_barplot(importance30, axis = TRUE), which = "row", width = unit(10, "cm"), 
                                rowname = anno_text(rownames(imp.30), just = "left", location = unit(0.05, 'npc')), 
                                show_annotation_name = TRUE)
draw(barplot.imp)

# in the code below, with cluster_cols = FALSE - the order is that specified in the feature table (.tsv) supplied. 
# So this may need to be edited

#are there Na's in the dataframe
sum(is.na(matrix))

heatmap.Na.imp<- ComplexHeatmap::pheatmap(matrix, scale= "row", cluster_rows = FALSE,
                                          treeheight_row = 100,
                                          col=viridis(365), cluster_cols = FALSE, cellheight = 10,
                                          cellwidth = 20,
                                                show_rownames = FALSE,
                                                right_annotation = barplot.imp,
                                                fontsize_row = 5)
draw(heatmap.Na.imp, heatmap_legend_side = "left")

pdf("/home/oded/Documents/post/articles/springs/random_forest/biofilm20/last5/heatmap/to_report/RanFor_biofilm20_Na_importance.pdf",  width = 15, height = 10, family="ArialMT")
draw(heatmap.Na.imp, heatmap_legend_side = "left")
dev.off()

png("/home/oded/Documents/post/articles/springs/random_forest/biofilm20/last5/heatmap/to_report/RanFor_biofilm20_Na_importance.png",  width = 14800, height = 9000, res = 1200)


heatmap.Na.imp.clust<- ComplexHeatmap::pheatmap(matrix, scale= "row", cluster_rows = TRUE, clustering_method = "average",
                         treeheight_row = 100,
                         col=viridis(365), cluster_cols = FALSE, cellheight = 10,
                         cellwidth = 20,
                         show_rownames = FALSE,
                         cutree_rows = 3,
                         right_annotation = barplot.imp,
                         fontsize_row = 5)

 

  # %IncMSE = percent Increase in Mean Squared Error (Variable Importance)"
