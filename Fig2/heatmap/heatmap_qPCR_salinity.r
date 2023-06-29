library(pheatmap)
library(RColorBrewer)
library(stats)
library(viridis)
library(stats)
library(dendsort)

# when importing raw name - use first column
metadata <- read.table(file = "/home/oded/Documents/post/articles/springs/qPCR/heatmap_qPCR/heatmap/statistics/complex.heatmap/heatmap-metadata.tsv", sep = "\t", header = T, row.names = 1)

str(metadata)

# transform dataframe to data.matrix
metadata_matrix <- data.matrix(metadata)

# for clustering all variables must be as numeric
metadata$Sep17<-as.numeric(metadata$Sep17)
metadata$Jan18<-as.numeric(metadata$Jan18)
metadata$Jun18<-as.numeric(metadata$Jun18)
metadata$Oct18<-as.numeric(metadata$Oct18)
metadata$Jan19<-as.numeric(metadata$Jan19)
str(metadata)


#omit NA's
#metadata <- metadata[complete.cases(metadata), ]

# first we calculate number of clusters.

# k-means clustering is very popular. However, it has some limitations: it requires the user to specify the number of clusters in advance
# and selects initial centroids randomly. The final k-means clustering solution is very
# sensitive to this initial random selection of cluster centers.

# see /home/oded/Documents/post/articles/springs/qPCR/heatmap_qPCR/complexheatmap/statistics/complex.heatmap
# page 39-40 for number of clusters to choose

# first we compute the number of clusters and plot
library(factoextra)
optimal_kmeans_cluster <- fviz_nbclust(metadata, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16, face  = 'bold'))

# now we see which samples were clustered together and some stats
set.seed(123)
km.res <- kmeans(metadata, 4, nstart = 25)
print(km.res)

# and we can visualize the clusters in a Principal Component Analysis (PCA) plot
fviz_cluster(km.res, data = metadata,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal())

# now we add  the point classifications of each sample to the original data 
metadata_matrix_cluster <- cbind(metadata_matrix, cluster = km.res$cluster)
head(metadata_matrix_cluster)
#now to see which culomn is the one of the clusters (the last one - this is the cbind function)
dim(metadata_matrix_cluster)

# now we will order the the last column
clust_o<- order(metadata_matrix_cluster[,6]) # order the last column

# order the matrix according to the order of the last column
metadata_matrix_cluster<- metadata_matrix_cluster[clust_o,]


# it is better to do Hierarchical K-Means Clustering (page 163) - this is done in pheatmap as well- 
res.hk <-hkmeans(metadata, 4)
names(res.hk)
fviz_dend(res.hk, cex = 0.6, palette = "jco",
          rect = TRUE, rect_border = "jco", rect_fill = TRUE)

fviz_cluster(res.hk, palette = "jco", repel = TRUE,
             ggtheme = theme_classic())

res.hk


             




# to add the clusters to the pheatmap

qPCRheatmap_all <- pheatmap(metadata_matrix_cluster[, 1:5], scale = "column", show_colnames = T, ,show_rownames = T,  cluster_cols = F, cluster_rows = T, cutree_rows = 4, color= inferno(256), cellwidth  = 50, cellheight = 6.8)

clusters <- cutree(qPCRheatmap_all$tree_row, k=4)[qPCRheatmap_all$tree_row[["order"]]]
annot_row <- data.frame(row.names = names(clusters),
                             cluster = as.factor(clusters))

qPCRheatmap_all2 <- pheatmap(metadata_matrix, scale = "column", show_colnames = T, ,show_rownames = T,  cluster_cols = F, cluster_rows = T, annotation_row = annot_row, cutree_rows = 4, color= inferno(256), cellwidth  = 50, cellheight = 8.5)

