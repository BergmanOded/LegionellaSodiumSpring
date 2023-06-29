library(tidyverse)
library(vegan)
library(devtools)
library(ggvenn)
library(VennDiagram)
#first reorder the feature table columns to match that of the matadata (Na order)

otu <- read.table(file = "/home/oded/Documents/post/articles/springs/feature_count/20/biofilm/last5/feature-table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 1, comment.char = "")

otu_Fuliya <- otu %>%
  select(contains('Fuliya')) %>% # 'select' function to choose the columns you want 
  mutate(Fuliya = rowSums(.))

otu_Tiberias <- otu %>%
  select(contains('Tiberias')) %>% # 'select' function to choose the columns you want 
  mutate('Tiberias Hot Springs' = rowSums(.))

otu_Tabgha <- otu %>%
  select(contains(c("Barbutim", "Maayan", "Gesher_hana", "Ein7", "Nur"))) %>% # 'select' function to choose the columns you want 
  mutate(Tabgha = rowSums(.))

otu_cbind <- cbind(otu_Tiberias, otu_Fuliya, otu_Tabgha)
str(otu_cbind)

otu_cluster <- otu_cbind %>%
  select("Tiberias Hot Springs", "Tabgha", "Fuliya")
str(otu_cluster)

write.csv(otu_cluster,"/home/oded/Documents/post/articles/springs/NGS/venn_diagram/biofilm/otu_cluster.csv", row.names = TRUE)

venn_tiberias <- rownames(otu_cluster)[otu_cluster[,"Tiberias Hot Springs"] > 0]
venn_tabgha <- rownames(otu_cluster)[otu_cluster[,"Tabgha"] > 0]
venn_fuliya <- rownames(otu_cluster)[otu_cluster[,"Fuliya"] > 0]

venn_springs_biofilm <- venn.diagram(
  x = list('Tiberias Hot Springs'=venn_tiberias, Tabgha=venn_tabgha, Fuliya=venn_fuliya),
  fill = c("#CC79A7", "#56B4E9", "#D55E00"),
  alpha = c(0.7, 0.7, 0.7), cex= 1.3, cat.cex = 1.5, filename=NULL, print.mode = c("raw", "percent"),
height = 4000, width = 8000, cat.just=list(c(0.55,-0.5), c(0.6,-28) , c(0.35,-0.5))
)
grid.draw(venn_springs_biofilm)

pdf("/home/oded/Documents/post/articles/springs/final R code/NGS/16S/venn_diagram/venn_biofilm_Leg.pdf",  width = 12, height = 6.5, family="ArialMT")
grid.draw(venn_springs_biofilm)
dev.off()



  