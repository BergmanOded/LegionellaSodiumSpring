# from QIIME2 to r - see qiime2_r_tutorial_script (latest date...)

library(tidyverse)
library(vegan)
library(devtools)
library(ggvenn)
library(VennDiagram)
#first reorder the feature table columns to match that of the matadata (Na order)

otu <- read.table(file = "/home/oded/Documents/post/articles/springs/18S/NGS/venn_diagram/water/feature-table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 1, comment.char = "")

otu_Fuliya <- otu %>%
  select(contains('Fuliya')) %>% # 'select' function to choose the columns you want 
  mutate(Fuliya = rowSums(.))

otu_Tiberias <- otu %>%
  select(contains('Tiberias')) %>% # 'select' function to choose the columns you want 
  mutate('Tiberias Hot Springs' = rowSums(.))

otu_Haon <- otu %>%
  select(contains('Haon')) %>% # 'select' function to choose the columns you want 
  mutate('Haon-Borehole' = rowSums(.))

otu_Tabgha <- otu %>%
  select(contains(c("Barbutim", "Maayan", "Gesher_hana", "Ein7", "Nur"))) %>% # 'select' function to choose the columns you want 
  mutate(Tabgha = rowSums(.))

otu_cbind <- cbind(otu_Tiberias, otu_Fuliya, otu_Haon, otu_Tabgha)
str(otu_cbind)

otu_cluster <- otu_cbind %>%
  select("Haon-Borehole", "Tiberias Hot Springs", "Tabgha", "Fuliya")
str(otu_cluster)

write.csv(otu_cluster,"/home/oded/Documents/post/articles/springs/18S/NGS/venn_diagram/water/otu_cluster.csv", row.names = TRUE)

venn_tiberias <- rownames(otu_cluster)[otu_cluster[,"Tiberias Hot Springs"] > 0]
venn_haon <- rownames(otu_cluster)[otu_cluster[,"Haon-Borehole"] > 0]
venn_tabgha <- rownames(otu_cluster)[otu_cluster[,"Tabgha"] > 0]
venn_fuliya <- rownames(otu_cluster)[otu_cluster[,"Fuliya"] > 0]

venn_springs_water <- venn.diagram(
  x = list('Tiberias Hot Springs'=venn_tiberias, 'Haon-Borehole'=venn_haon, Tabgha=venn_tabgha, Fuliya=venn_fuliya),
  fill = c("#CC79A7", "#009E73", "#56B4E9", "#D55E00"),
  alpha = c(0.7, 0.7, 0.7, 0.7), cex= 1.3, cat.cex = 1.5, filename=NULL, print.mode = c("raw", "percent"),
height = 4000, width = 8000, cat.just=list(c(0.37,0) , c(1,0) , c(1,0) , c(0.3,0))
)
grid.draw(venn_springs_water)

pdf("/home/oded/Documents/post/articles/springs/18S/NGS/venn_diagram/water/venn_water_springs20.pdf",  width = 15, height = 10, family="ArialMT")
grid.draw(venn_springs_water)
dev.off()

png("/home/oded/Documents/post/articles/springs/18S/NGS/venn_diagram/water/venn_water_springs20.png",  width = 14800, height = 9000, res = 1200)