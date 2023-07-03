library(tidyverse)
library(ggtree)
library(ggtreeExtra)
# BiocManager::install("ggtreeExtra")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(ape)
library(tidytree)
library(treeio)


#import the unrooted tree from mega (in newick )
unrooted.tree <- read.tree("/home/oded/Documents/post/articles/springs/tree/for_tree/percent/biofilm/0.5/last_5/tree_ML_springs_legionella_biofilm_0.5reads_HKY_short.nwk")
unrooted.tree

#total number of nodes
getNodeNum(unrooted.tree) #or
Nnode2(unrooted.tree)

# now just to have a better look at the tree (nodes and branches). 
#We can write it as a tidy data frame (a tbl_tree object)
unrooted.tree.x <- as_tibble(unrooted.tree)
unrooted.tree.x

# Adding tree information here: http://yulab-smu.top/treedata-book/chapter2.html#accesor-tidytree

# To add bootstrap values to the tree, I want to group them
# (i.e. below 50% and above 50%). This cannot be done easily, 
# the bootstrap is inter node data and is located in the lebel column of 
# the tree (see unrooted.tree.x for example). It is with the tip label 
# (i.e. the names of branches) - thus is not numerical.
# So we need to convert the tree to a tree data frame, do the grouping
# and convert back to a tree format. Then root and if we want, convert to phyloseq object.

#first we extract the node number and label (of bootstrap values only)
bootsv <-data.frame(unrooted.tree[["node.label"]])
ntip <- length(unrooted.tree[["tip.label"]])
nnode <- unrooted.tree[["Nnode"]]
Bootstrap <- data.frame(lapply(bootsv, function(x) as.numeric(sub("/.*", "",x))))
Bootstrap$node <- c((ntip + 1) : (ntip + nnode))
colnames(Bootstrap) <- c("label","node")
Bootstrap

# now we can make groups as we like, using cut
Bootstrap$label <- cut(as.numeric(Bootstrap$label), breaks = c(0,0.5,1),right = T,include.lowest = TRUE)
#Bootstrap$label <- cut(as.numeric(Bootstrap$label), breaks = c(0,0.5,0.75,0.9,1),right = T,include.lowest = TRUE)

# rename levels by the grouping
levels(Bootstrap$label) <- c("0<=Bootstrap<0.5","0.5<=Bootstrap<=1")
#levels(Bootstrap$label) <- c("0<=Bootstrap<0.5","0.5<=Bootstrap<0.75", "0.75<=Bootstrap<0.9", "0.9<=Bootstrap<=1")
# remove null levels if present
Bootstrap$label <- Bootstrap$label[, drop=T]

#Put the results back into the tree.
unrooted.tree[["node.label"]] <- Bootstrap$label

# now the unrooted tree will have the groups and not the original values

# root the tree to the outgroup (ape package)
rooted.tree <- root(unrooted.tree, outgroup = "C.burnetii_HM208383", resolve.root = TRUE)
rooted.tree

rooted.tree.x <- as_tibble(unrooted.tree)
rooted.tree.x


feature_table <- read.table(file = "/home/oded/Documents/post/articles/springs/tree/for_tree/percent/biofilm/0.5/last_5/feature-table.tsv", sep = ",", header = T, row.names = 1, 
                        skip = 1, comment.char = "")
str(feature_table)
# In ggtree and phyloseq, the trees are presentred in a way that the group is presentred for
# every sample. and as each ASV can appear in many samples, it will be annotated with
# multiple identical group identifires. The way I tried to go around it is to sum the
# stations of each group (spring cluster - in this example) 
# so we have a final ASV table with only groups instead of sample.
# as the metadata file must contain the samples in the otu_table. They must match
# So I added them to the metadata file.

# first I created df with the relevant samples and sum of them (in new column named as the cluster)
# Then bind them and finaly created a new df with only the sumed spring cluster column


feature_table_Fuliya <- feature_table %>%
  select(contains('Fuliya')) %>% # 'select' function to choose the columns you want 
  mutate(Fuliya = rowSums(.))

feature_table_Tiberias <- feature_table %>%
  select(contains('Tiberias')) %>% # 'select' function to choose the columns you want 
  mutate(Tiberias_Hot_Springs = rowSums(.))

feature_table_Tabgha <- feature_table %>%
  select(contains(c("Barbutim", "Maayan", "Gesher_hana", "Ein7", "Nur"))) %>% # 'select' function to choose the columns you want 
  mutate(Tabgha = rowSums(.))

#feature_table_reference_becteria <- feature_table %>%
# select(contains(c("reference"))) %>% # 'select' function to choose the columns you want 
#mutate(Reference_Bacteria = rowSums(.))

feature_table_outgroup <- feature_table %>%
  select(contains(c("Outgroup"))) %>% # 'select' function to choose the columns you want 
  mutate(Outgroup = rowSums(.))


feature_table_cbind <- cbind(feature_table_Tiberias, feature_table_Fuliya, feature_table_Tabgha, feature_table_outgroup)
str(feature_table_cbind)

feature_table_cluster <- feature_table_cbind %>%
  select("Tiberias_Hot_Springs", "Tabgha", "Fuliya", "Outgroup")
str(feature_table_cluster)

metadata <- read.table(file = "/home/oded/Documents/post/articles/springs/tree/for_tree/percent/biofilm/0.5/last_5/sample-metadata.tsv", sep = "\t", header = T, row.names = 1)

feature_table.tree = otu_table(as.matrix(feature_table_cluster), taxa_are_rows = TRUE)
taxa_names(feature_table.tree)
metadata.tree <- sample_data(metadata)

phyloseq.unrooted.tree <- phyloseq(feature_table.tree, metadata.tree, unrooted.tree)
phyloseq.rooted.tree <- phyloseq(feature_table.tree, metadata.tree, rooted.tree)

#match the colors to rest of article: #CC79A7, #009E73, #56B4E9, #D55E00
# Outgroup gray

phyloseq.rooted.tree.circular <- ggtree(phyloseq.rooted.tree, layout="circular") +  geom_tippoint(aes(x=x+1.1*hjust, color=Spring_Cluster), size = 2, position="identity") +
  scale_color_manual(values = c("#D55E00","gray","#56B4E9", "#CC79A7"), na.translate = F) + xlim(-0.00005, 0.5) +
  geom_tiplab(align=TRUE, linesize=0.25, offset = 0.04)+
  geom_point2(aes(subset = !isTip, shape  = label), size=1, show_guide = FALSE)+
  scale_shape_manual("Bootstrap Values", #legend name
                     values = c(32, 16),
                     na.translate = F,#remove na
                     labels = levels(Bootstrap$label))+
  geom_treescale(x = 0.19, y = 0, fontsize=5, linesize=2) + 
  theme(legend.position=c(0.9, 0.9), legend.title=element_text(size=10),
        legend.text=element_text(size=8)) +   
  hexpand(.3) +theme(plot.margin = unit(c(0,0,0,0), "mm"))


  phyloseq.rooted.tree.circular 
  

pdf("/home/oded/Documents/post/articles/springs/tree/for_tree/percent/biofilm/0.5/last_5/tree_ML_Leg_springs_biofilm_0.5reads.pdf", width = 15, height = 15, family="ArialMT")
phyloseq.rooted.tree.circular
dev.off()