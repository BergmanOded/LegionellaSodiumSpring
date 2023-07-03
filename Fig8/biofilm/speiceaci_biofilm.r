# first the general tototial, with the commant to run cross domain analysis:
# https://github.com/zdk123/SpiecEasi#cross-domain-interactions

# next we will use the following links to generate some extra analysis of the networks (with weights)
# https://biovcnet.github.io/_pages/NetworkScience_igraphviz
# https://biovcnet.github.io/_pages/NetworkScience_igraphcluster.html

# this link is a general explanation on networks... 
# https://bookdown.org/omarlizardo/_main/1-intro.html#intro
library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(phyloseq)
library(igraph)
library(ggplot2)
library(qiime2R)
library(viridis)
  # Cross domain interactions
  # the input for spiec.easi will be two OTU (ASV level for us) tables.
  # samples are columns, taxa are rows
  # identical name and order of samples
  # relative abundance is best (count data also possible).
  # and if needed a taxonomy table. Here we dont have as we work at the ASV level.
  # we took the tiberias OTU tables, and filtered them by the  random forest machine learning 
  # top 30 most important features. The resulting table was used to run spiec.easi. 


#our data:
  
otu_table_host <- read.table("/home/oded/Documents/post/articles/springs/speiceaci/biofilm20/THS/hosts_features/last5/feature-table.tsv", sep = "\t", header = T, row.names = 1, skip = 1, comment.char = "")  
otu_table_Leg <- read.table("/home/oded/Documents/post/articles/springs/speiceaci/biofilm20/THS/Legionella_features/last5/feature-table_speiceaci.tsv", sep = "\t", header = T, row.names = 1, skip = 1, comment.char = "")  
metadata <- read.table(file = "/home/oded/Documents/post/articles/springs/speiceaci/biofilm20/THS/sample-metadata.tsv", sep = "\t", header = T, row.names = 1)
taxonomy_Leg <- read.table("/home/oded/Documents/post/articles/springs/speiceaci/taxonomy_Leg.tsv", sep = "\t", header = T, row.names = 1, comment.char = "")  
taxonomy_hohst <- read.table("/home/oded/Documents/post/articles/springs/speiceaci/taxonomy_host.tsv", sep = "\t", header = T, row.names = 1, comment.char = "")  

# order ASVs files in the same column order
col_ord_host <- match(rownames(metadata), colnames(otu_table_host))
col_ord_host
otu_table_host_ordered  <- otu_table_host[ , col_ord_host]

col_ord_Leg <- match(rownames(metadata), colnames(otu_table_Leg))
col_ord_Leg
otu_table_Leg_ordered  <- otu_table_Leg[ , col_ord_Leg]

# convert to phyloseq (not sure in this is obligatiry, but if we have taxonomy we will use phyloseq anyway)

OTU_host = otu_table(as.matrix(otu_table_host_ordered), taxa_are_rows = TRUE)
OTU_Leg = otu_table(as.matrix(otu_table_Leg_ordered), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata)
taxasums_host<- tax_table(as.matrix(taxonomy_hohst))
taxasums_Leg<- tax_table(as.matrix(taxonomy_Leg))



ps_host <- phyloseq(OTU_host, SAMPLE, taxasums_host)
ps_Leg <- phyloseq(OTU_Leg, SAMPLE, taxasums_Leg)

# run the spiec.easi, to get a network
se.Leg.host <- spiec.easi(list(ps_host, ps_Leg), method='glasso', nlambda=40,
                      lambda.min.ratio=3e-3, pulsar.params = list(thresh = 0.1))


# This is the general presentation in the spiec.easi tutotial
dtype <- c(rep(1,ntaxa(ps_Leg)), rep(2,ntaxa(ps_host)))
plot(adj2igraph(getRefit(se.Leg.host)), vertex.color=dtype+1, vertex.attr=list(name=taxa_names(dtype)), vertex.size=9, add.rownames="TRUE", rmEmptyNodes = TRUE) 




# The adjacency matrix asks if two cases share a relationship or not.
# If two actors have a relationship, they are next to each other (adjacent) in the network
# We use the getRefit() function to extract the adjacency matrix from the spieceasi output. This is a square matrix with species 
# on rows and columns and a one if two species are connected (“adjacent”) in the network and a zero otherwise.
adj.mat <- getRefit(se.Leg.host)
table(as.numeric(adj.mat))


# we want to creat weighted networks
se.cor  <- cov2cor(as.matrix(getOptCov(se.Leg.host)))
weighted.adj.mat <- se.cor*getRefit(se.Leg.host)

#Let's take a loot at the weighted adjacency matrix
heatmap(as.matrix(weighted.adj.mat))

# so we want to plot the weighted adjacency matrix 
grph <- adj2igraph(weighted.adj.mat)
# lets look at the  vertices (nodes) 
V(grph)
# and at the edges of the network
E(grph)

# we want the width of the edges to represents the weights 
E(grph)$width <- abs(E(grph)$weight)*50


# and to plot it, as it is
plot(grph,vertex.size=1,
     vertex.label.cex=0.5,      
     vertex.label.color="black",
     layout=layout.circle(grph))

#Remove edges with very low weight and plot again
weight_threshold <- 0.01
grph.weight_threshold <- delete.edges(grph,which(abs(E(grph)$weight)<weight_threshold))

plot(grph.weight_threshold,vertex.size=1,
     vertex.label.cex=0.5,      
     vertex.label.color="black",
     layout=layout.circle(grph))

#Remove negative edges 

grph.pos <- delete.edges(grph,which(E(grph)$weight<0))
plot(grph.pos,
     vertex.label.cex=0.5,      
     vertex.label.color="black",
     edge.color="black",
     layout=layout_with_fr(grph.pos))

#Remove unconnected vertices
grph.pos.con <- delete.vertices(grph.pos,which(degree(grph.pos)<1))
plot(grph.pos.con, vertex.label.cex=0.5, vertex.label.color="black", 
     edge.color="black",
     layout=layout_with_fr(grph.pos.con)) 

# lets look at the  vertices (nodes) 
V(grph.pos.con)
# and at the edges of the network
E(grph.pos.con)

# Remove vertices not representing cross domain connections
# I only found a way if I remove them one at a time
# The number is the vertex "location" in  V(grph.pos.con.cross)
# The first omition is from grph.pos.con, the subsequent ones from grph.pos.con.cross

grph.pos.con.cross <- delete_vertices(grph.pos.con, 34) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 33) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 32) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 31) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 29) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 28) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 27) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 26) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 25) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 24) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 23) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 22) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 23) 
grph.pos.con.cross <- delete_vertices(grph.pos.con.cross, 22) 

# lets look at the  vertices (nodes) 
V(grph.pos.con.cross)
# and at the edges of the network
E(grph.pos.con.cross)



plot(grph.pos.con.cross, vertex.label.cex=0.5, vertex.label.color="black", 
     edge.color="black",
     layout=layout_with_fr(grph.pos.con.cross))

pdf("/home/oded/Documents/post/articles/springs/speiceaci/biofilm20/THS/to_report/biofilm_speiceaci.pdf")
dev.off()




# We can plot the degree distribution.
# Degree is the number of edges that a node (vertex) has. We see that most nodes are connected to few other nodes when consideringly only positive interactions above our weight threshold.
dd.grph.pos.con <- degree.distribution(grph.pos.con)
plot(0:(length(dd.grph.pos.con)-1), dd.grph.pos.con, type='b',
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

grph.pos.con_deg<-degree(grph.pos.con, v=V(grph.pos.con), mode="all")

fine = 500 # this will adjust the resolving power.

#this gives you the colors you want for every point
graphCol = viridis(fine)[as.numeric(cut(grph.pos.con_deg,breaks = fine))]

# now plot
plot(grph.pos.con, vertex.color=graphCol,
     edge.color="black",
     vertex.label.cex=0.5, nvertex.label.color="black",
     layout=layout_with_fr(grph.pos.con))


png(file = "/home/oded/Documents/post/articles/springs/speiceaci/biofilm20/THS/NETWORK_FIG/biofilm_speiceaci.png", width = 4000, height = 8000)

dev.off()

# check how far we are from the target stability threshold
getOptInd(se.Leg.host)
#and
sum(getRefit(se.Leg.host))/2

getStability(se.Leg.host)





ig.mb     <- adj2igraph(getRefit(se.Leg.host))
#vsize    <- rowMeans(clr(amgut1.filt, 1))+6
#am.coord <- layout.fruchterman.reingold(ig.mb)


spiec.deg <- degree(ig.mb)
hist(spiec.deg)
fit <- fit_power_law(spiec.deg)
fit



# A common issue that comes up with when running spiec.easi is coming up with an empty network after running StARS.
# look if there is a warning

# Looking at the stability along the lambda path:
se.Leg.host$select$stars$summary
