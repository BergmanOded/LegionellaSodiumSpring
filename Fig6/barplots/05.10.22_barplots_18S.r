library(tidyverse)
library(qiime2R)
library(dplyr)
library(tidytext)
library(phyloseq)
library(ggplot2)
library(ggh4x)
# barplot for ASVs in water + biofilm 20 reads in  2 samples min----
# we are interested in all freatures (to look at overall diversity)
metadata20 <- read.table(file = "/home/oded/Documents/post/articles/springs/18S/NGS/barplots/sample-metadata.tsv", sep = "\t", header = T, 
                         comment.char = "", row.names = 1, check.names = FALSE)


# make the levels to facet later
  metadata20$origin_new = factor(metadata20$origin_new, levels=c('Tiberias Hot Springs','Haon','Tabgha','Fuliya'))
metadata20$environment = factor(metadata20$environment, levels=c('Water','Biofilm'))

ASVs20 <- read.table(file = "/home/oded/Documents/post/articles/springs/18S/NGS/barplots/feature-table.tsv", sep = "\t", header = T, 
                     comment.char = "", row.names = 1, check.names = FALSE)

ASVs500 <- read.table(file = "/home/oded/Documents/post/articles/springs/18S/NGS/barplots/500/feature-table.tsv", sep = "\t", header = T, 
                     comment.char = "", row.names = 1, check.names = FALSE)
# I want the taxonomy at phylum, order and genus only for 20 reads filtering. 
taxonomy<-read_qza("/home/oded/Documents/post/qiime2/16S+18S_10.01.20/10.01.2020/raw_data_merged/18S/100-120/merge-artifacts/filter/cluster100/analysis/total/taxonomy/taxonomy-18S-total-138-99.qza")$data %>% parse_taxonomy()

taxonomy_kingdom <- read.table(file = "/home/oded/Documents/post/articles/springs/18S/NGS/barplots/taxonomy_ASVs.csv", sep = "\t", header = T, 
                     comment.char = "", row.names = 1, check.names = FALSE)


taxasums20.phylum<-summarize_taxa(ASVs20, taxonomy)$Phylum
taxasums20.order<-summarize_taxa(ASVs20, taxonomy)$Order
taxasums20.genus<-summarize_taxa(ASVs20, taxonomy)$Genus
taxasums500.genus<-summarize_taxa(ASVs500, taxonomy)$Genus

ASVs20 = otu_table(as.matrix(ASVs20), taxa_are_rows = TRUE)
ASVs500 = otu_table(as.matrix(ASVs500), taxa_are_rows = TRUE)

SAMPLE <- sample_data(metadata20)
taxasums_kingdom20 <- tax_table(as.matrix(taxonomy_kingdom))
taxasums20 <- tax_table(as.matrix(taxonomy))

# for phylum barplot - we need the taxonomy
BP20_phylum <- phyloseq(ASVs20, SAMPLE, taxasums20)
# for the ASVs20 and ASVs500 barplot we imput the ASV as 'Kingdom"
BP20 <- phyloseq(ASVs20, SAMPLE, taxasums_kingdom20)
BP500 <- phyloseq(ASVs500, SAMPLE, taxasums_kingdom20)

# convert phyloseq to relative abundance
BP20.rel.phylum = transform_sample_counts(BP20_phylum, function(x) x / sum(x) )
BP20.rel = transform_sample_counts(BP20, function(x) x / sum(x) )
BP500.rel = transform_sample_counts(BP500, function(x) x / sum(x) )

#define colors for ASVs20
library(RColorBrewer)
nb.cols20 <- 501
mycolors20 <- colorRampPalette(brewer.pal(11, "Paired"))(nb.cols20)

#define colors for ASVs500

#define colors
library(RColorBrewer)
nb.cols500.genus <- 96
mycolors500 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols500.genus)

#all samples, by facet spring cluster (origin), not ordered
barplot20.phylum <- plot_bar(BP20.rel.phylum, x= "Station.date", y = "Abundance", fill="Phylum") + facet_nested(environment~origin_new+station, scales = "free_x", space = "free_x") + 
  ylab("Abundance (proportion)") + xlab ("Sample") +
  scale_x_discrete( labels = c("Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019")) +
  theme(axis.title.y = element_text(color="black", size=8, vjust = 3),
        axis.text.y = element_text(color="black", size=8),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank()) +
  theme(legend.position = "top")+
  theme(strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"))+
  geom_bar(stat="identity", position="stack", color="black", size=0.02) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(panel.background = element_blank()) +  scale_fill_manual(values = mycolors20, guide="none") +

  theme(axis.line = element_line(colour = "black")) +
scale_fill_manual(values=c("mediumpurple1", "darksalmon" ))
barplot20.phylum


  #all samples, ASVs20 by facet spring cluster (origin)
barplot20 <- plot_bar(BP20.rel, x= "Station.date", y = "Abundance", fill="Kingdom") + facet_nested(environment~origin_new+station, scales = "free_x", space = "free_x") + 
  theme(legend.position = "none") + ylab("Abundance (proportion)") +
  scale_x_discrete( labels = c("Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019")) +
  theme(axis.title.y = element_text(color="black", size=8, vjust = 3),
        axis.text.y = element_text(color="black", size=8),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank()) +
  geom_bar(stat="identity", position="stack", color="black", size=0.02) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = mycolors20)
barplot20



#all samples, ASVs500 by facet spring cluster (origin)
barplot500 <- plot_bar(BP500.rel, x= "Station.date", y = "Abundance", fill="Kingdom") + facet_nested(environment~origin_new+station, scales = "free_x", space = "free_x") + 
  theme(legend.position = "none") + ylab("Abundance (proportion)") +
  scale_x_discrete( labels = c("Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019")) +
  theme(axis.title.y = element_text(color="black", size=8, vjust = 3),
        axis.text.y = element_text(color="black", size=8),
        axis.title.x = element_blank())+ 
  geom_bar(stat="identity", position="stack", color="black", size=0.02) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = mycolors500)
barplot500



barplot_combines_phylum_ASVs_final <- cowplot::plot_grid(barplot20.phylum, barplot20, barplot500, align = "v", ncol = 1, rel_heights = c(0.25, 0.22, 0.3), scale = 0.9)
barplot_combines_phylum_ASVs_final  
  
pdf("/home/oded/Documents/post/articles/springs/18S/NGS/barplots/12.20.21_18S_barplots.pdf",  width = 15, height = 10, family="ArialMT")
barplot_combines_phylum_ASVs_final
dev.off()