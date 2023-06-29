library(tidyverse)
library(qiime2R)
library(dplyr)
library(tidytext)
    library(phyloseq)
library(ggplot2)
library(ggh4x)
  # barplot for ASVs in water + biofilm 20 reads in  2 samples min----
# we are interested in all freatures (to look at overall diversity)
metadata20 <- read.table(file = "/home/oded/Documents/post/qiime2/16S+18S_10.01.20/10.01.2020/raw_data_merged/16S-Legionella/metadata/merge/spring/bar-plots-r/reverse_metatada_name/new/omit_06.2018/sample-metadata.tsv", sep = ",", header = T, 
                               comment.char = "", row.names = 1, check.names = FALSE)
  
  # make the levels to facet later
  metadata20$origin_new = factor(metadata20$origin_new, levels=c('Tiberias Hot Springs','Haon','Tabgha','Fuliya'))
metadata20$environment = factor(metadata20$environment, levels=c('Water','Biofilm'))

  ASVs20 <- read.table(file = "/home/oded/Documents/post/articles/springs/NGS/barplots/final_10.10.21/feature-table.tsv", sep = "\t", header = T, 
                           comment.char = "", row.names = 1, check.names = FALSE)

taxonomy<-read_qza("/home/oded/Documents/post/articles/springs/NGS/barplots/Final_corrected/ASVs/water+biofilm/feature_table_20/final4/taxonomy.qza")$data %>% parse_taxonomy()

ASVs = otu_table(as.matrix(ASVs20), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata20)
taxasums20 <- tax_table(as.matrix(taxonomy))
BP20 <- phyloseq(ASVs, SAMPLE, taxasums20)
# convert phyloseq to relative abundance
BP20.rel = transform_sample_counts(BP20, function(x) x / sum(x) )
# to get the taxonomy.qza of ASVs, I exported the original QIIME2 taxonomy.qza to TSV
#and copied the 'freature ID' to the 'taxon' colomn. Then imported back to QIIME2 as .qza
# once read with read_qza, it is persed with the feature ID in Kingdom. All else is NA.
# We latter use the Kingdom as identifier. 


#define colors
library(RColorBrewer)
nb.cols20 <- 1074
mycolors20 <- colorRampPalette(brewer.pal(11, "Paired"))(nb.cols20)
# geom_bar makes the borders within the bars

#all samples, by facet spring cluster (origin), not ordered
barplot20 <- plot_bar(BP20.rel, x= "Station.date", y = "Abundance", fill="Kingdom") + facet_nested(environment~origin_new+station, scales = "free_x", space = "free_x", nest_line = element_line(linetype = 1)) + 
ylab("Abundance (proportion)") + xlab ("Sample") +
scale_x_discrete( labels = c("Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019")) +
theme(axis.title.y = element_text(color="black", size=8, vjust = 3),
      axis.text.y = element_text(color="black", size=8),
                                              axis.title.x = element_blank(), 
                                              axis.text.x = element_blank()) +
  theme(strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"))+
                              scale_fill_manual(values = mycolors20, guide="none") +
 geom_bar(stat="identity", position="stack", color="black", size=0.02) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black"))
  barplot20 

# barplot for ASVs in water + biofilm filter for 2000 reads in  2 samples min (more abundant ASVs)----

  # Now we want to look at this only for more abundant ASVs, 
  # So we will do will filter water and biofilm samples seperatly 
  # (as filtering together meens a freature with less then 200o reads in each environment will be in).
  # So we have two matadata, freature_table and taxomony files for each.
  # Then we will merge them. All this in QIIME2.
  # From there it is the same in r.



ASVs2000 <- read_qza("/home/oded/Documents/post/qiime2/16S+18S_10.01.20/10.01.2020/raw_data_merged/16S-Legionella/merge-artifacts/filter/cluster100/analysis/springs/filter/frequency/20/taxonomy/barplots/Legionella-barplot_frequency2000-filter_table.qza")$data

# to get the taxonomy.qza of ASVs, I exported the original QIIME2 taxonomy.qza to TSV
#and copied the 'freature ID' to the 'taxon' colomn. Then imported back to QIIME2 as .qza
# once read with read_qza, it is persed with the feature ID in Kingdom. All else is NA.
# We latter use the Kingdom as identifier. 

taxasums2000<-summarize_taxa(ASVs2000, taxonomy)$Kingdom

ASVs2000 = otu_table(as.matrix(ASVs2000), taxa_are_rows = TRUE)
SAMPLE <- sample_data(metadata20)
taxasums20 <- tax_table(as.matrix(taxonomy))
BP2000 <- phyloseq(ASVs2000, SAMPLE, taxasums20)
# convert phyloseq to relative abundance
BP2000.rel = transform_sample_counts(BP2000, function(x) x / sum(x) )



#define colors
library(RColorBrewer)
nb.cols2000 <- 117
mycolors2000 <- colorRampPalette(brewer.pal(11, "Paired"))(nb.cols2000)


barplot2000 <- plot_bar(BP2000.rel, x= "Station.date", y = "Abundance", fill="Kingdom") + facet_nested(environment~origin_new+station, scales = "free_x", space = "free_x") + 
ylab("Abundance (proportion)") +
  scale_x_discrete( labels = c("Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan2019", "Sep-2017", "Jan-2018", "Oct-2018", "Jan-2019")) +
  theme(axis.title.y = element_text(color="black", size=8, vjust = 3),
        axis.text.y = element_text(color="black", size=8),
        axis.title.x = element_blank()) + 
  geom_bar(stat="identity", position="stack", color="black", size=0.02) +
theme(strip.background = element_blank(), strip.text.x = element_blank())+
  scale_fill_manual(values = mycolors2000, guide="none") +
  geom_bar(stat="identity", position="stack", color="black", size=0.02) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black"))
barplot2000




  
ASVs_20_2000_final <- cowplot::plot_grid(barplot20, barplot2000, align = "v", ncol = 1, rel_heights = c(0.4, 0.5), scale = 0.9)
ASVs_20_2000_final

  pdf("/home/oded/Documents/post/articles/springs/final R code/NGS/16S/barplots/Leg_barplots.pdf",  width = 15, height = 10, family="ArialMT")
ASVs_20_2000_final
dev.off()






