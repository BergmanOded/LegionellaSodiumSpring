metadata <- read.table(file = "/home/oded/Documents/post/articles/springs/18S/NGS/alpha_diversity/metadata.csv", sep = "\t", header = T)
str(metadata)

library(purrr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(rstatix)
library(ez)
library(ggpubr)


metadata_water <- dplyr::filter(metadata, grepl('water', environment))
metadata_biofilm <- dplyr::filter(metadata, grepl('biofilm', environment))


#statistics water----

metadata_water %>%
  group_by(collection_date) %>%
  get_summary_stats(shannon_entropy, faith_pd, pielou_evenness, observed_OTUs, type = "median")

metadata_water %>%
  group_by(cluster) %>%
  get_summary_stats(shannon_entropy, faith_pd, pielou_evenness, observed_OTUs, type = "median")

metadata_biofilm %>%
  group_by(collection_date) %>%
  get_summary_stats(shannon_entropy, faith_pd, pielou_evenness, observed_OTUs, type = "median")

metadata_biofilm %>%
  group_by(cluster) %>%
  get_summary_stats(shannon_entropy, faith_pd, pielou_evenness, observed_OTUs, type = "median")

metadata_biofilm %>%
  group_by(cluster) %>%
  get_summary_stats(shannon_entropy, type = "median")
# As in tests we have small sample sizes (4 and under), as well as uneven groups, I chose a non-parametric one way anova

kruskal.test(shannon_entropy ~ cluster,
             data = metadata_water)
kruskal.test(faith_pd ~ cluster,
             data = metadata_water)
kruskal.test(pielou_evenness ~ cluster,
             data = metadata_water)
kruskal.test(observed_OTUs ~ cluster,
             data = metadata_water)

kruskal.test(shannon_entropy ~ collection_date,
             data = metadata_water)
kruskal.test(faith_pd ~ collection_date,
             data = metadata_water)
kruskal.test(pielou_evenness ~ collection_date,
             data = metadata_water)
kruskal.test(observed_OTUs ~ collection_date,
             data = metadata_water) 

kruskal.test(shannon_entropy ~ cluster,
             data = metadata_biofilm)
kruskal.test(faith_pd ~ cluster,
             data = metadata_biofilm)
kruskal.test(pielou_evenness ~ cluster,
             data = metadata_biofilm)
kruskal.test(observed_OTUs ~ cluster,
             data = metadata_biofilm)

kruskal.test(shannon_entropy ~ collection_date,
             data = metadata_biofilm)
kruskal.test(faith_pd ~ collection_date,
             data = metadata_biofilm)
kruskal.test(pielou_evenness ~ collection_date,
             data = metadata_biofilm)
kruskal.test(observed_OTUs ~ collection_date,
             data = metadata_biofilm)

# cluster was significant in all metrics for both metadata_water and metadata_biofilm
# post hock

# Pairwise comparisons (p.adjust.method = Bonferroni-Holm) - only for origin (significant )

post.hoc.water.cluster.shannon <- metadata_water %>%
  wilcox_test(shannon_entropy ~ cluster, p.adjust.method = "BH")
post.hoc.water.cluster.shannon

post.hoc.water.cluster.faith_pd <- metadata_water %>%
  wilcox_test(faith_pd ~ cluster, p.adjust.method = "BH")
post.hoc.water.cluster.faith_pd

post.hoc.water.cluster.evenness <- metadata_water %>%
  wilcox_test(pielou_evenness ~ cluster, p.adjust.method = "BH")
post.hoc.water.cluster.evenness

post.hoc.water.cluster.OTU <- metadata_water %>%
  wilcox_test(observed_OTUs ~ cluster, p.adjust.method = "BH")
post.hoc.water.cluster.OTU


post.hoc.biofilm.cluster.shannon <- metadata_biofilm %>%
  wilcox_test(shannon_entropy ~ cluster, p.adjust.method = "BH")
post.hoc.biofilm.cluster.shannon

post.hoc.biofilm.cluster.faith_pd <- metadata_biofilm %>%
  wilcox_test(faith_pd ~ cluster, p.adjust.method = "BH")
post.hoc.biofilm.cluster.faith_pd

post.hoc.biofilm.cluster.evenness <- metadata_biofilm %>%
  wilcox_test(pielou_evenness ~ cluster, p.adjust.method = "BH")
post.hoc.biofilm.cluster.evenness

post.hoc.biofilm.cluster.OTU <- metadata_biofilm %>%
  wilcox_test(observed_OTUs ~ cluster, p.adjust.method = "BH")
post.hoc.biofilm.cluster.OTU


# diversity boxplot cluster----

str(metadata)

metadata_water$cluster <- factor(metadata_water$cluster, levels=c("Tiberias_Hot_Springs", "Haon", "Tabgha", "Fuliya"))
metadata_biofilm$cluster <- factor(metadata_biofilm$cluster, levels=c("Tiberias_Hot_Springs", "Haon", "Tabgha", "Fuliya"))

my.labels.water <- c("Tiberias Hot Springs", "Haon-Borehole", "Tabgha", "Fuliya")
my.labels.biofilm <- c("Tiberias Hot Springs", "Tabgha", "Fuliya")
# water

#shanon

p.water.cluster.shannon <- ggplot(metadata_water, aes(x=cluster, y=shannon_entropy, fill = cluster)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_y_continuous(limits=c(0, 10)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Shannon's entropy")))

  p.water.cluster.shannon.text <- p.water.cluster.shannon + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=8, vjust = 2.5),
  axis.text.x = element_blank(),
  axis.text.y = element_text(color="black", size=8),
  legend.title = element_text(color = "black", size=8, face ="bold"),
  legend.text = element_text(color = "black", size=8)) 

# match the colours to those of the PCA plot + statistics 
p.water.cluster.shannon.text.col <- p.water.cluster.shannon.text + scale_fill_manual(values = c("#CC79A7", "#009E73", "#56B4E9", "#D55E00", "#E69F00"))
post.hoc.water.cluster.shannon <- post.hoc.water.cluster.shannon %>% add_xy_position
p.water.stat.cluster.shannon.text.col <- p.water.cluster.shannon.text.col + stat_pvalue_manual(post.hoc.water.cluster.shannon, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.002, y.position = c(6, 7, 8, 9, 10))
p.water.stat.cluster.shannon.text.col


# faith_pd

p.water.cluster.faith_pd <- ggplot(metadata_water, aes(x=cluster, y=faith_pd, fill = cluster)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_y_continuous(breaks=seq(0,17.5,2.5), limits = c(0, 17.5)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Faith's PD")))

p.water.cluster.faith_pd.text <- p.water.cluster.faith_pd + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=8, vjust = 2.5),
  axis.text.x = element_blank(),
  axis.text.y = element_text(color="black", size=8),
  legend.title = element_text(color = "black", size=8, face ="bold"),
  legend.text = element_text(color = "black", size=8))  

# match the colours to those of the PCA plot  
p.water.cluster.faith_pd.text.col <- p.water.cluster.faith_pd.text + scale_fill_manual(values = c("#CC79A7", "#009E73", "#56B4E9", "#D55E00", "#E69F00"))
post.hoc.water.cluster.faith_pd <- post.hoc.water.cluster.faith_pd %>% add_xy_position
p.water.stat.cluster.faith_pd.text.col <- p.water.cluster.faith_pd.text.col + stat_pvalue_manual(post.hoc.water.cluster.faith_pd, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.002, y.position = c(13, 14.5, 16, 17.5))
p.water.stat.cluster.faith_pd.text.col


# pielou_evenness

p.water.cluster.evenness <- ggplot(metadata_water, aes(x=cluster, y=pielou_evenness, fill = cluster)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_y_continuous(limits=c(0, 2)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Pielou's evenness")))

p.water.cluster.evenness.text <- p.water.cluster.evenness + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=8, vjust = 2.5),
  axis.text.x = element_blank(),
  axis.text.y = element_text(color="black", size=8),
  legend.title = element_text(color = "black", size=8, face ="bold"),
  legend.text = element_text(color = "black", size=8))  

# match the colours to those of the PCA plot  
  p.water.cluster.evenness.text.col <- p.water.cluster.evenness.text + scale_fill_manual(values = c("#CC79A7", "#009E73", "#56B4E9", "#D55E00", "#E69F00"))
p.water.cluster.evenness.text.col
post.hoc.water.cluster.evenness <- post.hoc.water.cluster.evenness %>% add_xy_position
p.water.stat.cluster.evenness.text.col <- p.water.cluster.evenness.text.col + stat_pvalue_manual(post.hoc.water.cluster.evenness, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.002, y.position = c(1.3, 1.6))
p.water.stat.cluster.evenness.text.col
# observed_OTUs

p.water.cluster.OTU <- ggplot(metadata_water, aes(x=cluster, y=observed_OTUs, fill = cluster)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Spring Cluster", labels= my.labels.water) + scale_y_continuous(limits=c(0, 150)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Observed ASVs")))

p.water.cluster.OTU.text <- p.water.cluster.OTU + theme(
axis.title.x = element_blank(),
axis.title.y = element_text(color="black", size=8),
axis.text.y = element_text(color="black", size=8),
axis.text.x = element_text(color="black", size=8, vjust = -1),
legend.title = element_text(color = "black", size=8, face ="bold"),
legend.text = element_text(color = "black", size=8)) 



# match the colours to those of the PCA plot  
p.water.cluster.OTU.text.col <- p.water.cluster.OTU.text + scale_fill_manual(values = c("#CC79A7", "#009E73", "#56B4E9", "#D55E00", "#E69F00"))
post.hoc.water.cluster.OTU <- post.hoc.water.cluster.OTU %>% add_xy_position
p.water.stat.cluster.OTU.text.col <- p.water.cluster.OTU.text.col + stat_pvalue_manual(post.hoc.water.cluster.OTU, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.002, y.position = c(100, 115, 130, 145))
p.water.stat.cluster.OTU.text.col


#biofilm

#shanon

p.biofilm.cluster.shannon <- ggplot(metadata_biofilm, aes(x=cluster, y=shannon_entropy, fill = cluster)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_y_continuous(limits=c(0, 10)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Shannon's entropy")))

p.biofilm.cluster.shannon.text <- p.biofilm.cluster.shannon + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=8, vjust = 2.5),
  axis.text.x = element_blank(),
  axis.text.y = element_text(color="black", size=8),
  legend.title = element_text(color = "black", size=8, face ="bold"),
  legend.text = element_text(color = "black", size=8)) 

# match the colours to those of the PCA plot + statistics 
p.biofilm.cluster.shannon.text.col <- p.biofilm.cluster.shannon.text + scale_fill_manual(values = c("#CC79A7", "#56B4E9", "#D55E00"))
post.hoc.biofilm.cluster.shannon <- post.hoc.biofilm.cluster.shannon %>% add_y_position
p.biofilm.stat.cluster.shannon.text.col <- p.biofilm.cluster.shannon.text.col + stat_pvalue_manual(post.hoc.biofilm.cluster.shannon, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.002, y.position = c(8, 9.5))
p.biofilm.stat.cluster.shannon.text.col


# faith_pd

p.biofilm.cluster.faith_pd <- ggplot(metadata_biofilm, aes(x=cluster, y=faith_pd, fill = cluster)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_y_continuous(breaks=seq(0,17.5,2.5), limits = c(0, 17.5)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Faith's PD")))

p.biofilm.cluster.faith_pd.text <- p.biofilm.cluster.faith_pd + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=8, vjust = 2.5),
  axis.text.x = element_blank(),
  axis.text.y = element_text(color="black", size=8),
  legend.title = element_text(color = "black", size=8, face ="bold"),
  legend.text = element_text(color = "black", size=8))  

# match the colours to those of the PCA plot  
p.biofilm.cluster.faith_pd.text.col <- p.biofilm.cluster.faith_pd.text + scale_fill_manual(values = c("#CC79A7", "#56B4E9", "#D55E00"))
post.hoc.biofilm.cluster.faith_pd <- post.hoc.biofilm.cluster.faith_pd %>% add_y_position
p.biofilm.stat.cluster.faith_pd.text.col <- p.biofilm.cluster.faith_pd.text.col + stat_pvalue_manual(post.hoc.biofilm.cluster.faith_pd, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.002, y.position = c(15.5, 17.5))
p.biofilm.stat.cluster.faith_pd.text.col


# pielou_evenness

p.biofilm.cluster.evenness <- ggplot(metadata_biofilm, aes(x=cluster, y=pielou_evenness, fill = cluster)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_y_continuous(limits=c(0, 2)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Pielou's evenness")))

p.biofilm.cluster.evenness.text <- p.biofilm.cluster.evenness + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=8, vjust = 2.5),
  axis.text.x = element_blank(),
  axis.text.y = element_text(color="black", size=8),
  legend.title = element_text(color = "black", size=8, face ="bold"),
  legend.text = element_text(color = "black", size=8))  

# match the colours to those of the PCA plot  
p.biofilm.cluster.evenness.text.col <- p.biofilm.cluster.evenness.text + scale_fill_manual(values = c("#CC79A7", "#56B4E9", "#D55E00"))
p.biofilm.cluster.evenness.text.col
#post.hoc.biofilm.cluster.evenness <- post.hoc.biofilm.cluster.evenness %>% add_y_position
#p.biofilm.stat.cluster.evenness.text.col <- p.biofilm.cluster.evenness.text.col + stat_pvalue_manual(post.hoc.biofilm.cluster.evenness, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.002, y.position = c(1.5, 1.8))
#p.biofilm.stat.cluster.evenness.text.col

# observed_OTUs

p.biofilm.cluster.OTU <- ggplot(metadata_biofilm, aes(x=cluster, y=observed_OTUs, fill = cluster)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Spring Cluster", labels= my.labels.biofilm) + scale_y_continuous(limits=c(0, 150)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Observed ASVs")))

p.biofilm.cluster.OTU.text <- p.biofilm.cluster.OTU + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=8),
  axis.text.y = element_text(color="black", size=8),
  axis.text.x = element_text(color="black", size=8, vjust = -1),
  legend.title = element_text(color = "black", size=8, face ="bold"),
  legend.text = element_text(color = "black", size=8)) 

# match the colours to those of the PCA plot  
p.biofilm.cluster.OTU.text.col <- p.biofilm.cluster.OTU.text + scale_fill_manual(values = c("#CC79A7", "#56B4E9", "#D55E00"))
post.hoc.biofilm.cluster.OTU <- post.hoc.biofilm.cluster.OTU %>% add_y_position
p.biofilm.stat.cluster.OTU.text.col <- p.biofilm.cluster.OTU.text.col + stat_pvalue_manual(post.hoc.biofilm.cluster.OTU, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.002, y.position = c(130, 145))
p.biofilm.stat.cluster.OTU.text.col




# Collection date graph with the date ordered:

#date water

my.labels.date <- c("Sep17", "Jan18", "Jun18", "Oct18", "Jan19")


p.water.date.shannon <- ggplot(metadata_water, aes(x=collection_date, y=shannon_entropy)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Collection date", labels= my.labels.date) + scale_y_continuous(limits=c(0, 10)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Shannon's entropy"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=8, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black", size=8),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8))

p.water.date.shannon

p.water.date.faith_pd <- ggplot(metadata_water, aes(x=collection_date, y=faith_pd)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Collection date", labels= my.labels.date) + scale_y_continuous(breaks=seq(0,17.5,2.5), limits = c(0, 17.5)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Faith's PD"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=8, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black", size=8),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8))

p.water.date.faith_pd

p.water.date.evenness <- ggplot(metadata_water, aes(x=collection_date, y=pielou_evenness)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Collection date", labels= my.labels.date) + scale_y_continuous(limits=c(0, 2)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Pielou's evenness"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=8, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black", size=8),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8))

p.water.date.evenness

p.water.date.OTU <- ggplot(metadata_water, aes(x=collection_date, y=observed_OTUs)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Collection date", labels= my.labels.date) + scale_y_continuous(limits=c(0, 150)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Observed ASVs"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    axis.text.x = element_text(color="black", size=8, vjust = -1),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8)) 

p.water.date.OTU


#date biofilm

p.biofilm.date.shannon <- ggplot(metadata_biofilm, aes(x=collection_date, y=shannon_entropy)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Collection date", labels= my.labels.date) + scale_y_continuous(limits=c(0, 10)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Shannon's entropy"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=8, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black", size=8),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8))

p.biofilm.date.shannon

p.biofilm.date.faith_pd <- ggplot(metadata_biofilm, aes(x=collection_date, y=faith_pd)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Collection date", labels= my.labels.date) + scale_y_continuous(breaks=seq(0,17.5,2.5), limits = c(0, 17.5)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Faith's PD"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=8, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black", size=8),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8))

p.biofilm.date.faith_pd

p.biofilm.date.evenness <- ggplot(metadata_biofilm, aes(x=collection_date, y=pielou_evenness)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Collection date", labels= my.labels.date) + scale_y_continuous(limits=c(0, 2)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Pielou's evenness"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=8, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="black", size=8),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8))

p.biofilm.date.evenness

p.biofilm.date.OTU <- ggplot(metadata_biofilm, aes(x=collection_date, y=observed_OTUs)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Collection date", labels= my.labels.date) + scale_y_continuous(limits=c(0, 150)) + 
  theme(legend.position = "none") + labs(y=expression(paste("Observed ASVs"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    axis.text.x = element_text(color="black", size=8, vjust = -1),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8)) 

p.biofilm.date.OTU

#combine with cowplot

water_date_alpha_diversity <- cowplot::plot_grid(p.water.date.shannon, p.water.date.faith_pd, p.water.date.evenness, p.water.date.OTU, align = "v", ncol = 1, rel_heights = c(0.225, 0.225, 0.225, 0.325))

biofilm_date_alpha_diversity <- cowplot::plot_grid(p.biofilm.date.shannon, p.biofilm.date.faith_pd, p.biofilm.date.evenness, p.biofilm.date.OTU, align = "v", ncol = 1, rel_heights = c(0.225, 0.225, 0.225, 0.325))

water_biofilm_date_alpha_diversity <- cowplot::plot_grid(water_date_alpha_diversity, biofilm_date_alpha_diversity, align = "v", ncol = 2)


water_cluster_alpha_diversity <- cowplot::plot_grid(p.water.stat.cluster.shannon.text.col, p.water.stat.cluster.faith_pd.text.col, p.water.cluster.evenness.text.col, p.water.stat.cluster.OTU.text.col, align = "v", ncol = 1, rel_heights = c(0.225, 0.225, 0.225, 0.325))

biofilm_cluster_alpha_diversity <- cowplot::plot_grid(p.biofilm.stat.cluster.shannon.text.col, p.biofilm.stat.cluster.faith_pd.text.col, p.biofilm.cluster.evenness.text.col, p.biofilm.stat.cluster.OTU.text.col, align = "v", ncol = 1, rel_heights = c(0.225, 0.225, 0.225, 0.325))

water_biofilm_cluster_alpha_diversity <- cowplot::plot_grid(water_cluster_alpha_diversity, biofilm_cluster_alpha_diversity, align = "v", ncol = 2, rel_widths = c(0.55,0.45))

water_biofilm_cluster_date_alpha_diversity <- cowplot::plot_grid(water_biofilm_cluster_alpha_diversity, water_biofilm_date_alpha_diversity, align = "v", ncol = 2)
water_biofilm_cluster_date_alpha_diversity

pdf("/home/oded/Documents/post/articles/springs/18S/NGS/alpha_diversity/18S_Legionella_hosts_alpha_diversity.pdf",  width = 20, height = 10, family="ArialMT")

water_biofilm_cluster_date_alpha_diversity
dev.off()

png("/home/oded/Documents/post/articles/springs/18S/NGS/alpha_diversity/18S_Legionella_hosts_alpha_diversity.png",  width = 14800, height = 9000, res = 1200)

# alpha diversity of environment (water and biofilm)


kruskal.test(shannon_entropy ~ environment,
             data = metadata)
kruskal.test(faith_pd ~ environment,
             data = metadata)
kruskal.test(pielou_evenness ~ environment,
             data = metadata)
kruskal.test(observed_OTUs ~ environment,
             data = metadata)

my.labels.env <- c("Weter", "Biofilm")

#shannon

p.water.biofilm.env.shannon <- ggplot(metadata, aes(x=environment, y=shannon_entropy, fill = environment)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Environment", labels= my.labels.env) + scale_y_continuous(limits=c(0, 10)) + 
  labs(y=expression(paste("Shannon's entropy"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=8, vjust = 2.5),
    axis.text.x = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8))

# match the colours to those of the PCA plot  
p.water.biofilm.env.shannon.col <- p.water.biofilm.env.shannon + scale_fill_manual(values = c("#6b440b", "#7d80d8"))

#faith_pd

p.water.biofilm.env.faith <- ggplot(metadata, aes(x=environment, y=faith_pd, fill = environment)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Environment", labels= my.labels.env) + scale_y_continuous(limits=c(0, 12)) + 
  labs(y=expression(paste("Faith's PD"))) + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(color="black", size=8, vjust = 2.5),
    axis.text.x = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8))

# match the colours to those of the PCA plot  
p.water.biofilm.env.faith.col <- p.water.biofilm.env.faith + scale_fill_manual(values = c("#6b440b", "#7d80d8"))

#evenness

p.water.biofilm.env.evenness <- ggplot(metadata, aes(x=environment, y=pielou_evenness, fill = environment)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Environment", labels= my.labels.env) + scale_y_continuous(limits=c(0, 2)) + 
  labs(y=expression(paste("Pielou's evenness"))) + theme(
    axis.title.x = element_text(margin = margin(t = 22), color="black", size=8, vjust=5),
    axis.title.y = element_text(color="black", size=8, vjust = 2.5),
    axis.text.x = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8))

# match the colours to those of the PCA plot  
p.water.biofilm.env.evenness.col <- p.water.biofilm.env.evenness + scale_fill_manual(values = c("#6b440b", "#7d80d8"))


#observed OTUs

p.water.biofilm.env.OTU <- ggplot(metadata, aes(x=environment, y=observed_OTUs, fill = environment)) +   
  geom_boxplot() +
  theme_classic() + 
  scale_x_discrete(name ="Environment", labels= my.labels.env) + scale_y_continuous(limits=c(0, 150)) + 
  labs(y=expression(paste("Observed ASVs"))) + theme(
    axis.title.x = element_text(margin = margin(t = 22), color="black", size=8, vjust=5),
    axis.title.y = element_text(color="black", size=8, vjust = 2.5),
    axis.text.x = element_text(color="black", size=8),
    axis.text.y = element_text(color="black", size=8),
    legend.title = element_text(color = "black", size=8, face ="bold"),
    legend.text = element_text(color = "black", size=8))

# match the colours to those of the PCA plot  
p.water.biofilm.env.OTU.col <- p.water.biofilm.env.OTU + scale_fill_manual(values = c("#6b440b", "#7d80d8"))


water_biofilm_date_alpha_diversity <- cowplot::plot_grid(p.water.biofilm.env.shannon.col + theme(legend.position="none"), 
                                                         p.water.biofilm.env.faith.col+ theme(legend.position="none"), 
                                                         p.water.biofilm.env.evenness.col+ theme(legend.position="none"), 
                                                         p.water.biofilm.env.OTU.col+ theme(legend.position="none"), 
                                                         align = "v", ncol = 2, rel_heights = c(0.24, 0.28, 0.24, 0.28))
water_biofilm_date_alpha_diversity

legend_b <- get_legend( water_biofilm_date_alpha_diversity + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top"))

library(ggplot2)
library(cowplot)
library(rlang)

plot_grid(water_biofilm_date_alpha_diversity, legend_b, ncol = 1, rel_heights = c(1, .1))
