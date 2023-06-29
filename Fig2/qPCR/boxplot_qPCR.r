# Box plot qPCR springs cluster using Legionella primers

sample.metadatawater <- read.table(file = "/home/oded/Documents/post/articles/springs/qPCR/spring_cluster/water+biofilm/qPCR_Legionella_springs/statistics/water/sample-metadatawater.tsv", sep = ",", header = T)
sample.metadatabiofilm <- read.table(file = "/home/oded/Documents/post/articles/springs/qPCR/spring_cluster/water+biofilm/qPCR_Legionella_springs/statistics/biofilm/sample-metadatabiofilm.tsv", sep = ",", header = T)

str(sample.metadatawater)
str(sample.metadatabiofilm)

sample.metadatawater$Legionella.count<-as.numeric(sample.metadatawater$Legionella.count)
sample.metadatabiofilm$Legionella.count<-as.numeric(sample.metadatabiofilm$Legionella.count)


#add log transformations
sample.metadatawater$Legionella.count.Log <- log(1 +sample.metadatawater$Legionella.count)
sample.metadatabiofilm$Legionella.count.Log <- log(1 +sample.metadatabiofilm$Legionella.count)

str(sample.metadatawater)
str(sample.metadatabiofilm)


#outliers (we did not omit samples)
library(rstatix)

outliers_water_springs <- sample.metadatawater%>%
  group_by(origin) %>%
  identify_outliers(Legionella.count)

outliers_biofilm_springs <- sample.metadatabiofilm%>%
  group_by(origin) %>%
  identify_outliers(Legionella.count)

outliers_water_date <- sample.metadatawater%>%
  group_by(Date) %>%
  identify_outliers(Legionella.count)

outliers_biofilm_date <- sample.metadatabiofilm%>%
  group_by(Date) %>%
  identify_outliers(Legionella.count)

#check normality by group
sample.metadatawater %>%
  group_by(origin) %>%
  shapiro_test(Legionella.count)

sample.metadatabiofilm %>%
  group_by(origin) %>%
  shapiro_test(Legionella.count)

sample.metadatawater %>%
  group_by(origin) %>%
  shapiro_test(Legionella.count.Log)

sample.metadatabiofilm %>%
  group_by(origin) %>%
  shapiro_test(Legionella.count.Log)

# Homogneity of variance assumption
sample.metadatawater %>% levene_test(Legionella.count ~ origin)
sample.metadatabiofilm %>% levene_test(Legionella.count ~ origin)

sample.metadatawater %>% levene_test(Legionella.count.Log ~ origin)
sample.metadatabiofilm %>% levene_test(Legionella.count.Log ~ origin)

# Not all the groups distribute normally, levene's test significant for some groups
# we have extreme outliers
#and in addition the sample size of haon is n=5
# So I will run kruskal wallis on the Log trasformed data

#some stats first
origin.median.water.Leg <-sample.metadatawater %>%
  group_by(origin) %>%
  get_summary_stats(Legionella.count.Log, type = "median_iqr")


origin.median.biofilm.Leg <-sample.metadatabiofilm %>%
  group_by(origin) %>%
  get_summary_stats(Legionella.count.Log, type = "median_iqr")

date.median.water.Leg <- sample.metadatawater %>%
  group_by(Date) %>%
  get_summary_stats(Legionella.count.Log, type = "median_iqr")

date.median.biofilm.Leg <- sample.metadatabiofilm %>%
  group_by(Date) %>%
  get_summary_stats(Legionella.count.Log, type = "median_iqr")


# kruskal.test, for each environment separately (water and biofilm)

kruskal.water.spring <- kruskal.test(Legionella.count.Log ~ origin,
                                     data = sample.metadatawater)

kruskal.biofilm.spring <- kruskal.test(Legionella.count.Log ~ origin,
                                       data = sample.metadatabiofilm)

kruskal.water.date <- kruskal.test(Legionella.count.Log ~ Date,
                                   data = sample.metadatawater)

kruskal.biofilm.date <- kruskal.test(Legionella.count.Log ~ Date,
                                     data = sample.metadatabiofilm)

# first we need to reorder the X axis for spring cluster and date and only then do post hoc
# so that the order of the figures and the XY position of the post hoc bars will be the same

# indicate the levels to be presented (and ordered) for the x axis (water and biofilm - they are different)
sample.metadatawater$origin <- factor(sample.metadatawater$origin , levels=c("Tiberias Hot Springs", "Haon-Borehole", "Tabgha", "Fuliya"))
sample.metadatabiofilm$origin <- factor(sample.metadatabiofilm$origin , levels=c("Tiberias Hot Springs", "Tabgha", "Fuliya"))


# indicate the levels to be presented (and ordered) for the x axis (water and biofilm - they are different)
sample.metadatawater$Date <- factor(sample.metadatawater$Date , levels=c("Sep17", "Jan18", "Jun18", "Oct18", "Jan19"))
sample.metadatabiofilm$Date <- factor(sample.metadatabiofilm$Date , levels=c("Sep17", "Jan18", "Jun18", "Oct18", "Jan19"))
str(sample.metadatabiofilm)


# Pairwise comparisons (p.adjust.method = Bonferroni-Holm) - only for origin.1 (significant) grouped by environment 

post.hoc.water.sprins <- sample.metadatawater %>%
  wilcox_test(Legionella.count.Log ~ origin, p.adjust.method = "BH")

post.hoc.water.sprins

post.hoc.biofilm.sprins <- sample.metadatabiofilm %>%
  wilcox_test(Legionella.count ~ origin, p.adjust.method = "BH")

post.hoc.biofilm.sprins


post.hoc.biofilm.date <- sample.metadatabiofilm %>%
  wilcox_test(Legionella.count ~ Date, p.adjust.method = "BH")

post.hoc.biofilm.date

write.csv(post.hoc.water.sprins, "/home/oded/Documents/post/articles/springs/qPCR/spring_cluster/water+biofilm/qPCR_Legionella_springs/statistics/posthoc-water-sprins.csv", row.names = FALSE)
write.csv(post.hoc.biofilm.sprins,"/home/oded/Documents/post/articles/springs/qPCR/spring_cluster/water+biofilm/qPCR_Legionella_springs/statistics/posthoc-biofilm-sprins.csv", row.names = FALSE)
write.csv(post.hoc.biofilm.date,"/home/oded/Documents/post/articles/springs/qPCR/spring_cluster/water+biofilm/qPCR_Legionella_springs/statistics/posthoc-biofilm-date.csv", row.names = FALSE)

# add add_xy_position to add later to the graph

post.hoc.biofilm.sprins <- post.hoc.biofilm.sprins %>% add_xy_position
post.hoc.water.sprins <- post.hoc.water.sprins %>% add_xy_position
post.hoc.biofilm.date <- post.hoc.biofilm.date %>% add_xy_position


library(ggplot2)
library(ggpattern)
library(dplyr)
library(forcats)
library(ggpubr)
library(tidyverse)
library(extrafont)
font_import()
fonts()




# creat labeles for the x axis (water and biofilm - they are different)
my.labels <- c("THS", "Haon-Borehole", "Tabgha", "Fuliya")
my.labels1 <- c("THS", "Tabgha", "Fuliya")

# we will build two box plots (water, biofilm) and combine the using - grid.arrange
bxpp_water.spring <- ggboxplot(sample.metadatawater, x = "origin", y = "Legionella.count.Log", fill = "origin", width = 0.5) +
  scale_x_discrete(labels= my.labels) + scale_y_continuous(limits=c(0, 16)) +
  theme(legend.position = "none") + labs(y=expression(paste(italic("Legionella")," [Log10 (Cells/ml)]"))) +
  theme(plot.title = element_text(face="bold", size=10, hjust = 0.5))


bxpp_water.spring.text <- bxpp_water.spring + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=14, vjust = 2.5),
  axis.text.x = element_text(color="black", size=14, vjust = -2.5),
  axis.text.y = element_text(color="black", size=14),
  legend.title = element_text(color = "black", size = 14, face ="bold"),
  legend.text = element_text(color = "black", size = 14)) +
theme(text=element_text(family="ArialMT"))


bxpp_water.spring.text.col <- bxpp_water.spring.text + scale_fill_manual(values = c("#CC79A7", "#009E73", "#56B4E9", "#D55E00"))

bxpp_water.spring.text.col.stat <- bxpp_water.spring.text.col + stat_pvalue_manual(post.hoc.water.sprins, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.02, y.position = c(15, 16, 13, 14, 12))
bxpp_water.spring.text.col.stat.mar <- bxpp_water.spring.text.col.stat +
  theme(plot.margin = unit(c(25,4,25,8), "pt"))

bxpp_biofilm.spring <- ggboxplot(sample.metadatabiofilm, x = "origin", y = "Legionella.count.Log", fill = "origin", width = 0.35) +
  scale_x_discrete(labels= my.labels1) + scale_y_continuous(limits=c(0, 16)) +
  theme(legend.position = "none") + labs(y=expression(paste(italic("Legionella")," [Log10 (Cells/ml)]"))) +
  theme(plot.title = element_text(face="bold", size=10, hjust = 0.5))

bxpp_biofilm.spring.text <- bxpp_biofilm.spring + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=14, vjust = 2.5),
  axis.text.x = element_text(color="black", size=14, vjust = -2.5),
  axis.text.y = element_text(color="black", size=14),
  legend.title = element_text(color = "black", size = 14, face ="bold"),
  legend.text = element_text(color = "black", size = 14)) +
  theme(text=element_text(family="ArialMT"))

bxpp_biofilm.spring.text.col <- bxpp_biofilm.spring.text + scale_fill_manual(values = c("#CC79A7", "#56B4E9", "#D55E00"))

bxpp_biofilm.spring.text.col.stat <- bxpp_biofilm.spring.text.col + stat_pvalue_manual(post.hoc.biofilm.sprins, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.02, y.position = c(15, 16)) 
bxpp_biofilm.spring.text.col.stat

bxpp_biofilm.spring.text.col.stat.mar <- bxpp_biofilm.spring.text.col.stat +
  theme(plot.margin = unit(c(2,4,2,8), "pt"))

bxpp_springs_final <- cowplot::plot_grid(bxpp_water.spring.text.col.stat.mar, bxpp_biofilm.spring.text.col.stat, align = "h", ncol = 1, rel_heights = c(0.5, 0.5))
bxpp_springs_final


#a graph with the date ordered:

my.labels.date <- c("Sep17", "Jan18", "Jun18", "Oct18", "Jan19")

bxpp_water.date <- ggboxplot(sample.metadatawater, x = "Date", y = "Legionella.count.Log", width = 0.5) +
  scale_x_discrete(labels= my.labels.date) + scale_y_continuous(limits=c(0, 16)) +
  theme(legend.position = "none", axis.title.x = element_blank()) + labs(y=expression(paste(italic("Legionella")," [Log10 (Cells/ml)]"))) +
  theme(plot.title = element_text(face="bold", size=10, hjust = 0.5))

bxpp_water.date.text <- bxpp_water.date + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=14, vjust = 2.5),
  axis.text.x = element_text(color="black", size=14, vjust = -2.5),
  axis.text.y = element_text(color="black", size=14),
  legend.title = element_text(color = "black", size = 14, face ="bold"),
  legend.text = element_text(color = "black", size = 14)) +
  theme(text=element_text(family="ArialMT"))

bxpp_water.date.text

bxpp_water.date.text.mar <- bxpp_water.date.text +
  theme(plot.margin = unit(c(25,4,25,8), "pt"))

bxpp_biofilm.date <- ggboxplot(sample.metadatabiofilm, x = "Date", y = "Legionella.count.Log", width = 0.35) +
  scale_x_discrete(labels= my.labels.date) + scale_y_continuous(limits=c(0, 16)) +
  theme(legend.position = "none") + 
  labs(y=expression(paste(italic("Legionella")," [Log10 (Cells/ml)]"))) +
  theme(plot.title = element_text(face="bold", size=10, hjust = 0.5))

bxpp_biofilm.date.text <- bxpp_biofilm.date + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_text(color="black", size=14, vjust = 2.5),
  axis.text.x = element_text(color="black", size=14, vjust = -2.5),
  axis.text.y = element_text(color="black", size=14),
  legend.title = element_text(color = "black", size = 14, face ="bold"),
  legend.text = element_text(color = "black", size = 14)) +
  theme(text=element_text(family="ArialMT"))

bxpp_biofilm.date.text.stat <- bxpp_biofilm.date.text + stat_pvalue_manual(post.hoc.biofilm.date, label = "p.adj.signif", inherit.aes = FALSE, hide.ns=TRUE, tip.length = 0.02,  y.position = c(16, 15, 14, 13))
bxpp_biofilm.date.text.stat

bxpp_biofilm.date.text.stat.mar <- bxpp_biofilm.date.text.stat +
  theme(plot.margin = unit(c(2,4,2,8), "pt"))

bxpp_date_final <- cowplot::plot_grid(bxpp_water.date.text.mar, bxpp_biofilm.date.text.stat, align = "h", ncol = 1, rel_heights = c(0.5, 0.5))
bxpp_date_final


bxpp_final <- cowplot::plot_grid(bxpp_springs_final, bxpp_date_final, align = "v", ncol = 2, rel_heights = c(0.5, 0.5))
bxpp_final



pdf("/home/oded/Documents/post/articles/springs/final R code/qPCR/23.06.20_qPCR_springs_cluster_collection_date.pdf",  width = 15, height = 10, family="ArialMT")
bxpp_final
dev.off()



# create the legend for figures

# top legend

bxpp_water.spring.text.col.top <- bxpp_water.spring.text + scale_fill_manual(values = c("#CC79A7", "#009E73", "#56B4E9", "#D55E00")) + theme(legend.position = "top", legend.spacing.x = unit(0.3, 'cm')) + labs(fill='Spring Cluster:  ') 

leg_top <- get_legend(bxpp_water.spring.text.col.top)
as_ggplot(leg_top)

pdf("/home/oded/Documents/post/articles/springs/final R code/legend/Legend_springs_top.pdf",  width = 6, height = 1, family="ArialMTMT")

dev.off()  


# right legend

bxpp_water.spring.text.col.right <- bxpp_water.spring.text + scale_fill_manual(values = c("#CC79A7", "#009E73", "#56B4E9", "#D55E00")) + theme(legend.position = "right", legend.spacing.x = unit(0.3, 'cm'), legend.spacing.y = unit(0.2, 'cm')) + guides(fill = guide_legend(byrow = TRUE)) + labs(fill='Spring Cluster:  ') 

leg_right <- get_legend(bxpp_water.spring.text.col.right)

pdf("/home/oded/Documents/post/articles/springs/final R code/legend/Legend_springs_right.pdf",  width = 1.5, height = 2, family="ArialMT")
as_ggplot(leg_right)
dev.off()
