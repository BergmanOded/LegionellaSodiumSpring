#PCA total, biotic, abiotic - Legionella project.
#All variables numeric

#PCA water total data----
springs.metadata.water <- read.table("/home/oded/Documents/post/articles/springs/PCA/water/water-metadata.tsv", sep=",", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

str(springs.metadata.water)
springs.metadata.water$Na<-as.numeric(springs.metadata.water$Na)
springs.metadata.water$SO4<-as.numeric(springs.metadata.water$SO4)
springs.metadata.water$Conductivity<-as.numeric(springs.metadata.water$Conductivity)
springs.metadata.water$Legionella.count<-as.numeric(springs.metadata.water$Legionella.count)
springs.metadata.water$Origin.numeric<-as.numeric(springs.metadata.water$Origin.numeric)

# If we want to transform Leg cell number to three levels - low, medium, high
#library(dplyr)
#springs.metadata <- mutate(springs.metadata, Legionella.count.level = ifelse (Legionella.count %in% 0:1500, "low",
                                                                             # ifelse (Legionella.count %in% 1501:5000, "medium", "high")))

#log transformation (+1) - so no inf values+ combine the columns together to one file 
#(the numeric columns in the log file and the factor columns from the original file)
  log.springs.metadata <- log(1+springs.metadata.water[c(11:23)])
str(log.springs.metadata)
cbind.log.springs.metadata <-cbind(log.springs.metadata, springs.metadata.water[, c(1:10)])
shapiro.test(cbind.log.springs.metadata$Legionella.count.Log10)
str(cbind.log.springs.metadata)
# now we remove missing values and inf: 
cbind.log.springs.metadata.NA <- cbind.log.springs.metadata[complete.cases(cbind.log.springs.metadata), ]
str(cbind.log.springs.metadata.NA)


# first step is to identify the relevance of all the vectors
# 1. PCA initial vector analysis of all data.

# using prcomp on all the variables, we get the explained variencen for each PC
cbind.log.springs.metadata.NA.pca <- prcomp(cbind.log.springs.metadata.NA[,c(1:13)], center = TRUE,scale = TRUE)
summary(cbind.log.springs.metadata.NA.pca)

# initial graphs to see that the data spreads evenly throught
# the second graph is to look at the loadings. this will tell us what variables influence the PCA the most 
library(ggbiplot)
library(viridis)
plot(cbind.log.springs.metadata.NA.pca$x[,1], cbind.log.springs.metadata.NA.pca$x[,2])
loading.pcoa.water <- ggbiplot(cbind.log.springs.metadata.NA.pca, scale=0) + theme_classic()
# Now we look at the Rotation or the eigenvectors - these specify the orientation of the principal components relative to the original variables. 
#The elements of an eigenvector, that is, the values within a particular row of matrix A, are the weights aij. These values are called the loadings, 
#and they describe how much each variable contributes to a particular principal component.
#Large loadings (positive or negative) indicate that a particular variable has a strong relationship to a particular principal component
cbind.log.springs.metadata.NA.pca

# in addition we can see which principal componants are relevant. we square the sdev of each
# values > 1 are considered 
PC_sdev <- cbind.log.springs.metadata.NA.pca$sdev^2
PC_sdev

# plot all data
PCA_total <-cbind(cbind.log.springs.metadata.NA, cbind.log.springs.metadata.NA.pca$x[,1:3])
str(PCA_total)
summary(cbind.log.springs.metadata.NA.pca)
cbind.log.springs.metadata.NA.pca

library(ggplot2)
library(viridis)
theme<-theme(panel.background = element_blank(), axis.line = element_line(size = 0.5, colour = "black"), plot.margin=unit(c(1,1,1,1),"line"), legend.position="top", legend.box = "horizontal", legend.key = element_blank())
PcoA.water <- ggplot(PCA_total, aes(PC1,PC2, colour=origin, loadings = TRUE))+
  xlab("PC 1 (51.1%)") + ylab("PC 2 (18.7%)") +  
  #geom_text(label=rownames(cbind.log.PCA_station_total_03.08.20.NA), nudge_x = 0.5, nudge_y = 0.5, size=3) +
  geom_point(aes(size=Legionella.count.Log10, alpha=Legionella.count.Log10))+
  scale_size_continuous(range = c(-10,10))+
  scale_alpha_continuous(range=c(-0.2,1))+
  scale_colour_manual(values = c("#D55E00", "#009E73", "#56B4E9", "#CC79A7"))+
  #scale_shape_manual(values = c(0, 2, 5, 15, 17, 18))+
  geom_hline(yintercept=0,linetype="dashed", color = "red", size=0.5, alpha = 0.5) +
  geom_vline(xintercept=0,linetype="dashed", color = "red", size=0.5, alpha = 0.5) +
theme +
  theme(legend.key.size = unit(3,"line")) + labs(colour = "Spring cluster")

PcoA.water.nolegend <- PcoA.water + theme(legend.position="none")
PcoA.water.nolegend


#Setting aside variables known to be strongly correlated with others can substantially alter PCA results. 
#Therefore, there may be merit in discarding variables  beforehand, if you think they may be measuring the same underlying ecological aspect.
#If they only correlate, no need to discard

#calculate correlation matrix with p value - maybe only for relevant ones - at specific depth/year/stratification
# bsaed on normalization and n, choose pearson / spearman
library(Hmisc)

cor_rem <- cbind.log.springs.metadata.NA[, c(1:13)]
print(cor_rem)
cor_pca_Na_conduct <- rcorr(as.matrix(cor_rem), type=c("spearman"))
print(cor_pca_Na_conduct)

# make the presentation "flattened"
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
cor_total_PCA <- flattenCorrMatrix(cor_pca_Na_conduct$r, cor_pca_Na_conduct$P)

#plot the correlation between variables (include histogram)
library("PerformanceAnalytics")
library(GGally)
cor_pca_Na_conduct_plot <- cbind.log.springs.metadata.NA[, c(1:13)]
PcoA.water.correlation <- ggpairs(cor_pca_Na_conduct_plot[1:13], lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.1)))



#PCA water after omiting redundant chemical variables----

log.springs.metadata <- log(1+springs.metadata.water[c(11:23)])
str(log.springs.metadata)
shapiro.test(log.springs.metadata$Legionella.count)
cbind.log.springs.metadata <-cbind(log.springs.metadata, springs.metadata.water[, c(1:8, 10)])
str(cbind.log.springs.metadata)

#remove missing values and inf: 
cbind.log.springs.metadata.NA <- cbind.log.springs.metadata[complete.cases(cbind.log.springs.metadata), ]


# using prcomp of the variables we want

cbind.log.springs.metadata.rem.NA.pca <- prcomp(cbind.log.springs.metadata.NA[,c(1:3, 5, 8:11)], center = TRUE,scale = TRUE)
summary(cbind.log.springs.metadata.rem.NA.pca)

# initial graphs

plot(cbind.log.springs.metadata.rem.NA.pca$x[,1], cbind.log.springs.metadata.rem.NA.pca$x[,2])
loading.pcoa.water.red <- ggbiplot(cbind.log.springs.metadata.rem.NA.pca, scale=0) + theme_classic()
loading.pcoa.water.red
# Now for the PCA plot

  PCA_total.red <-cbind(cbind.log.springs.metadata.NA, cbind.log.springs.metadata.rem.NA.pca$x[,1:3])
head(PCA_total.red)
summary(cbind.log.springs.metadata.rem.NA.pca)
cbind.log.springs.metadata.rem.NA.pca
library(ggplot2)
library(viridis)
theme<-theme(panel.background = element_blank(), axis.line = element_line(size = 0.5, colour = "black"), plot.margin=unit(c(1,1,1,1),"line"), legend.position="top",  legend.key = element_blank())
PcoA.water.red <- ggplot(PCA_total.red, aes(PC1,PC2, colour = origin, loadings = TRUE, label = Station))+
  xlab("PC 1 (55.6%)") + ylab("PC 2 (22.1%)") +  
  #geom_text(label=rownames(cbind.log.PCA_station_total_03.08.20.NA), nudge_x = 0.5, nudge_y = 0.5, size=3) +
  geom_point(aes(size=Legionella.count.Log10, alpha=Legionella.count.Log10))+
  scale_size_continuous(range = c(-10,10))+
  scale_alpha_continuous(range=c(-0.2,1))+
  scale_colour_manual(values = c("#D55E00", "#009E73", "#56B4E9", "#CC79A7"))+
  #scale_shape_manual(values = c(0, 2, 5, 15, 17, 18))+
  geom_hline(yintercept=0,linetype="dashed", color = "red", size=0.5, alpha = 0.5) +
  geom_vline(xintercept=0,linetype="dashed", color = "red", size=0.5, alpha = 0.5) +
  theme +
  theme(legend.key.size = unit(3,"line")) + labs(colour = "Spring cluster")

PcoA.water.red

PcoA.water.red.nolegend <- PcoA.water.red + theme(legend.position="none")


 
#Setting aside variables known to be strongly correlated with others can substantially alter PCA results. 
#Therefore, there may be merit in discarding variables  beforehand, if you think they may be measuring the same underlying ecological aspect.
#If they only correlate, no need to discard
#calculate correlation matrix with p value - maybe only for relevant ones - at specific depth/year/stratification----
# bsaed on normalization and n, choose pearson / spearman
library(Hmisc)

cor_rem <- cbind.log.springs.metadata.NA[, c(1:3, 5, 8:11)]
print(cor_rem)
cor_pca_Na_conduct <- rcorr(as.matrix(cor_rem), type=c("spearman"))

# print the correlation and sig.
print(cor_pca_Na_conduct$r)
signif(cor_pca_Na_conduct$P,2)

# make the presentation "flattened"  in 4 column tables
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
correlation_var_water_PCA <- flattenCorrMatrix(cor_pca_Na_conduct$r, cor_pca_Na_conduct$P)

#plot the correlation between variables (include histogram)
library("PerformanceAnalytics")
chart.Correlation(cor_rem, histogram=TRUE, pch=19, method = "spearman")


#multiple regression (basic)

model <- lm(Legionella.count ~ Turbidity + Temp + Ph  + TP + Na + Nitrate + Fe + SO4, data = cbind.log.springs.metadata.NA)
summary(model)

model1 <- lm(Legionella.count ~ Temp + Ph  + TP+ Na + Fe, data = cbind.log.springs.metadata.NA)
summary(model1)

model2 <- lm(Legionella.count ~ Temp  + TP+ Na + Fe, data = cbind.log.springs.metadata.NA)
summary(model2)


# test multicollinearity (VIF)
library(car) 
car::vif(model)
car::vif(model1)
car::vif(model2)





#PCA biofilm total data----


springs.metadata.biofilm <- read.table("/home/oded/Documents/post/articles/springs/PCA/biofilm/biofilm-metadata.tsv", sep=",", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

str(springs.metadata.biofilm)
springs.metadata.biofilm$Na<-as.numeric(springs.metadata.biofilm$Na)
springs.metadata.biofilm$SO4<-as.numeric(springs.metadata.biofilm$SO4)
springs.metadata.biofilm$Conductivity<-as.numeric(springs.metadata.biofilm$Conductivity)
springs.metadata.biofilm$Year<-as.factor(springs.metadata.biofilm$Year)
springs.metadata.biofilm$Month<-as.factor(springs.metadata.biofilm$Month)
springs.metadata.biofilm$Legionella.count<-as.numeric(springs.metadata.biofilm$Legionella.count)
springs.metadata.biofilm$Origin.numeric<-as.numeric(springs.metadata.biofilm$Origin.numeric)
# transform Leg cell number to three levels - low, medium, high
#library(dplyr)
#springs.metadata <- mutate(springs.metadata, Legionella.count.level = ifelse (Legionella.count %in% 0:1500, "low",
#                                                                             ifelse (Legionella.count %in% 1501:5000, "medium", "high")))

#log transformation (+1) - so no inf values+ combine the columns together to one file 
#(the numeric columns in the log file and the factor columns from the original file)
log.springs.metadata.biofilm <- log(1+springs.metadata.biofilm[c(11:23)])
str(log.springs.metadata.biofilm)
cbind.log.springs.metadata.biofilm <-cbind(log.springs.metadata.biofilm, springs.metadata.biofilm[, c(1:10, 24)])
shapiro.test(cbind.log.springs.metadata.biofilm$Legionella.count.Log10)
str(cbind.log.springs.metadata.biofilm)
# now we remove missing values and inf: 
cbind.log.springs.metadata.NA.biofilm <- cbind.log.springs.metadata.biofilm[complete.cases(cbind.log.springs.metadata.biofilm), ]
str(cbind.log.springs.metadata.NA.biofilm)
#PCA total data----

# first step is to identify the relevance of all the vectors
# 1. PCA initial vector analysis of all data.

# using prcomp on all the variables, we get the explained variencen for each PC
cbind.log.springs.metadata.NA.pca.biofilm <- prcomp(cbind.log.springs.metadata.NA.biofilm[,c(1:13)], center = TRUE,scale = TRUE)
summary(cbind.log.springs.metadata.NA.pca.biofilm)

# initial graphs to see that the data spreads evenly throught
# the second graph is to look at the loadings. this will tell us what variables influence the PCA the most 
library(devtools)
library(ggbiplot)
plot(cbind.log.springs.metadata.NA.pca.biofilm$x[,1], cbind.log.springs.metadata.NA.pca.biofilm$x[,2])
loading.pcoa.biofilm <- ggbiplot(cbind.log.springs.metadata.NA.pca.biofilm, scale=0) + theme_classic()
loading.pcoa.biofilm
# Now we look at the Rotation or the eigenvectors - these specify the orientation of the principal components relative to the original variables. 
#The elements of an eigenvector, that is, the values within a particular row of matrix A, are the weights aij. These values are called the loadings, 
#and they describe how much each variable contributes to a particular principal component.
#Large loadings (positive or negative) indicate that a particular variable has a strong relationship to a particular principal component
cbind.log.springs.metadata.NA.pca.biofilm

# in addition we can see which principal componants are relevant. we square the sdev of each
# values > 1 are considered 
PC_sdev.biofilm <- cbind.log.springs.metadata.NA.pca.biofilm$sdev^2
PC_sdev.biofilm

# plot all data
PCA_total.biofilm <-cbind(cbind.log.springs.metadata.NA.biofilm, cbind.log.springs.metadata.NA.pca.biofilm$x[,1:3])
str(PCA_total.biofilm)
summary(cbind.log.springs.metadata.NA.pca.biofilm)
cbind.log.springs.metadata.NA.pca.biofilm

library(ggplot2)
library(viridis)
theme<-theme(panel.background = element_blank(), axis.line = element_line(size = 0.5, colour = "black"), plot.margin=unit(c(1,1,1,1),"line"), legend.position="top", legend.box = "horizontal",  legend.key = element_blank())
PcoA.biofilm <- ggplot(PCA_total.biofilm, aes(PC1,PC2, colour=origin, loadings = TRUE))+
  xlab("PC 1 (57.1%)") + ylab("PC 2 (13.1%)") +  
  #geom_text(label=rownames(cbind.log.PCA_station_total_03.08.20.NA), nudge_x = 0.5, nudge_y = 0.5, size=3) +
  geom_point(aes(size=Legionella.count.Log10, alpha=Legionella.count.Log10, ))+
  scale_size_continuous(range = c(-10,10))+
  scale_alpha_continuous(range=c(-0.2,1))+
  scale_colour_manual(values = c("#D55E00", "#56B4E9", "#CC79A7"))+
  #scale_shape_manual(values = c(0, 2, 5, 15, 17, 18))+
  geom_hline(yintercept=0,linetype="dashed", color = "red", size=0.5, alpha = 0.5) +
  geom_vline(xintercept=0,linetype="dashed", color = "red", size=0.5, alpha = 0.5) + 
   theme +
  theme(legend.key.size = unit(3,"line")) + labs(colour = "Spring cluster")

PcoA.biofilm <- PcoA.biofilm + theme(legend.title = element_text(size=10, face="bold"))


PcoA.biofilm.nolegend <- PcoA.biofilm + theme(legend.position="none") 
PcoA.biofilm.nolegend
#Setting aside variables known to be strongly correlated with others can substantially alter PCA results. 
#Therefore, there may be merit in discarding variables  beforehand, if you think they may be measuring the same underlying ecological aspect.
#If they only correlate, no need to discard

#calculate correlation matrix with p value - maybe only for relevant ones - at specific depth/year/stratification----
# bsaed on normalization and n, choose pearson / spearman
library(Hmisc)

cor_rem <- cbind.log.springs.metadata.NA.biofilm[, c(1:13)]
print(cor_rem)
cor_pca_Na_conduct <- rcorr(as.matrix(cor_rem), type=c("spearman"))
print(cor_pca_Na_conduct)

# make the presentation "flattened"
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
correlation_var_biofilm_PCA <- flattenCorrMatrix(cor_pca_Na_conduct$r, cor_pca_Na_conduct$P)

#plot the correlation between variables (include histogram)
library("PerformanceAnalytics")
library(GGally)
cor_pca_Na_conduct_plot <- cbind.log.springs.metadata.NA[, c(2:14)]
PcoA.biofilm.correlation <- ggpairs(cor_pca_Na_conduct_plot[1:13], lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.1)))

#PCA after omiting redundant chemical variables----

log.springs.metadata.biofilm <- log(1+springs.metadata.biofilm[c(11:23)])
str(log.springs.metadata.biofilm)
cbind.log.springs.metadata.biofilm <-cbind(log.springs.metadata.biofilm, springs.metadata.biofilm[, c(1:10)])
shapiro.test(cbind.log.springs.metadata.biofilm$Legionella.count.Log10)
str(cbind.log.springs.metadata.biofilm)

#remove missing values and inf: 
cbind.log.springs.metadata.NA.biofilm <- cbind.log.springs.metadata.biofilm[complete.cases(cbind.log.springs.metadata.biofilm), ]


# using prcomp of the variables we want

cbind.log.springs.metadata.rem.NA.pca.biofilm <- prcomp(cbind.log.springs.metadata.NA.biofilm[,c(1:3, 5, 8:11)], center = TRUE,scale = TRUE)
summary(cbind.log.springs.metadata.rem.NA.pca.biofilm)

  # initial graphs

plot(cbind.log.springs.metadata.rem.NA.pca.biofilm$x[,1], cbind.log.springs.metadata.rem.NA.pca.biofilm$x[,2])
loading.pcoa.biofilm.red <- ggbiplot(cbind.log.springs.metadata.rem.NA.pca.biofilm, scale=0) +theme_classic()
loading.pcoa.biofilm.red

# Now for the PCA plot

PCA_biofilm.red <-cbind(cbind.log.springs.metadata.NA.biofilm, cbind.log.springs.metadata.rem.NA.pca.biofilm$x[,1:3])
head(PCA_biofilm.red)
summary(cbind.log.springs.metadata.rem.NA.pca.biofilm)
cbind.log.springs.metadata.rem.NA.pca.biofilm
library(ggplot2)
library(viridis)
theme<-theme(panel.background = element_blank(), axis.line = element_line(size = 0.5, colour = "black"), plot.margin=unit(c(1,1,1,1),"line"), legend.position="top", legend.box = "horizontal",  legend.key = element_blank())
PcaA.biofilm.red <- ggplot(PCA_biofilm.red, aes(PC1,PC2, colour = origin, loadings = TRUE, label = Station))+
  xlab("PC 1 (64.1%)") + ylab("PC 2 (16%)") +  
  #geom_text(label=rownames(cbind.log.PCA_station_total_03.08.20.NA), nudge_x = 0.5, nudge_y = 0.5, size=3) +
  geom_point(aes(size=Legionella.count.Log10, alpha=Legionella.count.Log10))+
  scale_size_continuous(range = c(-10,10))+
  scale_alpha_continuous(range=c(-0.2,1))+
  scale_colour_manual(values = c("#D55E00", "#56B4E9", "#CC79A7"))+
  #geom_text()
  #stat_ellipse(type = "norm",   level = 0.95,   segments = 51) +
  #scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))+
  geom_hline(yintercept=0,linetype="dashed", color = "red", size=0.5, alpha = 0.5) +
  geom_vline(xintercept=0,linetype="dashed", color = "red", size=0.5, alpha = 0.5) +
 theme + theme(legend.key.size = unit(3,"line")) + labs(colour = "Spring cluster")

PcaA.biofilm.red

PcoA.biofilm.red.nolegend <- PcaA.biofilm.red + theme(legend.position="none")



#Setting aside variables known to be strongly correlated with others can substantially alter PCA results. 
#Therefore, there may be merit in discarding variables  beforehand, if you think they may be measuring the same underlying ecological aspect.
#If they only correlate, no need to discard
#calculate correlation matrix with p value - maybe only for relevant ones - at specific depth/year/stratification----
# bsaed on normalization and n, choose pearson / spearman
library(Hmisc)

cor_rem <- cbind.log.springs.metadata.NA.biofilm[, c(1:3, 5, 8:11)]
print(cor_rem)
cor_pca_Na_conduct <- rcorr(as.matrix(cor_rem), type=c("spearman"))
print(cor_pca_Na_conduct)

# make the presentation "flattened"
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
correlation_var_biofilm_PCA.red <- flattenCorrMatrix(cor_pca_Na_conduct$r, cor_pca_Na_conduct$P)

#plot the correlation between variables (include histogram)
library("PerformanceAnalytics")
cor_pca_Na_conduct_plot <- cbind.log.springs.metadata.NA[, c(2:4, 6, 9:12)]
chart.Correlation(cor_pca_Na_conduct_plot, histogram=TRUE, pch=19)


#multiple regression (basic)

model <- lm(Legionella.count ~ Turbidity + Temp + Ph  + TP + Na + Nitrate + Fe + SO4, data = cbind.log.springs.metadata.NA)
summary(model)

# test multicollinearity (VIF)
library(car) 
car::vif(model)


model1 <- lm(Legionella.count ~ Temp + TP + Na, data = cbind.log.springs.metadata.NA)
summary(model1)

model2 <- lm(Legionella.count ~ Temp + Na + Fe, data = cbind.log.springs.metadata.NA)
summary(model2)

car::vif(model1)
car::vif(model2)

# stepwise multiple regression
library(MASS)

# Stepwise regression model
step.model <- stepAIC(model, direction = "both", 
                      trace = FALSE)
summary(step.model)

# test multicollinearity (VIF)
library(car) 
car::vif(step.model)

# correlate Na and Legionella biofilm
cor_rem_Leg_Na <- springs.metadata[, c(9, 14, 19)]
print(cor_rem_Leg_Na)
cor_pca_Na_Legionella <- rcorr(as.matrix(cor_rem_Leg_Na), type=c("spearman"))
print(cor_pca_Na_Legionella)

# make the presentation "flattened"
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
flattenCorrMatrix(cor_pca_Na_Legionella$r, cor_pca_Na_Legionella$P)

#plot the correlation between variables (include histogram)
library("PerformanceAnalytics")
cor_pca_Na_conduct_plot <- cbind.log.springs.metadata.NA[, c(2:4, 6, 9:12)]
chart.Correlation(cor_pca_Na_conduct_plot, histogram=TRUE, pch=19)



# construct the total Figure (all data)


pcoa_combined <- plot_grid(PcoA.water.nolegend, PcoA.biofilm.nolegend, ncol = 1, labels=c("A", "C"), nrow = 2, rel_heights = c(15, 15))
loading_combined <- plot_grid(loading.pcoa.water, loading.pcoa.biofilm, ncol = 1, labels=c("B", "D"), nrow = 2)


library(cowplot)
  PCA_all_var <- ggdraw() +
  draw_plot(pcoa_combined, 0, 0, 0.5, 1) +
  draw_plot(loading_combined, 0.6, 0, 0.4, 1) 

  PCA_all_var

pdf("/home/oded/Documents/post/articles/springs/final R code/PCA/PCA_spring_chemistry_all.pdf",  width = 11, height = 8, family="ArialMT")
  
dev.off()
  
  
  
  
  
  # construct the reduced Figure 
  
  pcoa_combined.red <- plot_grid(PcoA.water.red.nolegend, PcoA.biofilm.red.nolegend, ncol = 1, nrow = 2, rel_heights = c(17, 17))
  loading_combined.red <- plot_grid(loading.pcoa.water.red, loading.pcoa.biofilm.red, ncol = 1, nrow = 2)
  
  
  library(cowplot)
  PCA_all_red <- ggdraw() +
    draw_plot(pcoa_combined.red, 0, 0, 0.5, 1) +
    draw_plot(loading_combined.red, 0.6, 0, 0.4, 1) 
  
  
  PCA_all_red
  
  pdf("/home/oded/Documents/post/articles/springs/final R code/PCA/PCA_spring_chemistry_red.pdf",  width = 11, height = 8, family="ArialMT")
  
  dev.off()
  
  
  
  
# save legend red
  
# top legend
  library(ggpubr)
  
  leg_top <- get_legend(PcaA.biofilm.red)
  as_ggplot(leg_top)
  
  
  # top legend
  library(ggpubr)
  
  PcoA.water.red.legend.top <- PcoA.water.red + theme(legend.position="top", legend.spacing.x = unit(0.3, 'cm'))
  
  leg_top <- get_legend(PcoA.water.red.legend.top) 
  as_ggplot(leg_top)
  
  pdf("/home/oded/Documents/post/articles/springs/final R code/legend/Legend_PCA_top.pdf",  width = 15, height = 1.5, family="ArialMT")
  
  dev.off()  
  
  
  PcoA.water.red.legend.right <- PcoA.water.red + theme(legend.position="right", legend.spacing.x = unit(0.3, 'cm'))
  leg_right <- get_legend(PcoA.water.red.legend.right)
  as_ggplot(leg_right)
  
  pdf("/home/oded/Documents/post/articles/springs/final R code/legend/Legend_PCA_right.pdf",  width = 2, height = 6, family="ArialMT")
  
  dev.off()  






