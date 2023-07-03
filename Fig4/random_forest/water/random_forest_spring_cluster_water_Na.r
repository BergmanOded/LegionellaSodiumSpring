# https://github.com/LangilleLab/microbiome_helper/wiki/Random-Forest-Tutorial
library(forcats)
library(randomForest)
library(plyr) # for the "arrange" function
library(rfUtilities) # to test model significance
library(caret) 
# to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function
library(dplyr)
library(reshape)
library(pheatmap)
library(pls)
# the caret package is not necessary for large sample size (50+ ?)

otu_table <- read.table("/home/oded/Documents/post/articles/springs/random_forest/water20/last5/feature-table.tsv", sep = "\t", header = T, row.names = 1, skip = 1, comment.char = "")  
metadata <- read.table("/home/oded/Documents/post/articles/springs/random_forest/water20/last5/Tabgha-Fuliya-combined/sample-metadata.tsv", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

# a basic explor of the data files (how many samples and ASVs/OTUs we have...)
dim(otu_table)
dim(metadata)
str(metadata)

# Pre-processing is one of the most important steps in machine learning.
# RF doesn't make any assumptions about how the data is distributed,
# so it often is not necessary to transform your data. 
# However, reducing noise in your input data will improve model performance. 
# An easy way to do this is to throw out features that are rare or have very low variance, across samples.
# Any cut-offs used at this step are slightly arbitrary and would depend on the dataset.

otu_nonzero_counts <- apply(otu_table, 1, function(y) sum(length(which(y > 0))))

# if margins are to large
par("mar")
par(mar=c(3,3,3,3))
# and plot
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

# For the spring article data we will use the filtered 2000 feature table. 
# We might also be interested in excluding OTUs that don't vary much across samples. 
# This usually isn't different from excluding OTUs based on their number of non-zero values.
# However, this can be useful for other types of feature tables. We can use nearZeroVar() function from the caret R package to perform this filtering.
# but in this tutotial we will filter features that are non-zero in less than a specified proportion:
# next two lines from (https://rstudio-pubs-static.s3.amazonaws.com/115631_7397b7cf67534479ae80f70546610eea.html#libraries)
# Random forests can handle sparse matrices, but we still want to prune out lots of our rare OTUs 
# which are just contributing noise. We will do this by eliminating OTUs with an average relative abundance below 0.0001


remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

# Based on the above histogram remove OTUs that have non-zero values in <= 10% of samples:
otu_table_rare_removed <- remove_rare(table=otu_table, cutoff_pro=0.0001)

# 254 OTUs are remaining,
dim(otu_table_rare_removed)
str(otu_table_rare_removed)
#To re-normalize_table so that each sample's column sums to 100 
#(you can use the ColSums function afterwards for a sanity check):

otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100


# Transforming Your Data (see more options at thhe end)

# One approach: subtracting each sample's mean (center) and then dividing by the sample's standard deviation (scale). 
# So centre-scale is like converting into a Z-score. 

otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE) 

# Running Model

# the metadata column will be the last column of each input dataframe
# Make a dataframe of training data with OTUs as column and samples as rows (transposed.)
# so we ca then add the state column to the end

otu_table_scaled_cluster <- data.frame(t(otu_table_scaled), check.names = FALSE)  
otu_table_scaled_cluster$origin_new <- metadata[rownames(otu_table_scaled_cluster), "origin_new"]  

# Now prep input tables for regression of inflammation score (Na):
# same thing with the Na

otu_table_scaled_Na <- data.frame(t(otu_table_scaled), check.names = FALSE)  
otu_table_scaled_Na$Na <- metadata[rownames(otu_table_scaled_Na), "Na"] 

# The seed number you choose is the starting point used in the generation 
# of a sequence of random numbers.
# So before running the models you should set the random seed so that they will be 
# reproducible.
set.seed(151)  

# The 2 parameters for a RF model are the number of trees in the forest (ntree) 
# and the number of features randomly sampled at each node in a tree (mtry). 
# Unless you have a reason to change mtry beforehand (or can optimize it with an independent partition of data) it's better to use the default values. The default mtry values differ for RF classification and regression.
# mtry in the below code is defult and thus does not appear
# The more trees you run in your forest the better the model will converge.
# here used 501 (computation time), which is similar to the default, 
# but in practice you may want to use something like 10,001 trees for a robust model 
# (depending on the computational time). 
# choose odd numbers of trees to ensure there are never any ties for binary classification models (have no idea what this means)




# Run RF to classify inflamed and control samples:
RF_cluster_classify <- randomForest( x=otu_table_scaled_cluster[,1:(ncol(otu_table_scaled_cluster)-1)] , y=otu_table_scaled_cluster[ , ncol(otu_table_scaled_cluster)] , ntree=1001, importance=TRUE, proximities=TRUE)
RF_cluster_classify


set.seed(998)
inTraining <- createDataPartition(otu_table_scaled_cluster$origin_new, p = .75, list = FALSE)
training <- otu_table_scaled_cluster[ inTraining,]
testing  <- otu_table_scaled_cluster[-inTraining,]
str(inTraining)
nrow(training)
nrow(testing)



fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

cluster_Fit <- train(
  origin_new ~ .,
  data = training,
  trControl = fitControl,
  method = "pls",
  ## Center and scale the predictors for the training
  ## set and all future samples.
  preProc = c("center", "scale"),
  tuneLength = 15
)



# Run RF to regress OTUs against inflammation score (IS):

RF_Na_regress <- randomForest( x=otu_table_scaled_Na[,1:(ncol(otu_table_scaled_Na)-1)] , y=otu_table_scaled_Na[ , ncol(otu_table_scaled_Na)] , ntree=1001, importance=TRUE, proximities=TRUE, ncomp=12)  
RF_Na_regress

# Assessing Model Fit
# Remember that if you partition your data into a training and test set beforehand 
# (which we didn't do in this case due to low sample sizes), you would be able to directly
# validate your trained model on independent data.
# Without this partitioning the simplest way to assess model fit is to take a 
# look at the model summaries by just typing the model name.

# The out-of-bag error rate is 0%, which is perfect! 
# The summary also gives us the mtry parameter, which was 25% in this case.

# We can also take a look at the regression RF summary
plot(RF_Na_regress)
# The mtry parameter is much higher (215/501) and that instead of out-of-bag error 
# and a confusion matrix, we can use the mean of squared residuals 
# and the % variance explained as performance metrics.
# This model is also performing extremely well but not as the the classification model.


# Permutation Test
# Although the above performance metrics provided by the models were encouraging, 
# often it isn't so clear that a model is working well. Testing whether your model's 
# performance metric is more extreme than expected by chance based on a permutation 
# test is one method to assess model significance.
# To run these significance tests with 1000 permutations:

  RF_cluster_classify_sig <- rf.significance( x=RF_cluster_classify ,  xdata=otu_table_scaled_cluster[,1:(ncol(otu_table_scaled_cluster)-1)] , nperm=1000 , ntree=1001 )  
  RF_cluster_classify_sig
  RF_Na_regress_sig <- rf.significance( x=RF_Na_regress ,  xdata=otu_table_scaled_Na[,1:(ncol(otu_table_scaled_Na)-1)] , nperm=1000 , ntree=1001 )  
  RF_Na_regress_sig
# Model OOB error = randomforest model error rate
# Random OOB error = random permutations  error rate 
# The % Var explained is a measure of how well out-of-bag predictions explain the target variance of the training set.

# Accuracy Estimated by Cross-validation
# One method to estimate model performance is to systematically partition the data 
# into training and test sets and repeatedly see how the model performs (cross-validation).
# The simplest cross-validation leave-one-out. The model is trained n times,
# where n is the number of samples, and one sample is left out each time for testing. 
# This provides estimates of model performance, which tend to be quite similar to the 
# internal measures (out-of-bag error and mean of squared residuals) described above. 
# We will run this cross-validation using the caret R package.

# The first step is defining the parameters we want to use for training:
# "LOOCV" = leave-one-out cross-validation
# caret can be used for much more sophisticated purposes.

fit_control <- trainControl( method = "LOOCV" )    

# These commands would run the leave-one-out cross-validation, 
# note the same ntree and mtry parameters are being used as above.

RF_cluster_classify_loocv <- train( otu_table_scaled_cluster[,1:(ncol(otu_table_scaled_cluster)-1)] , y=otu_table_scaled_cluster[, ncol(otu_table_scaled_cluster)] , method="rf", ntree=1001 , tuneGrid=data.frame( mtry=30) , trControl=fit_control )

RF_Na_regress_loocv <- train( otu_table_scaled_Na[,1:(ncol(otu_table_scaled_Na)-1)] , y=otu_table_scaled_Na[, ncol(otu_table_scaled_Na)] , method="rf", ntree=1001 , tuneGrid=data.frame( mtry=300) , trControl=fit_control )

# To see the performance metrics:

RF_cluster_classify_loocv$results    
RF_Na_regress_loocv$results  


#20 most important

imp <- importance(RF_Na_regress)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(X.IncMSE))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 20-50 predictors
imp.20 <- imp.sort[1:20, ]
imp.30 <- imp.sort[1:30, ]
imp.40 <- imp.sort[1:40, ]
imp.50 <- imp.sort[1:50, ]

imp.20$IncMSE = imp.20$X.IncMSE
imp.20 <- imp.20[, -2 ]

imp.30$IncMSE = imp.30$X.IncMSE
imp.30 <- imp.30[, -2 ]

imp.40$IncMSE = imp.40$X.IncMSE
imp.40 <- imp.40[, -2 ]

imp.50$IncMSE = imp.50$X.IncMSE
imp.50 <- imp.50[, -2 ]



write.csv(imp.50,"/home/oded/Documents/post/articles/springs/random_forest/water20/last5/heatmap/imp.50.csv", row.names = FALSE)
write.csv(imp.40,"/home/oded/Documents/post/articles/springs/random_forest/water20/last5/heatmap/imp.40.csv", row.names = FALSE)
write.csv(imp.30,"/home/oded/Documents/post/articles/springs/random_forest/water20/last5/heatmap/imp.30.csv", row.names = FALSE)
write.csv(imp.20,"/home/oded/Documents/post/articles/springs/random_forest/water20/last5/heatmap/imp.20.csv", row.names = FALSE)
  

