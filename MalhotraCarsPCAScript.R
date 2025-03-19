# Varenyam Malhotra
# April 15, 2023

# carspca.R
# calculates PCAs for 2004 car data
#
# set directory and clean up
setwd("/Users/vena/desktop/RFolder")
rm(list=ls())
#
# load the ggplot2 library
library(ggplot2)
#
# read in the data
data <- read.csv(file="cars-fixed04.dat.txt", stringsAsFactors=FALSE, header=TRUE)
head(data)
#
# filter the data
data <- data[, 8:ncol(data)]
head(data)
# Performing the PCA using prcomp
pca <- prcomp(data, scale=TRUE)
pca
summary(pca)
names(pca)
# Showing the x-values for all rows of the first column
pca$x[,1] 
# Showing the x-values for all rows of the second column
pca$x[,2] 
# Plotting these two values against each other
plot(pca$x[,1],pca$x[,2])
# Barplot based on variances
screeplot(pca, main="Scree Plot")
# Calculating actual PCAs based on variances and standard deviations, barplot distribution graph
pca.variance <- pca$sdev^2
pca.variance.per <- round(pca.variance/sum(pca.variance)*100,1)
pca.variance.per
barplot(pca.variance.per, main="Scree Plot", xlab="Principal Components", ylab="Percent Variable", ylim=c(0,100),names.arg=pca.variance.per)
# Now creating a new dataframe based on PCA values from two columns' rows
pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
pca.data
# And now using ggplot to plot the pca data from all rows of column 1 and column 2
ggplot(data=pca.data,aes(x=X,y=Y,label=Sample)) + geom_text() + xlab(paste("PC1 - ", pca.variance.per[1], "%", sep="")) + ylab(paste("PC2 - ", pca.variance.per[2], "%", sep="")) + theme_bw() + ggtitle("My PCA Plot")
# To figure out what is contributing to the pca values for column 1 and column 2 (loading scores)
loading_scores <- pca$rotation[,1]
loading_scores
# just talking absolute values to make all loading scores positive
car_scores <- abs(loading_scores)
car_scores
# the scores ranked from high to low (most to least contribution)
car_scores_ranked <- sort(car_scores, decreasing = TRUE)
car_scores_ranked
# top 10 contributors
top10 <- names(car_scores_ranked[1:10])
top10
# Now plotting pca as a 2D dataset
biplot(pca)

