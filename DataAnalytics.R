# Varenyam Malhotra
# March 25, 2023

#  Robert Gotwals
#  DataAnalytics.R
#  March, 2022
#  This code reads in a CSV (comma-separated values file)
# and does TWO things:
#  cleans up numerical and character missing data.  For numerical
# data, it replaces missing numerics with the mean of the available
# numbers
#
# For character data, representative values are imputed
# using a function found in "myfunctions.R", which must be sourced
# before running this script.
#
#  Setup up the working directory (folder) and cleanup.  Also, source in
# the myfunctions.R file
rm(list=ls())  #  this removes old files that were previously loaded
setwd("~/Desktop/RFolder")  #  all files must be in the 
#  designated folder
#
# read in the myfunctions.R file
source("myfunctions.R")
#
# Load libraries.  These must be installed before you can 
# use them, but you only need to install ONCE
install.packages("plotluck", "lubridate", "tidyverse")
library(plotluck)
library(lubridate)
library(tidyverse)
#
# Read in data.  Most datafiles will have a header row,
# which are column headings
filename <-"brainIQMissing.csv"
data <- read.csv(filename, header=TRUE)
#
# Run some code that will help you understand the data
#  You don't have to do all of these, but pick the ones
# you like best
dim(data)  #  shows the number of rows and columns.  You can also
# see this in the environment window on the right.
head(data)  # shows the first six rows
str(data)
summary(data)  # this shows the descriptive statistics, very helpful
# 
sapply(data, class)
# This next line counts how many columns are numerical
paste("There are ", sum(ifelse(sapply(data, class) == "numeric", 1, 0)), " numerical columns.")
#
#  CHECK:  is my data clean?   This code chunk will tell you 
# how many rows have one or more missing items.  You need to clean
# missing data before you can do any meaningful analytics
{
  clean <- ifelse(complete.cases(data)==TRUE,1,0)
paste("There are ",dim(data)[1]-sum(clean), " rows with missing data.")
}
#
# clean up missing numeric data.  Looking at my data, I see that
# columns 3 through 1 are numerical values.  Column 12 is a date.
# This code chunk replaces all missing numerical values with the
# mean of the existng data for that column.
indx <- which(sapply(data, is.numeric))
indx
for(i in indx) {                                
  data[, i][is.na(data[ , i])] <- mean(data[ , i], na.rm = TRUE)
}
#
# This dataset has data values, in the format mo/day/year, or mdy.
#  lubridate will convert these to a more standard YYYY-MO-DA format.
data$Date <- mdy(data$Date)
#
# Now we need to clean up missing character data.  In this example, only
# Column 2 ("Patient Type") is character data.  To use the random.impute
# function, you must have myfunctions.R loaded.
#  We impute the data, which gives a new column.  We replace the 
# original column with missing data with the newly-created column.
#  Then, we delete the imputed column. 
#data <- sample(data, length(data), replace  = TRUE)
#  random.impute.data.frame(data, c(2))
 #data$Import <- data$Import.imputed
#data$Import.imputed <- NULL
library(dplyr)

data <- data %>% mutate(c2 = ifelse(is.na(c2)), sample(In, Out), sum(is.na(c2)), replace = TRUE,c2)
#
#  Now we need to look at the statistics of our numerical data
#  For normal statistics to work, numerical values need to be 
# normally distributed, a "bell-shaped curve".  You can do analyses
# without them being normally distributed, but it's better if they are.
# First, we'll run plotluck, to see how the numerical values are 
#  distributed. 
#  Then, using the "rz.transform" function in myfunctions.r, we'll
# do a "rank-Z" transformation.  You could also take the logarithm
# base 10 instead of doing a rank-Z transform.  Finally, we'll run
# plotluck again to see how it went.
## IF YOU DON'T WANT TO TRANSFORM YOUR DATA, you comment it OUT
plotluck(data, .~1)
data[indx] <- lapply(data[indx], rz.transform)
plotluck(data, .~1)
#
# Once you have clean data, it's not a bad idea to save your file as a 
# clean file
write.csv(data, "protein_clean.csv")
#  Now we can do some more data analytics.  By running a
# pairwise plot, we can look at scatter plots, histograms,
# and correlation coefficients for numerical data
# It's best to see all of them, but that can create a 
# graph that is too big, so you can split them up
pairs(data[,3:7],upper.panel=panel.cor,diag.panel=panel.hist)
pairs(data[,8:11],upper.panel=panel.cor,diag.panel=panel.hist)
#
#  Now we can use built-in R functions to filter data.  For example,
# which countries have a value for red meat at or above its mean?

rm.filter <- data$RedMeat >= mean(data$RedMeat) 
rm.filter
data$Country[rm.filter]
# Which countries have red meat at or below its mean AND
# white meat at or below its mean?
{
  meat.filter <- data$RedMeat <= mean(data$RedMeat) & data$WhiteMeat <= mean(data$WhiteMeat)
print("Countries at or below means for red and white meat: ", quote=FALSE)
data$Country[meat.filter]
}
#
#
#  Now we cqn use tidyverse, the "grammar of data science", to
# select, filter, arrange, and summarize the data
#
# Select
#  I'll create a new dataset that only contains the columns of Country
# and RedMeat
redmeat <- select(data, Country, RedMeat)
head(redmeat)
# a new dataset that has country, RedMeat and WhiteMeat
meats <- select(data, Country, RedMeat:WhiteMeat)
head(meats)
#
# Filter
#  same as before, only using tidyverse
bigredmeat <- filter(data, RedMeat >= mean(RedMeat))
head(bigredmeat)
redwhitemeat <- filter(
  data,
  RedMeat >= mean(RedMeat),
  WhiteMeat >= mean(WhiteMeat)
)
head(redwhitemeat)
#
# Mutate:  this ADDS COLUMNS TO THE DATA FRAME, typically based on some
# function.  For example, I want to add two new columns, the ratio of white
# meat to red meat and the ratio of eggs to milk
data <- mutate(
  data,
  MeatRatio = WhiteMeat/RedMeat,
  DairyRatio = Eggs/Milk
)
#
# Arrange:  this will sort your data. In this example, I sort
# RedMeat from high to low (using a minus sign)
data <- arrange(data, RedMeat)
# Summarize:  I can use summarize to create general data results.  For
# example, I want to find the means of RedMeat and WhiteMeat
mean_meat <- summarize(
  data,
  mean_RedMeat = mean(RedMeat),
  mean_WhiteMeat = mean(WhiteMeat)
)
mean_meat
#
# Using group_by and a pipe (>%>), I can create another table
meat_summary <- data %>%
  group_by(Country) %>%
  summarize(
    mean_RedMeat = mean(RedMeat),
    mean_WhiteMeat = mean(WhiteMeat)
  )
head(meat_summary)
#
## PLOTTING
#  You can use standard R functions for plotting
#  Notice there is also a linear regression and
# a printing of the regression line
plot(data$WhiteMeat~data$RedMeat, main="Plot of RedMeat vs. WhiteMeat")
meatFit <- lm(data$WhiteMeat~data$RedMeat)
abline(meatFit)
plot(data$Eggs~data$Milk, main="Plot of Milk vs. Eggs")
dairyFit <- lm(data$Eggs~data$Milk)
abline(dairyFit)
#
#  ggplot2: ggplot2 is a plotting package
# that is part of tidyverse, and you can do
# more informative plots with it.  Here
# we plot Milk vs. Eggs, create a regression
# line, color-code the points by types of Import,
# with the size of the dots relative 
# to the MeatRation
lm <- ggplot(data=data, aes(x=Milk, y=Eggs)) + geom_point(aes(col=Import, size=MeatRatio)) + geom_smooth(method="lm", se=F) + labs(title= "Linear Model",  caption = "Data = Protein")
lm
#  We can also use a different regression
# method called "LOESS".  It fits a polynomial
# to the predictors
loess <- ggplot(data=data, aes(x=Milk, y=Eggs))+ 
  geom_point(aes(col=Import, size=MeatRatio)) + geom_smooth(method="loess", se=F) +labs(subtitle="Milk vs. Eggs", title= "LOESS Model",  caption = "Data = Protein")
loess
#
# EOF End of File
