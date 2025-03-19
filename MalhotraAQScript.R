# ---------------------------Heading--------------------------------------------
# Name: Varenyam Malhotra
# Date: September 11, 2022
# Program Name: The name of this program is MalhotraAQScript.R
# Program Objective: This program analyzes basic characteristics of air quality 
# data by creating multiple different plots to analyze relationships between 
# different variables of the data"
# ------------------------Setting Working Directory-----------------------------

# Setting the working directory 
setwd("/Users/vena/Desktop/RFolder/airquality")

# ----------------------------Lab 1 Part 1--------------------------------------

# Reading airquality.csv and loading it into a dataframe called airqualitydf
airqualitydf <- read.csv("airquality.csv")

# 1. Output the headers of the airqualitydf data to the names.txt file
sink("names.txt")
names(airqualitydf)
sink()

# 2. Printing the basic descriptive statistics into the summary.txt file
sink("summary.txt")
summary(airqualitydf)
sink()

# 3. Performing a multiple linear regression with wind and temperature as the x values 
# and ozone as the y value. 
# 4. The basic descriptive linear regression statistics are summarized to the multireg.txt file
sink("multireg.txt")
owtmultireg.fit <- lm(airqualitydf$Ozone~airqualitydf$Wind+airqualitydf$Temp)
summary(owtmultireg.fit)
sink()

# 5. Performing a linear regression between temperature (x value) and ozone (y value). 
# The line of best fit and points are plotted, and this data is 
# saved as a jpeg regressionplot.jpg
jpeg("regressionplot.jpg")
plot(airqualitydf$Ozone~airqualitydf$Temp, xlab = "Temp", ylab = "Ozone") 
ot.fit <- lm(airqualitydf$Ozone~airqualitydf$Temp)
abline(ot.fit)
dev.off()

# 6. Creating four horizontal boxplots of the ozone, wind, solar, and temperature 
# air quality data in a 2x2 matrix and then saving all four boxplots in 
# one jpeg file boxplot.jpg
jpeg("boxplot.jpg")
par(mfrow=c(2,2))
boxplot(airqualitydf$Ozone, main = "Ozone Boxplot", horizontal = TRUE)
boxplot(airqualitydf$Wind, main = "Wind Boxplot", horizontal = TRUE)
boxplot(airqualitydf$Solar.R, main = "Solar Boxplot", horizontal = TRUE)
boxplot(airqualitydf$Temp, main = "Temp Boxplot", horizontal = TRUE)
dev.off()

# ----------------------------Lab 1 Part 2--------------------------------------

# Getting the baseline BIC number for Ozone
bic.ozoneozone <- BIC(lm(airqualitydf$Ozone~1)) # 1148.801
# 1. Temp effect on Ozone
bic.tempozone <- BIC(lm(airqualitydf$Ozone~airqualitydf$Temp)) # 1075.967
# 2. Solar.R effect on Ozone
bic.solarrozone <- BIC(lm(airqualitydf$Ozone~airqualitydf$Solar.R)) # 1091.843

# Getting the baseline BIC number for Solar.R
bic.solarrsolarr <- BIC(lm(airqualitydf$Solar.R~1)) # 1737.428
# 3. Ozone effect on Solar.R
bic.ozonesolarr <- BIC(lm(airqualitydf$Solar.R~airqualitydf$Ozone)) # 1315.552

# Getting the baseline BIC number for Temp
bic.temptemp <- BIC(lm(airqualitydf$Temp~1)) # 1131.027
# 5. Ozone effect on Temp
bic.ozonetemp <- BIC(lm(airqualitydf$Temp~airqualitydf$Ozone)) # 786.8075

#-------------------------------BIC Model Analysis------------------------------
# 1. Temp effect on Ozone
# As the BIC number for temperature's effect on ozone (1075.967) is over 70 less
# than the baseline BIC number for ozone (1148.801), I can conclude that temperature
# has a strong effect on ozone levels.
# 2. Solar.R effect on Ozone
# The BIC number for Solar.R's effect on ozone (1091.843) is over 45 less than the
# baseline BIC number for ozone, so I can conclude that Solar.R has a strong effect
# on ozone levels, but not as strong as the effect that temperature has. 
#
# 3. Ozone effect on Solar.R
# As the BIC number for ozone's effect on Solar.R (1315.552) is over 7420 less
# than the baseline BIC number for Solar.R (1737.428), I can conclude that ozone 
# levels have a monumental effect on Solar.R. 
# 
# 5. Ozone effect on Temp
# As the BIC number for ozone's effect on temperature (786.8075) is over 340 less
# than the baseline BIC number for temperature (1131.027), I can conclude that
# ozone levels have a very strong effect on temperature. 

#-----------------------------BIC Models Output File----------------------------

sink("bicscores.txt")
sprintf("Temp --> Ozone %f", bic.tempozone)
sprintf("Solar.R --> Ozone %f", bic.solarrozone)
sprintf("Ozone --> Solar.R %f", bic.ozonesolarr)
sprintf("Ozone --> Temp %f", bic.ozonetemp)
sink()




