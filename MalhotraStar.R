# ---------------------------Heading--------------------------------------------

# Name: Varenyam Malhotra
# Date: September 12, 2022
# Program Name: The name of this program is MalhotraStar.R
# Program Objective: This program loads and analyzes the data in the dataset HIP_star.csv

# ------------------------Setting Working Directory-----------------------------

# Setting the working directory 
setwd("/Users/vena/Desktop/RFolder/starplotting")

# ---------------------------Star Plotting Lab----------------------------------

# 1. Reading the HIP_star.csv file and assigning it to the stars dataframe.
stars <- read.csv("HIP_star.csv")

# 2. Printing the names of all items in the stars dataframe and outputting the 
#    names to the starnames.txt output file.
sink("starnames.txt")
names(stars)
sink()

# 3. Plotting the RA data (x axis) and DE data (y axis) and creating the starchart.jpg
#    file containing this plot. 
jpeg("starchart.jpg")
plot(stars$DE~stars$RA, xlab = "RA", ylab = "DE", main = "Plot of all stars 
     with right acension (RA) and declination (DE)")
dev.off()

# 4. Plotting box plots of the Vmag, Plx, pmRA, and pmDE stars data in a 2x2 matrix
#    and then saving that matrix to the jpeg boxplots.jpg
jpeg("boxplots.jpg")
par(mfrow=c(2,2))
boxplot(stars$Vmag, main = "Vmag")
boxplot(stars$Plx, main = "Plx")
boxplot(stars$pmRA, main = "pmRA")
boxplot(stars$pmDE, main = "pmDE")
dev.off()

# 5.Filtering out all stars with Vmag values between (and including) 8 and 9. 
#   Then plotting the filtered Vmag values (x axis) vs. the BV values (y axis).
#   The plot is saved to the jpeg VmagvsBV.jpg
jpeg("VmagvsBV.jpg")
vmagfiltered <- stars[which(stars$Vmag >= 8 & stars$Vmag <= 9), ] 
plot(vmagfiltered$B.V~vmagfiltered$Vmag, xlab = "Vmag[vmagfiltered]", ylab = "B.V[vmagfiltered]", 
     main = "Plot of Filtered Vmag values and B.V]")
dev.off()

# 6. Calculating the log (base 10) of the Plx values and saving those calculations
#    in the variable logPlx.
logPlx <- log10(vmagfiltered$Plx)

# 7. Creating a new variable logL based on a specific equation involving the 
#    filtered Vmag values and the newly created logPlx variable.
logL <- (15 - vmagfiltered$Vmag - 5*logPlx)/2.5

# 8. Plotting the BV values vs. the logL values with asteriks as points on the graph.
#    A line of best fit is also plotted on the graph. 
#    The graph is saved in the jpeg file BVvslogL.jpg  
jpeg("BVvslogL.jpg")
plot(logL~vmagfiltered$B.V, pch = 8, xlab = "B.V[vmagfiltered]", ylab= "logL[vmagfiltered]", 
     main = "Plot of B-V and logL (luminosity)")
abline(lm(logL~vmagfiltered$B.V))
dev.off()

