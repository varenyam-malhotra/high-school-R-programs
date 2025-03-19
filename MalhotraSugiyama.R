#-------------------------------------------------------------------------------------------------
# R script for analyzing QTL data using the template provided by Mr Gotwals. 
# Varenyam Malhotra
# October 13, 2022
#-------------------------------------------------------------------------------------------------

# clean things up and set up working directory
rm(list=ls())
setwd("/Users/vena/Desktop/RFolder/hyper")

# load the QTL library after installation using install.packages("qtl")
library(qtl)

# Load the data!
cross <- read.cross("csv", file="sugiyamashort.csv", genotypes = c("B", "C", "H"), 
                    na.strings = "-", alleles = c("B", "C"))

# Sometimes the genetic markers are too close. Jittermap will move them apart slightly so my results are better.
jittermap(cross)

#  A summary of the cross gives me some basic data providing number of individuals, phenotypes, chromosomes, markers
summary(cross)
# I need to see what phenotypes are in the dataset
names(cross$pheno)
# you can see the genotypes in the data as well!
names(cross$geno)
# take a look at my data, make sure it's pretty clean
#  I should NOT see any really big red spots, especially in the bottom right corner under the diagonal
cross <- est.rf(cross)
plotRF(cross)
plot.map(cross)
plotMissing(cross)

# ------------Histograms ------------------------------------------------
# 3 (a) Histogram of BP phenotype 
bp <- cross$pheno$BP_final
hist(bp)
# 3 (b) Histogram of Resting Heart Rate phenotype
restingheartrate <- cross$pheno$HR_final
hist(restingheartrate)
# 3 (c) Histogram of Heart Weight
heartweight <- cross$pheno$heart_wt
hist(heartweight)
# -----------------------------------------------------------------------

# Pairwise recombination fractions and LOD scores.
plotRF(cross)
# -----------------------------------------------------------------------
# 4 (a) Genetic Map 
plot.map(cross)
# 4 (b)  Map of Missing Data 
plotMissing(cross)
# ------------------------------------------------------------------------

# ------------Blood Pressure ----------------------------------------------
#  Diagnostics for Blood Pressure. If data is relatively clean, we get a 45 degree diagonal line
qqnorm(bp)
qqline(bp)

#calculating the mainscan
# first you need to generate what the scan SHOULD look like
# Calculate a genetic probability map.
cross <- calc.genoprob(cross, step=2.0, off.end = 0.0, error.prob = 1.0e-4, 
                       map.function = "haldane", stepwidth = "fixed")

# Run a simulated geno probability calculations
cross <- sim.geno(cross, step=2.0, off.end = 0.0, error.prob = 1.0e-4, 
                  map.function = "haldane", stepwidth = "fixed")

# Perform the mainscan for the QTL
#  I'm only going to run this for 100 "permulations" -- typically you do 500-1000, but that takes a LONG time.
cross.scanBP <- scanone(cross, pheno.col = 3, model="normal", method = "em")
cross.scanBP.perm <- scanone(cross, pheno.col = 3, model="normal", method = "em", n.perm = 100)
#
# plot the mainscan for Blood Pressure
plot(cross.scanBP, main= "Blood Pressure Mainscan Plot")
#  I'm putting threshold likes at 63% confidence, 90% confidence, and 95% confidence.
threshBP <- summary(cross.scanBP.perm, alpha = c(0.37, 0.10, 0.05))
abline(h=threshBP[1], col="blue")
abline(h=threshBP[2], col="red")
abline(h=threshBP[3], col="green")

#  mainscan summary text-based output for Blood pressure
summary(cross.scanBP, perm=cross.scanBP.perm, alpha=0.05)

# ---------------------Resting Heart Rate ------------------------------
qqnorm(restingheartrate)
qqline(restingheartrate)

cross.scanHtR <- scanone(cross, pheno.col = 4, model="normal", method = "em")
cross.scanHtR.perm <- scanone(cross, pheno.col = 4, model="normal", method = "em", n.perm = 100)

# plot the mainscan for Heart Rate
plot(cross.scanHtR, main= "Heart Rate Mainscan Plot")
#  I'm putting threshold likes at 63% confidence, 90% confidence, and 95% confidence.
threshHtR <- summary(cross.scanHtR.perm, alpha = c(0.37, 0.10, 0.05))
abline(h=threshHtR[1], col="blue")
abline(h=threshHtR[2], col="red")
abline(h=threshHtR[3], col="green")

#  mainscan summary text for Heart Rate
summary(cross.scanHtR, perm=cross.scanHtR.perm, alpha=0.05)

# ---------------------- Heart Weight  ------------------------------
qqnorm(heartweight)
qqline(heartweight)

cross.scanHtWeight <- scanone(cross, pheno.col = 6, model="normal", method = "em")
cross.scanHtWeight.perm <- scanone(cross, pheno.col = 6, model="normal", method = "em", n.perm = 100)

# plot the mainscan for Heart Weight
plot(cross.scanHtWeight, main= "Heart Weight Mainscan Plot")
#  I'm putting threshold likes at 63% confidence, 90% confidence, and 95% confidence.
threshHtW <- summary(cross.scanHtWeight.perm, alpha = c(0.37, 0.10, 0.05))
abline(h=threshHtW[1], col="blue")
abline(h=threshHtW[2], col="red")
abline(h=threshHtW[3], col="green")

#  mainscan summary text for Heart Weight
summary(cross.scanHtWeight, perm=cross.scanHtWeight.perm, alpha=0.05)
#-------------------------------------------------------------------
# do an effect plot
# once you see an effect plot, you'll understand what it does!
# second effect plot
#  I'm going to do a confidence interval plot, try to zoom in on the scan
# All done

#BP Effect plot
#chromosome 7
firstbpchr <- find.marker(cross, chr=7, pos=44.71)
effectplot(cross, pheno.col=3, mname1=firstbpchr)
#chromosome 15
secondbpchr<- find.marker(cross, chr=15, pos=3.96)
effectplot(cross, pheno.col=3, mname1=secondbpchr)

#BP Confidence Intervals
#chromosome 7
CIchr7 <- bayesint(cross.scanBP, lodcolumn=1, chr=7, prob=0.95)
CIchr7
plot(cross.scanBP, chr=7, main="Confidence Interval for Chr7")
lines(x=CIchr7[c(1,3),2], y=c(0,0), type="l", col="green", lwd=4)
#chromosome 15
CIchr15 <- bayesint(cross.scanBP, lodcolumn=1, chr=15, prob=0.95)
CIchr15
plot(cross.scanBP, chr=15, main="Confidence Interval for Chr15")
lines(x=CIchr15[c(1,3),2], y=c(0,0), type="l", col="green", lwd=4)


#Resting Heart Rate Effect Plot
#chromosome 2
onlyheartratechr <- find.marker(cross, chr=2, pos=57.7)
effectplot(cross, pheno.col=4, mname1=onlyheartratechr)

#Resting Heart Rate Confidence Interval
#chromosome 2
CIchr2 <- bayesint(cross.scanHtR, lodcolumn=1, chr=2, prob=0.95)
CIchr2
plot(cross.scanBP, chr=2, main="Confidence Interval for Chr2")
lines(x=CIchr2[c(1,3),2], y=c(0,0), type="l", col="green", lwd=4)


#Heart Weight Effect Plot
#chromosome 12
onlyheartweightchr <- find.marker(cross, chr=12, pos=52.2)
effectplot(cross, pheno.col=6, mname1=onlyheartweightchr)

#Heart Weight Confidence Interval
#chromosome 12
CIchr12 <- bayesint(cross.scanHtWeight, lodcolumn=1, chr=12, prob=0.95)
CIchr12
plot(cross.scanBP, chr=12, main="Confidence Interval for Chr12")
lines(x=CIchr12[c(1,3),2], y=c(0,0), type="l", col="green", lwd=4)


#EOF