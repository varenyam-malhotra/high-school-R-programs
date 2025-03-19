#-------------------------------------------------------------------------------------------------
# R script for analyzing QTL data (Sugiyama2: https://phenome.jax.org/projects/Sugiyama2)
# using the template provided by Mr Gotwals. 
# Varenyam Malhotra
# October 17, 2022
#-------------------------------------------------------------------------------------------------

# clean things up and set up working directory
rm(list=ls())
setwd("/Users/vena/Desktop/RFolder/hyper")

# load the QTL library after installation using install.packages("qtl")
library(qtl)

# Load the data!
cross <- read.cross("csv", file="Sugiyama2001_B6xA_BC_B37_Data.csv", genotypes = c("B", "H"), na.strings = "-",
                    alleles = c("B", "C"))

# Sometimes the genetic markers are too close. Jittermap will move them apart slightly so my results are better.
jittermap(cross)

#  A summary of the cross gives me some basic data providing number of individuals, phenotypes, chromosomes, markers
summary(cross)
# I need to see what phenotypes are in the dataset
names(cross$pheno)

# take a look at my data, make sure it's pretty clean
#  I should NOT see any really big red spots, especially in the bottom right corner under the diagonal
cross <- est.rf(cross)

# ------------Histograms ------------------------------------------------
# 3 (a) Histogram of BP phenotype 
bp <- cross$pheno$BP_salt
hist(bp)
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
plot(cross.scanBP, main= "Blood Pressure Mainscan")
#  I'm putting threshold likes at 63% confidence, 90% confidence, and 95% confidence.
thresh <- summary(cross.scanBP.perm, alpha = c(0.37, 0.10, 0.05))
abline(h=thresh[1], col="blue")
abline(h=thresh[2], col="red")
abline(h=thresh[3], col="green")

#  mainscan summary text for Blood pressure
summary(cross.scanBP, perm=cross.scanBP.perm, alpha=0.05)

# do an effect plot
# once you see an effect plot, you'll understand what it does!
#BP Effect plots
#chromosome 1
firstbpchr <- find.marker(cross, chr=1, pos=41.4)
effectplot(cross, pheno.col=3, mname1=firstbpchr)
#chromosome 4
secondbpchr<- find.marker(cross, chr=4, pos=32.5)
effectplot(cross, pheno.col=3, mname1=secondbpchr)

#BP Confidence Intervals
#chromosome 1
CIchr1 <- bayesint(cross.scanBP, lodcolumn=1, chr=1, prob=0.95)
CIchr1
plot(cross.scanBP, chr=1, main="Confidence Interval for Chr1")
lines(x=CIchr1[c(1,3),2], y=c(0,0), type="l", col="green", lwd=4)
#chromosome 4
CIchr4 <- bayesint(cross.scanBP, lodcolumn=1, chr=4, prob=0.95)
CIchr4
plot(cross.scanBP, chr=4, main="Confidence Interval for Chr4")
lines(x=CIchr4[c(1,3),2], y=c(0,0), type="l", col="green", lwd=4)



# All done
#EOF