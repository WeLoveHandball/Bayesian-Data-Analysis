## Code Term Paper Chris Ferrat
# We use package BMS with Bodyfat-dataset

rm(list = ls())

WorkDir <- ("D:/Dropbox/Bearbeitete Dokumente Uni Basel/20Seminararbeit Bayesian Data Analysis - Topics in Statistics, Econometrics and Data Science/Code und Daten") ## change directory!
setwd(WorkDir)

#All cores of PC work now together
options(mc.cores = parallel::detectCores())

#	Chapter 4 - Loading and preparing dataset
#########################################

#install.packages("mfp")
library(mfp)
data(bodyfat)

colnames(bodyfat)[colnames(bodyfat)=="siri"] <- "bodyfat"

# Corrections due to package "mfp" #https://www.rdocumentation.org/packages/mfp/versions/1.5.2/topics/bodyfat acces on 30.11.2019
bodyfat$height[bodyfat$case==42] <- 69.5
bodyfat_cleaned <- (bodyfat[c(1:47,49:75,77:95,97:181,183:252),])

# Ordering columns that bodyfat in first column to be the dependent variable
bodyfat <- bodyfat_cleaned[,c(3,5:17)]

#	Descriptive Statistics
#########################################
str(bodyfat)
summary(bodyfat)

plot(bodyfat)
hist(bodyfat$bodyfat)
hist(bodyfat$abdomen)

# Chapter 4.1 - Creating a BMA object with dataset Bodyfat
############################################
set.seed(2019)
library(BMS)
bms_uniform_UIP <- bms(bodyfat, g="UIP", mprior="uniform", nmodel = 2000, mcmc="bd", burn = 2000, iter = 10000, user.int=FALSE)

# Table 1: See PIP and coefficients of the averaged model
estimates.bma(bms_uniform_UIP, exact = TRUE, order.by.pip = TRUE, include.constant = FALSE,incl.possign = TRUE, std.coefs = FALSE, condi.coef = FALSE)
#equivalent to coef(bms_uniform_UIP,...)
#exact = FALSE means coefficients PIP, coefficients is calculated MCMC sampled models

# Table 2: Additional Information of the output
info.bma(bms_uniform_UIP) #equivalent to summary(bms_uniform_UIP)

# Figure 1: PMP Correlation
plotConv(bms_uniform_UIP[1:100]) #convergence of the mcmc samples with exact weights

# All in 1 #(estimates.bma plus info.bma)
print(bms_uniform_UIP) #PIP sampled from MCMC

# Figure 2: Image of the best 400 models
image(bms_uniform_UIP[1:400])

# Table 3: Topmodels 1 and 2 - Variable inclusion yes=1/no=0
topmodels.bma(bms_uniform_UIP[1:2])

# Table 4: Betas of the 2 best models
beta.draws.bma(bms_uniform_UIP[1:2])#beta.draws.bma(bms_uniform_UIP[1:2])# or[, 1:2] at the end
beta.draws.bma(bms_uniform_UIP[1:2], stdev=TRUE) #If stdev=TRUE, it returns their posterior standard deviations.

# Table 5: Posterior model probability of the 2 best models
pmp.bma(bms_uniform_UIP[1:2])

# Table 6: Closer look at the first Topmodel
bms_uniform_UIP[1]$topmod

par(mfrow = c(1,1))

# Chapter 4.2 - Uniform, fixed and flexible priors
# Figure 3: Uniform case - UIP=N=248
bms_uniform_UIP_N_248 <- bms(bodyfat, mprior="uniform", g="UIP", nmodel = 2000, mcmc="bd", burn = 2000, iter = 10000, user.int=FALSE)
plotModelsize(bms_uniform_UIP_N_248, main="Posterior Model Size Distribution (Uniform_g_UIP=N=248) - Non-Informative", sub = "Mean: 4.4619", ylab = "Probability Mass", ylim=c(0,0.4))
#plot(bms_uniform_UIP) shows plot and mcmc convergence simultaneously

# Figure 4: Plot binominal fixed prior for Model Size 3
bms_fixed3_UIP <- bms(bodyfat, mcmc="bd", mprior="fixed", mprior.size = 3, g="UIP", nmodel = 2000, user.int=FALSE)
plotModelsize(bms_fixed3_UIP, main="Posterior Model Size Distribution (FixedModelSize3_UIP) - Informative", sub = "Mean: 3.2223", ylab = "Probability Mass")

# Figure 5: Plot binominal fixed prior for Model Size 6
bms_fixed6_UIP <- bms(bodyfat, mcmc="bd", mprior="fixed", mprior.size = 6, g="UIP", nmodel = 2000, user.int=FALSE)
plotModelsize(bms_fixed6_UIP, main="Posterior Model Size Distribution (FixedModelSize6_UIP) - Informative", sub = "Mean: 4.356", ylab = "Probability Mass")

# Figure 6: Plot binominal fixed prior for Model Size 10
bms_fixed10_UIP <- bms(bodyfat, mcmc="bd", mprior="fixed", mprior.size = 10, g="UIP", nmodel = 2000, user.int=FALSE)
plotModelsize(bms_fixed10_UIP, main="Posterior Model Size Distribution (FixedModelSize10_UIP) - Informative", sub = "Mean: 6.5063", ylim = c(0,0.4), ylab = "Probability Mass")

# Figure 14: All 3 fixed plots together
par(mfrow = c(3,1))
plotModelsize(bms_fixed3_UIP, main="Posterior Model Size Distribution (FixedModelSize3_UIP) - Informative", sub = "Mean: 3.2223", ylab = "Probability Mass")
plotModelsize(bms_fixed6_UIP, main="Posterior Model Size Distribution (FixedModelSize6_UIP) - Informative", sub = "Mean: 4.356", ylab = "Probability Mass")
plotModelsize(bms_fixed10_UIP, main="Posterior Model Size Distribution (FixedModelSize10_UIP) - Informative", sub = "Mean: 6.5063", ylim = c(0,0.4), ylab = "Probability Mass")

par(mfrow = c(1,1))

# Beta-Binominal for g equal to 0.2, 10, 50 and 248(UIP)
# Figure 7: Beta-Binominal for g equal to 0.2
bms_BB_0.2 <- bms(bodyfat, mprior="random",  g=0.2, nmodel = 2000, mcmc="bd", burn = 2000, iter = 10000, user.int=FALSE)
plotModelsize(bms_BB_0.2, main="Posterior Model Size Distribution (BB_g_0.2)", sub = "Mean: 9.2551", ylim = c(0,0.4), ylab = "Probability Mass")

# Figure 8: Beta-Binominal for g equal to 10
bms_BB_10 <- bms(bodyfat, mprior="random",  g=10, nmodel = 2000, mcmc="bd", burn = 2000, iter = 10000, user.int=FALSE)
plotModelsize(bms_BB_10, main="Posterior Model Size Distribution (BB_g_10)", sub = "Mean: 7.6234", ylim = c(0,0.4), ylab = "Probability Mass")

# Figure 9: Beta-Binominal for g equal to 50
bms_BB_50 <- bms(bodyfat, mprior="random",  g=50, nmodel = 2000, mcmc="bd", burn = 2000, iter = 10000, user.int=FALSE)
plotModelsize(bms_BB_50, main="Posterior Model Size Distribution (BB_g_50)", sub = "Mean: 5.1489", ylim = c(0,0.4), ylab = "Probability Mass")

# Figure 10: Beta-Binominal for g equal to 248(UIP)
bms_BB_UIP <- bms(bodyfat, mprior="random",  g="UIP", nmodel = 2000, mcmc="bd", burn = 2000, iter = 10000, user.int=FALSE)
plotModelsize(bms_BB_UIP, main="Posterior Model Size Distribution (BB_UIP_g_248)", sub = "Mean: 3.6417", ylim = c(0,0.4), ylab = "Probability Mass")

# Figure 15: All 4 Beta-Binomial plots together
par(mfrow = c(4,1))
plotModelsize(bms_BB_0.2, main="Posterior Model Size Distribution (BB_g_0.2)", sub = "Mean: 9.2551", ylim = c(0,0.4), ylab = "Probability Mass")
plotModelsize(bms_BB_10, main="Posterior Model Size Distribution (BB_g_10)", sub = "Mean: 7.6234", ylim = c(0,0.4), ylab = "Probability Mass")
plotModelsize(bms_BB_50, main="Posterior Model Size Distribution (BB_g_50)", sub = "Mean: 5.1489", ylim = c(0,0.4), ylab = "Probability Mass")
plotModelsize(bms_BB_UIP, main="Posterior Model Size Distribution (BB_UIP_g_248)", sub = "Mean: 3.6417", ylim = c(0,0.4), ylab = "Probability Mass")

# Chapter 4.3 - Posterior Density for different coefficients with large and small g
# uniform g=UIP=N=248 is in the code above, case g=0.2 is below
bms_uniform_g0.2 <- bms(bodyfat, mprior="uniform", g=0.2, nmodel = 2000, mcmc="bd", burn = 2000, iter = 10000, user.int=FALSE)
#plotModelsize(bms_uniform_g0.2, ylim = c(0,0.4), ylab = "Probability Mass")

par(mfrow = c(1,2))

# Figure 11: Final comparison Model space of uniform with large and small g
plotModelsize(bms_uniform_UIP, ylim=c(0,0.4), main="Posterior Model Size Distribution (Uniform_g=N=248)", ylab = "Probability Mass", sub = "g-Prior = 248, Mean: 4.4619")
plotModelsize(bms_uniform_g0.2, ylim=c(0,0.4), main="Posterior Model Size Distribution (Uniform_g=0.2)", ylab = "Probability Mass", sub = "g-Prior = 0.2, Mean: 7.0393")

# Figure 12: Final comparison Coefficients - Biceps (large VS small g)
density(bms_uniform_UIP, reg = "biceps", addons="esplzb", ylim=c(0,1.5), sub = "g-Prior = 248")
density(bms_uniform_g0.2, reg = "biceps", addons ="esplzb", ylim=c(0,1.5), sub = "g-Prior = 0.2", xlim=c(-0.4,0.9))

# Figure 13: Final comparison Coefficients - Abdomen (large VS small g)
density(bms_uniform_UIP, reg = "abdomen", addons="esplzb", ylim=c(0,9), sub = "g-Prior = 248")
density(bms_uniform_g0.2, reg = "abdomen", addons ="esplzb", ylim=c(0,9), sub = "g-Prior = 0.2")

######################### END #########################
