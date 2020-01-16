# NOAC AF Cost-effectiveness model
# Main script to run the model
# Howard Thom 14-03-2015
# Updated for Wei Fang and Zhenru Wang to use noac.net.benefit() 29-05-2017

# Optimized to calculate cohort vectors for treatments in parallel (increases speed by factor of 3.4)
# Updated to sample n.cycle transition matrices for each age (very slow and memory intensive!)

baseline.directory<-dirname(rstudioapi::getSourceEditorContext()$path)

library(expm)
library(nnet)
library(MASS) # mvrnorm
library(MVN)
library(mvtnorm)

n.samples<-1000
n.cycles<-120
initial.age<-70

# CE thresholds at which to generate the NB
lambdas<-c(1:50)*1000 
lambdas<-20000
# Run requisite scripts and load functions
source("utility.functions.1.R")
source("generate.transition.matrix.7.R")
source("NOAC.AF.net.benefit.1.R")
source("age.independent.generate.probabilities_2.R")
source("EVPI_std.R")
source("EVPI_std_p.R")

# Load outputs from WinBUGS for log hazard ratios, baseline log hazards, and no treamtent hazard ratios
bugs.loghr<-read.csv(file=paste(baseline.directory,"/data/bugs.loghr.csv",sep=""))
bugs.baseline<-read.csv(file=paste(baseline.directory,"/data/bugs.baseline.csv",sep=""))
hr.no.treatment<-read.csv(file=paste(baseline.directory,"/data/hr.no.treatment.csv",sep=""))
# Remove the empty first columns
bugs.loghr<-bugs.loghr[,-1]
bugs.baseline<-bugs.baseline[,-1]
hr.no.treatment<-hr.no.treatment[,-1]

# Multivariate normal distribution fit
test = bugs.loghr[,8:35]
sigma_loghr = cov(test)
mu_loghr = colMeans(test)

