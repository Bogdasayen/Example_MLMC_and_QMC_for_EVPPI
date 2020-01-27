
baseline.directory<-"/Desktop/DOAC_QMC"
setwd(paste(baseline.directory,"/code",sep=""))

pkg <- c("nnet","expm","Matrix","MASS","boot","tictoc","randtoolbox","foreach","doParallel","parallel","doRNG","reshape2")
lapply(pkg, library, character.only = TRUE)

n.samples<-100
n.cycles<-120
initial.age<-70

# CE thresholds at which to generate the NB
lambdas<-c(1:50)*1000 
lambdas<-20000

# Run requisite scripts and load functions
source("utility.functions.1.R")
source("generate.transition.matrix.7.R")
source("NOAC.AF.net.benefit.1.R")
source("recursive_sort_3.R")
source("age.independent.generate.probabilities_2.R")
source("sfunc.R")
source("sort_loghr.R")
source("sort_baseline.R")
source("sort_notreat.R")

# Load outputs from WinBUGS for log hazard ratios, baseline log hazards, and no treamtent hazard ratios
bugs.loghr<-read.csv(file=paste(baseline.directory,"/data/bugs.loghr.csv",sep=""))
bugs.baseline<-read.csv(file=paste(baseline.directory,"/data/bugs.baseline.csv",sep=""))
hr.no.treatment<-read.csv(file=paste(baseline.directory,"/data/hr.no.treatment.csv",sep=""))
# Remove the empty first columns
bugs.loghr<-bugs.loghr[,-1]
bugs.baseline<-bugs.baseline[,-1]
hr.no.treatment<-hr.no.treatment[,-1]

########################################

data_loghr <- bugs.loghr[,seq(8,35)]
data_baseline <- bugs.baseline
data_notreat <- hr.no.treatment

mean_loghr <- apply(data_loghr,2,mean)
mean_baseline <- apply(data_baseline,2,mean)
mean_notreat <- apply(data_notreat,2,mean)

cov_loghr <- cov(data_loghr)
cov_baseline <- cov(data_baseline)
cov_notreat <- cov(data_notreat)

L_loghr <- chol(cov_loghr)
L_baseline <- chol(cov_baseline)

cov_loghr.m <- melt(cov_loghr)
cov_baseline.m <- melt(cov_baseline)
cov_notreat.m <- melt(cov_notreat)

cor_loghr <- cor(data_loghr)
cor_baseline <- cor(data_baseline)
cor_notreat <- cor(data_notreat)
cor_loghr.m <- melt(cor_loghr)
cor_baseline.m <- melt(cor_baseline)
cor_notreat.m <- melt(cor_notreat)

########################################

bugs.loghr <- bugs.loghr[-seq(1,13615),]
bugs.baseline <- bugs.baseline[-seq(1,13615),]
hr.no.treatment <- rbind(hr.no.treatment,hr.no.treatment[-seq(1,54462),])

# Sorting MCMC samples
# sort_ind <- c(8,9) # sort index of loghr
sort_ind <- c(13,18) # sort index of loghr
X1 <- bugs.loghr[,sort_ind[1]]
X2 <- bugs.loghr[,sort_ind[2]]
X3 <- seq(1,dim(bugs.loghr)[1])
x <- cbind(X1,X2,X3)
x <- recursive_sort_3(0,x,1,2,14)
bugs.loghr <- bugs.loghr[x[,3],]

sort_ind <- c(1,2) # sort index of baseline
X1 <- bugs.baseline[,sort_ind[1]]
X2 <- bugs.baseline[,sort_ind[2]]
X3 <- seq(1,dim(bugs.baseline)[1])
x <- cbind(X1,X2,X3)
x <- recursive_sort_3(0,x,1,2,14)
bugs.baseline <- bugs.baseline[x[,3],]

# sort_ind <- c(1,2) # sort index of notreat
sort_ind <- c(2,3) # sort index of notreat
X1 <- hr.no.treatment[,sort_ind[1]]
X2 <- hr.no.treatment[,sort_ind[2]]
X3 <- seq(1,dim(hr.no.treatment)[1])
x <- cbind(X1,X2,X3)
x <- recursive_sort_3(0,x,1,2,16)
hr.no.treatment <- hr.no.treatment[x[,3],]

rm(X1,X2,X3,x) 


# if parallel
numCores <- detectCores()  
cl <- makeCluster(numCores)
registerDoParallel(cl)
# stopCluster(cl)

# EVPI & EVPPI source functions
source("EVPI_QMC_par.R")
source("EVPPI_All_Cost_QMC_p2.R")
source("EVPPI_All_Utility_QMC_p.R")
source("EVPPI_Complex_QMC_MVN.R")
source("EVPPI_HR_all_QMC_MVN.R")
source("EVPPI_HR_all_QMC_p2.R")
source("EVPPI_HR_baseline_QMC_MVN.R")
source("EVPPI_HR_loghr_QMC_p.R")
source("EVPPI_HR_loghr_QMC_MVN.R")
source("EVPPI_HR_notreat_QMC_p.R")
source("EVPPI_Switch_Prob_QMC_p.R")

# Examples to run EVPPI functions above
EVPPI_HR_loghr_QMC_MVN(512,256,8)

EVPPI_HR_baseline_QMC_MVN(512,128,8)

EVPPI_Switch_Prob_QMC_p(256,128,8)

EVPPI_All_Utility_QMC_p(256,256,8)






