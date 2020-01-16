library(foreach)
library(doParallel)
library(parallel)
library(tictoc)
require(doRNG)
numCores <- detectCores()  
cl <- makeCluster(numCores)
registerDoParallel(cl)
# the above code prepare the parallel cores and clusters
## need to run the above code to get the parallel computing work if you want to use MLMC with parallel
stopCluster(cl)
# clear the cluster after each use

## EVPI

N = 80

# for each calculation of each core, we only calculate 1000 samples for efficiency
NN = 10
inputs <- 1:(N/NN)
set.seed(666)
tic()
# This loop does the parallel calculation
results <- foreach(i=inputs)  %dorng%  {
  EVPI_std_p(NN)
}
toc()
# Due to some function unavailable in the parallel loop, we do statitics outside 
sum1 = rep(0,2)
for(i in inputs){
  NetB = matrix(unlist(results[i]),NN,5)
  Ind = which.is.max(colMeans(NetB))
  EVPI_sample = apply(NetB,1,max)-NetB[,Ind]
  
  sum1[1] = sum1[1]+sum(EVPI_sample)
  sum1[2] = sum1[2]+sum(EVPI_sample^2)
}

EVPI = sum1[1]/N
EVPI_std = sqrt(sum1[2]/N-EVPI^2)/sqrt(N)
# 416.7916   0.771

##############################################################################################
## EVPPI
N = 64
M = 8

NN = 64
inputs <- 1:(N*M/NN)

source("EVPPI_Cost_EventCost_std_p.R")
source("EVPPI_Cost_TreatmentCost_std_p.R")
source("EVPPI_Cost_StateCost_std_p.R")
source("EVPPI_Prob_Switch_std_p.R")
source("EVPPI_Prob_Effect_std_p.R")
source("EVPPI_HR_all_std_p.R")
source("EVPPI_HR_baseline_std_p.R")
source("EVPPI_HR_loghr_std_p.R")
source("EVPPI_HR_NT_std_p.R")
source("EVPPI_Utility_Age_std_p.R")
source("EVPPI_Utility_state_std_p.R")
source("EVPPI_Utility_event_std_p.R")
source("EVPPI_All_Cost_std_p.R")
source("EVPPI_All_Utility_std_p.R")
source("EVPPI_NOAC_std_p.R")
source("EVPPI_NOAC_Complex_std_p.R")
source("age.independent.generate.probabilities_2.R")
source("EVPPI_NOAC_std_p_t.R")


# This loop does the parallel calculation
# each loop corresponds to EVPPI of one set of parameters 
# just need to run one and then jump to line 140
results1 <- foreach(i=inputs)  %dorng%  {
  EVPPI_Cost_EventCost_std_p(M,NN/M)
}
results2 <- foreach(i=inputs)  %dorng%  {
  EVPPI_Cost_TreatmentCost_std_p(M,NN/M)
}
results3 <- foreach(i=inputs)  %dorng%  {
  EVPPI_Cost_StateCost_std_p(M,NN/M)
}
results4 <- foreach(i=inputs)  %dorng%  {
  EVPPI_Prob_Switch_std_p(M,NN/M)
}
results5 <- foreach(i=inputs)  %dorng%  {
  EVPPI_Prob_Effect_std_p(M,NN/M)
}
results6 <- foreach(i=inputs)  %dorng%  {
  EVPPI_HR_all_std_p(M,NN/M)
}
results7 <- foreach(i=inputs)  %dorng%  {
  EVPPI_HR_baseline_std_p(M,NN/M)
}
results8 <- foreach(i=inputs)  %dorng%  {
  EVPPI_HR_loghr_std_p(M,NN/M)
}
results9 <- foreach(i=inputs)  %dorng%  {
  EVPPI_HR_NT_std_p(M,NN/M)
}
results10 <- foreach(i=inputs)  %dorng%  {
  EVPPI_Utility_Age_std_p(M,NN/M)
}
results11 <- foreach(i=inputs)  %dorng%  {
  EVPPI_Utility_state_std_p(M,NN/M)
}
results12 <- foreach(i=inputs)  %dorng%  {
  EVPPI_Utility_event_std_p(M,NN/M)
}

results13 <- foreach(i=inputs)  %dorng%  {
  EVPPI_All_Cost_std_p(M,NN/M)
}
results14 <- foreach(i=inputs)  %dorng% {
  EVPPI_All_Utility_std_p(M,NN/M)
}

# Use multivariate normal distribution to approximate the covariance structure between MCMC samples
# firts calculate the mean and variance of the Multivariate normal 
test = bugs.loghr[,8:35]
sigma_loghr = cov(test)
mu_loghr = colMeans(test)

results15 <- foreach(i=inputs,.packages='MASS') %dopar%{
  EVPPI_NOAC_std_p(M,NN/M)
}
results16 <- foreach(i=inputs,.packages='MASS') %dopar%{
  EVPPI_NOAC_Complex_std_p(M,NN/M)
}
# eatimate the parameters of the Multivariate t distribution 
test_t = cov.trob(test)
sigma_loghr_t = test_t$cov
mu_loghr_t = test_t$center

results_t <- foreach(i=inputs,.packages='mvtnorm')  %dorng%  {
  EVPPI_NOAC_std_p_t(M,NN/M)
}

##############################################################
### change results1 to the one you just get (reusltN) to do the analysis below
results = results1 
##############################################################


# The following code can calculate the the EVPPI value for any parameters after the parallel computing
# First summerize the parallel results to a new matrix
NetB = matrix(NA,M*N,5)
for(i in inputs){
  NetB[((i-1)*NN+1):(i*NN),] = matrix(unlist(results[i]),NN,5)
}

# Calculate net benefit of the optimal decision for each sample 
NetB_max = apply(NetB,1,max)
NetB_max_sample = apply(matrix(NetB_max,N,M,byrow=TRUE),1,mean)

# find the optimal decision without any information and calculate the Net benefit
Ind = which.is.max(colMeans(NetB))
NetB_low_sample = apply(matrix(NetB[,Ind],N,M,byrow=TRUE),1,mean)

# find the optimal decision with partial information and calculate the net benefit
NetB_mid = matrix(NA,N,5)
for(i in 1:5){
  NetB_mid[,i] = apply(matrix(NetB[,i],N,M,byrow=TRUE),1,mean)
}
NetB_mid_sample = apply(NetB_mid,1,max)

# Construct the samples for both DIFF and EVPPI itself
DIFF_sample = NetB_max_sample- NetB_mid_sample
EVPPI_sample = NetB_mid_sample- NetB_low_sample

# Calculate the estimated value and Monte Carlo standard deviation
DIFF = mean(DIFF_sample)
DIFF_std = sqrt(var(DIFF_sample)/N)
EVPPI = mean(EVPPI_sample)
EVPPI_std = sqrt(var(EVPPI_sample)/N)

