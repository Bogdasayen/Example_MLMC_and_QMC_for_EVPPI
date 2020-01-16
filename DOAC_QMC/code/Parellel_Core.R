library(foreach)
library(doParallel)
library(parallel)
library(tictoc)
require(doRNG)
numCores <- detectCores()  
cl <- makeCluster(numCores)
registerDoParallel(cl)
# the above code prepare the parallel cores and clusters
stopCluster(cl)

## EVPI

K = 10
N = 2000

for(j in seq(1,K)){

# for each calculation of each core, we only calculate 1000 samples for efficiency
NN = 1000
inputs <- 1:(N/NN)

U <- sobol(N, dim = 57, init = TRUE, scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)
U <- sobol(N, dim = 57, init = FALSE, scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)

# tic()
# This loop does the parallel calculation
results <- foreach(i=inputs)  %dorng%  {
    index <- seq((i-1)*NN+1,i*NN)
    EVPI_qmc_p(NN,U,index)
}
# toc()
# Due to some function unavailable in the parallel loop, we do statitics outside 
sum1 = rep(0,2)
for(i in inputs){
  NetB = matrix(unlist(results[i]),NN,5)
  Ind = which.is.max(colMeans(NetB))
  EVPI_sample = apply(NetB,1,max)-NetB[,Ind]
  
  sum1[1] = sum1[1]+sum(EVPI_sample)
  sum1[2] = sum1[2]+sum(EVPI_sample^2)
}

Sum[j] = sum1[1]/N
Sum2[j] = Sum[j]^2
# EVPI_std = sqrt(sum1[2]/N-EVPI^2)/sqrt(N)
}

EVPI = sum(Sum)/K
EVPI_std = sqrt(sum(Sum2)/K - EVPI^2)/sqrt(K)

# 411.047   6.61

##############################################################################################
## EVPPI
N = 1024
M = 1024

NN = 1024
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
# This loop does the parallel calculation
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

# The following code can calculate the the EVPPI value for any parameters
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

