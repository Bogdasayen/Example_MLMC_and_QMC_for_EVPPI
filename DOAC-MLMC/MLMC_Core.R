require(ggplot2)
require(grid)
require(Rcpp)
require(tictoc)
require(doRNG)

source("mlmc.R") #MLMC code to calculate the value to some accuracy
source("mlmc.test.R") #test function of the convergence rates and output the MLMC results
source("plot.mlmc.test.R") #function to plot the result
source("multiplot.R") #this is the function from ggplot

source("EVPI_l_p.R")
source("EVPI_std_p.R")
source("EVPPI_l_p.R")

source('EVPPI_Cost_EventCost_std_p.R')

set.seed(666) # Set random seed to ensure the same output 
# without this MLMC will give different output each time
tic()
tst <- mlmc.test(EVPPI_l_p, M=2, N=128,
                 L=5, N0=128,
                 eps.v=c(60),
                 Lmin=2, Lmax=10)
toc()

plot(tst,which=c("var", "mean", "Nl", "cost"),cols=2) 





