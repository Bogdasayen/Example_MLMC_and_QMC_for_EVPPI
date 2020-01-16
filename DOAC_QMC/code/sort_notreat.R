
# ind in seq(1,7), d=7
# MCMC.noTreatment.samples = hr.no.treatment[MCMC.selection,]

sort_notreat<-function(U)
{
    # source("recursive_sort_3.R")
    # library(randtoolbox)
    
    n <- 16
    
    # X1 <- hr.no.treatment[,ind[1]]
    # X2 <- hr.no.treatment[,ind[2]]
    # X3 <- seq(1,dim(hr.no.treatment)[1])
    # x <- cbind(X1,X2,X3)
    # 
    # MCMC.selection <- sample(1:59999, 5537, replace = FALSE)
    # x <- rbind(x,x[MCMC.selection,])  # Now we have 2^16 samples
    # x <- recursive_sort_3(0,x,1,2,n)
    
    N <- dim(U)[1]
    Y <- rep(0,N)
    # U <- sobol(N, dim = 2, init = TRUE, scrambling = 1, seed = sample(1:10^6, 1), normal = FALSE)
    # U <- sobol(N, dim = 2, init = FALSE, scrambling = 1, seed = sample(1:10^6, 1), normal = FALSE)
    
    for(k in seq(1,N)){
        i <- floor(2^(n/2)*U[k,1])
        j <- floor(2^(n/2)*U[k,2])
        ibits <- as.numeric(intToBits(i)[1:n])
        jbits <- as.numeric(intToBits(j)[1:n])
        Y[k] <- 4^seq(0,n-1)%*%(jbits+2*ibits)+1
        # Y[k] <- x[ind+1,3]
    }
    
    return(list(QMCsample_notreat=Y))
}
