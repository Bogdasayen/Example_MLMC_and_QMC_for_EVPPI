
mc_EVPPI_lor2 <- function(N,n){
    
    source("net_benefit.R")
    pkg <- c("Matrix","MASS","boot","tictoc","randtoolbox")
    lapply(pkg, library, character.only = TRUE)
    
    tic()
    set.seed(666)
    
    EVPPI_sum <- matrix(0,3,2)
    
    mean_rec <- c(0.99,1.33)
    var_rec <- matrix(c(0.22,0.15,0.15,0.2),2,2)
    mean_rel <- c(-1.48,-0.4)
    var_rel <- matrix(c(0.14,0.05,0.05,0.11),2,2)
    
    for(N1 in seq(1,N,by=1e4)){
        N2 <- min(1e4, N-N1+1)
    
        lor_rec2 <- rnorm(N2,mean_rec[1],sqrt(var_rec[1,1]))
        lor_rel2 <- rnorm(N2,mean_rel[1],sqrt(var_rel[1,1]))
        
        lor_rec2 <- c(t(replicate(n,lor_rec2)))
        lor_rel2 <- c(t(replicate(n,lor_rel2)))
        
        mean_rec3 <- mean_rec[2] + var_rec[1,2]/var_rec[1,1]*(lor_rec2 - mean_rec[1])
        var_rec3 <- var_rec[2,2] - var_rec[1,2]*var_rec[2,1]/var_rec[1,1]
        mean_rel3 <- mean_rel[2] + var_rel[1,2]/var_rel[1,1]*(lor_rel2 - mean_rel[1])
        var_rel3 <- var_rel[2,2] - var_rel[1,2]*var_rel[2,1]/var_rel[1,1]
        
        
        C_t <- c(0,300,30)
        lambda <- 2e4
            
        C_rec <- matrix(rnorm(n*N2,1000,50),n,N2)
        C_rel <- matrix(rnorm(n*N2,2000,100),n,N2)
        C_norec <- matrix(rnorm(n*N2,2500,125),n,N2)
        Q_rec <- matrix(rnorm(n*N2,26,2),n,N2)
        Q_rel  <- matrix(rnorm(n*N2,23,3),n,N2)
        Q_norec <- matrix(rnorm(n*N2,20,4),n,N2)
        
        lor_rec3 <-  rnorm(n*N2,0,sqrt(var_rec3))+mean_rec3
        lor_rel3 <-  rnorm(n*N2,0,sqrt(var_rel3))+mean_rel3
        
        P_rec <- matrix(0,N2*n,3)
        P_rel <- matrix(0,N2*n,3)
        P_rec[,1] <- rbeta(N2*n,6,200)
        P_rec[,2] <- inv.logit(logit(P_rec[,1])+lor_rec2)
        P_rec[,3] <- inv.logit(logit(P_rec[,1])+lor_rec3)
        
        P_rel[,1] <- rbeta(N2*n,2,100)
        P_rel[,2] <- inv.logit(logit(P_rel[,1])+lor_rel2)
        P_rel[,3] <- inv.logit(logit(P_rel[,1])+lor_rel3)
        
        for(i in seq(1,N2)){
            ind <- seq((i-1)*n+1,i*n)
            
            Result <- net_benefit(lambda,P_rec[ind,],P_rel[ind,],C_rec[,i],C_rel[,i],C_norec[,i],
                                  Q_rec[,i],Q_rel[,i],Q_norec[,i],C_t)
            NB <- colMeans(Result$NB)
            NB <- max(NB) - NB;
            EVPPI_sum[,1] <- EVPPI_sum[,1] + NB
            EVPPI_sum[,2] <- EVPPI_sum[,2] + NB^2
        }
    }
    
    EVPPI <- min(EVPPI_sum[,1]/N)
    ind <- which.min(EVPPI_sum[,1])
    var <- (EVPPI_sum[ind,2]/N-EVPPI^2)/N
    toc()
    cat(sprintf("EVPPI (lor2) = %.4f +/- %.4f, var = %.4f, N=%.1e. \n\n ", EVPPI,3*sqrt(var),var,N))
    return(list(EVPPI=EVPPI,var=var))
}