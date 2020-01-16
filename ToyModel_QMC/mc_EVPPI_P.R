
mc_EVPPI_P <- function(N,n){

    pkg <- c("Matrix","MASS","boot","tictoc","randtoolbox")
    lapply(pkg, library, character.only = TRUE)
    
    tic()
    set.seed(666)
    
    EVPPI_sum <- matrix(0,3,2)
    
    mean_rec <- c(0.99,1.33)
    var_rec <- matrix(c(0.22,0.15,0.15,0.2),2,2)
    mean_rel <- c(-1.48,-0.4)
    var_rel <- matrix(c(0.14,0.05,0.05,0.11),2,2)
    
    C_t <- c(0,300,30)
    lambda <- 2e4
    
    for(N1 in seq(1,N,by=1e4)){
        N2 <- min(1e4, N-N1+1)
        
        lor_rec <- mvrnorm(N2,mean_rec,var_rec)
        lor_rel <- mvrnorm(N2,mean_rel,var_rel)
        P_rec <- matrix(0,N2,3)
        P_rel <- matrix(0,N2,3)
        P_rec[,1] <- rbeta(N2,6,200)
        P_rec[,2] <- inv.logit(logit(P_rec[,1])+lor_rec[,1])
        P_rec[,3] <- inv.logit(logit(P_rec[,1])+lor_rec[,2])
        
        P_rel[,1] <- rbeta(N2,2,100)
        P_rel[,2] <- inv.logit(logit(P_rel[,1])+lor_rel[,1])
        P_rel[,3] <- inv.logit(logit(P_rel[,1])+lor_rel[,2])
        
        C_rec <- matrix(rnorm(n*N2,1000,50),n,N2)
        C_rel <- matrix(rnorm(n*N2,2000,100),n,N2)
        C_norec <- matrix(rnorm(n*N2,2500,125),n,N2)
        
        Q_rec <- matrix(rnorm(n*N2,26,2),n,N2)
        Q_rel  <- matrix(rnorm(n*N2,23,3),n,N2)
        Q_norec <- matrix(rnorm(n*N2,20,4),n,N2)
        
        QALY = matrix(0,n,3)
        Cost = matrix(0,n,3)
        
        for(i in seq(1,N2)){
            QALY[,1] <- P_rec[i,1]*(1-P_rel[i,1])*Q_rec[,i] + P_rec[i,1]*P_rel[i,1]*Q_rel[,i] + (1-P_rec[i,1])*Q_norec[,i]
            QALY[,2] <- P_rec[i,2]*(1-P_rel[i,2])*Q_rec[,i] + P_rec[i,2]*P_rel[i,2]*Q_rel[,i] + (1-P_rec[i,2])*Q_norec[,i]
            QALY[,3] <- P_rec[i,3]*(1-P_rel[i,3])*Q_rec[,i] + P_rec[i,3]*P_rel[i,3]*Q_rel[,i] + (1-P_rec[i,3])*Q_norec[,i]
            
            Cost[,1] <- P_rec[i,1]*(1-P_rel[i,1])*C_rec[,i] + P_rec[i,1]*P_rel[i,1]*C_rel[,i] + (1-P_rec[i,1])*C_norec[,i] 
            Cost[,2] <- P_rec[i,2]*(1-P_rel[i,2])*C_rec[,i] + P_rec[i,2]*P_rel[i,2]*C_rel[,i] + (1-P_rec[i,2])*C_norec[,i] + C_t[2]
            Cost[,3] <- P_rec[i,3]*(1-P_rel[i,3])*C_rec[,i] + P_rec[i,3]*P_rel[i,3]*C_rel[,i] + (1-P_rec[i,3])*C_norec[,i] + C_t[3] 
            
            NB <- colMeans(lambda*QALY - Cost)
            NB <- max(NB) - NB;
            EVPPI_sum[,1] <- EVPPI_sum[,1] + NB
            EVPPI_sum[,2] <- EVPPI_sum[,2] + NB^2
        }
    }

    EVPPI <- min(EVPPI_sum[,1]/N)
    ind <- which.min(EVPPI_sum[,1])
    var <- (EVPPI_sum[ind,2]/N-EVPPI^2)/N
    toc()
    cat(sprintf("EVPPI (P) = %.4f +/- %.4f, var = %.4f, N=%.1e, n=%d. \n\n ", EVPPI,3*sqrt(var),var,N,n))
    return(list(EVPPI=EVPPI,var=var))
}