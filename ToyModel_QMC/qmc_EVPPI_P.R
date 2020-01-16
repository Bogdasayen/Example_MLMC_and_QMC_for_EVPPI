
qmc_EVPPI_P <- function(N,n,M){
    
    source("sfunc.R")
    pkg <- c("Matrix","MASS","boot","tictoc","randtoolbox")
    lapply(pkg, library, character.only = TRUE)
    
    tic()
    set.seed(sample(seq(1,10^6),1))
    
    mean_rec <- c(0.99,1.33)
    var_rec <- matrix(c(0.22,0.15,0.15,0.2),2,2)
    L_rec = chol(var_rec)
    mean_rel <- c(-1.48,-0.4)
    var_rel <- matrix(c(0.14,0.05,0.05,0.11),2,2)
    L_rel = chol(var_rel)
    
    C_t <- c(0,300,30)
    lambda <- 2e4
    
    EVPPI <- matrix(0,M,2)
    
    for(m in seq(1,M)){
        EVPPI_sum <- rep(0,3)
        
        for(N1 in seq(1,N,by=1e4)){
            N2 <- min(1e4, N-N1+1)
            
            U <- sobol(N2, dim = 6, init = TRUE, scrambling = 1, seed = 666, normal = FALSE)
            
            P_rec <- matrix(0,N2,3)
            P_rel <- matrix(0,N2,3)
            
            X_rec = qnorm(U[,1:2])
            lor_rec <- t(replicate(N2,mean_rec)) + X_rec%*%L_rec
            X_rel = qnorm(U[,3:4])
            lor_rel <- t(replicate(N2,mean_rel)) + X_rel%*%L_rel
            
            P_rec[,1] <- Rbeta.inv(U[,5], 6, 200, log=FALSE)
            P_rec[,2] <- inv.logit(logit(P_rec[,1])+lor_rec[,1])
            P_rec[,3] <- inv.logit(logit(P_rec[,1])+lor_rec[,2])
            P_rel[,1] <- Rbeta.inv(U[,6], 2, 100, log=FALSE)
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
                EVPPI_sum <- EVPPI_sum + NB
            }
        }
        EVPPI[m,1] <- min(EVPPI_sum/N)
        EVPPI[m,2] <- EVPPI[m,1]^2
    }
    
    
    EVPPI_QMC <- sum(EVPPI[,1])/M
    EVPPI_var_QMC <- sum(EVPPI[,2])/M - EVPPI_QMC^2
    var <- EVPPI_var_QMC/M
    toc()
    cat(sprintf("QMC: EVPPI (P) = %.4f +/- %.4f, var = %.4f, N=%.1e, n=%d, M=%d. \n\n ", EVPPI_QMC,3*sqrt(var),var,N,n,M))
    return(list(EVPPI_QMC=EVPPI_QMC,var=var))
}