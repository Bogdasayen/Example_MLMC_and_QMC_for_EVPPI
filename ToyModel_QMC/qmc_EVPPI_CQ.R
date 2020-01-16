
qmc_EVPPI_CQ <- function(N,n,M){
    
    source("sfunc.R")
    pkg <- c("Matrix","MASS","boot","tictoc","randtoolbox")
    lapply(pkg, library, character.only = TRUE)
    
    tic()
    set.seed(666)
    
    mean_rec <- c(0.99,1.33)
    var_rec <- matrix(c(0.22,0.15,0.15,0.2),2,2)
    mean_rel <- c(-1.48,-0.4)
    var_rel <- matrix(c(0.14,0.05,0.05,0.11),2,2)
    
    C_t <- c(0,300,30)
    lambda <- 2e4
    
    EVPPI <- matrix(0,M,2)
    
    for(m in seq(1,M)){
        EVPPI_sum <- rep(0,3)
        for(N1 in seq(1,N,by=1e4)){
            N2 <- min(1e4, N-N1+1)
            
            U <- sobol(N, dim = 6, init = TRUE, scrambling = 1, seed = 666, normal = FALSE)
            C_rec <- qnorm(U[,1],1000,50)
            C_rel <- qnorm(U[,2],2000,100)
            C_norec <- qnorm(U[,3],2500,125)
            Q_rec <- qnorm(U[,4],26,2)
            Q_rel  <- qnorm(U[,5],23,3)
            Q_norec <- qnorm(U[,6],20,4)
            
            P_rec <- matrix(0,N2*n,3)
            P_rel <- matrix(0,N2*n,3)
            
            lor_rec <- mvrnorm(N2*n,mean_rec,var_rec)
            lor_rel <- mvrnorm(N2*n,mean_rel,var_rel)
            
            P_rec[,1] <- rbeta(N2*n,6,200)
            P_rec[,2] <- inv.logit(logit(P_rec[,1])+lor_rec[,1])
            P_rec[,3] <- inv.logit(logit(P_rec[,1])+lor_rec[,2])
            
            P_rel[,1] <- rbeta(N2*n,2,100)
            P_rel[,2] <- inv.logit(logit(P_rel[,1])+lor_rel[,1])
            P_rel[,3] <- inv.logit(logit(P_rel[,1])+lor_rel[,2])
            
            for(i in seq(1,N2)){
                ind <- seq((i-1)*n+1,i*n)
                
                QALY <- P_rec[ind,]*(1-P_rel[ind,])*Q_rec[i] + P_rec[ind,]*P_rel[ind,]*Q_rel[i] + (1-P_rec[ind,])*Q_norec[i]
                Cost <- P_rec[ind,]*(1-P_rel[ind,])*C_rec[i] + P_rec[ind,]*P_rel[ind,]*C_rel[i] + (1-P_rec[ind,])*C_norec[i] 
                Cost[,2] <- Cost[,2] + C_t[2]
                Cost[,3] <- Cost[,3] + C_t[3]
                
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
    cat(sprintf("QMC: EVPPI (CQ) = %.4f +/- %.4f, var = %.4f, N=%.1e, n=%d, M=%d. \n\n ", EVPPI_QMC,3*sqrt(var),var,N,n,M))
    return(list(EVPPI_QMC=EVPPI_QMC,var=var))
}