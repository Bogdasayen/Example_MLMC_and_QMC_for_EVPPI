qmc_EVPI <- function(N,M){
    
    source("net_benefit.R")
    source("sfunc.R")
    
    require(Matrix)
    require(MASS) # mvrnorm
    require(boot) # logit inv.logit
    require(randtoolbox)
    require(tictoc)
    
    tic()
    set.seed(666)
    
    EVPI <- matrix(0,M,2)
    QALY <- matrix(0,M,3)
    COST <- matrix(0,M,3)
    
    mean_rec <- c(0.99,1.33)
    var_rec <- matrix(c(0.22,0.15,0.15,0.2),2,2)
    L_rec = chol(var_rec)
    mean_rel <- c(-1.48,-0.4)
    var_rel <- matrix(c(0.14,0.05,0.05,0.11),2,2)
    L_rel = chol(var_rel)
    
    U <- sobol(N, dim = 12, init = TRUE, scrambling = 1, seed = 666, normal = FALSE)
    
    for(m in seq(1,M)){
        
        U <- sobol(N, dim = 12, init = FALSE, scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)
        
        P_rec <- matrix(0,N,3)
        P_rel <- matrix(0,N,3)
        
        X_rec = qnorm(U[,1:2])
        lor_rec <- t(replicate(N,mean_rec)) + X_rec%*%L_rec
        X_rel = qnorm(U[,3:4])
        lor_rel <- t(replicate(N,mean_rel)) + X_rel%*%L_rel
        
        P_rec[,1] <- Rbeta.inv(U[,5], 6, 200, log=FALSE)
        P_rec[,2] <- inv.logit(logit(P_rec[,1])+lor_rec[,1])
        P_rec[,3] <- inv.logit(logit(P_rec[,1])+lor_rec[,2])
        
        P_rel[,1] <- Rbeta.inv(U[,6], 2, 100, log=FALSE)
        P_rel[,2] <- inv.logit(logit(P_rel[,1])+lor_rel[,1])
        P_rel[,3] <- inv.logit(logit(P_rel[,1])+lor_rel[,2])
        
        C_rec <- qnorm(U[,7],1000,50)
        C_rel <- qnorm(U[,8],2000,100)
        C_norec <- qnorm(U[,9],2500,125)
        
        Q_rec <- qnorm(U[,10],26,2)
        Q_rel  <- qnorm(U[,11],23,3)
        Q_norec <- qnorm(U[,12],20,4)
        
        C_t <- c(0,300,30)
        lambda <- 2e4
        
        Result <- net_benefit(lambda,P_rec,P_rel,C_rec,C_rel,C_norec,Q_rec,Q_rel,Q_norec,C_t)
        
        NB <- Result$NB
        
        EVPI[m,1] <- mean(apply(NB,1,max)) - max(apply(NB,2,mean))
        EVPI[m,2] <- EVPI[m,1]^2
        QALY[m,] = colSums(Result$QALY)/N;
        COST[m,] = colSums(Result$COST)/N;
         
    }
    
    EVPI_QMC <- sum(EVPI[,1])/M
    EVPI_var_QMC <- sum(EVPI[,2])/M - EVPI_QMC^2;   
    var <- EVPI_var_QMC/M
    QALY_mean <- colSums(QALY)/M
    COST_mean <- colSums(COST)/M
    toc()
    
    cat(sprintf("QMC: EVPI = %.4f +/- %.4f, var = %.4f, N=%.1e. \n\n ", EVPI_QMC,3*sqrt(var),var,N))
    return(list(EVPI=EVPI_QMC,var=var,QALY_mean=QALY_mean,COST_mean=COST_mean))
}