
mc_EVPI <- function(N){
    
    pkg <- c("Matrix","MASS","boot","tictoc")
    lapply(pkg, library, character.only = TRUE)
    
    tic()
    set.seed(666)
    
    EVPI_sum <- matrix(0,3,2)
    QALY_mean <- rep(0,3)
    COST_mean <- rep(0,3)
    
    mean_rec <- c(0.99,1.33)
    var_rec <- matrix(c(0.22,0.15,0.15,0.2),2,2)
    mean_rel <- c(-1.48,-0.4)
    var_rel <- matrix(c(0.14,0.05,0.05,0.11),2,2)
    
    for(N1 in seq(1,N,by=1e4)){
        N2 <- min(1e4, N-N1+1)
        
        P_rec <- matrix(0,N2,3)
        P_rel <- matrix(0,N2,3)
        lor_rec <- mvrnorm(N2,mean_rec,var_rec)
        lor_rel <- mvrnorm(N2,mean_rel,var_rel)
        
        P_rec[,1] <- rbeta(N2,6,200)
        P_rec[,2] <- inv.logit(logit(P_rec[,1])+lor_rec[,1])
        P_rec[,3] <- inv.logit(logit(P_rec[,1])+lor_rec[,2])
        
        P_rel[,1] <- rbeta(N2,2,100)
        P_rel[,2] <- inv.logit(logit(P_rel[,1])+lor_rel[,1])
        P_rel[,3] <- inv.logit(logit(P_rel[,1])+lor_rel[,2])
        
        C_rec <- rnorm(N2,1000,50)
        C_rel <- rnorm(N2,2000,100)
        C_norec <- rnorm(N2,2500,125)
        
        Q_rec <- rnorm(N2,26,2)
        Q_rel  <- rnorm(N2,23,3)
        Q_norec <- rnorm(N2,20,4)
        
        C_t <- c(0,300,30)
        lambda <- 2e4
        
        QALY <- P_rec*(1-P_rel)*Q_rec + P_rec*P_rel*Q_rel + (1-P_rec)*Q_norec
        Cost <- P_rec*(1-P_rel)*C_rec + P_rec*P_rel*C_rel + (1-P_rec)*C_norec 
        
        Cost[,2] <- Cost[,2] +C_t[2]
        Cost[,3] <- Cost[,3] +C_t[3]
        
        NB <- lambda*QALY - Cost
        QALY_mean <- QALY_mean + colSums(QALY)
        COST_mean <- COST_mean + colSums(Cost)
        
        NB <- replicate(3,apply(NB,1,max)) -NB
        EVPI_sum[,1] <- EVPI_sum[,1] + colSums(NB)
        EVPI_sum[,2] <- EVPI_sum[,2] + colSums(NB^2)
    }
    
    EVPI <- min(EVPI_sum[,1]/N)
    ind <- which.min(EVPI_sum[,1])
    EVPI_var <- EVPI_sum[ind,2]/N-EVPI^2    
    var <- EVPI_var/N
    QALY_mean <- QALY_mean/N
    COST_mean <- COST_mean/N
    toc()
    
    cat(sprintf("EVPI (P) = %.4f +/- %.4f, var = %.4f, N=%.1e. \n\n ", EVPI,3*sqrt(var),var,N))
    return(list(EVPI=EVPI,var=var,QALY_mean=QALY_mean,COST_mean=COST_mean))
}