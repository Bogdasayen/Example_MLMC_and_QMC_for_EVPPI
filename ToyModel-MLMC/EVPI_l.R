EVPI_l <- function(l, N) {
    require(Matrix)
    require(MASS) # mvrnorm
    require(boot) # logit
    
    lamda <- 20000
    C_t1 <- 300
    C_t2 <- 30
    M <- 2^(l+1)
  
    mu_rec <- c(0.99,1.33)
    sigma_rec <- matrix(c(0.22,0.15,0.15,0.20),2,2,byrow = TRUE)
    mu_rel <- c(-1.48,-0.4)
    sigma_rel <- matrix(c(0.14,0.05,0.05,0.11),2,2,byrow = TRUE)
    
    sum1 <- rep(0, 7)
    
    for(N1 in seq(1, N, by=10000)) {
      N2 <- min(10000, N-N1+1)
      
      P_nt_rec <- matrix(rbeta(M*N2,6,200),M,N2,byrow = TRUE)
      P_nt_rel <- matrix(rbeta(M*N2,2,100),M,N2,byrow = TRUE)
      
      lor_rec <- mvrnorm(M*N2, mu_rec, sigma_rec)
      lor_rel <- mvrnorm(M*N2, mu_rel, sigma_rel)
      
      lor_t1_rec <- matrix(lor_rec[,1],M,N2,byrow = TRUE)
      lor_t1_rel <- matrix(lor_rel[,1],M,N2,byrow = TRUE)
      P_t1_rec <- inv.logit(logit(P_nt_rec)+lor_t1_rec)
      P_t1_rel <- inv.logit(logit(P_nt_rel)+lor_t1_rel)
      
      lor_t2_rec <- matrix(lor_rec[,2],M,N2,byrow = TRUE)
      lor_t2_rel <- matrix(lor_rel[,2],M,N2,byrow = TRUE)
      P_t2_rec <- inv.logit(logit(P_nt_rec)+lor_t2_rec)
      P_t2_rel <- inv.logit(logit(P_nt_rel)+lor_t2_rel)
      
      C_rec = matrix(rnorm(M*N2,1000,50),M,N2,byrow = TRUE)
      C_rel = matrix(rnorm(M*N2,2000,100),M,N2,byrow = TRUE)
      C_no_rec = matrix(rnorm(M*N2,2500,125),M,N2,byrow = TRUE)
      
      Q_rec = matrix(rnorm(M*N2,26,2),M,N2,byrow = TRUE)
      Q_rel = matrix(rnorm(M*N2,23,3),M,N2,byrow = TRUE)
      Q_no_rec = matrix(rnorm(M*N2,20,4),M,N2,byrow = TRUE)
      
      NB_nt = (lamda*(P_nt_rec*(1-P_nt_rel)*Q_rec + P_nt_rec*P_nt_rel*Q_rel
                     +(1-P_nt_rec)*Q_no_rec)
      - (P_nt_rec*(1-P_nt_rel)*C_rec + P_nt_rec*P_nt_rel*C_rel
         +(1-P_nt_rec)*C_no_rec ))
      NB_t1 = (lamda*(P_t1_rec*(1-P_t1_rel)*Q_rec + P_t1_rec*P_t1_rel*Q_rel
                     +(1-P_t1_rec)*Q_no_rec)
      - (C_t1 + P_t1_rec*(1-P_t1_rel)*C_rec + P_t1_rec*P_t1_rel*C_rel
         +(1-P_t1_rec)*C_no_rec ))
      NB_t2 = (lamda*(P_t2_rec*(1-P_t2_rel)*Q_rec + P_t2_rec*P_t2_rel*Q_rel
                     +(1-P_t2_rec)*Q_no_rec)
      - (C_t2 + P_t2_rec*(1-P_t2_rel)*C_rec + P_t2_rec*P_t2_rel*C_rel
         +(1-P_t2_rec)*C_no_rec ))
      
      Pb = colMeans(pmax(pmax(NB_nt,NB_t1),NB_t2))
      Pf = Pb-pmax(pmax(colMeans(NB_nt), colMeans(NB_t1)), colMeans(NB_t2))
      if(M==2)
      {
        Pc = Pb-0.5*(pmax(pmax(NB_nt[1,], NB_t1[1,]),NB_t2[1,])
                     +pmax(pmax(NB_nt[2,], NB_t1[2,]), NB_t2[2,]))
        
      } else {
      Pc = Pb-0.5*(pmax(pmax(colMeans(NB_nt[1:(M/2),]), colMeans(NB_t1[1:(M/2),])),colMeans(NB_t2[1:(M/2),]))
                   +pmax(pmax(colMeans(NB_nt[((M/2)+1):M,]), colMeans(NB_t1[((M/2)+1):M,])), colMeans(NB_t2[((M/2)+1):M,])))
      }
      sum1[1] = sum1[1] + sum(Pf-Pc);
      sum1[2] = sum1[2] + sum((Pf-Pc)^2);
      sum1[3] = sum1[3] + sum((Pf-Pc)^3);
      sum1[4] = sum1[4] + sum((Pf-Pc)^4);
      sum1[5] = sum1[5] + sum(Pf);
      sum1[6] = sum1[6] + sum(Pf^2);
      sum1[7] = sum1[7] + M*N2;
    }
    sum1
}