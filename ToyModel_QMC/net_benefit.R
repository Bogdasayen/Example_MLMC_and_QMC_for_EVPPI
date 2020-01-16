
net_benefit <- function(lambda,P_rec,P_rel,C_rec,C_rel,C_norec,Q_rec,Q_rel,Q_norec,C_t){
  QALY <- P_rec*(1-P_rel)*Q_rec + P_rec*P_rel*Q_rel + (1-P_rec)*Q_norec
  Cost <- P_rec*(1-P_rel)*C_rec + P_rec*P_rel*C_rel + (1-P_rec)*C_norec 
  
  Cost[,2] <- Cost[,2] +C_t[2]
  Cost[,3] <- Cost[,3] +C_t[3]
  
  NB <- lambda*QALY - Cost
  
  return( list(NB=NB,QALY=QALY,COST=Cost))
}