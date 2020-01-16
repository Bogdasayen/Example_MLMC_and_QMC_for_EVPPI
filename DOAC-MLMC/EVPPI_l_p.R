EVPPI_l_p<-function(l,N)
{
  # this function generate all the random samples for the EVPPI calculation
  # 34 normal, 8 uniform, 15 beta, MCMC 7, MCMC 28, MCMC 7 
  # then pass all the samples to new age.independent.generate.probabilities_2 function to get
  # the samples of parameters and then calculate the net benefit function
  # N is the number of outer samples, M=2^l is the number of inner samples
  sum1 <- rep(0, 7)
  M = 2^(l+1)
  Np = max(M,128)
   
  inputs = 1:ceiling(M*N/Np)

  
  results <- foreach(i=inputs,.export=c('bugs.baseline', 'bugs.loghr','hr.no.treatment','age.independent.generate.probabilities_2','n.events',
                                        'event.names','event.state.codes','hr.future.stroke','hr.future.death',
                                        'hr.future.bleed','hr.future.ich','hr.future.tiase','lifetables','n.treatments',
                                        'treatment.names','state.names','n.health.states','noac.net.benefit','n.states',
                                        'initial.age','n.cycles','generate.probabilities','t.names','hr.death.age',
                                        'nondeath.health.states','event.codes','next.state.name','treatment.switch.indices',
                                        'ich.treatment.switch.indices','mi.treatment.switch.indices','lambdas',
                                        'EVPPI_Cost_EventCost_std_p','sigma_loghr','mu_loghr'),.packages='MASS') %dorng% {
  NN=min(Np, N*M-(i-1)*Np)
  #####################################################
  EVPPI_Cost_EventCost_std_p(M,NN/M)  # if you want to calculate EVPPI for other parameters, you need to change this line using other corresponding function
  # and change the function name in the last line (line 22) of the condition of foreach to the same one.
  ###########################################33########
  }                        
  NetB = matrix(NA,M*N,5)
  for(i in inputs){
    NN <- min(Np, N*M-(i-1)*Np)
    nn = min(i*Np,N*M)
    NetB[(nn-NN+1):nn,] = matrix(unlist(results[i]),NN,5)
  }
  NetB_max = apply(NetB,1,max)
  NetB_max_sample = apply(matrix(NetB_max,N,M,byrow=TRUE),1,mean)
  NetB_low_c_1 = matrix(NA,N,5)
  NetB_low_c_2 = matrix(NA,N,5)
  NetB_low_f = matrix(NA,N,5)
  for(n in 1:5){
    temp = matrix(NetB[,n],N,M,byrow=TRUE)
    NetB_low_f[,n] = apply(temp,1,mean)
    if(l==0){
      NetB_low_c_1[,n] = temp[,1]
      NetB_low_c_2[,n] = temp[,2]
    }else{
      NetB_low_c_1[,n] = apply(temp[,1:max(M/2,2)],1,mean)
      NetB_low_c_2[,n] = apply(temp[,min(M/2+1,M-1):M],1,mean)
    }
  }
  NetB_low_f_sample = apply(NetB_low_f,1,max)
  NetB_low_c_sample = (apply(NetB_low_c_1,1,max)+apply(NetB_low_c_2,1,max))/2
  Pf = NetB_max_sample - NetB_low_f_sample
  Pc = NetB_max_sample - NetB_low_c_sample
  
  sum1[1] = sum1[1] + sum(Pf-Pc);
  sum1[2] = sum1[2] + sum((Pf-Pc)^2);
  sum1[3] = sum1[3] + sum((Pf-Pc)^3);
  sum1[4] = sum1[4] + sum((Pf-Pc)^4);
  sum1[5] = sum1[5] + sum(Pf);
  sum1[6] = sum1[6] + sum(Pf^2);
  sum1[7] = sum1[7] + M*N;
  
  return(sum1)
}
