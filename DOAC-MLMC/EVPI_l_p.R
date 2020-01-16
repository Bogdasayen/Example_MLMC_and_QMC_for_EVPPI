EVPI_l_p<-function(l,N)
{
  # this function generate all the random samples for the EVPPI calculation
  # 34 normal, 8 uniform, 15 beta, MCMC 7, MCMC 28, MCMC 7 
  # then pass all the samples to new age.independent.generate.probabilities_2 function to get
  # the samples of parameters and then calculate the net benefit function
  # N is the number of outer samples, M=2^l is the number of inner samples
  sum1 <- rep(0, 7)
  M = 2^(l+1)
  inputs = 1:ceiling(M*N/100)
  results <- foreach(i=inputs,.export=c('bugs.baseline', 'bugs.loghr','hr.no.treatment','age.independent.generate.probabilities_2','n.events',
                                        'event.names','event.state.codes','hr.future.stroke','hr.future.death',
                                        'hr.future.bleed','hr.future.ich','hr.future.tiase','lifetables','n.treatments',
                                        'treatment.names','state.names','n.health.states','noac.net.benefit','n.states',
                                        'initial.age','n.cycles','generate.probabilities','t.names','hr.death.age',
                                        'nondeath.health.states','event.codes','next.state.name','treatment.switch.indices',
                                        'ich.treatment.switch.indices','mi.treatment.switch.indices','lambdas')) %dopar% {
    NN=min(100, N*M-(i-1)*100)
    Event.cost.samples = matrix(NA,NN,6)
    for(m in 1:4){
      Event.cost.samples[,m] = runif(NN)
    }
    Event.cost.samples[,5] = rnorm(NN)
    Event.cost.samples[,6] = rnorm(NN)
    
    # Treatment cost
    # Dimension 1, uniform
    Treatment.cost.samples = matrix(NA,NN,1)
    Treatment.cost.samples[,1] = runif(NN)
    
    # Healthstate cost
    # Dimension 2, normal
    Health.cost.samples = matrix(NA,NN,2)
    Health.cost.samples[,1] = rnorm(NN)
    Health.cost.samples[,2] = rnorm(NN)
    
    
    # Switch Probability
    # Dimension 4, first 3 beta (0.1,0.9), last 1 beta (0.3,0.7)
    Switch.probability.samples = matrix(NA,NN,4)
    Switch.probability.samples[,1] = rbeta(NN,0.1,0.9)
    Switch.probability.samples[,2] = rbeta(NN,0.1,0.9)
    Switch.probability.samples[,3] = rbeta(NN,0.1,0.9)
    Switch.probability.samples[,4] = rbeta(NN,0.3,0.7)
    
    
    # Effects of previous events
    # Dimension 24, all normal
    Effect.history.samples = matrix(NA,NN,24)
    for(m in 1:24){
      Effect.history.samples[,m] = rnorm(NN)
    }
    
    # Random selection of MCMC samples
    # Dimension 3, selcet 3 set of MCMC samples
    MCMC.selection = sample(1:29999, NN, replace = TRUE)
    MCMC.baseline.samples = bugs.baseline[MCMC.selection,]
    MCMC.selection = sample(1:29999, NN, replace = TRUE)
    MCMC.loghr.samples = bugs.loghr[MCMC.selection,]
    MCMC.selection = sample(1:59999, NN, replace = TRUE)
    MCMC.noTreatment.samples = hr.no.treatment[MCMC.selection,]
    
    # Utility factor for different age
    # Dimension 4, only for 65 and 75, all beta
    Utility.age.samples = matrix(NA,NN,4)
    Utility.age.samples[,1] = rbeta(NN,388.47,109.57)
    Utility.age.samples[,2] = rbeta(NN,551.74,155.62)
    Utility.age.samples[,3] = rbeta(NN,191.17,63.72)
    Utility.age.samples[,4] = rbeta(NN,406.37,165.98)
    
    # Health state utilities
    # Dimension 4, first 3 normal, last 1 beta
    Utility.state.samples = matrix(NA,NN,4)
    Utility.state.samples[,1] = rnorm(NN)
    Utility.state.samples[,2] = rnorm(NN)
    Utility.state.samples[,3] = rnorm(NN)
    Utility.state.samples[,4] = rbeta(NN,3.941,1.385)
    
    # Events utilities
    # Dimension 6, first 3 uniform, last 3 normal
    Utility.event.samples = matrix(NA,NN,6)
    Utility.event.samples[,1] = runif(NN)
    Utility.event.samples[,2] = runif(NN)
    Utility.event.samples[,3] = runif(NN)
    Utility.event.samples[,4] = rnorm(NN)
    Utility.event.samples[,5] = rnorm(NN)
    Utility.event.samples[,6] = rnorm(NN)
    
    age.independent.samples<-age.independent.generate.probabilities_2(NN,
                                                                      Event.cost.samples,Treatment.cost.samples,Health.cost.samples,
                                                                      Switch.probability.samples,Effect.history.samples,
                                                                      MCMC.baseline.samples,MCMC.loghr.samples,MCMC.noTreatment.samples,
                                                                      Utility.age.samples,Utility.state.samples,Utility.event.samples
    )
    
    model.outputs<-noac.net.benefit(n.samples=NN,n.cycles=n.cycles,initial.age=initial.age,lambdas=lambdas,age.independent.samples=age.independent.samples)
    
    NetB = model.outputs$NB
  }
  
    NetB = matrix(NA,M*N,5)
  for(i in inputs){
    NN <- min(100, N*M-(i-1)*100)
    nn = min(i*100,N*M)
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
      if(M==2){
        NetB_low_c_1[,n] = temp[,1]
        NetB_low_c_2[,n] = temp[,2]
      }else{
        NetB_low_c_1[,n] = apply(temp[,1:(M/2)],1,mean)
        NetB_low_c_2[,n] = apply(temp[,(M/2+1):M],1,mean)
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










