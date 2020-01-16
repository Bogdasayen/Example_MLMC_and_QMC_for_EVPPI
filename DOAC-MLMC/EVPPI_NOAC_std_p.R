EVPPI_NOAC_std_p<-function(M,N)
{
  # this function generate all the random samples for the EVPPI calculation
  # 34 normal, 8 uniform, 15 beta, MCMC 7, MCMC 28, MCMC 7 
  # then pass all the samples to new age.independent.generate.probabilities_2 function to get
  # the samples of parameters and then calculate the net benefit function
  # N is the number of outer samples, M=2^l is the number of inner samples
  NN <- M*N
  # Event cost (Outer samples, only N samples but M repeats)
  # Dimension 6, first 4 are normal and last 2 are uniform, for the facility of QMC version
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
  
  # Random selection of MCMC samples  (Outer samples, only N samples but M repeats)
  # Dimension 3, selcet 3 set of MCMC samples
  MCMC.selection = sample(1:29999, NN, replace = TRUE)
  MCMC.baseline.samples = bugs.baseline[MCMC.selection,]
  MCMC.selection = sample(1:59999, NN, replace = TRUE)
  MCMC.noTreatment.samples = hr.no.treatment[MCMC.selection,]
  
  # NOAC samples
  temp = sample(1:29999, N, replace=TRUE)
  MCMC.selection = rep(temp, times = 1, each = M)
  MCMC.loghr.samples = bugs.loghr[MCMC.selection,]
  
  index <- (1:14)
  nindex <- setdiff(1:28,index)
  sigma_coef <- sigma_loghr[nindex,index]%*%solve(sigma_loghr[index,index])
  sigma_cond <- sigma_loghr[nindex,nindex] - sigma_coef%*%sigma_loghr[index,nindex]
  
  mu_con <- matrix(mu_loghr[nindex],length(nindex),NN) + sigma_coef%*%(t(MCMC.loghr.samples[,index+7])-matrix(mu_loghr[index],length(index),NN))
  
  MCMC.loghr.samples[,nindex+7] <-  t(mu_con)+ mvrnorm(NN, rep(0,length(nindex)), sigma_cond)
  
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
  return(NetB)
}