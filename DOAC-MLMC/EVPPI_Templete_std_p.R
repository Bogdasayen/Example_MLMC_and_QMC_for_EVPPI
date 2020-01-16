EVPPI_Templete_std_p<-function(M,N)
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
    temp = runif(N)
    Event.cost.samples[,m] = rep(temp, times = 1, each = M)
  }
  temp = rnorm(N)
  Event.cost.samples[,5] = rep(temp, times = 1, each = M)
  temp = rnorm(N)
  Event.cost.samples[,6] = rep(temp, times = 1, each = M)
  
  # Treatment cost (Outer samples, only N samples but M repeats) 
  # Dimension 1, uniform
  Treatment.cost.samples = matrix(NA,NN,1)
  temp = runif(N)
  Treatment.cost.samples[,1] = rep(temp, times = 1, each = M)
  
  # Healthstate cost
  # Dimension 2, normal (Outer samples, only N samples but M repeats)
  Health.cost.samples = matrix(NA,NN,2)
  temp = rnorm(N)
  Health.cost.samples[,1] = rep(temp, times = 1, each = M)
  temp = rnorm(N)
  Health.cost.samples[,2] = rep(temp, times = 1, each = M)
  
  
  # Switch Probability  (Outer samples, only N samples but M repeats)
  # Dimension 4, first 3 beta (0.1,0.9), last 1 beta (0.3,0.7)
  Switch.probability.samples = matrix(NA,NN,4)
  temp = rbeta(N,0.1,0.9)
  Switch.probability.samples[,1] = rep(temp, times = 1, each = M)
  temp = rbeta(N,0.1,0.9)
  Switch.probability.samples[,2] = rep(temp, times = 1, each = M)
  temp = rbeta(N,0.1,0.9)
  Switch.probability.samples[,3] = rep(temp, times = 1, each = M)
  temp = rbeta(N,0.3,0.7)
  Switch.probability.samples[,4] = rep(temp, times = 1, each = M)
  
  
  # Effects of previous events (Outer samples, only N samples but M repeats)
  # Dimension 24, all normal
  Effect.history.samples = matrix(NA,NN,24)
  for(m in 1:24){
    temp = rnorm(N)
    Effect.history.samples[,m] = rep(temp, times = 1, each = M)
  }
  
  # Random selection of MCMC samples  (Outer samples, only N samples but M repeats)
  # Dimension 3, selcet 3 set of MCMC samples
  temp = sample(1:29999, N, replace=TRUE)
  MCMC.selection = rep(temp, times = 1, each = M)
  MCMC.baseline.samples = bugs.baseline[MCMC.selection,]
  temp = sample(1:29999, N, replace=TRUE)
  MCMC.selection = rep(temp, times = 1, each = M)
  MCMC.loghr.samples = bugs.loghr[MCMC.selection,]
  temp = sample(1:59999, N, replace=TRUE)
  MCMC.selection = rep(temp, times = 1, each = M)
  MCMC.noTreatment.samples = hr.no.treatment[MCMC.selection,]
  
  # Utility factor for different age  (Outer samples, only N samples but M repeats)
  # Dimension 4, only for 65 and 75, all beta
  Utility.age.samples = matrix(NA,NN,4)
  temp = rbeta(N, 388.47, 109.57)
  Utility.age.samples[,1] = rep(temp, times = 1, each = M)
  temp = rbeta(N, 551.74, 155.62)
  Utility.age.samples[,2] = rep(temp, times = 1, each = M)
  temp = rbeta(N, 191.17, 63.72)
  Utility.age.samples[,3] = rep(temp, times = 1, each = M)
  temp = rbeta(N, 406.37, 165.98)
  Utility.age.samples[,4] = rep(temp, times = 1, each = M)
  
  # Health state utilities  (Outer samples, only N samples but M repeats)
  # Dimension 4, first 3 normal, last 1 beta
  Utility.state.samples = matrix(NA,NN,4)
  temp = rnorm(N)
  Utility.state.samples[,1] = rep(temp, times = 1, each = M)
  temp = rnorm(N)
  Utility.state.samples[,2] = rep(temp, times = 1, each = M)
  temp = rnorm(N)
  Utility.state.samples[,3] = rep(temp, times = 1, each = M)
  temp = rbeta(N, 3.941, 1.385)
  Utility.state.samples[,4] = rep(temp, times = 1, each = M)
  
  # Events utilities  (Outer samples, only N samples but M repeats)
  # Dimension 6, first 3 uniform, last 3 normal
  Utility.event.samples = matrix(NA,NN,6)
  temp = runif(N)
  Utility.event.samples[,1] = rep(temp, times = 1, each = M)
  temp = runif(N)
  Utility.event.samples[,2] = rep(temp, times = 1, each = M)
  temp = runif(N)
  Utility.event.samples[,3] = rep(temp, times = 1, each = M)
  temp = rnorm(N)
  Utility.event.samples[,4] = rep(temp, times = 1, each = M)
  temp = rnorm(N)
  Utility.event.samples[,5] = rep(temp, times = 1, each = M)
  temp = rnorm(N)
  Utility.event.samples[,6] = rep(temp, times = 1, each = M)
  
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