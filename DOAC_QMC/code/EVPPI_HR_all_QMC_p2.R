
EVPPI_HR_all_QMC_p2<-function(N,M,K)
{
    # this function generate all the random samples for the EVPPI calculation
    # 34 normal, 8 uniform, 15 beta, MCMC 7, MCMC 28, MCMC 7 
    # then pass all the samples to new age.independent.generate.probabilities_2 function to get
    # the samples of parameters and then calculate the net benefit function
    # N is the number of outer samples, M=2^l is the number of inner samples
    
    # source("age.independent.generate.probabilities_2.R")
    inputs <- seq(1,K)
    
    tic()
    results <- foreach(i=inputs,.export=c('bugs.baseline', 'bugs.loghr','hr.no.treatment','age.independent.generate.probabilities_2',
                                          'n.events','event.names','event.state.codes','hr.future.stroke','hr.future.death',
                                          'hr.future.bleed','hr.future.ich','hr.future.tiase','lifetables','n.treatments','treatment.names',
                                          'state.names','n.health.states','noac.net.benefit','n.states','initial.age','n.cycles',
                                          'generate.probabilities','t.names','hr.death.age','nondeath.health.states','event.codes',
                                          'next.state.name','treatment.switch.indices','ich.treatment.switch.indices',
                                          'mi.treatment.switch.indices','lambdas','sobol','which.is.max','Rbeta.inv','sort_loghr',
                                          'sort_baseline','sort_notreat'))  %dorng%  {
                                              
                                              U <- sobol(N, dim = 6, init = TRUE,  scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)
                                              U <- sobol(N, dim = 6, init = FALSE, scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)
                                              
                                              N_path <- 1024
                                              NetB <- matrix(NA,M*N,5)
                                              inputs <- 1:ceiling(N*M/N_path)
                                              
                                              for(j in inputs) {
                                                  NN <- min(N_path, M*N-(j-1)*N_path)
                                                  index <- seq(1+(j-1)*N_path/M, min(j*N_path/M,N))
                                                  
                                                  U2 <- sobol(M, dim = 51, init = TRUE,  scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)
                                                  U2 <- matrix(0,NN,51)
                                                  
                                                  for(k in seq(1,NN/M)) {
                                                      U2[seq((k-1)*M+1,k*M),] <- sobol(M, dim = 51, init = FALSE,  scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)
                                                  }
                                                  
                                                  # Event cost 
                                                  # Dimension 6, first 4 are normal and last 2 are uniform, for the facility of QMC version
                                                  Event.cost.samples = matrix(NA,NN,6)
                                                  
                                                  for(m in 1:4){
                                                      Event.cost.samples[,m] = U2[,m]
                                                  }
                                                  Event.cost.samples[,5] = qnorm(U2[,5])
                                                  Event.cost.samples[,6] = qnorm(U2[,6])
                                                  
                                                  # Treatment cost (Inner samples, only N samples but M repeats) 
                                                  # Dimension 1, uniform
                                                  Treatment.cost.samples = matrix(NA,NN,1)
                                                  Treatment.cost.samples[,1] = U2[,7]
                                                  
                                                  # Healthstate cost
                                                  # Dimension 2, normal (Outer samples, only N samples but M repeats)
                                                  Health.cost.samples = matrix(NA,NN,2)
                                                  Health.cost.samples[,1] = qnorm(U2[,8])
                                                  Health.cost.samples[,2] = qnorm(U2[,9])
                                                  
                                                  # Switch Probability
                                                  # Dimension 4, first 3 beta (0.1,0.9), last 1 beta (0.3,0.7)
                                                  Switch.probability.samples = matrix(NA,NN,4)
                                                  Switch.probability.samples[,1] = Rbeta.inv(U2[,10], 0.1, 0.9, log=FALSE)
                                                  Switch.probability.samples[,2] = Rbeta.inv(U2[,11], 0.1, 0.9, log=FALSE) 
                                                  Switch.probability.samples[,3] = Rbeta.inv(U2[,12], 0.1, 0.9, log=FALSE) 
                                                  Switch.probability.samples[,4] = Rbeta.inv(U2[,13], 0.3, 0.7, log=FALSE)
                                                  
                                                  # Effects of previous events
                                                  # Dimension 24, all normal
                                                  Effect.history.samples = matrix(NA,NN,24)
                                                  for(m in 1:24){
                                                      Effect.history.samples[,m] = qnorm(U2[,13+m])
                                                  }
                                                  
                                                  # Random selection of MCMC samples
                                                  # Dimension 3, selcet 3 set of MCMC samples
                                                  temp = sort_baseline(U[index,c(3,4)])$QMCsample_baseline
                                                  MCMC.selection = rep(temp, times = 1, each = M)
                                                  MCMC.baseline.samples = bugs.baseline[MCMC.selection,]
                                                  
                                                  temp = sort_loghr(U[index,c(1,2)])$QMCsample_loghr
                                                  MCMC.selection = rep(temp, times = 1, each = M)
                                                  MCMC.loghr.samples = bugs.loghr[MCMC.selection,]
                                                  
                                                  temp = sort_notreat(U[index,c(5,6)])$QMCsample_notreat
                                                  MCMC.selection = rep(temp, times = 1, each = M)
                                                  MCMC.noTreatment.samples = hr.no.treatment[MCMC.selection,]
                                                  
                                                  # Utility factor for different age
                                                  # Dimension 4, only for 65 and 75, all beta
                                                  Utility.age.samples = matrix(NA,NN,4)
                                                  Utility.age.samples[,1] = Rbeta.inv(U2[,38], 388.47,109.57, log=FALSE) #rbeta(NN,388.47,109.57)
                                                  Utility.age.samples[,2] = Rbeta.inv(U2[,39], 551.74,155.62, log=FALSE) #rbeta(NN,551.74,155.62)
                                                  Utility.age.samples[,3] = Rbeta.inv(U2[,40], 191.17,63.72, log=FALSE) #rbeta(NN,191.17,63.72)
                                                  Utility.age.samples[,4] = Rbeta.inv(U2[,41], 406.37,165.98, log=FALSE) #rbeta(NN,406.37,165.98)
                                                  
                                                  # Health state utilities
                                                  # Dimension 4, first 3 normal, last 1 beta
                                                  Utility.state.samples = matrix(NA,NN,4)
                                                  Utility.state.samples[,1] = qnorm(U2[,42])
                                                  Utility.state.samples[,2] = qnorm(U2[,43])
                                                  Utility.state.samples[,3] = qnorm(U2[,44])
                                                  Utility.state.samples[,4] = Rbeta.inv(U2[,45], 3.941,1.385, log=FALSE) #rbeta(NN,3.941,1.385)
                                                  
                                                  
                                                  # Events utilities
                                                  # Dimension 6, first 3 uniform, last 3 normal
                                                  Utility.event.samples = matrix(NA,NN,6)
                                                  Utility.event.samples[,1] = U2[,46]
                                                  Utility.event.samples[,2] = U2[,47]
                                                  Utility.event.samples[,3] = U2[,48]
                                                  Utility.event.samples[,4] = qnorm(U2[,49])
                                                  Utility.event.samples[,5] = qnorm(U2[,50])
                                                  Utility.event.samples[,6] = qnorm(U2[,51])
                                                  
                                                  age.independent.samples<-age.independent.generate.probabilities_2(NN,
                                                                                                                    Event.cost.samples,Treatment.cost.samples,Health.cost.samples,
                                                                                                                    Switch.probability.samples,Effect.history.samples,
                                                                                                                    MCMC.baseline.samples,MCMC.loghr.samples,MCMC.noTreatment.samples,
                                                                                                                    Utility.age.samples,Utility.state.samples,Utility.event.samples
                                                  )
                                                  
                                                  model.outputs<-noac.net.benefit(n.samples=NN,n.cycles=n.cycles,initial.age=initial.age,lambdas=lambdas,age.independent.samples=age.independent.samples)
                                                  
                                                  index2 <- seq(1+(j-1)*N_path, min(j*N_path,M*N))
                                                  NetB[index2,] <- matrix(unlist(model.outputs$NB),NN,5)
                                              }
                                              
                                              # Calculate net benefit of the optimal decision for each sample 
                                              NetB_max = apply(NetB,1,max)
                                              NetB_max_sample = apply(matrix(NetB_max,N,M,byrow=TRUE),1,mean)
                                              
                                              # find the optimal decision without any information and calculate the Net benefit
                                              Ind = which.is.max(colMeans(NetB))
                                              NetB_low_sample = apply(matrix(NetB[,Ind],N,M,byrow=TRUE),1,mean)
                                              
                                              # find the optimal decision with partial information and calculate the net benefit
                                              NetB_mid = matrix(NA,N,5)
                                              for(i in 1:5){
                                                  NetB_mid[,i] = apply(matrix(NetB[,i],N,M,byrow=TRUE),1,mean)
                                              }
                                              NetB_mid_sample = apply(NetB_mid,1,max)
                                              
                                              # Construct the samples for both DIFF and EVPPI itself
                                              DIFF_sample = NetB_max_sample- NetB_mid_sample
                                              EVPPI_sample = NetB_mid_sample- NetB_low_sample
                                              
                                              # Calculate the estimated value and Monte Carlo standard deviation
                                              DIFF = mean(DIFF_sample)
                                              EVPPI = mean(EVPPI_sample)
                                              list(DIFF_p = DIFF, EVPPI_p = EVPPI)
                                          }
    
    Res <- matrix(unlist(results),2,K)
    
    DIFF <- mean(Res[1,])
    std_DIFF <- sqrt((mean(Res[1,]*Res[1,]) - DIFF^2)/K)
    EVPPI <- mean(Res[2,])
    std_EVPPI <- sqrt((mean(Res[2,]*Res[2,]) - EVPPI^2)/K)
    
    toc()
    cat(sprintf("QMC (EVPPI_HR_all): \nDIFF  = %.4f +/- %.4f,  std_DIFF  = %.4f.\nEVPPI = %.4f +/- %.4f, std_EVPPI = %.4f.\nN=%d, M=%d, K=%d. \n\n ", DIFF,3*std_DIFF,std_DIFF,EVPPI,3*std_EVPPI,std_EVPPI,N,M,K))
    return(list(DIFF = c(DIFF, std_DIFF), EVPPI = c(EVPPI, std_EVPPI)))
}