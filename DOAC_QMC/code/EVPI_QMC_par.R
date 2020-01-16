
EVPI_QMC_par<-function(N,M)
{
    # N <= 1e5, M<=16
    tic()
    inputs = seq(1,M)
    results <- foreach(i=inputs,.export=c('bugs.baseline', 'bugs.loghr','hr.no.treatment','age.independent.generate.probabilities_2',
                                          'n.events','event.names','event.state.codes','hr.future.stroke','hr.future.death',
                                          'hr.future.bleed','hr.future.ich','hr.future.tiase','lifetables','n.treatments','treatment.names',
                                          'state.names','n.health.states','noac.net.benefit','n.states','initial.age','n.cycles',
                                          'generate.probabilities','t.names','hr.death.age','nondeath.health.states','event.codes',
                                          'next.state.name','treatment.switch.indices','ich.treatment.switch.indices',
                                          'mi.treatment.switch.indices','lambdas','sobol','which.is.max','Rbeta.inv','sort_loghr',
                                          'sort_baseline','sort_notreat')) %dopar% {
              
                  U <- sobol(N, dim = 57, init = TRUE, scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)
                  U <- sobol(N, dim = 57, init = FALSE, scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)
                  
                  U <- U[seq(1,N),]
                  
                  sum1 <- 0
          
                  for(N1 in seq(1, N, by=1000)) {
                      NN <- min(1000, N-N1+1)
                      index <- seq(N1,N1+NN-1)
                      # Event cost 
                      # Dimension 6, first 4 are normal and last 2 are uniform, for the facility of QMC version
                      Event.cost.samples = matrix(NA,NN,6)
                      for(m in 1:4){
                          Event.cost.samples[,m] = U[index,m]
                      }
                      Event.cost.samples[,5] = qnorm(U[index,5])
                      Event.cost.samples[,6] = qnorm(U[index,6])
                      
                      # Treatment cost
                      # Dimension 1, uniform
                      Treatment.cost.samples = matrix(NA,NN,1)
                      Treatment.cost.samples[,1] = U[index,7]
                      
                      # Healthstate cost
                      # Dimension 2, normal
                      Health.cost.samples = matrix(NA,NN,2)
                      Health.cost.samples[,1] = qnorm(U[index,8])
                      Health.cost.samples[,2] = qnorm(U[index,9])
                      
                      
                      # Switch Probability
                      # Dimension 4, first 3 beta (0.1,0.9), last 1 beta (0.3,0.7)
                      Switch.probability.samples = matrix(NA,NN,4)
                      Switch.probability.samples[,1] = Rbeta.inv(U[index,16], 0.1, 0.9, log=FALSE)
                      Switch.probability.samples[,2] = Rbeta.inv(U[index,17], 0.1, 0.9, log=FALSE) 
                      Switch.probability.samples[,3] = Rbeta.inv(U[index,18], 0.1, 0.9, log=FALSE) 
                      Switch.probability.samples[,4] = Rbeta.inv(U[index,19], 0.3, 0.7, log=FALSE) 
                      
                      
                      # Effects of previous events
                      # Dimension 24, all normal
                      Effect.history.samples = matrix(NA,NN,24)
                      for(m in 1:24){
                          Effect.history.samples[,m] = qnorm(U[index,19+m])
                      }
                      
                      # Random selection of MCMC samples
                      # Dimension 3, selcet 3 set of MCMC samples
                      
                      MCMC.selection = sort_loghr(U[index,c(10,11)])$QMCsample_loghr
                      MCMC.loghr.samples = bugs.loghr[MCMC.selection,]
                      MCMC.selection = sort_baseline(U[index,c(12,13)])$QMCsample_baseline
                      MCMC.baseline.samples = bugs.baseline[MCMC.selection,]
                      MCMC.selection = sort_notreat(U[index,c(14,15)])$QMCsample_notreat
                      MCMC.noTreatment.samples = hr.no.treatment[MCMC.selection,]
                      
                      # MCMC.selection = sample(1:16384, NN, replace = TRUE)
                      # MCMC.loghr.samples = bugs.loghr[MCMC.selection,]
                      # MCMC.selection = sample(1:16384, NN, replace = TRUE)
                      # MCMC.baseline.samples = bugs.baseline[MCMC.selection,]
                      # MCMC.selection = sample(1:65536, NN, replace = TRUE)
                      # MCMC.noTreatment.samples = hr.no.treatment[MCMC.selection,]
                      
                      # Utility factor for different age
                      # Dimension 4, only for 65 and 75, all beta
                      Utility.age.samples = matrix(NA,NN,4)
                      Utility.age.samples[,1] = Rbeta.inv(U[index,44], 388.47,109.57, log=FALSE) #rbeta(NN,388.47,109.57)
                      Utility.age.samples[,2] = Rbeta.inv(U[index,45], 551.74,155.62, log=FALSE) #rbeta(NN,551.74,155.62)
                      Utility.age.samples[,3] = Rbeta.inv(U[index,46], 191.17,63.72, log=FALSE) #rbeta(NN,191.17,63.72)
                      Utility.age.samples[,4] = Rbeta.inv(U[index,47], 406.37,165.98, log=FALSE) #rbeta(NN,406.37,165.98)
                      
                      # Health state utilities
                      # Dimension 4, first 3 normal, last 1 beta
                      Utility.state.samples = matrix(NA,NN,4)
                      Utility.state.samples[,1] = qnorm(U[index,48])
                      Utility.state.samples[,2] = qnorm(U[index,49])
                      Utility.state.samples[,3] = qnorm(U[index,50])
                      Utility.state.samples[,4] = Rbeta.inv(U[index,51], 3.941,1.385, log=FALSE) #rbeta(NN,3.941,1.385)
                      
                      # Events utilities
                      # Dimension 6, first 3 uniform, last 3 normal
                      Utility.event.samples = matrix(NA,NN,6)
                      Utility.event.samples[,1] = U[index,52]
                      Utility.event.samples[,2] = U[index,53]
                      Utility.event.samples[,3] = U[index,54]
                      Utility.event.samples[,4] = qnorm(U[index,55])
                      Utility.event.samples[,5] = qnorm(U[index,56])
                      Utility.event.samples[,6] = qnorm(U[index,57])
                      
                      # New function to generate the probabilities in the Markov
                      age.independent.samples<-age.independent.generate.probabilities_2(NN,
                                                                                        Event.cost.samples,Treatment.cost.samples,Health.cost.samples,
                                                                                        Switch.probability.samples,Effect.history.samples,
                                                                                        MCMC.baseline.samples,MCMC.loghr.samples,MCMC.noTreatment.samples,
                                                                                        Utility.age.samples,Utility.state.samples,Utility.event.samples
                      )
                      # Calculate the net Benefits of each treatment for each sample
                      model.outputs<-noac.net.benefit(n.samples=NN,n.cycles=n.cycles,initial.age=initial.age,lambdas=lambdas,age.independent.samples=age.independent.samples)
                      NetB <- model.outputs$NB
                      
                      Ind <- which.is.max(colMeans(NetB)) # find the optimal treatment without any information 2
                      EVPI_sample <- apply(NetB,1,max)-NetB[,Ind,]
                      
                      sum1 <- sum1+sum(EVPI_sample)
                  }
                  EVPI = sum1/N
              }
    toc()
    
    EVPI <- unlist(results)
    EVPI_2 <- EVPI*EVPI
    
    EVPI_QMC <- sum(EVPI)/M
    EVPI_var_QMC <- sum(EVPI_2)/M - EVPI_QMC^2
    std_EVPI <- sqrt(EVPI_var_QMC/M)
    cat(sprintf("QMC: EVPI = %.4f +/- %.4f, std = %.4f, N=%.4e, M=%d. \n\n ", EVPI_QMC,3*std_EVPI,std_EVPI,N,M))
    
    # stopCluster(cl)
    # stopImplicitCluster()
    return(list(EVPI_QMC=EVPI_QMC,std_EVPI=std_EVPI))
    
}
