

EVPPI_Complex_QMC_MVN<-function(N,M,K)
{
    inputs <- seq(1,K)
    
    tic()
    results <- foreach(i=inputs,.export=c('bugs.baseline', 'bugs.loghr','hr.no.treatment','age.independent.generate.probabilities_2',
                                          'n.events','event.names','event.state.codes','hr.future.stroke','hr.future.death',
                                          'hr.future.bleed','hr.future.ich','hr.future.tiase','lifetables','n.treatments','treatment.names',
                                          'state.names','n.health.states','noac.net.benefit','n.states','initial.age','n.cycles',
                                          'generate.probabilities','t.names','hr.death.age','nondeath.health.states','event.codes',
                                          'next.state.name','treatment.switch.indices','ich.treatment.switch.indices',
                                          'mi.treatment.switch.indices','lambdas','sobol','which.is.max','Rbeta.inv','sort_loghr',
                                          'sort_baseline','sort_notreat','mean_baseline','L_baseline','mean_loghr','L_loghr'))  %dorng%  {
                                              
                                              U <- sobol(N, dim = 53, init = TRUE,  scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)
                                              U <- sobol(N, dim = 53, init = FALSE, scrambling = 1, seed = sample(seq(1,10^6),1), normal = FALSE)
                                              
                                              N_path <- 1024
                                              NetB <- matrix(NA,M*N,5)
                                              inputs <- 1:ceiling(N*M/N_path)
                                              
                                              for(j in inputs) {
                                                  NN <- min(N_path, M*N-(j-1)*N_path)
                                                  index <- seq(1+(j-1)*N_path/M, min(j*N_path/M,N))
                                                  
                                                  # Event cost 
                                                  # Dimension 6, first 4 are normal and last 2 are uniform, for the facility of QMC version
                                                  Event.cost.samples = matrix(NA,NN,6)
                                                  
                                                  for(m in 1:4){
                                                      temp = U[,m+21]
                                                      Event.cost.samples[,m] = rep(temp, times = 1, each = M)
                                                  }
                                                  temp = qnorm(U[,26])
                                                  Event.cost.samples[,5] = rep(temp, times = 1, each = M)
                                                  temp = qnorm(U[,27])
                                                  Event.cost.samples[,6] = rep(temp, times = 1, each = M)
                                                  
                                                  # Treatment cost (Inner samples, only N samples but M repeats) 
                                                  # Dimension 1, uniform
                                                  Treatment.cost.samples = matrix(NA,NN,1)
                                                  Treatment.cost.samples[,1] = runif(NN)
                                                  
                                                  # Healthstate cost
                                                  # Dimension 2, normal (Outer samples, only N samples but M repeats)
                                                  Health.cost.samples = matrix(NA,NN,2)
                                                  temp = qnorm(U[,28])
                                                  Health.cost.samples[,1] = rep(temp, times = 1, each = M)
                                                  temp = qnorm(U[,29])
                                                  Health.cost.samples[,2] = rep(temp, times = 1, each = M)
                                                  
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
                                                  
                                                  MCMC.baseline.samples = bugs.baseline[seq(1,length(index)*M),]
                                                  X_baseline <- qnorm(U[index,15:21])
                                                  temp <- rep(1,length(index))%*% t(mean_baseline) + X_baseline%*%L_baseline
                                                  MCMC.baseline.samples[,seq(1,7)] <- matrix(rep(temp,times=1,each =M),NN,7)
                                                  
                                                  
                                                  MCMC.loghr.samples = bugs.loghr[seq(1,length(index)*M),]
                                                  X_loghr <- qnorm(U[index,c(1:14,40:53)])
                                                  temp <- rep(1,length(index))%*% t(mean_loghr) + X_loghr%*%L_loghr
                                                  MCMC.loghr.samples[,seq(8,35)] <- matrix(rep(temp,times=1,each =M),NN,28)
                                                  
                                                  # index_con <- (1:14)
                                                  # nindex <- setdiff(1:28,index_con)
                                                  # sigma_coef <- cov_loghr[nindex,index_con]%*%solve(cov_loghr[index_con,index_con])
                                                  # sigma_cond <- cov_loghr[nindex,nindex] - sigma_coef%*%cov_loghr[index_con,nindex]
                                                  # 
                                                  # mu_con <- matrix(mean_loghr[nindex],length(nindex),NN) + sigma_coef%*%(t(MCMC.loghr.samples[,index_con+7])-matrix(mean_loghr[index_con],length(index_con),NN))
                                                  # 
                                                  # MCMC.loghr.samples[,nindex+7] <-  t(mu_con)+ mvrnorm(NN, rep(0,length(nindex)), sigma_cond)
                                                  
                                                  
                                                  MCMC.selection = sample(1:65536, NN, replace = TRUE)
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
                                                  temp = qnorm(U[index,30])
                                                  Utility.state.samples[,1] = rep(temp, times = 1, each = M)
                                                  temp = qnorm(U[index,31])
                                                  Utility.state.samples[,2] = rep(temp, times = 1, each = M)
                                                  temp = qnorm(U[index,32])
                                                  Utility.state.samples[,3] = rep(temp, times = 1, each = M)
                                                  temp = Rbeta.inv(U[index,33], 3.941,1.385, log=FALSE)
                                                  Utility.state.samples[,4] = rep(temp, times = 1, each = M)
                                                  
                                                  # Events utilities
                                                  # Dimension 6, first 3 uniform, last 3 normal
                                                  Utility.event.samples = matrix(NA,NN,6)
                                                  temp = U[index,34]
                                                  Utility.event.samples[,1] = rep(temp, times = 1, each = M)
                                                  temp = U[index,35]
                                                  Utility.event.samples[,2] = rep(temp, times = 1, each = M)
                                                  temp = U[index,36]
                                                  Utility.event.samples[,3] = rep(temp, times = 1, each = M)
                                                  temp = U[index,37]
                                                  Utility.event.samples[,4] = rep(temp, times = 1, each = M)
                                                  temp = U[index,38]
                                                  Utility.event.samples[,5] = rep(temp, times = 1, each = M)
                                                  temp = U[index,39]
                                                  Utility.event.samples[,6] = rep(temp, times = 1, each = M)
                                                  
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
    cat(sprintf("QMC (EVPPI_HR_loghr): \nDIFF  = %.4f +/- %.4f,  std_DIFF  = %.4f.\nEVPPI = %.4f +/- %.4f, std_EVPPI = %.4f.\nN=%d, M=%d, K=%d. \n\n ", DIFF,3*std_DIFF,std_DIFF,EVPPI,3*std_EVPPI,std_EVPPI,N,M,K))
    return(list(DIFF = c(DIFF, std_DIFF), EVPPI = c(EVPPI, std_EVPPI)))
}