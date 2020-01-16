age.independent.generate.probabilities_2 <- function(n.samples,
  Event.cost.samples,Treatment.cost.samples,Health.cost.samples,
  Switch.probability.samples,Effect.history.samples,
  MCMC.baseline.samples,MCMC.loghr.samples,MCMC.noTreatment.samples,
  Utility.age.samples,Utility.state.samples,Utility.event.samples
){
  # This is a new version revised by Wei to facilate the calculation of EVPPI
  # bugs.object.baseline$sims.array[i.sample,1,j.outcome] is the baseline
  # log hazard for the j.outcome event
  # bugs.object.fixed$sims.array[i.sample,1,(k.treat-2)*17+j.outcome] is
  # the log hazard ratio of treatment k.treat relative to Coumarin (INR 2-3)
  
  # Mean costs for now
  event.costs<-t(matrix(c(2956,10844,0,1064,1751.5,2373,24234,0),ncol=n.samples,nrow=n.events))
  # SE acute costs
  event.costs[,6]<- 1186.5+Event.cost.samples[,1]*(3559.5-1186.5)
  # TIA acute costs
  event.costs[,4]<- 532+Event.cost.samples[,2]*(1596-532)
  # Bleed acute costs
  event.costs[,5]<- 875.75+Event.cost.samples[,3]*(2627.25-875.75)
  # MI acute costs
  event.costs[,1]<- 2415.24+Event.cost.samples[,4]*(7245.72-2415.24)
  # Ischemic stroke acute costs
  S.acute.cost.mean=11626
  S.acute.cost.SD=16868
  S.acute.cost.SE=S.acute.cost.SD/sqrt(162)
  event.costs[,2]<-S.acute.cost.mean+S.acute.cost.SE*Event.cost.samples[,5]
  # ICH acute costs
  I.acute.cost.mean=11453
  I.acute.cost.SD=13815
  I.acute.cost.SE=I.acute.cost.SD/sqrt(17)
  event.costs[,7]<-I.acute.cost.mean+I.acute.cost.SE*Event.cost.samples[,6]
  
  colnames(event.costs)<-event.names
  
  # Probability patient will switch given they have each type of event
  event.switch.probs<-matrix(0,nrow=n.samples,ncol=length(event.names))
  colnames(event.switch.probs)<-event.names
  # Only switch after MI if on Dabigatran
  event.switch.probs[,"MI"]<-rep(1,n.samples)
  is.disc.params<-list("alpha"=0.1,"beta"=0.9) #beta.parameters(0.3,0.7) # (0.5,0.25)
  event.switch.probs[,"Ischemic stroke"]<-Switch.probability.samples[,1]
  b.disc.params<-list("alpha"=0.3,"beta"=0.7) #beta.parameters(0.3,0.7) #(0.5,0.4)
  event.switch.probs[,"Clinically relevant bleeding"]<-Switch.probability.samples[,4]
  tia.disc.params<-list("alpha"=0.1,"beta"=0.9) #beta.parameters(0.1,0.9) #(0.05,0.1)
  event.switch.probs[,"Transient ischemic attack (TIA)"]<-Switch.probability.samples[,2]
  se.disc.params<-list("alpha"=0.1,"beta"=0.9) #beta.parameters(0.1,0.9) #(0.05,0.1)
  event.switch.probs[,"SE"]<-Switch.probability.samples[,3]
  event.switch.probs[,"ICH"]<-rep(1,n.samples)
  event.switch.probs[,"Stay"]<-rep(0,n.samples)
  event.switch.probs[,"Death (all causes)"]<-rep(0,n.samples)
  
  # Hazard ratios for effect of previous events on future events
  hr.event.history<-array(1,dim=c(n.samples,n.events,n.events-1))
  # Effect of prior bleeds on future events
  hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="MI"]<- 1 # No evidence
  hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="Ischemic stroke"]<-exp(hr.future.stroke["Bleed","Logmean"]+hr.future.stroke["Bleed","logSD"]*Effect.history.samples[,1])
  # Had to asssume effect on death same as that of stroke
  hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="Death (all causes)"]<-exp(hr.future.death["Stroke","Logmean"]+hr.future.death["Stroke","logSD"]*Effect.history.samples[,2])
  hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="Transient ischemic attach (TIA)"]<-exp(hr.future.tiase["Bleed","Logmean"]+hr.future.tiase["Bleed","logSD"]*Effect.history.samples[,3])
  hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="Clinically relevant bleeding"]<-exp(hr.future.bleed["Bleed","Logmean"]+hr.future.bleed["Bleed","logSD"]*Effect.history.samples[,4])
  hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="SE"]<-exp(hr.future.tiase["Bleed","Logmean"]+hr.future.tiase["Bleed","logSD"]*Effect.history.samples[,5])
  hr.event.history[,event.state.codes==" B ",event.names[-n.events]=="ICH"]<-exp(hr.future.ich["Bleed","Logmean"]+hr.future.ich["Bleed","logSD"]*Effect.history.samples[,6])
  
  # Effect of prior Intracranial hemorrhage of future events
  hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="MI"]<-1
  hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="Ischemic stroke"]<-exp(hr.future.stroke["ICH","Logmean"]+hr.future.stroke["ICH","logSD"]*Effect.history.samples[,7])
  # Had to assume effect on death the same as that of stroke
  hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="Death (all causes)"]<-exp(hr.future.death["Stroke","Logmean"]+hr.future.death["Stroke","logSD"]*Effect.history.samples[,8])
  hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="Transient ischemic attach (TIA)"]<-exp(hr.future.tiase["ICH","Logmean"]+hr.future.tiase["ICH","logSD"]*Effect.history.samples[,9])
  hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="Clinically relevant bleeding"]<-exp(hr.future.bleed["ICH","Logmean"]+hr.future.bleed["ICH","logSD"]*Effect.history.samples[,10])
  hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="SE"]<-exp(hr.future.tiase["ICH","Logmean"]+hr.future.tiase["ICH","logSD"]*Effect.history.samples[,11])
  hr.event.history[,event.state.codes==" I ",event.names[-n.events]=="ICH"]<-exp(hr.future.ich["ICH","Logmean"]+hr.future.ich["ICH","logSD"]*Effect.history.samples[,12])
  # Effect of prior MI on future events
  hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="MI"]<-1
  hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="Ischemic stroke"]<-exp(hr.future.stroke["MI","Logmean"]+hr.future.stroke["MI","logSD"]*Effect.history.samples[,13])
  hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="Death (all causes)"]<-exp(hr.future.death["MI","Logmean"]+hr.future.death["MI","logSD"]*Effect.history.samples[,14])
  hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="Transient ischemic attach (TIA)"]<-exp(hr.future.tiase["MI","Logmean"]+hr.future.tiase["MI","logSD"]*Effect.history.samples[,15])
  hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="Clinically relevant bleeding"]<-exp(hr.future.bleed["MI","Logmean"]+hr.future.bleed["MI","logSD"]*Effect.history.samples[,16])
  hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="SE"]<-exp(hr.future.tiase["MI","Logmean"]+hr.future.tiase["MI","logSD"]*Effect.history.samples[,17])
  hr.event.history[,event.state.codes==" M ",event.names[-n.events]=="ICH"]<-exp(hr.future.ich["MI","Logmean"]+hr.future.ich["MI","logSD"]*Effect.history.samples[,18])
  # Effect of prior strokes on future events
  hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="MI"]<-1
  hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="Ischemic stroke"]<-exp(hr.future.stroke["Stroke","Logmean"]+hr.future.stroke["Stroke","logSD"]*Effect.history.samples[,19])
  hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="Death (all causes)"]<-exp(hr.future.death["Stroke","Logmean"]+hr.future.death["Stroke","logSD"]*Effect.history.samples[,20])
  hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="Transient ischemic attach (TIA)"]<-exp(hr.future.tiase["Stroke","Logmean"]+hr.future.tiase["Stroke","logSD"]*Effect.history.samples[,21])
  hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="Clinically relevant bleeding"]<-exp(hr.future.bleed["Stroke","Logmean"]+hr.future.bleed["Stroke","logSD"]*Effect.history.samples[,22])
  hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="SE"]<-exp(hr.future.tiase["Stroke","Logmean"]+hr.future.tiase["Stroke","logSD"]*Effect.history.samples[,23])
  hr.event.history[,event.state.codes==" S ",event.names[-n.events]=="ICH"]<-exp(hr.future.ich["Stroke","Logmean"]+hr.future.ich["Stroke","logSD"]*Effect.history.samples[,24])	
  
  
  # Proportional utility decrements from Kind et al. 1999
  # Beta distributions for each age range estimated by Pete Bryden
  # Use weighted average of males (60%) and females (40%) to represent AF pop
  # Ratio of each category to 70 year olds is the proportional decrement/increment
  
  # Matrix containing absolute utilities
  kind.age.utility<-matrix(NA,nrow=n.samples,ncol=8)
  colnames(kind.age.utility)<-c(35,45,55,65,75,85,95,105)
  
  kind.age.utility[,"35"]<-0.6*rbeta(n.samples,656.7,64.95)+0.4*rbeta(n.samples,1006.6,99.5)
  kind.age.utility[,"45"]<-0.6*rbeta(n.samples,341.41,65.03)+0.4*rbeta(n.samples,544.1,96.02)
  kind.age.utility[,"55"]<-0.6*rbeta(n.samples,330.43,93.2)+0.4*rbeta(n.samples,526.59,123.52)
  kind.age.utility[,"65"]<-0.6*Utility.age.samples[,1]+0.4*Utility.age.samples[,2]
  kind.age.utility[,"75"]<-0.6*Utility.age.samples[,3]+0.4*Utility.age.samples[,4]
  
  # Assume absolute utilities for ages 90, 100 and 110 are the same as 80 (not true)
  kind.age.utility[,"85"]<-kind.age.utility[,"95"]<-kind.age.utility[,"105"]<-kind.age.utility[,"75"]
  
  age.utility.factor<-kind.age.utility/kind.age.utility[,4]
  
  
  # Treatment.costs costs (all but Warfarin costs are fixed)
  treatment.costs<-matrix(0,nrow=n.samples,ncol=n.treatments)
  # Uniform distribution on Warfarin costs for now
  treatment.costs[,1]<-52.57+Treatment.cost.samples[,1]*(157.70-52.57)
  treatment.costs[,2]<-200.42
  treatment.costs[,3]<-200.44
  treatment.costs[,4]<-200.44
  treatment.costs[,5]<-191.63
  colnames(treatment.costs)<-treatment.names
  
  # Health state costs (divided by four to go from annual to quarterly costs)
  # Stroke
  S.cost.mean=3613
  S.cost.SD=4235
  S.cost.SE=S.cost.SD/sqrt(136)
  S.cost<-S.cost.mean+S.cost.SE*Health.cost.samples[,1]
  # ICH (Assume it is similar to stroke)
  I.cost.mean=3613
  I.cost.SD=4235
  I.cost.SE=I.cost.SD/sqrt(136)
  I.cost<- I.cost.mean+I.cost.SE*Health.cost.samples[,2]
  # MI adds only an instant cost, so this post-state has 0 management cost
  M.cost<-rep(0,n.samples)
  
  # Health state utilities (from sources identified in Bayer Table 49)
  # These are combined proportionally
  # All utilities are later divided by 4 to make them 3-monthly
  AF.utility<-0.779+0.0045*Utility.state.samples[,1]
  S.utility<-0.69+0.025*Utility.state.samples[,2]/AF.utility
  M.utility<-0.718+0.0163*Utility.state.samples[,3]/AF.utility
  I.utility<-Utility.state.samples[,4]/AF.utility
  # Need evidence for post bleed utility (for now assume same as Stroke)
  B.utility<-S.utility
  
  # Cap the utilities at 1
  AF.utility[AF.utility>1]<-1
  S.utility[S.utility>1]<-1
  M.utility[M.utility>1]<-1
  I.utility[I.utility>1]<-1
  B.utility[B.utility>1]<-1
  
  
  state.utilities<-state.costs<-matrix(0,nrow=n.samples,ncol=length(state.names))
  colnames(state.utilities)<-colnames(state.costs)<-state.names
  for(i.treatment in 1:n.treatments){
    for(i.health.state in 1:n.health.states){
      state.costs[,(i.treatment-1)*n.health.states+i.health.state]<-treatment.costs[,i.treatment]
    }
  }
  
  # Use a for loop to add effects of previous events (should use something faster)
  for(i.sample in 1:n.samples)
  {
    # History of 1 event
    # Add MI costs 
    state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)==-1]<-state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)==-1]+M.cost[i.sample]
    # Add ICH costs
    state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)==-1]<-state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)==-1]+I.cost[i.sample]
    # Add stroke costs 
    state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)!=-1]<-state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)!=-1]+S.cost[i.sample]
    # History of 2 events
    # Add bleed costs (excluding those that are both ICH)
    # Add max of stroke or ICH cost to states with history of both events but not Bleed
    state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)!=-1]<-state.costs[i.sample,gregexpr(" M ",state.names)==-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)!=-1]+max(I.cost[i.sample],S.cost[i.sample])
    # Add max of MI or ICH cost to states with history of both events but not Stroke
    state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)==-1]<-state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)==-1]+max(M.cost[i.sample],I.cost[i.sample])
    # Add max of MI or stroke cost to states with history of both events but not ICH
    state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)!=-1]<-state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)==-1 & gregexpr(" S ",state.names)!=-1]+max(M.cost[i.sample],S.cost[i.sample])
    # History of 3 events
    # Add max of MI, stroke, or ICH costs
    state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)!=-1]<-state.costs[i.sample,gregexpr(" M ",state.names)!=-1 & gregexpr(" I ",state.names)!=-1 & gregexpr(" S ",state.names)!=-1]+max(M.cost[i.sample],I.cost[i.sample],S.cost[i.sample])
    
    # Utilities are multiplicative		
    state.utilities[i.sample,]<-AF.utility[i.sample]
    state.utilities[i.sample,gregexpr(" B ",state.names)!=-1]<-state.utilities[i.sample,gregexpr(" B ",state.names)!=-1]*B.utility[i.sample]
    state.utilities[i.sample,gregexpr(" I ",state.names)!=-1]<-state.utilities[i.sample,gregexpr(" I ",state.names)!=-1]*I.utility[i.sample]
    state.utilities[i.sample,gregexpr(" M ",state.names)!=-1]<-state.utilities[i.sample,gregexpr(" M ",state.names)!=-1]*M.utility[i.sample]
    state.utilities[i.sample,gregexpr(" S ",state.names)!=-1]<-state.utilities[i.sample,gregexpr(" S ",state.names)!=-1]*S.utility[i.sample]
  }
  
  
  # Costs and utilities for the death state should be zero
  state.utilities[,"Dead"]<-0
  state.costs[,"Dead"]<-0
  
  # Disutilities
  # Old dummy disutilities
  #event.utilities=t(matrix(c(-0.1,-0.1,0,-0.01,-0.03,-0.01,-0.1,0),ncol=n.samples,nrow=n.events))
  event.utilities<-matrix(0,nrow=n.samples,ncol=n.events)
  colnames(event.utilities)<-event.names
  event.utilities[,"MI"]<-0.683+0.0156*Utility.event.samples[,4]-AF.utility
  event.utilities[,"Ischemic stroke"]<-1.5*(-0.59)+Utility.event.samples[,1]*2*0.59
  event.utilities[,"Transient ischemic attack (TIA)"]<-1.5*(-0.131)+Utility.event.samples[,2]*2*0.131
  event.utilities[,"Clinically relevant bleeding"]<--0.03+0.001531*Utility.event.samples[,5]
  event.utilities[,"SE"]<-1.5*(-0.131)+Utility.event.samples[,3]*2*0.131
  event.utilities[,"ICH"]<-0.6+0.064*Utility.event.samples[,6]-AF.utility
  
  # Finally, scale the state utilities to be 3-monthly
  state.utilities<-state.utilities/4
  
  
  return(list("event.costs"=event.costs,"event.utilities"=event.utilities,
              "event.switch.probs"=event.switch.probs,"hr.event.history"=hr.event.history,
              "bugs.loghr"=MCMC.loghr.samples,"bugs.baseline"=MCMC.baseline.samples,
              "hr.no.treatment"=MCMC.noTreatment.samples,
              "age.utility.factor"=age.utility.factor,
              "state.costs"=state.costs,"state.utilities"=state.utilities))
}