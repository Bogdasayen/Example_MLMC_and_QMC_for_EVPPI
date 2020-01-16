# Function to convert life-tables into hazard ratios for age

generate.hr.death.age<-function(base.age,lifetables)
{
		# Use life tables for hazard of death
		# Life tables give probability for male/female death within one year
		# Assume a 60/40 split male/female
		# Calculate a hazard from this using Prob=1-exp(-hazard)
		# So annual hazard=-log(1-Prob)
		hazard.death<-rep(NA,max(lifetables[,"Age"])-base.age+1)
		for(i in 1:(max(lifetables[,"Age"])-base.age+1)){
			hazard.death[i]<--log(1-((0.6*lifetables[lifetables[,"Age"]==(base.age+i-1),"male.qx"]+0.4*lifetables[lifetables[,"Age"]==(base.age+i-1),"female.qx"])))
		}
		hr.death<-hazard.death/hazard.death[1]
		names(hr.death)<-c(base.age:(max(lifetables[,"Age"])))
		return(hr.death)
}