# Short script to convert the bugs outputs to simple csv files for
# use in the NOACs model

# Fixed event codes
event.codes<-c(5,1,6,7,17,15,16)
n.events<-length(event.names)

# bugs.object.baseline$sims.array[i.sample,1,j.outcome] is the baseline
# log hazard for the j.outcome event
bugs.baseline<-bugs.object.baseline$sims.array[,1,event.codes]

# bugs.object.fixed$sims.array[i.sample,1,(k.treat-2)*17+j.outcome] is
# the log hazard ratio of treatment k.treat relative to Coumarin (INR 2-3)
bugs.loghr<-matrix(nrow=dim(bugs.object.fixed$sims.array)[1],ncol=n.treatments*length(event.codes))
# Start at 2 as we don't want/have hazard ratios for reference (Coumarin (INR 2-3))
# Don't go to final treatment "No treatment" 
for(i.treatment in 2:(n.treatments-1))
{
	k.treat<-which(t.names==treatment.names[i.treatment])
	bugs.loghr[,(i.treatment-1)*length(event.codes)+c(1:length(event.codes))]<-bugs.object.fixed$sims.array[,1,(k.treat-2)*17+event.codes]
}

# Compare to make sure it is working correctly
bugs.loghr[1,(i.treatment-1)*length(event.codes)+c(1:length(event.codes))]
bugs.object.fixed$sims.array[1,1,(k.treat-2)*17+event.codes]

# No treatment meta-analysis results
load(file=paste(baseline.directory,"/data/hr.no.treatment.rda",sep=""))
# No evidence comparing MI rates in Warfarin and placebo
hr.no.treatment[,"MI"]<-rep(1,dim(hr.no.treatment)[1])


# Save as csv files for use in the NOACs net benefit function
write.csv(bugs.loghr,file=paste(baseline.directory,"/data/bugs.loghr.csv",sep=""))
write.csv(bugs.baseline,file=paste(baseline.directory,"/data/bugs.baseline.csv",sep=""))
write.csv(hr.no.treatment,file=paste(baseline.directory,"/data/hr.no.treatment.csv",sep=""))







