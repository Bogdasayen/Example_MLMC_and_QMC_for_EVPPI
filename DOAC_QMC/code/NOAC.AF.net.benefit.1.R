# Function to generate the costs, qalys, and net benefits of the NOACs model


noac.net.benefit<-function(n.samples,n.cycles,initial.age,lambdas,age.independent.samples)
{
	# Vectors of total costs and total qalys output from the model
	total.qalys<-total.costs<-matrix(0,nrow=n.samples,ncol=(n.treatments-1))
	colnames(total.costs)<-colnames(total.qalys)<-treatment.names[1:(n.treatments-1)]


	# Set up cohort vector
	cohort.vector<-array(0,dim=c(n.samples,(n.treatments-1),n.states))
	# Initialise the cohort vectors for each treatment
	for(i.treatment in 1:(n.treatments-1))cohort.vector[,i.treatment,(n.health.states)*(i.treatment-1)+1]<-1

	# Run model for n.cycles (eg. n.cycles=120 means years=30)
	age<-initial.age
	old.age<-age-1

	for(i.cycle in 1:n.cycles)
	{
		# Only resample every four cycles
		if((floor(age)-old.age)==1)
		{
			sampled.probabilities<-generate.probabilities(n.samples,ages=floor(age),event.costs=age.independent.samples$event.costs,
							event.utilities=age.independent.samples$event.utilities,event.switch.probs=age.independent.samples$event.switch.probs,
							hr.event.history=age.independent.samples$hr.event.history,bugs.loghr=age.independent.samples$bugs.loghr,
							bugs.baseline=age.independent.samples$bugs.baseline,hr.no.treatment=age.independent.samples$hr.no.treatment)
			old.age<-age
		}

		# Half cycle correction (subtract half of initial QALYs and costs)

		for(i.sample in 1:n.samples)
		{
		  if(i.cycle==1)
		  {
		    total.costs[i.sample,]<--0.5*((1/1.035)^(1/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.costs+sampled.probabilities$transient.costs[[floor(age)]])[i.sample,]
		    total.qalys[i.sample,]<--age.independent.samples$age.utility.factor[1,toString(round(age/10)*10-5)]*0.5*((1/1.035)^(1/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.utilities+sampled.probabilities$transient.utilities[[floor(age)]])[i.sample,]
		  }
			# Add costs and QALYs for current cohort vector
			total.costs[i.sample,]<-total.costs[i.sample,]+((1/1.035)^(i.cycle/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.costs+sampled.probabilities$transient.costs[[floor(age)]])[i.sample,]
			total.qalys[i.sample,]<-total.qalys[i.sample,]+age.independent.samples$age.utility.factor[1,toString(round(age/10)*10-5)]*((1/1.035)^(i.cycle/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.utilities+sampled.probabilities$transient.utilities[[floor(age)]])[i.sample,]
			# Matrix multiplication to update cohort vector
			cohort.vector[i.sample,,]<-cohort.vector[i.sample,,]%*%sampled.probabilities$transition.matrix[[floor(age)]][i.sample,,]		
		} # End sample loop

		age<-age+0.25

	} # End cycle loop

	# Half cycle correction (add half of final QALYs and costs)
	age<-age+0.25
	# Only resample every four cycles
	if((floor(age)-old.age)==1)
	{
		sampled.probabilities<-generate.probabilities(n.samples,ages=floor(age),event.costs=age.independent.samples$event.costs,
						event.utilities=age.independent.samples$event.utilities,event.switch.probs=age.independent.samples$event.switch.probs,
						hr.event.history=age.independent.samples$hr.event.history,bugs.loghr=age.independent.samples$bugs.loghr,
						bugs.baseline=age.independent.samples$bugs.baseline,hr.no.treatment=age.independent.samples$hr.no.treatment)
	old.age<-age
	}
	# Sum the costs and QALYs over all n.cycles, discounting at 3.5% per year.
	for(i.sample in 1:n.samples){
	total.costs[i.sample,]<-total.costs[i.sample,]+0.5*((1/1.035)^((n.cycles+1)/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.costs+sampled.probabilities$transient.costs[[floor(age)]])[i.sample,]
	total.qalys[i.sample,]<-total.qalys[i.sample,]+age.independent.samples$age.utility.factor[1,toString(round(age/10)*10-5)]*0.5*((1/1.035)^((n.cycles+1)/4))*cohort.vector[i.sample,,]%*%(age.independent.samples$state.utilities+sampled.probabilities$transient.utilities[[floor(age)]])[i.sample,]
  }
	# Generate the net benefit and incremental net benefit for output
	INB<-array(NA, dim=c(n.samples,n.treatments-1,length(lambdas)))
	NB<-array(NA, dim=c(n.samples,n.treatments-1,length(lambdas)))
	# Calculate the incremental net benefit for each treatment
	for(i.treatment in 2:(n.treatments-1))
	{
	for(i.lambda in 1:length(lambdas))
	{
		INB[,i.treatment-1,i.lambda]<-lambdas[i.lambda]*(total.qalys[,i.treatment]-total.qalys[,1])-(total.costs[,i.treatment]-total.costs[,1])
	}
	}
	for(i.treatment in 1:(n.treatments-1))for(i.lambda in 1:length(lambdas)){NB[,i.treatment,i.lambda]<-lambdas[i.lambda]*total.qalys[,i.treatment]-total.costs[,i.treatment]}

	return(list("total.costs"=total.costs,"total.qalys"=total.qalys,"NB"=NB,"INB"=INB))

} # End net benefit function

