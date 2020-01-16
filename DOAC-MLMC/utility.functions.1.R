# NOAC AF Cost-effectiveness model
# Script with utility functions

gamma.parameters<-function(gamma.mean,gamma.sd)
{
	gamma.shape=(gamma.mean^2)/(gamma.sd^2)
	gamma.scale=(gamma.sd^2)/gamma.mean
	return(list("shape"=gamma.shape,"scale"=gamma.scale))
}

beta.parameters<-function(beta.mean,beta.sd)
{
	beta.var<-beta.sd^2
  alpha <- ((1 - beta.mean) / beta.var - 1 / beta.mean) * beta.mean ^ 2
  beta <- alpha * (1 / beta.mean - 1)
  return(list(alpha = alpha, beta = beta))
}
