# functions to simulate density-dependent trees
# magnus
# modified 6.2.14

library("ape")
library("TESS")
library("parallel")
source("simPlot.R")
if(.Platform$pkgType == "mac.binary")	dyn.load("../cpp/Rfunc_mac.so")
if(.Platform$pkgType == "source")		dyn.load("../cpp/Rfunc.so")

cores	<- detectCores()

# random ultrametric tree
randUMT	<- function(nt, lambda=1, mu=0)
{
	tree 	<- sim.globalBiDe.taxa(n=1, nTaxa=nt, lambda=lambda, mu=mu)	
	return( tree[[1]] )
}

# simulate tip trait data: 2 traits
genTree	<- function(tree, a=0, sigma=1, dt=1, nTraits=1, kernel="CE", lim=0) 
{
	dt 	<- 0.01 * dt

	# Require tree$edge to be a matrix of integers not doubles. TESS gives doubles.
	tree$edge <- apply(tree$edge, c(1,2), function(x) x <- as.integer(x))	

	tip_count      <- length( tree$tip.label )
	st             <- tree$edge[,1]
	end            <- tree$edge[,2]
	length         <- tree$edge.length
	edge_count     <- length(st)

	if(kernel=="BM")	k=0
	if(kernel=="CE")	k=1
	if(kernel=="LIM")	k=2

	result	<- .C 	("genTree", ec=edge_count, nc = tip_count,
			nt = as.integer(nTraits), kernel=as.integer(k),
			a=as.double(a), start=st, end, len=as.double(length),
			sigma=as.double(sigma), dt=as.double(dt), lim=as.double(lim),
			tip=rep(0.0, edge_count*nTraits))


	traits <- matrix(result$tip, ncol=nTraits)
	traits <- traits[ which(end <= tip_count) ,]
	
	return( data.frame(traits) )
}

# generate vcv-matrix of simulated trees. a=0 is BM
simVCV	<- function(tree, a=0, sigma=1, reps=1e3, dt=1) 
{
	x	<- t( replicate(reps, genTree(tree, a, sigma, dt)$trait1))
	return( cov(x) )
}

# obtain difference between data and a single simulated tree (for ABC not user)
get_dif	<- function(tree, data, a, sigma, force=FALSE, dt=1) 
{
	ntips	<- length(data[,1])
	nTraits	<- length(data[1,])
	new		<- genTree(tree, a, sigma, dt, nTraits=nTraits)		# simulate dataset

	if(a != 0 | force == TRUE)
	{
		difs					<- as.matrix(dist(data))	# euclidian distance
		difs[which(difs==0)]	<- 1e100					# ignore matrix diagonal
		Dgap					<- apply(difs, 1, min)	

		difs					<- as.matrix(dist(new))		# euclidian distance
		difs[which(difs==0)]	<- 1e100					# ignore matrix diagonal
		Ngap					<- apply(difs, 1, min)

		# Use summary statistics: mean and sd of gaps between neighbours
		return( abs(mean(Dgap) - mean(Ngap)) + abs(sd(Dgap) - sd(Ngap)))
	} else {
		Dtrait_var		<- var(data)	
		Ntrait_var		<- var(new)	
		return( sum( abs(Dtrait_var - Ntrait_var) ) )
	}	

}

# Fit model to tree and tip data
ABC	<- function	(tree, data, min=0, max=20, reps=1e3, e=NA, a=NA, sigma=NA, dt=1)
{
	use 	<- rep(FALSE, reps)
	if(is.na(sigma)) 	sig	<- runif(reps, min, max)
		else		sig	<- rep(sigma, reps)
	if(is.na(a)) 		atry	<- runif(reps, min, max)
		else		atry	<- rep(a, reps)

	# Setting e empirically
	treps	<- reps/1e1
	if(is.na(a) & is.na(sigma)) 
		test 	<- replicate(treps, get_dif(tree, data, runif(1, min, max), runif(1, min, max), dt=dt) )
	else if(is.double(a))
		test 	<- replicate(treps, get_dif(tree, data, a, runif(1, min, max), dt=dt) )
	else if(is.double(sigma))
		test 	<- replicate(treps, get_dif(tree, data, runif(1, min, max), sigma, dt=dt) )

	if(is.na(e))	ep	<- min(test)
	else		ep	<- e * min(test)

	use <- mcmapply(function(use, atry, sig)	
		{
			dif <- get_dif(tree, data, atry, sig, dt=dt)
			if(dif < ep)	use	<- TRUE
			else		use	<- FALSE
		}, use, atry, sig, mc.cores=cores)

	sigmaEst 	<- mean( sig[which(use == TRUE)] )
	aEst 		<- mean( atry[which(use == TRUE)] )

	hits 	<- 100*( length(which(use == TRUE)) / length(use) )	
	hits 	<- paste(hits, "%")

	if(!is.na(aEst))
	{
		if(is.na(a) & is.na(sigma))	
			return( data.frame(aEst, sigmaEst, hits, ep) )
		else if(is.double(a))	
			return( data.frame(sigmaEst, hits, ep) )
		else if(is.double(sigma))
			return( data.frame(aEst, hits, ep) )
	} else {
		x <- ABC(tree, data, min=min, max=max, reps=reps, e=e, a=a, sigma=sigma, dt=dt)
		return(x)
	}
}

dir_dif	<- function(tree, data, a, sigma, force=FALSE, dt=1, trait_sd=1) 
{
	ntips	<- length(data[,1])
	nTraits	<- length(data[1,])
	new		<- genTree(tree, a, sigma, dt, nTraits=nTraits)	

	tip_prob <- 1:ntips
	for (i in 1:ntips)
	{
		# get probability of a single tip trait value
		tip_prob[i]	<- dnorm(data[i,], mean=new[i,], sd=trait_sd)
	}

	# full probability is produict of all tip probabilities
	return( prod(tip_prob) )
}

dirABC	<- function	(tree, data, min=0, max=20, reps=1e3, e=NA, dt=1, factor=1e12)
{
	use 	<- rep(FALSE, reps)
	sig		<- runif(reps, min, max)
	atry	<- runif(reps, min, max)

	use <- mcmapply(function(use, atry, sig)	
		{
			dif <- dir_dif(tree, data, atry, sig, dt=dt, trait_sd=trait_sd)
			rn 	<- runif(1)						# random number between 0 and 1
			if((factor*dif) > rn)	use	<- TRUE	# so if dif=0.1/factor, probability of use=TRUE is 0.1
			else			use	<- FALSE
		}, use, atry, sig, mc.cores=cores)

	sigmaEst 	<- mean( sig[which(use == TRUE)] )
	aEst 		<- mean( atry[which(use == TRUE)] )
	hits 		<- 100*( length(which(use == TRUE)) / length(use) )	

	return( data.frame(aEst, sigmaEst, hits) )
}


# Get liklihood for given a, sigma. Called by AIC
model_lik   <- function(tree, data, reps=1e3, e, a=0, sigma, min=0, max=20, dt=1)
{
	x	<- 1:reps
	x <- mclapply(x, function(x)	{
			dif 		<- get_dif(tree, data, a, sigma, force=TRUE, dt=dt)
			if(dif < e)	x <- 1
			else		x <- 0	
			}, mc.cores=cores) 	

	count <- sum(as.integer(x))
	lik 	<- count / reps
	return(lik)
}

# Test AIC values for density model vs BM model
LRT	<- function(tree, data, reps=1e3, min=0, max=20, dt=1) 
{
	brown	<- ABC(tree, data, min=min, max=max, reps=reps, e=NA, a=0, dt=dt)
	full	<- ABC(tree, data, min=min, max=max, reps=reps, e=NA, dt=dt)
	Ea	<- full$aEst
	Ea_s	<- full$sigmaEst
	Es	<- brown$sigmaEst

	e	<- full$e

	lik_BM  <- model_lik(tree, data, e=e, reps=10*reps, a=0, sigma=Es, dt=dt)
	lik_den	<- model_lik(tree, data, e=e, reps=10*reps, a=Ea, sigma=Ea_s, dt=dt)

	D	<- -2 * log( lik_BM / lik_den )

	#AIC_BM		<- 2*1 - 2 * log(lik_BM)
	#AIC_density	<- 2*2 - 2 * log(lik_den)
	#relative_lik	<- exp( (AIC_BM - AIC_density) / 2)
	return( data.frame(lik_BM, lik_den, D) )
}
