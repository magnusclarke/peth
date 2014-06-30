# functions to simulate density-dependent trees
# magnus
# modified 6.2.14

library("ape")
library("TESS")
library("parallel")
library("ks")
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
	if(kernel=="NF")	k=3

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
get_dif	<- function(tree, data, a, sigma, force=FALSE, dt=1, kernel="CE", lim=0) 
{
	ntips	<- length(data[,1])
	nTraits	<- length(data[1,])
	new		<- genTree(tree, a, sigma, dt, nTraits=nTraits, kernel=kernel, lim=lim)		# simulate dataset

	if(a != 0 | force == TRUE)
	{
		difs					<- as.matrix(dist(data))	# euclidian distance
		difs[which(difs==0)]	<- NA						# ignore matrix diagonal
		Dgap					<- apply(difs, 1, min, na.rm=T)	

		difs					<- as.matrix(dist(new))		# euclidian distance
		difs[which(difs==0)]	<- NA						# ignore matrix diagonal
		Ngap					<- apply(difs, 1, min, na.rm=T)

		# Use summary statistics: mean and sd of gaps between neighbours
		return( abs(mean(Dgap) - mean(Ngap)) + abs(sd(Dgap) - sd(Ngap)))
	} else {
		Dtrait_var		<- var(data)	
		Ntrait_var		<- var(new)	
		return( sum( abs(Dtrait_var - Ntrait_var) ) )
	}	

}

# Fit model to tree and tip data
ABC	<- function	(tree, data, min=0, max=10, reps=1e3, e=NA, a=NA, sigma=NA, dt=1, kernel="CE", lim=0, plot=FALSE)
{
	use 	<- rep(FALSE, reps)
	if(is.na(sigma)) 	sig	<- runif(reps, min, max)
		else		sig	<- rep(sigma, reps)
	if(is.na(a)) 		atry	<- runif(reps, min, max)
		else		atry	<- rep(a, reps)

	# Setting e empirically
	treps	<- as.integer( sqrt(reps) )		#reps/1e1
	if(is.na(a) & is.na(sigma)) 
		test 	<- replicate(treps, get_dif(tree, data, runif(1, min, max), runif(1, min, max), dt=dt, kernel=kernel, lim=lim) )
	else if(is.double(a))
		test 	<- replicate(treps, get_dif(tree, data, a, runif(1, min, max), dt=dt, kernel=kernel, lim=lim) )
	else if(is.double(sigma))
		test 	<- replicate(treps, get_dif(tree, data, runif(1, min, max), sigma, dt=dt, kernel=kernel, lim=lim) )

	if(is.na(e))	ep	<- min(test)
	else		ep	<- e * min(test)

	use <- mcmapply(function(use, atry, sig)	
		{
			dif <- get_dif(tree, data, atry, sig, dt=dt)
			if(dif < ep)	use	<- TRUE
			else			use	<- FALSE
		}, use, atry, sig, mc.cores=cores)

	Usig	<- sig[which(use == TRUE)]
	Uatry	<- atry[which(use == TRUE)]
	dat 	<- matrix(ncol=2, nrow=length(Usig))
	dat[,1]	<- Usig
	dat[,2]	<- Uatry

	return(dat)
}


manualLRT   <- function(tree, data, min=0, max=10, reps=1e3, e=NA, a=NA, sigma=NA, dt=1, kernel1="LIM", kernel2="BM", lim=5, plot=F)
{
	sig 	<- runif(reps, min, max)
    atry    <- runif(reps, min, max)
   	H0_dist <- 1:reps 					# H0 distances to truth
   	H1_dist <- 1:reps 					# H1 distances to truth

   	# Simulate H0 (BM model)
    H0_dist <- mcmapply(function(H0_dist, atry, sig)
    {
        H0_dist <- get_dif(tree, data, atry, sig, dt=dt, kernel=kernel2, lim=lim)
    }, H0_dist, atry, sig, mc.cores=cores)

   	# Simulate H1 (competition with limits)
    H1_dist <- mcmapply(function(H1_dist, atry, sig)
    {
        H1_dist <- get_dif(tree, data, atry, sig, dt=dt, kernel=kernel1, lim=lim)
    }, H1_dist, atry, sig, mc.cores=cores)

    # Get H0 simulation 1% of way from smallest
    cutoff	<- quantile(H0_dist, 0.005)

    # How many simulations smaller than cutoff?
    H0_lik 	<- length( which(H0_dist < cutoff) )
    H1_lik 	<- length( which(H1_dist < cutoff) )
    LRT				<- -2 * log( H0_lik / H1_lik )

	return( data.frame(H0_lik, H1_lik, LRT) )
}

manualLRTold2   <- function(tree, data, min=0, max=10, reps=1e3, e=NA, a=NA, sigma=NA, dt=1, kernel1="LIM", kernel2="BM", lim=5, plot=F)
{
    use1    <- rep(FALSE, reps)
    use2    <- rep(FALSE, reps)
    sig <- runif(reps, min, max)
    atry    <- runif(reps, min, max)

    # Setting e empirically for kernel1
    treps   <- as.integer (1e3)			#( sqrt(reps) )     #reps/1e1
    if(is.na(a) & is.na(sigma))
    {
        test    <- replicate(treps, get_dif(tree, data, runif(1, min, max), runif(1, min, max), dt=dt, kernel=kernel1, lim=lim) )
    }    else if(is.double(a))
    {
        test    <- replicate(treps, get_dif(tree, data, a, runif(1, min, max), dt=dt, kernel=kernel1, lim=lim) )
    }    else if(is.double(sigma))
    { 
        test    <- replicate(treps, get_dif(tree, data, runif(1, min, max), sigma, dt=dt, kernel=kernel1, lim=lim) )
    }
    if(is.na(e))    ep1  <- min(test)    else            ep1  <- e * min(test)


    # Setting e empirically for kernel2
    treps   <- as.integer( sqrt(reps) )     #reps/1e1
    if(is.na(a) & is.na(sigma))
    {
        test    <- replicate(treps, get_dif(tree, data, runif(1, min, max), runif(1, min, max), dt=dt, kernel=kernel2, lim=lim) )
    }    else if(is.double(a))
    {
        test    <- replicate(treps, get_dif(tree, data, a, runif(1, min, max), dt=dt, kernel=kernel2, lim=lim) )
    }    else if(is.double(sigma))
    { 
        test    <- replicate(treps, get_dif(tree, data, runif(1, min, max), sigma, dt=dt, kernel=kernel2, lim=lim) )
    }
    if(is.na(e))    ep2  <- min(test)    else            ep2  <- e * min(test)

    use1 <- mcmapply(function(use1, atry, sig)
        {
            dif <- get_dif(tree, data, atry, sig, dt=dt, kernel=kernel1, lim=lim)
            if(dif < ep1)    use1    <- TRUE
            else            use1    <- FALSE
        }, use1, atry, sig, mc.cores=cores)

    use2 <- mcmapply(function(use2, atry, sig)
        {
            dif <- get_dif(tree, data, atry, sig, dt=dt, kernel=kernel2, lim=lim)
            if(dif < ep2)    use2    <- TRUE
            else            use2    <- FALSE
        }, use2, atry, sig, mc.cores=cores)

    a_accepted		<- atry[ which(use1 == T) ]
    sig_accepted1	<- sig[ which(use1 == T) ]
    sig_accepted2	<- sig[ which(use2 == T) ]

    sig_est1	<- mean(sig_accepted1)
	sig_est2	<- mean(sig_accepted2)
	a_est	<- mean(a_accepted)			# goes with H1 & sig_est1

    # Setting e empirically for LRT
    treps   <- as.integer( sqrt(reps) )
    test    <- replicate(treps, get_dif(tree, data, a_est, sig_est1, dt=dt, kernel=kernel1, lim=lim) )
    if(is.na(e))    ep  <- min(test)    else	ep  <- e * min(test)

    use1 <- mcmapply(function(use)
        {
            dif <- get_dif(tree, data, a_est, sig_est1, dt=dt, kernel=kernel1, lim=lim)
            if(dif < ep)	use1    <- TRUE
            else			use1    <- FALSE
        }, use1, mc.cores=cores)

    lik1	<- length(which(use1==TRUE))

    use2 <- mcmapply(function(use)
        {
            dif <- get_dif(tree, data, 0, sig_est2, dt=dt, kernel=kernel2, lim=lim)
            if(dif < ep)	use2    <- TRUE
            else			use2    <- FALSE
        }, use2, mc.cores=cores)

	lik2	<- length(which(use2==TRUE))

    LRT				<- -2 * log( lik2 / lik1 )

	return( data.frame(sig_est1, a_est, sig_est2, lik1, lik2, LRT) )
}



manualLRTold	<- function(tree, data, min=0, max=10, reps=1e3, e=NA, a=NA, sigma=NA, dt=1, kernel1="LIM", kernel2="BM", lim=0, plot=F)
{
	use1 	<- rep(FALSE, reps)
	use2 	<- rep(FALSE, reps)
	if(is.na(sigma)) 	sig	<- runif(reps, min, max)
		else			sig	<- rep(sigma, reps)
	if(is.na(a)) 		atry	<- runif(reps, min, max)
		else			atry	<- rep(a, reps)

	# Setting e empirically
	treps	<- as.integer( sqrt(reps) )		#reps/1e1
	if(is.na(a) & is.na(sigma)) 
		test 	<- replicate(treps, get_dif(tree, data, runif(1, min, max), runif(1, min, max), dt=dt, kernel=kernel1, lim=lim) )
	else if(is.double(a))
		test 	<- replicate(treps, get_dif(tree, data, a, runif(1, min, max), dt=dt, kernel=kernel1, lim=lim) )
	else if(is.double(sigma))
		test 	<- replicate(treps, get_dif(tree, data, runif(1, min, max), sigma, dt=dt, kernel=kernel1, lim=lim) )

	if(is.na(e))	ep	<- min(test)
	else			ep	<- e * min(test)

	use1 <- mcmapply(function(use, atry, sig)	
		{
			dif <- get_dif(tree, data, atry, sig, dt=dt, kernel=kernel1)
			if(dif < ep)	use1	<- TRUE
			else			use1	<- FALSE
		}, use1, atry, sig, mc.cores=cores)

	use2 <- mcmapply(function(use, atry, sig)	
		{
			dif <- get_dif(tree, data, atry, sig, dt=dt, kernel=kernel2)
			if(dif < ep)	use2	<- TRUE
			else			use2	<- FALSE
		}, use2, atry, sig, mc.cores=cores)

	Usig1	<- sig[which(use1 == TRUE)]
	Uatry1	<- atry[which(use1 == TRUE)]
	dat1 	<- matrix(ncol=2, nrow=length(Usig1))
	dat1[,1]	<- Usig1
	dat1[,2]	<- Uatry1

	Usig2	<- sig[which(use2 == TRUE)]
	Uatry2	<- atry[which(use2 == TRUE)]
	dat2 	<- matrix(ncol=2, nrow=length(Usig2))
	dat2[,1]	<- Usig2
	dat2[,2]	<- Uatry2

	library("ks")

	k1	<- kde(dat1, xmin=c(0, 0), xmax=c(max,max))
	k2	<- kde(dat2, xmin=c(0, 0), xmax=c(max,max))

	k_max_index 	<- which(k1$estimate == max(k1$estimate), arr.ind = TRUE)
	H1_lik 			<- k$estimate[k_max_index[1], k_max_index[2]]
	H1_est 			<- c(unlist(k$eval.points)[k_max_index[1]], unlist(k$eval.points)[length(k$estimate[,1]) + k_max_index[2]])

	k0_max_index	<- which(k2$estimate == max(k2$estimate), arr.ind = TRUE)
	H0_lik 			<- k0$estimate[k0_max_index[1], k0_max_index[2]]
	H0_est 			<- c(unlist(k0$eval.points)[k0_max_index[1]], unlist(k0$eval.points)[length(k0$estimate[,1]) + k0_max_index[2]])

	LRT				<- -2 * log( H0_lik / H1_lik )

	if(plot==TRUE)	plot(k, "persp")

	return( data.frame(H0_est, H0_lik, H1_est, H1_lik, LRT) )
}


LRT	<- function(dat, min=0, max=10, plot=FALSE) 
{
	if(length(dat)<2){
		H0_est 	<- 0
		H0_lik 	<- 0
		H1_est 	<- 0
		H1_lik 	<- 0
		LRT 	<- 0
	} else {
		k 	<- kde(dat, xmin=c(0, 0), xmax=c(max,max))
		k0	<- kde(dat, xmin=c(0, 0), xmax=c(max,0))		# sigma to max, a to 0.

		k_max_index 	<- which(k$estimate == max(k$estimate), arr.ind = TRUE)
		H1_lik 			<- k$estimate[k_max_index[1], k_max_index[2]]
		H1_est 			<- c(unlist(k$eval.points)[k_max_index[1]], unlist(k$eval.points)[length(k$estimate[,1]) + k_max_index[2]])


		k0_max_index	<- which(k0$estimate == max(k0$estimate), arr.ind = TRUE)
		H0_lik 			<- k0$estimate[k0_max_index[1], k0_max_index[2]]
		H0_est 			<- c(unlist(k0$eval.points)[k0_max_index[1]], unlist(k0$eval.points)[length(k0$estimate[,1]) + k0_max_index[2]])

		LRT				<- -2 * log( H0_lik / H1_lik )

		if(plot==TRUE)	plot(k, "persp")
	}
	return( data.frame(H0_est, H0_lik, H1_est, H1_lik, LRT) )
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

# This still needs to be switched to use ks package; or, to be integrated into ABC as an option.
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
# LRT	<- function(tree, data, reps=1e3, min=0, max=20, dt=1) 
# {
# 	brown	<- ABC(tree, data, min=min, max=max, reps=reps, e=NA, a=0, dt=dt)
# 	full	<- ABC(tree, data, min=min, max=max, reps=reps, e=NA, dt=dt)
# 	Ea		<- full$aEst
# 	Ea_s	<- full$sigmaEst
# 	Es		<- brown$sigmaEst

# 	e		<- full$e

# 	lik_BM  <- model_lik(tree, data, e=e, reps=10*reps, a=0, sigma=Es, dt=dt)
# 	lik_den	<- model_lik(tree, data, e=e, reps=10*reps, a=Ea, sigma=Ea_s, dt=dt)

# 	D		<- -2 * log( lik_BM / lik_den )

# 	#AIC_BM		<- 2*1 - 2 * log(lik_BM)
# 	#AIC_density	<- 2*2 - 2 * log(lik_den)
# 	#relative_lik	<- exp( (AIC_BM - AIC_density) / 2)
# 	return( data.frame(lik_BM, lik_den, D) )
# }
