# functions to simulate density-dependent phylogenetic data
# magnusclarke@gmail.com
# modified 2016

require('ks')
require('picante')
source('tree.R')
if(.Platform$pkgType == "mac.binary")	dyn.load("../cpp/Rfunc_mac.so")
if(.Platform$pkgType == "source")		dyn.load("../cpp/Rfunc.so")

#--------------------------------------------------------------------------------------#
#--- Get a dataset simulated under BM + competition, for a given tree. ----------------#
#--- Returns trait values for tips in order corresponding to ape tree tips. -----------#
#--------------------------------------------------------------------------------------#
sim = function(tree, dt=0.001, sigma=1, a=0, ntraits=1, symp=NA, allo=NA)
{
	if(class(tree)=="phylo")
	{
		tree = ape2peth(tree)
	} else if(class(tree)!="pethtree") {
		stop("Tree incorrectly formatted.")
	}

	num_tips = length(tree$data_order)
	splitting_nodes = tree$splitting_nodes - 1 	# R counts from 1; c counts from 0.
	times = tree$times
	tval = rep(0, num_tips*ntraits)	

	if(all(is.na(symp)))
	{
		symp = replicate(num_tips^2, 0)
		allo = replicate(num_tips^2, 9e9)
	} else if(class(symp)=="matrix") {
		symp = as.vector(symp)
		allo = as.vector(allo)
	}

	result = .C ("pathsim", ntip=as.integer(num_tips), dt=as.double(dt), 
				rate = as.double(sigma^2), a=as.double(a), r_intervals=as.double(times), 
				splitters=as.integer(splitting_nodes), tval = as.double(tval),
				ntraits=as.integer(ntraits), symp=as.double(symp), allo=as.double(allo)
				)
	
	tval = result$tval

	ape_tval = matrix(ncol=ntraits, nrow=num_tips)
	for (i in 1:ntraits) 
	{
		ape_tval[,i] = tval[seq(i, length(tval), by=ntraits)]
		ape_tval[,i] = ape_tval[,i][tree$data_order]
	}

	result$tval = ape_tval

	result$symp = NULL
	result$allo = NULL

	return(result)
}

#--------------------------------------------------------------------------------------#
#--- Generate a matrix of the times at which lineages come into sympatry. -------------#
#--------------------------------------------------------------------------------------#
symp_matrix = function(tree, delay=0)
{
	if(class(tree)=="phylo")
	{
		t = ape2peth(tree)
	} else if(class(tree)!="pethtree") {
		stop("Tree incorrectly formatted.")
	} else {
		t = tree
	}

	ntip = length(t$data_order)
	s = matrix(nrow=ntip, ncol=ntip, 0)
	if(class(tree)=="phylo")
	{
		rownames(s) = colnames(s) = tree$tip.label
	}

	# find ages (time from root) of tip lineages
	age = 1:ntip
	age[1]=age[2]=0
	for(i in 1:ntip)
	{
		age[i+2] = sum(t$times[1:i])
	}

	# delay measured in terms of mean time between speciation events.
	time_between_splits = mean(abs(diff(t$times)))
	delay = delay * time_between_splits

	# apply starttimes and delay to symp-matrix s
	for (i in 1:ntip) 
	{
		for (j in 1:ntip)
		{
			s[i,j] = max(c(age[i], age[j])) + delay
		}
	}
	return(s)
}

#--------------------------------------------------------------------------------------#
#---------- Generate vcv-matrix of simulated trees. -----------------------------------#
#--------------------------------------------------------------------------------------#
as_vcv	= function(tree, sigma=1, a=0, reps=1e4, dt=0.01) 
{
	x	= t( replicate(reps, sim(tree, sigma=sigma, dt=dt, a=a)$tval))
	return( cov(x) )
}

#--------------------------------------------------------------------------------------#
#---- Get difference between given data and a single simulated dataset. ---------------#
#--------------------------------------------------------------------------------------#
get_dif	= function(tree, data, a, sigma, dt=1, nTraits=1, use_K=FALSE, symp) 
{
	ntips	= length(data[,1])
	nTraits	= length(data[1,])

	new		= sim(tree=tree, a=a, sigma=sigma, dt=dt, ntraits=nTraits, symp=symp)$tval
	
    difs					= as.matrix(dist(data))			# Euclidian distance
	difs[which(difs==0)]	= NA							# Ignore matrix diagonal
	Dgap					= apply(difs, 1, min, na.rm=T)	

	difs					= as.matrix(dist(new))			
	difs[which(difs==0)]	= NA						
	Ngap					= apply(difs, 1, min, na.rm=T)

	# Summary statistics: mean and sd of gaps between neighbours. Plus Blomberg's K optionally.
    if(use_K)
    {
		dataK 	= tryCatch(Kcalc(data[,1], tree, F), error=function(err){return(1)})
		newK 	= tryCatch(Kcalc(new[,1], tree, F), error=function(err){return(1)})
		return( abs(mean(Dgap) - mean(Ngap)) * abs(sd(Dgap) - sd(Ngap)) * abs(dataK - newK))
    } else {
		return( abs(mean(Dgap) - mean(Ngap)) * abs(sd(Dgap) - sd(Ngap)))
	}
}

#--------------------------------------------------------------------------------------#
#-------- Get summary statistics for a given tree and dataset -------------------------#
#--------------------------------------------------------------------------------------#
summary_stats = function(tree, data, use_K=FALSE)
{
	difs = as.matrix(dist(data))			# Euclidian distance
	difs[which(difs==0)] = NA				# Ignore matrix diagonal
	gap	= apply(difs, 1, min, na.rm=T)	

	# Summary statistics: mean and sd of gaps between neighbours. Plus Blomberg's K optionally.
	if(use_K)
	{
		k 	= tryCatch(Kcalc(data, tree, F), error=function(err){return(1)})
		stats = c(mean(gap), sd(gap), k)
	} else {
		stats = c(mean(gap), sd(gap) )
	}
	return( stats )
}

#--------------------------------------------------------------------------------------#
#-------- Generate file with many simulation parameters and summary statistics --------#
#--------------------------------------------------------------------------------------#
param_stats = function(tree, file='param_stats.out', reps=1e3, max_sigma=8, max_a=8, 
				symp=NA, allo = NA, dt=0.001, use_K=FALSE)
{
	sig = runif(reps, 0, max_sigma)
	atry = runif(reps, 0, max_a)

	for(i in 1:reps)
	{
		d = sim(tree=tree, a=atry[i], sigma=sig[i], dt=dt, ntraits=1, symp=symp, allo=allo)$tval
		stats = summary_stats(tree=tree, data=d, use_K=use_K)
	   	write(c(sig[i], atry[i], stats), file=file, append=TRUE, sep=",")
	}
}

#--------------------------------------------------------------------------------------#
#---------- Likelihood ratio: BM versus competition -----------------------------------#
#--------------------------------------------------------------------------------------#
lrt	= function(tree, data, file=NA, posteriorSize=500, use_K=FALSE, dt=0.001, max_sigma=8, max_a=8, sim_dat=NA)
{
   	# Read simulations into R
	if(is.na(sim_dat))
	{
		x 	= read.csv(file)
	} else {
		x = sim_dat
	}
   	sig = x[,1]
   	atry = x[,2]
   	stat1 = x[,3]
   	stat2 = x[,4]

	# Get summary stats for the true data
	tstat = summary_stats(tree=tree, data=data, use_K=use_K)
	tstat1 = tstat[1]
	tstat2 = tstat[2]

	if(use_K)
	{
		stat3 = x[,5]
		tstat3 = tstat[3]
		diff = abs(stat1-tstat1)  * abs(stat2-tstat2) * abs(stat3 - tstat3)
	} else {
		diff = abs(stat1-tstat1)  * abs(stat2-tstat2) 
	}

    # Get simulations from nth smallest to smallest
    H1_post		= order(diff)[1:posteriorSize]

    Usig		= sig[H1_post]
	Uatry		= atry[H1_post]
	H1_post		= matrix(ncol=2, nrow=length(Usig))
	H1_post[,1]	= Usig
	H1_post[,2]	= Uatry

	k 	= kde(H1_post, xmin=c(0, 0), xmax=c(max_sigma,max_a))
	k0	= kde(H1_post, xmin=c(0, 0), xmax=c(max_sigma,0))		# sigma to max, a to 0.

	# Use kernel smoothing to estimate likelihood maxima with and without competition.
	k_max_index 	= which(k$estimate == max(k$estimate), arr.ind = TRUE)
	H1_lik 			= k$estimate[k_max_index[1], k_max_index[2]]
	H1_est 			= c(unlist(k$eval.points)[k_max_index[1]], unlist(k$eval.points)[length(k$estimate[,1]) + k_max_index[2]])

	k0_max_index	= which(k0$estimate == max(k0$estimate), arr.ind = TRUE)
	H0_lik 			= k0$estimate[k0_max_index[1], k0_max_index[2]]
	H0_est 			= c(unlist(k0$eval.points)[k0_max_index[1]], unlist(k0$eval.points)[length(k0$estimate[,1]) + k0_max_index[2]])

	LRT				= -2 * log( H0_lik / H1_lik )

	return( data.frame(H0_est, H0_lik, H1_est, H1_lik, LRT) )
}


lrt_unnested   = function(tree, data, file=NA, posteriorSize=500, use_K=FALSE, reps=1e4,
						dt=0.001, max_sigma=8, max_a=8, min_param=0, max_param=8, prior=NA, nullh='EB', post2=100)
{
   	# Read simulations into R
	if(is.na(prior))
	{
		x 	= read.csv(file)
	} else {
		x = prior
	}
   	sig = x[,1]
   	param = x[,2]
   	stat1 = x[,3]
   	stat2 = x[,4]
   	model = x[,5]

	# Get summary stats for the true data
	tstat = summary_stats(tree=tree, data=data, use_K=use_K)
	tstat1 = tstat[1]
	tstat2 = tstat[2]

	if(use_K)
	{
		stat3 = x[,6]
		tstat3 = tstat[3]
		diff = abs(stat1-tstat1)  * abs(stat2-tstat2) * abs(stat3 - tstat3)
	} else {
		diff = abs(stat1-tstat1)  * abs(stat2-tstat2) 
	}

    #-------- Now, generate separate posteriors for H0 and H1, and estimate params
    sig0 = sig[which(model==nullh)]
   	param0 = param[which(model==nullh)]
   	stat10 = stat1[which(model==nullh)]
   	stat20 = stat2[which(model==nullh)]
   	model0 = model[which(model==nullh)]
   	diff0 = diff[which(model==nullh)]
    H0_post = order(diff0)[1:posteriorSize]

    Usig		= sig[H0_post]
	Uparam		= param[H0_post]
	H0_post		= matrix(ncol=2, nrow=length(Usig))
	H0_post[,1]	= Usig
	H0_post[,2]	= Uparam

	k0 	= kde(H0_post, xmin=c(0, min_param), xmax=c(max_sigma,max_param))

	# Use kernel smoothing to estimate likelihood maxima.
	k0_max_index 	= which(k0$estimate == max(k0$estimate), arr.ind = TRUE)
	H0_est 			= c(unlist(k0$eval.points)[k0_max_index[1]], unlist(k0$eval.points)[length(k0$estimate[,1]) + k0_max_index[2]])

	#------------------------------

    sig1 = sig[which(model=='comp')]
   	param1 = param[which(model=='comp')]
   	stat11 = stat1[which(model=='comp')]
   	stat21 = stat2[which(model=='comp')]
   	model1 = model[which(model=='comp')]
   	diff1 = diff[which(model=='comp')]
    H1_post = order(diff1)[1:posteriorSize]

    Usig1		= sig[H1_post]
	Uparam1		= param[H1_post]
	H1_post		= matrix(ncol=2, nrow=length(Usig1))
	H1_post[,1]	= Usig1
	H1_post[,2]	= Uparam1

	k1 	= kde(H1_post, xmin=c(0, 0), xmax=c(max_sigma,max_a))

	# Use kernel smoothing to estimate likelihood maxima.
	k1_max_index 	= which(k1$estimate == max(k1$estimate), arr.ind = TRUE)
	H1_est 			= c(unlist(k1$eval.points)[k1_max_index[1]], unlist(k1$eval.points)[length(k1$estimate[,1]) + k1_max_index[2]])

	#----- Generate new datasets using fitted models

	if(use_K)
	{
		d1 = matrix(nrow=reps, ncol=4)
	} else {
		d1 = matrix(nrow=reps, ncol=3)
	}

	for(i in 1:reps)
	{
		d = sim(tree=tree, a=H1_est[2], sigma=H1_est[1])$tval
		s = summary_stats(tree=tree, data=d, use_K=use_K)
		d1[i,1:2] = s[1:2]
		d1[i,3] = 'comp'
		if(use_K)
		{
			d1[i,4] = s[3]
		}
	}

	if(use_K)
	{
		d0 = matrix(nrow=reps, ncol=4)
	} else {
		d0 = matrix(nrow=reps, ncol=3)
	}
	for(i in 1:reps)
	{
		if(nullh=='EB')
		{
			transformed_tree = transf.branch.lengths(tree, model='EB', parameters=list(rate=H0_est[2]))$tree
			d = rTraitCont(transformed_tree, model='BM', sigma=H0_est[1])
		} else if(nullh=='OU')
		{
			d = rTraitCont(tree, model='OU', sigma=H0_est[1], alpha=H0_est[2])
		}
		s = summary_stats(tree=tree, data=d, use_K=use_K)
		d0[i,1:2] = s[1:2]
		d0[i,3] = nullh
		if(use_K)
		{
			d0[i,4] = s[3]
		}
	}

   	newstat1 = as.numeric( c( d0[,1], d1[,1] ) )
   	newstat2 = as.numeric( c( d0[,2], d1[,2] ) )
   	newmodel = c( d0[,3], d1[,3] )

	if(use_K)
	{
		newstat3 = as.numeric( c( d0[,4], d1[,4] ) )
		newdiff = abs(newstat1-tstat1)  * abs(newstat2-tstat2) * abs(newstat3 - tstat3)
	} else {
		newdiff = abs(newstat1-tstat1)  * abs(newstat2-tstat2) 
	}

    # Get simulations from nth smallest distance to smallest
    newposterior = order(newdiff)[1:post2]

    # number of null sims in posterior
    H0_lik = length(which(newmodel[newposterior]==nullh))

    # number of alternate sims in posterior
    H1_lik = length(which(newmodel[newposterior]=='comp'))

    LRTstat = -2 * log(H0_lik / H1_lik)

	return( data.frame(H0_est, H1_est, LRTstat) )
}











#--------------------------------------------------------------------------------------#
#---------- Likelihood ratio: BM versus competition -----------------------------------#
#--------------------------------------------------------------------------------------#
oldlrt	= function(tree, data, min=0, max_sigma=10, max_a=5, reps=1e3, dt=0.001, 
	file="sample.out", posteriorSize=500, use_K=FALSE, symp=NA)
{
	if(file.exists(file))
	{
		print(paste(file, 'already exists!'))
		return(c())
	}

	# Simulate and write to file as we go. Single threaded.
	for(i in 1:reps)
   	{
 		sig 	= runif(1, min, max_sigma)
   		atry 	= runif(1, min, max_a)
  		dist 	= get_dif(tree, data, atry, sig, dt=dt, use_K=use_K, symp=symp)
	   	write(c(sig, atry, dist), file=file, append=TRUE, sep=",")
   	}

   	# Read simulations back into R
   	x 	= read.csv(file)
   	sig = x[,1]
   	atry= x[,2]
   	dist= x[,3]

    # Get simulations from nth smallest to smallest
    H1_post		= order(dist)[1:posteriorSize]

    Usig		= sig[H1_post]
	Uatry		= atry[H1_post]
	H1_post		= matrix(ncol=2, nrow=length(Usig))
	H1_post[,1]	= Usig
	H1_post[,2]	= Uatry

	error_bar_sig = sd(Usig)
	error_bar_a	 = sd(Uatry)
	error_bar = c(error_bar_sig, error_bar_a)

	k 	= kde(H1_post, xmin=c(0, 0), xmax=c(max_sigma,max_a))
	k0	= kde(H1_post, xmin=c(0, 0), xmax=c(max_sigma,0))		# sigma to max, a to 0.

	# Use kernel smoothing to estimate likelihood maxima with and without competition.
	k_max_index 	= which(k$estimate == max(k$estimate), arr.ind = TRUE)
	H1_lik 			= k$estimate[k_max_index[1], k_max_index[2]]
	H1_est 			= c(unlist(k$eval.points)[k_max_index[1]], unlist(k$eval.points)[length(k$estimate[,1]) + k_max_index[2]])

	k0_max_index	= which(k0$estimate == max(k0$estimate), arr.ind = TRUE)
	H0_lik 			= k0$estimate[k0_max_index[1], k0_max_index[2]]
	H0_est 			= c(unlist(k0$eval.points)[k0_max_index[1]], unlist(k0$eval.points)[length(k0$estimate[,1]) + k0_max_index[2]])

	LRT				= -2 * log( H0_lik / H1_lik )

	file.remove(file)

	return( data.frame(H0_est, H0_lik, H1_est, H1_lik, LRT) )
}

#--------------------------------------------------------------------------------------#
#---------------------- Legacy functions. Don't use these. ----------------------------#
#--------------------------------------------------------------------------------------#

genTree	= function(tree, a=0, sigma=1, sigma2 = 1, dt=1, nTraits=1, kernel="CE", lim=0) 
{
	x = data.frame( sim(tree=tree, dt=0.01*dt, sigma=sigma, a=a, ntraits=nTraits)$tval )
	names(x) = 'traits'
	return(x)
}

asVCV=as_vcv
rUMT=rand_umt

LRT = function(tree, data, a, sigma, dt=0.01, nTraits, kernel, lim, sstat, 
			   reps, posteriorSize, max_sigma=5, max_a=5, min=0, file='sample.out')
{
	return(lrt(tree=tree, data=data, min=min, max_sigma=max_sigma, max_a=max_a, 
			   reps=reps, dt=dt, file=file, posteriorSize=posteriorSize, use_K=FALSE))
}

#--------------------------------------------------------------------------------------#



