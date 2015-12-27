# functions to simulate density-dependent phylogenetic data
# magnusclarke@gmail.com
# modified 2015

require('ks')
source('tree.R')
if(.Platform$pkgType == "mac.binary")	dyn.load("../cpp/Rfunc_mac.so")
if(.Platform$pkgType == "source")		dyn.load("../cpp/Rfunc.so")

#--------------------------------------------------------------------------------------#
#--- Get a dataset simulated under BM + competition, for a given tree. ----------------#
#--- Returns trait values for tips in order corresponding to ape tree tips. -----------#
#--------------------------------------------------------------------------------------#
sim = function(tree, dt=0.01, sigma=1, a=0)
{
	num_tips = length(tree$data_order)
	splitting_nodes = tree$splitting_nodes - 1 		# R counts from 1; c counts from 0.
	times = tree$times
	tval = rep(0, num_tips)	

	result = .C ("pathsim", ntip=as.integer(num_tips), dt=as.double(dt), 
				rate = as.double(sigma^2), a=as.double(a), r_intervals=as.double(times), 
				splitters=as.integer(splitting_nodes), tval = as.double(tval)
				)

	result$tval = result$tval[tree$data_order]

	return(result)
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
get_dif	= function(tree, data, a, sigma, dt=1, nTraits=1, use_K=FALSE) 
{
	ntips	= length(data[,1])
	nTraits	= length(data[1,])

	new		= genTree(tree=tree, a=a, sigma=sigma, dt=dt, nTraits=nTraits)
	
    difs					= as.matrix(dist(data))			# Euclidian distance
	difs[which(difs==0)]	= NA							# Ignore matrix diagonal
	Dgap					= apply(difs, 1, min, na.rm=T)	

	difs					= as.matrix(dist(new))			
	difs[which(difs==0)]	= NA						
	Ngap					= apply(difs, 1, min, na.rm=T)

	# Summary statistics: mean and sd of gaps between neighbours. Plus Blomberg's K if true.
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
#---------- Likelihood ratio: BM versus competition -----------------------------------#
#--------------------------------------------------------------------------------------#
lrt	= function(tree, data, min=0, max_sigma=10, max_a=5, reps=1e3, dt=0.01, file="sample.out", posteriorSize=500, use_K=FALSE)
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
  		dist 	= get_dif(tree, data, atry, sig, dt=dt, use_K=use_K)
	   	write(c(sig, atry, dist), file=file, append=TRUE, sep=",")
   	}

   	# Read simulations back into R
   	x 	= read.csv(file)
   	sig = x[,1]
   	atry= x[,2]
   	dist= x[,3]

    # Get simulations from 500th smallest to smallest
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
#---------------------- Legacy functions ----------------------------------------------#
#--------------------------------------------------------------------------------------#

# tree must be a pethtree class object, or will throw bad_alloc!
genTree	= function(tree, a=0, sigma=1, sigma2 = 1, dt=1, nTraits=1, kernel="CE", lim=0) 
{
	x = data.frame( sim(tree=tree, dt=0.01*dt, sigma=sigma, a=a)$tval )
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
