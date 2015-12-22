source('a2p.R')
if(.Platform$pkgType == "mac.binary")	dyn.load("cpp/Rfunc_mac.so")
if(.Platform$pkgType == "source")		dyn.load("cpp/Rfunc.so")

tre = rUMT(5)
t = ape2peth(tre)

require('ks')

run = function(tree, dt=0.01, sigma=1, a=0)
{
	num_tips = length(tree$data_order)
	splitting_nodes = tree$splitting_nodes - 1 		# c counts from 0!!
	times = tree$times
	tval = rep(0, num_tips)	

	result = .C ("pathsim", ntip=as.integer(num_tips), dt=as.double(dt), 
				rate = as.double(sigma), a=as.double(a), r_intervals=as.double(times), 
				splitters=as.integer(splitting_nodes), tval = as.double(tval)
				)

	result$tval = result$tval[tree$data_order]

	return(result)
}

#-------------------- Old peth code: -----------------------------------#

# tree must be a pethtree class object, or will throw bad_alloc!
genTree	= function(tree, a=0, sigma=1, sigma2 = 1, dt=1, nTraits=1, kernel="CE", lim=0) 
{
	x = data.frame( run(tree=tree, dt=0.01*dt, sigma=sigma, a=a)$tval )
	names(x) = 'traits'
	return(x)
}

# generate vcv-matrix of simulated trees.
asVCV	= function(tree, sigma=1, a=0, reps=1e4, dt=0.01) 
{
	x	= t( replicate(reps, run(tree, sigma=sigma, dt=dt, a=a)$tval))
	return( cov(x) )
}

# obtain difference between data and a single simulated tree (for ABC not user)
get_dif	= function(tree, data, a, sigma, sigma2="NA", force=FALSE, dt=1, kernel="CE", lim=0, nTraits=1, sstat="std") 
{
	if(sigma2=="NA") 	sigma2 = sigma 		# For if we aren't using a 2-rate model. Shouldn't matter really.

	ntips	= length(data[,1])
	nTraits	= length(data[1,])
	new		= genTree(tree=tree, a=a, sigma=sigma, sigma2=sigma2, dt=dt, nTraits=nTraits, kernel=kernel, lim=lim)		# simulate dataset

	# Get Blomberg's K for first trait
    if(sstat=="K")
    {
    dataK 	= tryCatch(Kcalc(data[,1], tree, F), error=function(err){return(1)})  #dataK + Kcalc(data[,i], tree, F)
    newK 	= tryCatch(Kcalc(new[,1], tree, F), error=function(err){return(1)})  #newK  + Kcalc(new[,i], tree, F)
    }
	
    difs					= as.matrix(dist(data))			# euclidian distance
	difs[which(difs==0)]	= NA							# ignore matrix diagonal
	Dgap					= apply(difs, 1, min, na.rm=T)	

	difs					= as.matrix(dist(new))			# euclidian distance
	difs[which(difs==0)]	= NA							# ignore matrix diagonal
	Ngap					= apply(difs, 1, min, na.rm=T)

	# Use summary statistics: mean and sd of gaps between neighbours
    if(sstat=="std") 	return( abs(mean(Dgap) - mean(Ngap)) * abs(sd(Dgap) - sd(Ngap)))
    if(sstat=="K")		return( abs(mean(Dgap) - mean(Ngap)) * abs(sd(Dgap) - sd(Ngap)) * abs(dataK - newK))
}

LRT	= function(tree, data, min=0, max_sigma=10, max_a=5, reps=1e3, e=NA, a=NA, sigma=NA, dt=1, lim=5, file="sample.out", posteriorSize=500, sstat="std", kernel='CE')
{
	if(file.exists(file))
	{
		print(paste(file, 'already exists!'))
		return(c())
	}

    lim = max(abs(data))
	# Simulate and write to file as we go. Single threaded.
	for(i in 1:reps)
   	{
 		sig 	= runif(1, min, max_sigma)
   		atry 	= runif(1, min, max_a)
  		dist 	= get_dif(tree, data, atry, sig, sigma2="NA", dt=dt, kernel=kernel, lim=lim, sstat=sstat)

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



