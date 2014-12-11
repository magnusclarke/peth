# FUNCTIONs to simulate density-dependent trees
# magnus
# modified 6.2.14

library("ape")
library("TESS")
library("parallel")
library("ks")
library("picante")

if(.Platform$pkgType == "mac.binary")	dyn.load("../cpp/Rfunc_mac.so")
if(.Platform$pkgType == "source")		dyn.load("../cpp/Rfunc.so")

cores	= detectCores()

# random ultrametric tree - RENAME
rUMT	= function(nt, lambda=1, mu=0)
{
	tree 	= sim.globalBiDe.taxa(n=1, nTaxa=nt, lambda=lambda, mu=mu)	
	return( tree[[1]] )
}

# simulate tip trait data: 2 traits
genTree	= function(tree, a=0, sigma=1, sigma2 = 1, dt=1, nTraits=1, kernel="CE", lim=0) 
{
	dt 	= 0.01 * dt

	if(is.na(sigma2)) 	sigma2=sigma

	# Require tree$edge to be a matrix of integers not doubles. TESS gives doubles.
	tree$edge = apply(tree$edge, c(1,2), function(x) x = as.integer(x))	
	rcut = length(tree$edge) / 2		# Change BM rate at midpoint

	tip_count      = length( tree$tip.label )
	st             = tree$edge[,1]
	end            = tree$edge[,2]
	length         = tree$edge.length
	edge_count     = length(st)

	if(kernel=="BM")	k=0
	if(kernel=="CE")	k=1
	if(kernel=="LIM")	k=2
	if(kernel=="NF")	k=3
	if(kernel=="RC")	k=4

	result	= .C 	("genTree", ec=edge_count, nc = tip_count,
			nt = as.integer(nTraits), kernel=as.integer(k), ratecut=as.integer(rcut),
			a=as.double(a), start=st, end, len=as.double(length),
			sigma=as.double(sigma), sigma2=as.double(sigma2), dt=as.double(dt), lim=as.double(lim),
			tip=rep(0.0, edge_count*nTraits))


	traits = matrix(result$tip, ncol=nTraits)
	traits = traits[ which(end <= tip_count) ,]
	
	return( data.frame(traits) )
}

# generate vcv-matrix of simulated trees. a=0 is BM
asVCV	= function(tree, a=0, sigma=1, reps=1e3, dt=1) 
{
	x	= t( replicate(reps, genTree(tree, a, sigma, dt=dt)$traits))
	return( cov(x) )
}

# obtain difference between data and a single simulated tree (for ABC not user)
get_dif	= function(tree, data, a, sigma, sigma2="NA", force=FALSE, dt=1, kernel="CE", lim=0, nTraits=1, sstat="std") 
{
	if(sigma2=="NA") 	sigma2 = sigma 		# For if we aren't using a 2-rate model. Shouldn't matter really.

	ntips	= length(data[,1])
	nTraits	= length(data[1,])
	new		= genTree(tree=tree, a=a, sigma=sigma, sigma2=sigma2, dt=dt, nTraits=nTraits, kernel=kernel, lim=lim)		# simulate dataset

	# Get Blomberg's K averaged over traits
	if(sstat=="K")
	{
		dataK 	= 0
		newK 	= 0
		for(i in nTraits){
			dataK 	= dataK + Kcalc(data[,i], tree, F)
			newK 	= newK  + Kcalc(new[,i], tree, F)
		}
		dataK = dataK / nTraits
		newK = newK / nTraits
	}

	difs					= as.matrix(dist(data))			# euclidian distance
	difs[which(difs==0)]	= NA							# ignore matrix diagonal
	Dgap					= apply(difs, 1, min, na.rm=T)	

	difs					= as.matrix(dist(new))			# euclidian distance
	difs[which(difs==0)]	= NA							# ignore matrix diagonal
	Ngap					= apply(difs, 1, min, na.rm=T)

	# Use summary statistics: mean and sd of gaps between neighbours
    if(sstat=="std") 	return( abs(mean(Dgap) - mean(Ngap)) + abs(sd(Dgap) - sd(Ngap)))
    if(sstat=="K")		return( abs(mean(Dgap) - mean(Ngap))^2 + abs(sd(Dgap) - sd(Ngap))^2 + abs(dataK - newK)^2)
    if(sstat=="Kutsuk")	return( sum(abs(new - data)) )    # Kutsukake method: compare absolute values (slow)
}

LRT 	= function(tree, data, min=0, max=10, reps=1e3, e=NA, a=NA, sigma=NA, sigma2=NA, dt=1, kernel1="CE", kernel2="BM", lim=5, plot=F, file="sample.out", sstat="std")
{
	if(kernel1=="CE" && kernel2=="BM")
	{
		return( nestedLRT(tree=tree, data=data, min=min, max=max, reps=reps, e=NA, a=NA, sigma=NA, dt=dt, plot=plot, file=file, sstat=sstat) )
	} else {
		return( unnestedLRT(tree=tree, data=data, min=min, max=max, reps=reps, e=NA, a=NA, sigma=NA, sigma2=NA, dt=dt, kernel1=kernel1, kernel2=kernel2, lim=lim, plot=plot) )
	}
}

nestedLRT	= function(tree, data, min=0, max=10, reps=1e3, e=NA, a=NA, sigma=NA, dt=1, lim=5, plot=FALSE, file="sample.out", posteriorSize=500, sstat="std")
{
	# Simulate and write to file as we go. Single threaded.
	for(i in 1:reps)
   	{
   		sig 	= runif(1, min, max)
   		atry 	= runif(1, min, max)
  		dist 	= get_dif(tree, data, atry, sig, sigma2="NA", dt=dt, kernel="CE", lim=lim, sstat=sstat)
	   	write(c(sig, atry, dist), file=file, append=TRUE, sep=",")
   	}

   	# Read simulations back into R
   	x 	= read.csv(file)
   	sig = x[,1]
   	atry= x[,2]
   	dist= x[,3]

	#----------- OLD MULTITHREADED BT MEMORY INTENSIVE METHOD---------------#
	
	# sig 	= runif(reps, min, max)
	# atry    = runif(reps, min, max)
	# H1_dist = 1:reps 				

	# # Simulate and get distances to truth
	# H1_dist = mcmapply(function(H1_dist, atry, sig) {
	#    H1_dist = get_dif(tree, data, atry, sig, sigma2="NA", dt=dt, kernel="CE", lim=lim)
	# }, H1_dist, atry, sig, mc.cores=cores)

	#-----------------------------------------------------------------------#

    # Get simulation 0.5% of way from smallest distance to truth
    # cutoff	= quantile(H1_dist, 0.005)
    # H1_post		= which(H1_dist < cutoff)

    # Get simulations from 500th smallest to smallest
    H1_post		= order(dist)[1:posteriorSize]

    Usig		= sig[H1_post]
	Uatry		= atry[H1_post]
	H1_post		= matrix(ncol=2, nrow=length(Usig))
	H1_post[,1]	= Usig
	H1_post[,2]	= Uatry

	k 	= kde(H1_post, xmin=c(0, 0), xmax=c(max,max))
	k0	= kde(H1_post, xmin=c(0, 0), xmax=c(max,0))		# sigma to max, a to 0.

	k_max_index 	= which(k$estimate == max(k$estimate), arr.ind = TRUE)
	H1_lik 			= k$estimate[k_max_index[1], k_max_index[2]]
	H1_est 			= c(unlist(k$eval.points)[k_max_index[1]], unlist(k$eval.points)[length(k$estimate[,1]) + k_max_index[2]])


	k0_max_index	= which(k0$estimate == max(k0$estimate), arr.ind = TRUE)
	H0_lik 			= k0$estimate[k0_max_index[1], k0_max_index[2]]
	H0_est 			= c(unlist(k0$eval.points)[k0_max_index[1]], unlist(k0$eval.points)[length(k0$estimate[,1]) + k0_max_index[2]])

	LRT				= -2 * log( H0_lik / H1_lik )

	if(plot==TRUE)	plot(k, "persp")

	return( data.frame(H0_est, H0_lik, H1_est, H1_lik, LRT) )
}

unnestedLRT   = function(tree, data, min=0, max=10, reps=1e3, e=NA, a=NA, sigma=NA, sigma2=NA, dt=1, kernel1="LIM", kernel2="BM", lim=5, plot=F)
{
	sig 	= runif(reps, min, max)
	sig2 	= runif(reps, min, max)
    atry    = runif(reps, min, max)
   	H0_dist = 1:reps 					# H0 distances to truth
   	H1_dist = 1:reps 					# H1 distances to truth

   	#-------------------- Simulate H0  -------------------------------------------#
    H0_dist = mcmapply(function(H0_dist, atry, sig, sig2) {
        H0_dist = get_dif(tree=tree, data=data, a=atry, sigma=sig, sigma2=sig2, dt=dt, kernel=kernel2, lim=lim)
    }, H0_dist, atry, sig, sig2, mc.cores=cores)

    # Get H0 simulation 0.5% of way from smallest
    cutoffH0	= quantile(H0_dist, 0.005)

    # Which simulations smaller than cutoff?
    H0_post		= which(H0_dist < cutoffH0)
	Usig		= sig[H0_post]
	Usig2		= sig2[H0_post]
	Uatry		= atry[H0_post]
	H0_post		= matrix(ncol=2, nrow=length(Usig))
	H0_post[,1]	= Usig
	H0_post[,2]	= Usig2
	# H0_post[,2]	= Uatry

    # Get H0 MLEs with kernel density function to find likelihood peak
	k 	= kde(H0_post, xmin=c(0, 0), xmax=c(max,max))
	k_max_index 	= which(k$estimate == max(k$estimate), arr.ind = TRUE)
	H0_lik 			= k$estimate[k_max_index[1], k_max_index[2]]
	H0_est 			= c(unlist(k$eval.points)[k_max_index[1]], unlist(k$eval.points)[length(k$estimate[,1]) + k_max_index[2]])

	#-------------------- Simulate H1 -------------------------------------------#
    H1_dist = mcmapply(function(H1_dist, atry, sig, sig2) {
        H1_dist = get_dif(tree, data, atry, sig, sig2, dt=dt, kernel=kernel1, lim=lim)
    }, H1_dist, atry, sig, sig2, mc.cores=cores)

    # Get H1 simulation 0.5% of way from smallest
    cutoff	= quantile(H1_dist, 0.005)

    # Which simulations smaller than cutoff?
    H1_post	= which(H1_dist < cutoff)
    Usig		= sig[H1_post]
	Uatry		= atry[H1_post]
	H1_post		= matrix(ncol=2, nrow=length(Usig))
	H1_post[,1]	= Usig
	H1_post[,2]	= Uatry

    # Get H0 MLEs with kernel density function to find likelihood peak
	k1 	= kde(H1_post, xmin=c(0, 0), xmax=c(max,max))
	k1_max_index 	= which(k1$estimate == max(k1$estimate), arr.ind = TRUE)
	H1_lik 			= k1$estimate[k1_max_index[1], k1_max_index[2]]
	H1_est 			= c(unlist(k1$eval.points)[k1_max_index[1]], unlist(k1$eval.points)[length(k1$estimate[,1]) + k1_max_index[2]])

	#-------------------------- resimulate H0 and H1 with MLEs and compare likelihoods ----------#

    H0_dist = mcmapply(function(H0_dist)
    {
        H0_dist = get_dif(tree, data, a=0, sigma=H0_est[1], sigma2=H0_est[2], dt=dt, kernel=kernel2, lim=lim, force=TRUE)
    }, H0_dist, mc.cores=cores)

    H1_dist = mcmapply(function(H1_dist)
    {
        H1_dist = get_dif(tree, data, sigma=H1_est[1], a=H1_est[2], dt=dt, kernel=kernel1, lim=lim)
    }, H1_dist, mc.cores=cores)

    # Get H0 simulation 0.5% of way from smallest
    cutoff	= quantile(H0_dist, 0.005)

    # How many simulations smaller than cutoff?
    H0_lik 	= length( which(H0_dist < cutoff) )
    H1_lik 	= length( which(H1_dist < cutoff) )
    LRT		= -2 * log( H0_lik / H1_lik )

	# return( data.frame(H0_lik, H1_lik, LRT) )
	return( data.frame(H0_est, H0_lik, H1_est, H1_lik, LRT) )
}
