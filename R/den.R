# Functions to simulate traits on trees, with interspecific interactions
# magnus

library("ape")
library("TESS")
library("ks")
library("picante")

if(.Platform$pkgType == "mac.binary")	dyn.load("../cpp/Rfunc_mac.so")
if(.Platform$pkgType == "source")		dyn.load("../cpp/Rfunc.so")

# random ultrametric tree - RENAME
rUMT	= function(nt, lambda=1, mu=0)
{
	tree 	= sim.globalBiDe.taxa(n=1, nTaxa=nt, lambda=lambda, mu=mu)	
	return( tree[[1]] )
}

# Simulate tip trait data
genTree	= function(tree, a=0, sigma=1, dt=1, nTraits=1, kernel="CE", lim=0) 
{
	dt 	= 0.01 * dt

	# Require tree$edge to be a matrix of integers not doubles. TESS gives doubles.
	tree$edge = apply(tree$edge, c(1,2), function(x) x = as.integer(x))	

	tip_count      = length( tree$tip.label )
	st             = tree$edge[,1]
	end            = tree$edge[,2]
	length         = tree$edge.length
	edge_count     = length(st)

	if(kernel=="BM")	k=0
	if(kernel=="CE")	k=1
	if(kernel=="LIM")	k=2
	if(kernel=="NF")	k=3

	result	= .C("genTree", ec=edge_count, nc = tip_count,
				nt = as.integer(nTraits), kernel=as.integer(k), ratecut=as.integer(1),
				a=as.double(a), start=st, end, len=as.double(length),
				sigma=as.double(sigma), sigma2=as.double(1), dt=as.double(dt), lim=as.double(lim),
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
get_dif	= function(tree, data, a, sigma, force=FALSE, dt=1, kernel="CE", lim=0, nTraits=1, sstat="std") 
{
	ntips	= length(data[,1])
	nTraits	= length(data[1,])
	new		= genTree(tree=tree, a=a, sigma=sigma, dt=dt, nTraits=nTraits, kernel=kernel, lim=lim)		# simulate dataset

	# Get Blomberg's K for first trait
    if(sstat=="K")
    {
	    dataK 	= Kcalc(data[,1], tree, F)
	    newK 	= Kcalc(new[,1], tree, F)
    }
	
    difs					= as.matrix(dist(data))			# euclidian distance
	difs[which(difs==0)]	= NA							# ignore matrix diagonal
	Dgap					= apply(difs, 1, min, na.rm=T)	

	difs					= as.matrix(dist(new))			# euclidian distance
	difs[which(difs==0)]	= NA							# ignore matrix diagonal
	Ngap					= apply(difs, 1, min, na.rm=T)

	# Use summary statistics: mean and sd of gaps between neighbours
	if(sstat=="std") 	return( abs(mean(Dgap) - mean(Ngap)) *	abs(sd(Dgap) - sd(Ngap)))
    if(sstat=="K")		return( abs(mean(Dgap) - mean(Ngap)) * abs(sd(Dgap) - sd(Ngap)) * abs(dataK - newK))
    if(sstat=="Kutsuk")	return( sum(abs(new - data)) )    # Kutsukake method: compare absolute values (slow)
}

LRT 	= function(tree, data, min=0, max=10, reps=1e3, e=NA, a=NA, sigma=NA, dt=1, kernel1="CE", kernel2="BM", lim=5, plot=F, file="sample.out", posteriorSize=500, sstat="std")
{
	if(kernel1=="CE" && kernel2=="BM")
	{
		return( nestedLRT(tree=tree, data=data, min=min, max=max, reps=reps, e=NA, a=NA, sigma=NA, dt=dt, plot=plot, file=file, posteriorSize=posteriorSize, sstat=sstat) )
	} else {
		return( unnestedLRT(tree=tree, data=data, min=min, max=max, reps=reps, e=NA, a=NA, sigma=NA, dt=dt, kernel1=kernel1, kernel2=kernel2, lim=lim, plot=plot) )
	}
}

nestedLRT	= function(tree, data, min=0, max=10, reps=1e3, e=NA, a=NA, sigma=NA, dt=1, lim=5, plot=FALSE, useFile=FALSE, file="sample.out", posteriorSize=500, sstat="std")
{
	# Simulate and write to file as we go. Single threaded.
	# For mutithreaded performance, run with same outfile many times, and 
	# then run with useFile=TRUE.
	if(useFile==FALSE)
	{
		for(i in 1:reps)
	   	{
	   		sig 	= runif(1, min, max)
	   		atry 	= runif(1, min, max)
	  		dist 	= get_dif(tree, data, atry, sig, dt=dt, kernel="CE", lim=lim, sstat=sstat)
		   	write(c(sig, atry, dist), file=file, append=TRUE, sep=",")
	   	}
	}

   	# Read simulations back into R
   	x 	= read.csv(file)
   	sig = x[,1]
   	atry= x[,2]
   	dist= x[,3]

    # Get simulations from posteriorSize'th smallest to smallest
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
