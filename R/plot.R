# Script to plot evolution on trees with character displacement.
# Slow and messy  -- just for plotting single simulations.
# Mainly a transcription of sim.cpp into R.
# Magnus
# 2016

source('tree.R')

##---------  Colour-blind palette -------------------------------#

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
				"#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
				"#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##---------------------------------------------------------------#

plotsim = function(tre, deltat=0.0025, sigma=1, a=0, ntraits=1, symp=NA, plot=TRUE, plotlim=7, limit=NA, alpha=0)
{
	if(class(tre)=="phylo")
	{
		tree = ape2peth(tre)
	} else if(class(tre)!="pethtree") {
		stop("Tree incorrectly formatted.")
	}

	num_tips = length(tree$data_order)
	splitting_nodes = tree$splitting_nodes 
	times = tree$times

	subclades=FALSE

	myalpha <<- alpha

	if(all(is.na(symp)))
	{
		symp = matrix(0, nrow=num_tips, ncol=num_tips)
	} else if(class(symp)=="matrix") {
		symp = symp
	} else if(symp=='subclades')
	{
		subclades = TRUE

		symp = matrix(0, nrow=num_tips, ncol=num_tips)

		# Find two monophyletic subclades
		sub1 = 2
		sub1 = c(sub1, which(splitting_nodes == 2)+1 )
		for (n in sub1) 
		{
			sub1 = c( sub1 , which(splitting_nodes == n)+1)
		}

		sub2 = which(! 1:num_tips %in% sub1)

		symp[sub1 , sub2] = 9e9
		symp[sub2 , sub1] = 9e9
	}

	result = pathsim(tre=tre, ntip=as.integer(num_tips), deltat=as.double(deltat), 
				rate = as.double(sigma^2), a=as.double(a), intervals=as.double(times), 
				splitters=as.integer(splitting_nodes),
				ntraits=as.integer(ntraits), symp=symp, limit=limit
				)
	
	plotlim = max(abs(pdat), na.rm=TRUE)

	if(plot==TRUE)
	{
		plot(pdat[,1], pch=20, cex=0.25, cxy=0.25, type="l", ylim=c(-plotlim, plotlim),
			, xlab="", ylab="", xaxt="n", yaxt="n", bty='n', col=cbPalette[1])

		if(subclades==TRUE)
		{
			for(i in sub2[-1])
			{
				lines(pdat[,i], pch=20, cex=0.25, cxy=0.25, type="l", col=cbPalette[1])
			}

			for(i in sub1)
			{
				lines(pdat[,i], pch=20, cex=0.25, cxy=0.25, type="l", col=cbPalette[2])
			}
		} else {
			for(i in 2:num_tips)
			{
				lines(pdat[,i], pch=20, cex=0.25, cxy=0.25, type="l", col=cbPalette[i%%7+1])
			}			
		}
	}


	ape_tval = matrix(ncol=ntraits, nrow=num_tips)
	for (i in 1:ntraits) 
	{
		ape_tval[,i] = result[seq(i, length(result), by=ntraits)]
		ape_tval[,i] = ape_tval[,i][tree$data_order]
	}


	return(ape_tval)

}

#--- Generate a matrix of the times at which lineages come into sympatry. -------------#
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
	# age = age[t$data_order]

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


# Random Number Generator  
d = function()
{
	rnorm(1)
}

# get simulation
pathsim = function(tre, ntip, deltat, rate, a, intervals, splitters, ntraits, symp, limit)
{
	num_segment  <<-ntip-1
	mydt <<- deltat

	mysymp <<- symp / mydt 	# mysymp measured in timestep units

	tree_speciators <<- splitters

	lim <<- limit

	mya <<- a

	mytime <<- 0
	x <<- replicate(ntip, 0)
	rate <<- rate * sqrt(mydt)
	total_time <<- max(node.depth.edgelength(tre))
	total_time_steps <<- total_time / mydt
	segment_steps <<- intervals / mydt


	pdat <<- matrix(NA, ncol=ntip, nrow=total_time_steps)


	myroot <<- 0
	tval <<- myroot

	sumDist <<- 0
	sumSqDist <<- 0
	dists <<- 0

	path()

	return( tval )
}

# Simulate dataset
path = function()
{
	for(i in 1:num_segment)
	{
		xx <<- as.integer( tree_speciators[i] )
		tval <<- c(tval, tval[xx])
		evolve_segment(segment_steps[i])
	}
} 

# 	Modify trait values for one segment (i.e. between speciation
#	events). Needs parameter: number of time steps within that segment. 
evolve_segment = function(nsteps)
{
	for (step in 1:nsteps)
	{
		step_segment()
		pdat[mytime , 1:length(tval)] <<- tval
	}
}

#	Update trait values for one time step.
step_segment = function()
{
	num_species <<- length(tval)
	for (species in 1:num_species)
	{
		step_species(species)
	}

	mytime <<- mytime + 1
}

# Do one evolutionary step on one species. 
step_species = function(species)
{
	tval[species] <<- tval[species] + d() * rate - myalpha * tval[species] * mydt		# OU
	
	# Loop over species i that are interacting with this species.
	for (i in 1:num_species)
	{	
		if(mytime > mysymp[species,i])
		{
			interaction(species, i)
		}
	}

	if(!is.na(lim))
	{
		if(tval[species] > lim)
		{
			tval[species] <<- lim 
		}

		if(tval[species] < -lim)
		{
			tval[species] <<- -lim 
		}
	}
}

# Interaction between two species; updates trait values. 
interaction = function(s1, s2)
{
	update_distance(s1, s2)
	
	# Compute resource distribution overlap g (0 to 1) between two species.
	q 	= 0.5 * sqrt(sumSqDist)
	pn 	= pnorm(q, 0, 1, FALSE, FALSE)
	g 	= 2 * mydt * mya * pn
	
	# Update trait values.
	if (sumDist != 0)
	{
		tval[s1] <<- tval[s1] + g*(dists/sumDist)
		tval[s2] <<- tval[s2] - g*(dists/sumDist)
	}
}

update_distance = function(s1, s2)
{
	sumDist <<- 0
	sumSqDist <<- 0

	# Loop over traits and get trait distances between species s1 and s2.
	dists			<<- tval[s1] - tval[s2]
	sqDists			<<- dists*dists
	sumDist			<<- sumDist + abs(dists)
	sumSqDist		<<- sumSqDist + sqDists;
}


