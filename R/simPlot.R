
# script to plot density dependent trees
# sourced by den.R
# slow and messy  -- just for plotting single simulations

P_runSim	<- function(node_count, seg, time, a, sigma, pdat, ptime, stime) 
{
	dt	<- 0.002
	# Time loop; j labels time
	for (j in 1:(time/dt)) 
	{
		# i labels species; seg is a vector of trait values
		for (i in 1:length(seg)) 
		{
			seg[i]	<- seg[i] + sigma*rnorm(1, 0, sqrt(dt)) 	# do BM evolution on seg
		}

		# i labels species ; SURELY WE SHOULD HAVE i AND j FOR EVERY SPECIES!!
		# for(i in 1:(length(seg)-1))
		# {
		# 	q <- abs(seg[i] - seg[i+1])			#distance from species i to species i+1
		# 	if(q <= 2.2) {
		# 		pn = 0.1 * q * (4.4 - q)
		# 	} else if (q>2.2 && q<2.6) {
		# 		pn = 0.49
		# 	} else if (q > 2.6) {
		# 		pn = 0.50
		# 	} else {
		# 		pn = 0.50
		# 	}
		# 	g = 0.5 * dt * a * (0.5 - pn)		# interspecific effect, scaled with dt and a
		# 	seg[i] = seg[i] - g					# species i gets pushed 'down' away from species i+1
		# 	seg[i+1] = seg[i+1] + g				# species i+1 gets pushed 'up'
		# }


		for(i in 1:length(seg))
		{
			for (k in 1:length(seg))
			{
				q <- abs(seg[i] - seg[k])			# distance from species i to species k
				if(q <= 2.2) {
					pn = 0.1 * q * (4.4 - q)
				} else if (q>2.2 && q<2.6) {
					pn = 0.49
				} else if (q > 2.6) {
					pn = 0.50
				} else {
					pn = 0.50
				}
				g = 0.5 * dt * a * (0.5 - pn)		# interspecific effect, scaled with dt and a
				if(seg[i]>seg[k]){
					seg[i] = seg[i] + g				# species i gets pushed away from species k
					seg[k] = seg[k] - g				# species k gets pushed away from species i
				} else if(seg[k]>=seg[i]){
					seg[k] = seg[k] + g	
					seg[i] = seg[i] - g	
				}
			}
		}


		stime <- stime + dt	
		for(i in 1:length(seg))
		{		
			pdat[[i]] <- c(pdat[[i]], seg[i])		
			#ptime <- c(ptime, stime)
		}
		if(length(seg) < node_count)
		{
			for(i in (length(seg)+1):node_count)
			{
				pdat[[i]] <- c(pdat[[i]], NA)
			}
		}
		ptime <- c(ptime, stime)
	}
	return(list(seg, pdat, ptime, stime))
}


simPlot	<- function(tree, a=0, sigma=1) 
{
	root_val<- 0

	root_node	<- tree$edge[1,1]
	node_count	<- length( tree$tip.label )		# really tip count
	st	<- tree$edge[,1]
	end	<- tree$edge[,2]
	length	<- tree$edge.length
	value	<- rep(0, length(st))
	running	<- rep(F, length(st))
	running[ which(st == root_node) ]	<- T

	pdat <- as.list( rep(NA, node_count) )
	ptime <- c(NA)	
	stime <- 0

	while(any(running==T)){
		time	<- min(length[ which(running == T) ]) 
		runResult	<- P_runSim(node_count, value[which(running==T)], time, a, sigma, pdat, ptime, stime)	
		pdat <- runResult[[2]]		
		ptime <- runResult[[3]]	 
		stime <- runResult[[4]]	

		value[ which(running==T) ] <- runResult[[1]]

		length[ which(running == T) ]	<- length[ which(running == T) ] - time
		x<- which(length == 0 & running == T)
		node	<- end[x]
		val	<- value[x]
		new_row	<- which(st %in% node)
		value[new_row]	<- val
		running[new_row]<- T		
		running[ which(length == 0) ]	<- F
	}

	plot(pdat[[1]]~ptime, pch=20, cex=0.1, type="p", ylim=c(-5, 5),
		, xlab="", ylab="", xaxt="n", yaxt="n")#, bty='n')
	for(i in 2:node_count)
	{
		lines(pdat[[i]]~ptime, pch=20, cex=0.1, type="p")#, col=i)
	}
	return(value[which(end %in% 1:node_count)])
}
