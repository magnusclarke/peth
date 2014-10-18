
# script to plot density dependent trees
# sourced by den.R
# slow and messy  -- just for plotting single simulations

P_runSim	<- function(node_count, seg, time, a, sigma, pdat, ptime, stime) 
{
	dt	<- 0.001
	for (j in 1:(time/dt)) 
	{
		for (i in 1:length(seg)) 
		{
			seg[i]	<- seg[i] + sigma*rnorm(1, 0, sqrt(dt)) 
		}

		for(i in 1:(length(seg)-1))
		{
			q <- 0.5 * abs(seg[i] - seg[i+1])
			if(q <= 2.2)		pn = 0.1 * q * (4.4 - q)
			else if (q>2.2 && q<2.6)	pn = 0.49
			else if (q > 2.6)	pn = 0.50
			else			pn = 0.50
			g = dt * a * (0.5 - pn)
			seg[i] = seg[i] - g	
			seg[i+1] = seg[i+1] + g	
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


	# pplot <- data.frame(ptime, pdat) 	
	# plot(pplot, pch=20, cex=0.25)		
	plot(pdat[[1]]~ptime, pch=20, cex=0.1, type="p", ylim=c(-5*sigma*(1), 5*sigma*(1)),
		, xlab="", ylab="", xaxt="n", yaxt="n")
	for(i in 2:node_count)
	{
		lines(pdat[[i]]~ptime, pch=20, cex=0.1, type="p")#, col=i)
	}
	return(value[which(end %in% 1:node_count)])

	# return( length(pdat[[5]]) )
	#return( length(ptime) )
}
