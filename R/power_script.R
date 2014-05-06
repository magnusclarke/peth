require("ks")
source("den.R")

tree <- randUMT(30)
data <- genTree(tree, sigma=1, a=1)
min=0
max=3
reps=1e4
a=NA
sigma=NA
dt=1

use 	<- rep(FALSE, reps)
sig		<- runif(reps, min, max)
atry	<- runif(reps, min, max)

treps	<- reps/1e1
test 	<- replicate(treps, get_dif(tree, data, runif(1, min, max), runif(1, min, max), dt=dt) )
ep	<- 2*min(test)

use <- mcmapply(function(use, atry, sig)	
	{
		dif <- get_dif(tree, data, atry+1, sig, dt=dt)
		if(dif < ep)	use	<- TRUE
		else		use	<- FALSE
	}, use, atry, sig, mc.cores=cores)

Usig	<- sig[which(use == TRUE)]
Uatry	<- atry[which(use == TRUE)]
dat 	<- matrix(ncol=2, nrow=length(Usig))
dat[,1]	<- Usig
dat[,2]	<- Uatry

k 	<- kde(dat)
k0

