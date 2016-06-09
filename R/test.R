source('sim.R')

tree = rand_umt(20)

true_data = sim(tree)$tval

max_sigma = 5
max_a = 5
min_param = 0
max_param = 2

reps = 1e3

sig = runif(reps, 0, max_sigma)
atry = runif(reps, 0, max_a)
param = runif(reps, min_param, max_param)

file = 'test.out'

# -------- Generate alternative datasets under competition
for(i in 1:reps)
{
	d = sim(tree=tree, a=atry[i], sigma=sig[i])$tval
	stats = summary_stats(tree=tree, data=d, use_K=FALSE)
	model = 'comp'
   	write(c(sig[i], atry[i], stats, model), file=file, ncolumns=5, append=TRUE, sep=",")
}

# -------- Generate null datasets under EB or OU
for(i in 1:reps)
{
	d = rTraitCont(tree, model='OU', sigma=sig[i], alpha=param[i])

	stats = summary_stats(tree=tree, data=d, use_K=FALSE)
	model = 'OU'
   	write(c(sig[i], atry[i], stats, model), file=file, ncolumns=5, append=TRUE, sep=",")
}

source('sim.R'); lrt_unnested(tree, data=true_data, file=file, reps=1e3, max_sigma=5, 
				max_a=5, min_param=0, max_param=2, posteriorSize=250, nullh='OU')
