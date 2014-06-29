library("ape")
library("ks")

dat <- read.csv("morphodata.csv", sep="\t")

tre	<- read.tree("skinktree183.txt")

dat_label	<- paste(dat[,10], "_", dat[,11], sep="")

myDat 		<- dat$toe

# remove data for species not in tree
#unwanted	<- which(! dat_label %in% tre$tip.label)
wanted		<- match(tre$tip.label, dat_label)
myDat 		<- myDat[wanted]
dat_label 	<- dat_label[wanted]

# normalise myDat around 0
myDat 		<- myDat - mean(myDat)
myDat 		<- data.frame(myDat)

# LRT
#post_depth_pect	<- ABC(tre, myDat, reps=1e3, max=5)
#post_head_width		<- ABC(tre, myDat, reps=1e3, max=3)
#post_head_width2		<- ABC(tre, myDat, reps=2e3, max=3)
#post_shank				<- ABC(tre, myDat, reps=2e3, max=3)
post_toe				<- ABC(tre, myDat, reps=2e3, max=3)
LRT(post_head_width2)


# plot likelihood surface
head_width_surf 	<- kde(post_head_width)
depth_pect_surf 	<- kde(post_depth_pect)
shank_surf 			<- kde(post_shank)
toe_surf 			<- kde(post_toe)

png("scinc_plots.png", width=1000, height=1000, res=130)
	par(mfrow=c(2,2))
	plot(tre, show.tip.label=F)
	plot(head_width_surf, xlab="BM rate", ylab="competition strength")#xlim=c(0.1, 0.4), ylim=c(0, 1.5), 
	plot(shank_surf, xlab="BM rate", ylab="competition strength")
	plot(toe_surf, xlab="BM rate", ylab="competition strength")
dev.off()






#----------------------------------------------------------------------------




### Repeat for BM vs lim
dir	<- getwd()
source("~/Documents/phd/Iceberg/peth/R/den.R", chdir=T)
setwd(dir)
library("ape")
library("ks")

dat <- read.csv("morphodata.csv", sep="\t")

tre	<- read.tree("skinktree183.txt")

dat_label	<- paste(dat[,10], "_", dat[,11], sep="")

myDat 		<- dat$toe

# remove data for species not in tree
#unwanted	<- which(! dat_label %in% tre$tip.label)
wanted		<- match(tre$tip.label, dat_label)
myDat 		<- myDat[wanted]
dat_label 	<- dat_label[wanted]

# normalise myDat around 0
myDat 		<- myDat - mean(myDat)
myDat 		<- data.frame(myDat)

# Set hard limits at the most extreme trait value of the sample
sim_lim		<- c(max(myDat), -min(myDat))
sim_lim 	<- max(sim_lim)

lrt_toe	<- manualLRT(tre, dat=myDat, min=0, max=2, reps=1e3, e=NA, a=NA, sigma=NA, dt=1, kernel1="BM", kernel2="BM", lim=sim_lim, plot=F)




# LRT
#post_depth_pect	<- ABC(tre, myDat, reps=1e3, max=5)
#post_head_width		<- ABC(tre, myDat, reps=1e3, max=3)
#post_head_width2		<- ABC(tre, myDat, reps=2e3, max=3)
#post_shank				<- ABC(tre, myDat, reps=2e3, max=3)
post_toe				<- ABC(tre, myDat, reps=2e3, max=3)
LRT(post_head_width2)


# plot likelihood surface
head_width_surf 	<- kde(post_head_width)
depth_pect_surf 	<- kde(post_depth_pect)
shank_surf 			<- kde(post_shank)
toe_surf 			<- kde(post_toe)

png("scinc_plots.png", width=1000, height=1000, res=130)
	par(mfrow=c(2,2))
	plot(tre, show.tip.label=F)
	plot(head_width_surf, xlab="BM rate", ylab="competition strength")#xlim=c(0.1, 0.4), ylim=c(0, 1.5), 
	plot(shank_surf, xlab="BM rate", ylab="competition strength")
	plot(toe_surf, xlab="BM rate", ylab="competition strength")
dev.off()


