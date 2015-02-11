source('den.R')

#tre = rUMT(150, mu=0.0)
tre = read.tree('~/Desktop/Harmon2010/geospiza.phy')
dat = genTree(tre, nTraits=2, sigma=0.3, a=1, cov=0.99, dt=0.01)

#ddat = read.table('~/Desktop/Harmon2010/geospiza.dat', header=T)
#tre = drop.tip(tre, 'olivacea')

#dat = matrix(nrow = 13, ncol=2)
#dat[,1] = ddat$culmenL
#dat[,2] = ddat$beakD

plot(dat)

pic1 = pic(dat[,1], tre, scaled=F)
pic2 = pic(dat[,2], tre, scaled=F)

points(pic1, pic2, col=2)

print( cor(dat[,1], dat[,2]) )
print( cor(pic1, pic2) )
