source('den.R')

tre = rUMT(150, mu=0)
dat = genTree(tre, nTraits=2, sigma=5, a=25, cov=0.9, dt=0.2)

plot(dat)

pic1 = pic(dat[,1], tre)
pic2 = pic(dat[,2], tre)

points(pic1, pic2, col=2)

print( cor(dat[,1], dat[,2]) )
print( cor(pic1, pic2) )
