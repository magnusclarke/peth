source('den.R')
t = rUMT(20)
d = genTree(t, sigma=2, a=0)
print(sd(d$traits))
d = genTree(t, sigma=1, a=1)
print(sd(d$traits))

