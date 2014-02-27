dir	<- getwd()
source("R/den.R", chdir=T)
setwd(dir)

# Create a random ultrametric tree with 4 tips
tree	<- randUMT(4)

# Plot simulated values over the tree with BM
simPlot(tree)
# And with density-dependence
simPlot(tree, a=5)

# Print VCV-matrix representation of tree
vcv(tree)
# Print empirical VCV matrices from simulations
simVCV(tree)
simVCV(tree, a=5)

# Simulate dataset of trait values
tree	<- randUMT(30)
dat	<- genTree(tree, a=0, sigma=10)

# Find MLE of sigma with ABC
ABC(tree, dat, a=0)

# Find MLE for density-dependence parameter, given sigma
dat	<- genTree(tree, a=5, sigma=8)
ABC(tree, dat, sigma=8)

# Jointly estimate a and sigma (takes a minute)
dat		<- genTree(tree, a=5, sigma=10)
ABC(tree, dat, reps=1e4)

# Estimate likelihood of density model relative to BM (takes a few minutes)
dat		<- genTree(tree, a=2, sigma=2)
AIC(tree, dat) 
