options(warn=-1)        # Disable warning messages for whole script

cat('Initialising peth test...')

dir = getwd()
suppressMessages(source('../R/sim.R', chdir=T))

setwd(dir)
cat('\n')

#----------- Check sensible distributions from simulations--------#
check_dist = function(dat, kernel){
    cat(paste('Checking distributions of simulated traits (', kernel, ')... ', sep=''))
    if(mean(dat) < 4 && mean(dat) > -4){
        mean_test= 'OK'
    } else {
        mean_test='nope'
    }
    if(sd(dat) < 10 && sd(dat) > 0.5){
        sd_test= 'OK'
    } else{
        sd_test='nope'
    }
    if(mean_test=='OK' && sd_test=='OK')
    {
        cat('OK\n')
    } else {
        cat('Error: distribution too non-normal!\n')
    }
}
tre = rUMT(1000)
ptre = ape2peth(tre)
datBM = genTree(ptre)$traits
datCOMP = genTree(ptre, a=0, kernel='CE')$traits
datLIM = genTree(ptre, a=2, kernel='LIM')$traits
check_dist(datBM, 'BM')
check_dist(datCOMP, 'Competition')
check_dist(datLIM, 'Bounded competition')
#-----------------------------------------------------------------#

#----------- Check correct number of traits simulated ------------#
cat('Checking trait data dimensionality... ')
tre = rUMT(5)
ptre = ape2peth(tre)
dat = genTree(ptre, nTraits=3)
dat2 = genTree(ptre, a=5, nTraits=6)
if( length(dat[1,])==3 && length(dat2[1,])==6)
{
    cat('OK\n')
} else {
    cat('Error: incorrect number of trait dimensions generated!\n')
}
#-----------------------------------------------------------------#

#----------- Compare VCV matrices for test dataset ---------------#
cat('Checking VCV matrices for test data... ')
tre=0;dat=0
load('test.tre')
load('test.dat')
ptre = ape2peth(tre)
testvcv = asVCV(ptre)
distance = sum(abs(vcv - testvcv))
if(distance < 0){
    cat('Error: negative difference between vcv matrices!\n')
} else if(distance < 4){
    cat('OK\n')
} else {
    cat('Error: vcv matrix too inaccurate!\n')
}

cat('Checking VCV matrix with competition... ')
compvcv = asVCV(ptre, a=2)
sumvcv = sum(abs(vcv))
sumcompvcv = sum(abs(compvcv))
if(sumcompvcv - sumvcv > 6) 
{
    cat('OK\n')
} else {
    cat('Error: competition vcv not significantly larger than BM vcv!\n')
}
#-----------------------------------------------------------------#

#----------- See if LRT results are sensible ---------------------#
tre = rUMT(20)
ptre = ape2peth(tre)
dat = genTree(ptre)
test_lrt = function(lrt1){
    if(any(is.na(lrt1)))
    {
        cat('Error: LRT produced NA values!\n')
    } else if(lrt1$LRT[1] < -10 || lrt1$LRT[1] > 10){
        cat('Error: Extreme likelihood ratios produced!\n')
    } else {
        cat('OK\n')
    }
}
cat('Checking whether std LRT results are sensible... ')
lrt1 = LRT(ptre, dat, sstat='std', reps=1e2, posteriorSize=50); file.remove('sample.out')
test_lrt(lrt1)
cat('Checking whether LRT results with signal (K) are sensible... ')
lrt2 = LRT(ptre, dat, sstat='K', reps=1e2, posteriorSize=50); file.remove('sample.out')
test_lrt(lrt2)
#-----------------------------------------------------------------#

#----------- Simulation speed timer ------------------------------#
tre = 0
load('test.tre')
ptre = ape2peth(tre)
cat('\nTiming a thousand simulations - small tree (par time 0.57s)...\n\t')
t = system.time(replicate(1e3, genTree(ptre, a=1)))
cat(t[1])
cat(' seconds\n')
load('test_large.tre')
ptre = ape2peth(tre)
cat('\nTiming a thousand simulations - large tree (par time 1.65s)...\n\t')
t = system.time(replicate(1e3, genTree(ptre, a=1)))
cat(t[1])
cat(' seconds\n')
#-----------------------------------------------------------------#
