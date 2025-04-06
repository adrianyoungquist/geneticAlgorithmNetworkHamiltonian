# This main function carries out automated discovery of network Hamiltonian models
# that can self-assemble into the 2-ribbon amyloid fibril topological structure.
# The code implements the genetic algorithm described in the article:
#
# "A Genetic Algorithm for Automated Parameterization of Network Hamiltonian 
#                    Models of Amyloid Fibril Formation"
#
# by: Gianmarc Grazioli, Andy Tao, Inika Bhatia, and Patrick Regan 
#
# currently under review at Digital Discovery - Royal Society of Chemistry
# preprint: https://chemrxiv.org/engage/chemrxiv/article-details/64a3d048ba3e99daef7c671a
#
# The example used in this code is very similar to the code used to generate 
# the figure in the aforementioned article showing the evolution of different 
# generations of models as the genetic algorithm converges on a region of parameter 
# space that produces maximal fibril yield for 2-ribbon type amyloid fibril structures.
# To give an idea of run time, this main function runs in under 10 minutes on a 
# 2019 MacBook Pro with a 2.4 GHz 8-Core Intel Core i9 processor and 64 GB of RAM.
#
# Code written by Gianmarc Grazioli 

source("fibril_assay.R")
source("R/evolution.R")
source("R/breeding.R")
source("R/simulate_from_param_list.R")
source("sample_hyperspheroid.R")
library(ergm)
library(sna)
library(parallel)

begin=proc.time()[[3]]
seedVal = 301656
set.seed(seedVal) # set random seed if consistent results desired

nCount=48 #number of nodes 
stats = "edges+kstar(2)+nsp(1:2)"
statCount = 4
fib_type = "2-ribbon"

init_param = c(100.0, -31.13533, 10.70222, 26.31755) #initial best parameter
spheroid.pt.count = 50 #number of random parameters in spheroid surrounding init_param
edgeVal = 100 #fix edge value to this
radius = 2 #radius of sphere to ultimately distort to spheroid 
gen_zero_matrix = sample_hyperspheroid(length(init_param), spheroid.pt.count, radius, init_param, 
                                scaleFac = .15, fixEdge = T, fixEdgeValue=100)
gen_zero_params = list() #create list of parameters from matrix
for (i in 1:nrow(gen_zero_matrix)) {gen_zero_params[[i]] <- gen_zero_matrix[i,]}

cat("Simulating from random hyperspheroid...\n")
spheroid_sims = simulate_from_param_list(gen_zero_params, stats, fib_type, nCount, 
                                         np=16, burnTime=sum(2^(1:20)))

cat("Take the highest fibril fractions from the spheroid points.\n")
gen_zero = spheroid_sims$params[which(spheroid_sims$yields > 0 )]
if(length(gen_zero) > 8){gen_zero = gen_zero[1:8]} 

cat("Begin evolution from generation zero of size ", length(gen_zero), "...\n")
myGenerations_new = evolution(gen_zero, 3, 3, 2, "2-ribbon", stats, nCount, 
                              reps = 16, noiseVarInit = .10, topFrac = .25, maxSurvivors=4, 
                              fixEdge=TRUE, fixEdgeValue=100, smartVar = T,
                              smartPts = T, useLineDensity=T, minLineDensity=3.3, 
                              printParams=FALSE, burnTime = 2^20, childMax = 20)

finalReps = 32 #number of times to repeat final simulations of the best of each generation
reps_of_the_best = list(list(),list(),list(),list(),list()) 
for (i in 1:finalReps) {reps_of_the_best[[1]][[i]] = init_param}
for (i in 1:finalReps) {reps_of_the_best[[2]][[i]] = spheroid_sims$params[which.max(spheroid_sims$yields)][[1]]} 
for (gen in 1:length(myGenerations_new)) {
  for (i in 1:finalReps) {reps_of_the_best[[2+gen]][[i]] = myGenerations_new[[gen]]$params[[1]]}
}

finalSims = list()
for (i in 1:length(reps_of_the_best)) {
  cat("Running final simulations", i, "of", length(reps_of_the_best),"\n")
  finalSims[[i]] <- simulate_from_param_list(reps_of_the_best[[i]], stats, fib_type, nCount, 
                                             np=16, burnTime=sum(2^(1:20)))
}

mean_fibril_fractions = c(); for (i in 1:length(finalSims)) {
  mean_fibril_fractions = c(mean_fibril_fractions, mean(finalSims[[i]]$yields))
}

plot(mean_fibril_fractions, type = 'l')
points(mean_fibril_fractions)
save(myGenerations_new, file = "myGenerations_new.RData")

end=proc.time()[[3]]
timeSec = end - begin; minutes = floor(timeSec/60); seconds = round(timeSec - minutes*60);
cat("Calculation took ", minutes, "minutes and ", seconds, "seconds.\n")
