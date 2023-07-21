# This function is a convenient way to run a large number of ergm simulations in 
# parallel, given a list of model parameters.
#
# written by Gianmarc Grazioli 

source("fibril_assay.R")
library(ergm)
library(sna)
library(parallel)

simulate_from_param_list <- function(paramList, stats, fibType, nodeCount, 
                                     np=16, burnTime=sum(2^(1:20))){
  getYield <- function(netIn){
    counts<-fibril_assay(netIn)
    return(sum(counts[[fibType]])/as.double(nrow(counts)))
  }
  
  run_sim_for_parallel <- function(param){
    initGraph <- network.initialize(nodeCount, directed=F)
    formula <- as.formula(initGraph~stats)
    simResult = simulate(as.formula(paste("initGraph", stats, sep="~")), coef=param, 
                         control=control.simulate.formula(MCMC.burn=burnTime,MCMC.interval=1),
                         constraint=~bd(maxout=12),mc.cores=np)
    yield = getYield(simResult)
    return(list("simResult" = simResult, "yield" = yield, "param" = param))
  }
  gen_zero <- mclapply(paramList, run_sim_for_parallel, mc.cores = np)
  yields = c(); params = list(); simResults = list();
  for (i in 1:length(gen_zero)) {
    params[[i]] = gen_zero[[i]]$param
    simResults[[i]] = gen_zero[[i]]$simResult
    yields = c(yields, gen_zero[[i]]$yield)
  }
  return(list("params" = params, "stats" = stats, "yields" = yields))
}