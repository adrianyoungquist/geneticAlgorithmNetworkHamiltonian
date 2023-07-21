# This code handles the entire evolution process for the genetic algorithm, with
# the exception of the breeding, i.e. sharing of parameter information from 
# parent models to child models. 
#
# written by Gianmarc Grazioli 

source("fibril_assay.R")
source("breeding.R")
library(ergm)
library(sna)
library(parallel)

runOneSim <- function(param, suffStats, nodeCount, burnTime=2^20, np=16){
  initGraph <- network.initialize(nodeCount,directed=F)
  formula <- as.formula(initGraph~suffStats)
  out<-simulate(as.formula(paste("initGraph", suffStats, sep="~")), coef=param, 
                control=control.simulate.formula(MCMC.burn=burnTime,MCMC.interval=1),
                constraint=~bd(maxout=12),mc.cores=np)
  return(out)
}

simAllparams <- function(params, suffStats, nodeCount, burnTime=2^20, np=16){
  funcForParallel <- function(p){
    sim<-runOneSim(p, suffStats=suffStats, nodeCount=nodeCount, burnTime=burnTime, np=np)
    return(list("sim" = sim, "params" = p, "stats" = suffStats))
  }
  sims <- mclapply(params, funcForParallel, mc.cores = np) 
  return(sims)
}

thinTheHerd <- function(sims, fibrilType, topFrac=.5, maxSurvivors=20){
  fibFracsWithParams <- list()
  for (par in 1:length(sims[[1]])) {#loop over each parameter in the rep
    fibFracTot = 0
    for (rep in 1:length(sims)) {#loop over reps containing one of each parameter
      nets <- sims[[rep]][[par]]$"sim"
      params <- sims[[rep]][[par]]$"params"
      stats <- sims[[rep]][[par]]$"stats"
      fibFracData <- get(fibrilType, fibril_assay(nets))
      fibFracTot <- fibFracTot + sum(fibFracData)/as.double(length(fibFracData))
    }
    fibFracsWithParams[[par]] <- list("fibFrac" = fibFracTot/length(sims), "params" = params)
  }
  topN = floor(topFrac*length(fibFracsWithParams))
  if(topN < 2){cat("topN in thinTheHerd() dropped below 2\n"); topN=2}
  if(topN > maxSurvivors){topN = maxSurvivors}
  
  fibFracsOnly <- c()
  paramsOnly <- list()
  for (i in 1:length(fibFracsWithParams)) {
    fibFracsOnly <- c(fibFracsOnly, fibFracsWithParams[[i]]$fibFrac)
    paramsOnly[[i]] <- fibFracsWithParams[[i]]$params
  }
  
  theFittest = order(fibFracsOnly, decreasing = TRUE)[1:topN]
  bestYields = sort(fibFracsOnly, decreasing = TRUE)[1:topN]
  return(list("params"=paramsOnly[theFittest], "stats"=stats, "yields"=bestYields))
}

evolution <- function(origParents, numGenerations, numLinKids, numMutantsInit, 
                      fibrilType, suffStats, nodeCount, topFrac=.5, noiseVarInit=.1, 
                      reps=1, burnTime=2^20, np=16, maxSurvivors=20, fixEdge=F, 
                      fixEdgeValue=100, initYieldCut=0.0, smartVar = T, 
                      fudge=0.0, smartPts = T, useLineDensity=TRUE, 
                      minLineDensity=.5, printParams=FALSE, pairMax = 10000, childMax = 100000){
  generations <- list()
  parents <- origParents
  yieldCut <- initYieldCut
  noiseVar <- noiseVarInit
  numMutants <- numMutantsInit
  if(printParams){
    sink("all_params.csv")
    sink()
  }
  for (g in 1:numGenerations) {
    cat("Generation ", g, " is born!\n", sep = '')
    newGeneration <- breedGeneration(parents, numLinKids, numMutants, noiseVar = noiseVar, 
                                     useLineDensity=useLineDensity, minLineDensity=minLineDensity,
                                     pairMax = pairMax, childMax = childMax)
    if (fixEdge){
      for (p in 1:length(newGeneration)) {
        newGeneration[[p]][1] <- fixEdgeValue 
      }
    }
    if(printParams){
      sink("all_params.csv", append = TRUE)
      for(p in newGeneration){
        #print(c(p,g))
        lineToPrint = c(unlist(p), g)        
        cat(lineToPrint, sep = ', ')
        #cat(c(p,g), sep = ', ')
        cat('\n')
      }
    sink()
    }
    simReps <- list()
    cat("Simulating",length(newGeneration)*reps,"parameters...\n")
    for (r in 1:reps) {
      simReps[[r]] <- simAllparams(newGeneration, suffStats, nodeCount, burnTime=burnTime, np=np)
    } 
    survivors <- thinTheHerd(simReps, fibrilType, topFrac = topFrac, maxSurvivors=maxSurvivors)
    cat(length(survivors$params), "survivors remaining\n")
    cat("Best yields =\n", survivors$yields,"\n")
    generations[[g]] <- survivors
    # Prevent backtracking to weaker parameters with these lines:
    if(g > 1){
      if( (survivors$yields[1] >= (generations[[g-1]]$yields[1])*(1.0 - fudge)) & (survivors$yields[1] >= yieldCut) ){
        parents <- survivors$params
        yieldCut <- survivors$yields[1]
        noiseVar <- noiseVarInit
        numMutants <- numMutantsInit
      } else {
          cat("Offspring did not surpass the parents, rerun the parents.\n")
          if(smartVar){
            noiseVar = noiseVar + noiseVar*.05
            cat("Increasing noise variance to ", noiseVar, ".\n", sep = '')
            if(smartPts){
              numMutants = numMutants + 1
              cat("Increasing number of mutants to ", numMutants, ".\n", sep = '')
            }
          }
        }  
    } else {
        if(survivors$yields[1] >= yieldCut){
          parents <- survivors$params  
        } else {cat("First generation did not surpass the initYieldCut, rerun the originals.\n")}
      } 
  }
  return(generations)
}

testIt = FALSE
if(testIt){
  nCount=48 #number of nodes
  stats = "edges+kstar(2)+nsp(1:2)+esp(0)"
  theOriginals = list(c(157.912 -log(nCount),-24,-3.3,-2.7,-20), 
                      c(157.912 -log(nCount),-24,-3.3,-2.7,-20))
  myGenerations = evolution(theOriginals, 3, 2, 2, "1,2 2-ribbon", stats, nCount)
}
