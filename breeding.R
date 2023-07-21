# Given two parent ERGM parameters, this returns a new set of parameters that
# includes both parents, children points along the line connecting the parents, 
# and mutant children that are noised up versions of the other children.
#
# written by Gianmarc Grazioli 

get.all.idx.pairs <- function(numList){
  len = length(numList)
  out <- vector("list", len*(len-1)/2)
  ctr=1
  for (i in 1:len) {
    for (j in 1:len) {
      if(i < j){
        out[[ctr]] = c(i,j)
        ctr = ctr + 1
      }
    }
  }
  return(out)
}

breedTwoParents <- function(vec1, vec2, numPointsInit, useLineDensity=TRUE, minLineDensity=.5){
  basis = vec2 - vec1
  if(useLineDensity){
    distance <- sqrt(sum(basis^2))
    lineDensity <- numPointsInit/distance
    if(lineDensity < minLineDensity){
      numPoints <- round(minLineDensity*distance)
    }else{numPoints <- numPointsInit}
  }else{numPoints <- numPointsInit}
  steps = seq(from = 0, to = 1, length.out = numPoints)
  out <- list()
  for(i in 1:length(steps)){
    out[[i]] <- vec1 + steps[[i]]*basis
  }
  if(length(out) > 2){
    out = out[2:(length(out)-1)]
    #parents returned later, so only return the children
  }
  return(out)
}

addNoiseToPoint <-function(vec, noiseVar=.05){
  out = c()
  for (v in vec){
    out = c(out, rnorm(1, v, noiseVar))
  }
  return(out)
}

breedGeneration <- function(parents, numLinKids, numMutants, noiseVar=.1, 
                            useLineDensity=TRUE, minLineDensity=.5, pairMax=10000, childMax = 100000){
  # create list of all possible parent pairs, then use breedTwoParents() to create 
  # an entire generation of new parameters.
  if(numLinKids < 3){print("If numLinKids < 3, then no points between will be produced.")}
  myPairs <- get.all.idx.pairs(parents)
  if(length(myPairs) > pairMax){
    print("Max pairs reached, proliferation may be too high.")
    myPairs = sample(myPairs, pairMax)
  }
  out = list()
  for (p in myPairs){
    par1 = parents[[p[[1]]]]
    par2 = parents[[p[[2]]]]
    if(sum(abs(par1 - par2)) > .000000001){# avoid underflow errors
      if(sqrt(sum((par1 - par2)^2))*minLineDensity > 1){
        # if pair not too closely related so that minLineDensity*distance would 
        #produce 1 or more children in the space between, breed them:
        out <- c(out, breedTwoParents(par1, par2, numLinKids, 
                                      useLineDensity=useLineDensity, 
                                      minLineDensity=minLineDensity))  
      }
    }
  }
  if(length(out) < 2){
    print("No breeding pairs formed, returning parents, try increasing minLineDensity")
    out = parents
  }
  if(length(out) >= 2 && numMutants > 0){
    mutants <- list()
    ctr = 1
    for (i in 1:length(out)){
      for (j in 1:numMutants){
        mutants[[ctr]] <- addNoiseToPoint(out[[i]], noiseVar = noiseVar)
        ctr = ctr + 1
      }
    }
    out = c(out, mutants)
    if(length(out) > childMax){
      out = sample(out, childMax)
    }
    parentMutants <- list()
    ctr = 1
    for (i in 1:length(parents)){
      for (j in 1:numMutants){
        parentMutants[[ctr]] <- addNoiseToPoint(parents[[i]], noiseVar = noiseVar)
        ctr = ctr + 1
      }
    }
    out = c(out, parentMutants, parents)
  }
  return(out)
}

testIt = F
if(testIt){

gen0 = list(c(1,2), c(2,7), c(3, 13), c(10,2), c(11,7), c(12, 13))
gen0_x = lapply(gen0, '[[', 1)
gen0_y = lapply(gen0, '[[', 2)
plot(gen0_x, gen0_y, col=rgb(.9, 0, .7,.5), pch=19)

gen1 = breedGeneration(gen0, 3, 1, noiseVar = .05, minLineDensity = 2.50)#was 2.5
print("gen 1 breeding complete")
gen1_x = lapply(gen1, '[[', 1)
gen1_y = lapply(gen1, '[[', 2)
plot(gen1_x, gen1_y, col=rgb(.9, 0, .7,.5), pch=19)

gen2 = breedGeneration(gen1, 3, 1, noiseVar = .15, minLineDensity = .4, childMax = 4000)
print("gen 2 breeding complete")
gen2_x = lapply(gen2, '[[', 1)
gen2_y = lapply(gen2, '[[', 2)

plot(gen2_x, gen2_y, col=rgb(0, .8, .8, .2), pch=19)
points(gen1_x, gen1_y, col=rgb(0, 0, .7, .5), pch=19)
points(gen0_x, gen0_y, col=rgb(.88, 0, .88, 1), pch=19)
cat("gen2 has",length(gen2),"points.")
}
