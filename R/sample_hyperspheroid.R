# This function samples randomly distributed points within a
# spheroid or hyperspheroid of any number of dimensions. The
# points are centered around a user-defined center, and the
# spread in each dimension is a function of the magnitude of 
# position in each dimension (i.e. if the center is at the point
# {5, 100}, the spread will be wider about y = 100 than x = 5).
#
# written by Gianmarc Grazioli 

normIt<-function(vec){
  norm<-sqrt(sum(vec^2))
  vec/norm
}

getHyperSphereFromOneCenter<-function(d, n, radius, center){
  coords<-list()
  for(i in 1:n){
    rawPts<-rnorm(d+1)
    coords[[i]]<-radius*((normIt(rawPts))[1:d])
  }
  out<-list()
  for(i in 1:d){
    out[[i]]<-unlist(lapply(coords,'[[',i))
  }
  out<-as.matrix(do.call(cbind, out))
  out<-sweep(out, MARGIN = 2, center, '+')
  out
}


sample_hyperspheroid <- function(d, n, radius, center, scaleFac = .15, fixEdge = FALSE, fixEdgeValue=100){
  # scaleFac is the max spread in each dimension from center. For example:
  # if c(2,3), scaleFac=.1 would go from 1.8 to 2.2 for 2 and from 2.7 to
  # 3.3 for 3.
  sphere = getHyperSphereFromOneCenter(d, n, radius, center)
  sphere = sweep(sphere, MARGIN = 2, center, '-')
  for (c in 1:ncol(sphere)) {
    cMax = max(abs(sphere[,c]))
    sphere[,c] = (sphere[,c]/cMax)*(center[c]*scaleFac)
  }
  sphere = sweep(sphere, MARGIN = 2, center, '+')
  if(fixEdge){
    sphere = cbind(rep(fixEdgeValue, nrow(sphere)), sphere[,2:ncol(sphere)])
  }
  sphere = rbind(sphere, center) #include original center
  return(sphere)
}

testIt = FALSE
if(testIt){
  points2D = sample_hyperspheroid(2, 1000, 5, c(30, 15))
  plot(points2D, asp = 1)
}
