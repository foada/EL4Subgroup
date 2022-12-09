#### gFunc ####
gFunc = function(mu, beta){
  theta = rbind(mu, beta) #vec(r,1)
  g = sapply(1:length(yVec), 
             function(i){zMat[,i] * (yVec[i] - t(zMat[,i]) %*% theta)})
  return(g) #mat(r,n)
}

gFuncFir = function(k){
  gFir = sapply(1:length(yVec), 
                function(i){- zMat[,i] * zMat[k,i]})
  return(gFir) #mat(r,n)
}


#### log*(z) ####
logFunc = function(z){
  if(z >= c) return(log(z))
  return(
    log(c) - 1.5 + 2 * z / c - z^2 / (2*c^2)
  )
}
logFuncFir = function(z){
  if(z >= c) return(1/z)
  return( 2 / c - z / (c^2) )
}
logFuncSec = function(z){
  if(z >= c) return( - 1/(z^2) )
  return( - 1 / (c^2) )
}

#### penalFunc ####
penalFunc = function(x,yita){
  gamma = 3.7
  f=function(a){min(1, max((gamma - a/yita),0) / (gamma-1))} 
  return(
    yita * integrate(
      Vectorize(f),
      0, x )$value)
}

penalFuncFir = function(x,yita){
  b = 3.7
  return(yita * (
    (x <= yita) + 
      max((b*yita -x),0) / ((b-1)*yita) * (x > yita)
  ))
}

penalFuncSec = function(x0,yita){
  return( penalFuncFir(abs(x0), yita) / abs(x0) )
}

penalApprFunc = function(x0,x,yita){
  return(
    penalFunc(x0,yita) + 
      1/2 * penalFuncFir(x0,yita)/x0 * (x^2 - x0^2)
  )
}

#### lFunc ####
lFunc = function(params,yitaVec){
  mu = params[[1]]; beta = params[[2]];
  lambda = params[[3]]; kesi = params[[4]]
    
  fir = sum(sapply(1 + t(lambda) %*% gFunc(mu,beta), logFunc))
  
  sec = n * sum(sapply(abs(kesi), penalFunc, yitaVec[1]))
  thr = n * sum(sapply(abs(lambda), penalFunc,yitaVec[2]))
  
  return(fir + sec - thr)
}