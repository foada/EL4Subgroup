#### gFunc ####
gFunc = function(beta){
  g = sapply(1:length(yVec), 
             function(i){zMat[,i] * (yVec[i] - t(zMat[,i]) %*% beta)})
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
#计算中已经带有绝对值
penalFunc = function(x,yita){
  ## where yita is tuning parameter
  a <- 3.7 ## a = 3.7 suggested by Fan and Li(2001)
  abs.x <- abs(x)
  temp.1 <- yita * abs.x  ## for the case abs.beta <= yita
  temp.2 <- (2*a*yita*abs.x - abs.x^2 - yita^2)/(2*(a-1.0)) ## for the case abs.x <= a*yita
  temp.3 <- yita^2 * (a+1.0)/2
  
  scad.val <- temp.1*(abs.x <= yita) + temp.2*(abs.x > yita)*(abs.x <= a * yita) + temp.3*(abs.x > a * yita)
  return(scad.val)
}

#已经是带绝对值的函数the abs has already exists
penalFuncFir <- function(x,yita){
  ## where lambda is tuning parameter
  a <- 3.7 ## a = 3.7 suggested by Fan and Li(2001)
  abs.x <- abs(x)
  
  #temp.1 <- lambda  ## for the case abs.beta <= lambda
  temp.1 <- yita*sign(x) #改
  #temp.2 <- (a*lambda - abs.beta)/(a-1.0) ## for the case abs.beta <= a*lambda
  temp.2 <- max( a*yita - abs.x,0)/(a-1.0)*sign(x) #改
  
  grad.scad <- temp.1*(abs.x <= yita) + temp.2*(abs.x > yita)*(abs.x <= a * yita)
  
  return(grad.scad)
}


penalFuncSec = function(x,yita){
  thresh <- yita/4 #?
  
  init.x <- x
  index.1 <- (abs(init.x) < 1e-4)
  index.2 <- (abs(init.x) < thresh)
  index.3 <- as.logical((1-index.1)*index.2)
  init.x[index.1] <- thresh #?
  init.x[index.3] <- sign(x[index.3])*thresh
  
  temp <- penalFuncFir(abs(init.x),yita)
  
  return(temp/abs(init.x))
}

penalApprFunc = function(x0,x,yita){
  return(
    penalFunc(x0,yita) + 
      1/2 * penalFuncFir(x0,yita)/x0 * (x^2 - x0^2)
  )
}


#### lFunc ####
lFunc = function(params,yitaVec){
   beta = params[[1]]
  lambda = params[[2]] 
  
  fir = sum(sapply(1 + t(lambda) %*% gFunc(beta), logFunc))
  sec = n * sum(sapply(beta, penalFunc, yitaVec[1]))
  thr = n * sum(sapply(lambda, penalFunc,yitaVec[2])) 
  return(fir + sec - thr)
}