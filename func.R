#### gFunc ####
gFunc = function(kesi, beta){
  theta = rbind(kesi, beta) #vec(r,1)
  yVec_pair=deltaMat %*% matrix(yVec,ncol=1)
  pair=function(i){
    z=yVec_pair[i,] - t(zMat[,i]) %*% theta
    if(z==0){
      return(epsilon)
    }
    return(z)
  }
  g = sapply(1:length(yVec_pair), 
             function(i){zMat[,i] * pair(i)})
  return(g) #mat(r,m)
}

gFuncFir = function(k){
  gFir = sapply(1:m, 
                function(i){- zMat[,i] * zMat[k,i]})
  return(gFir) #mat(r,m)
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
  kesi = params[[1]]; beta = params[[2]];
  lambda = params[[3]]
    
  fir = sum(sapply(1 + t(lambda) %*% gFunc(kesi,beta), logFunc))
  sec = m * sum(sapply(abs(kesi), penalFunc, yitaVec[1]))
  thr = m * sum(sapply(abs(lambda)[1:m], penalFunc,yitaVec[2])) #只对lambda的前m个分量做惩罚
  fth = m * sum(DMat %*% kesi)
  return(fir + sec - thr + fth)
}