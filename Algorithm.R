#--------------------------- Iteration -------------------------------
computeBIC.function = function(yitaVec, initParams){
  # r = computeR.function(initParams)
  r = 1
  j=0
  oldParams = initParams
  
  while(r >= epsilon & j<10000){
    #使用j的原因是防止迭代过多次仍然不能收敛，用j控制一下让能跑出来结果
    newParams = iterate.function(oldParams,yitaVec)
    r = computeR.function(newParams)
    oldParams = newParams
    j=j+1
  }
  
  newMu = newParams[[1]]
  k = length(unique(newMu))
  BIC = 2 * lFunc(newParams, yitaVec) + log(n) * (k + p)
  
  return(list(BIC,newParams))
} 


computeR.function = function(params){
  mu = params[[1]]; kesi = params[[4]]
  r = sqrt( sum( (deltaMat %*% mu - kesi) ^ 2 ) )
  return(r)
}


iterate.function = function(oldParams,yitaVec){
  oldMu = oldParams[[1]]; oldBeta = oldParams[[2]]
  oldLambda = oldParams[[3]]; oldKesi = oldParams[[4]]
  oldNu = oldParams[[5]]
  
  newMu = iterateMu.function(oldMu, oldBeta, oldLambda, oldKesi, oldNu)
  for(k in 1:p){
    temp=iterateBeta.function(k, oldBeta, newMu, oldLambda, yitaVec[2])
    #update beta_k
    oldBeta=temp[[1]]
    #the update of lambda while updating beta_k
    oldLambda=temp[[2]]
  }
  newBeta = oldBeta
  newLambda = oldLambda
  newDelta = iterateDelta.function(newMu,oldNu)
  newKesi = iterateKesi.function(newDelta,yitaVec[1])
  newNu = iterateNu.function(oldNu,newMu,newKesi)
  
  newParams = list(newMu, newBeta, newLambda, newKesi, newNu)
  return(newParams)
}


iterateMu.function = function(oldMu, oldBeta, oldLambda, oldKesi, oldNu){
  gMat = gFunc(oldMu, oldBeta)
  myfunc = function(gi){logFuncFir(1 + t(oldLambda) %*% gi)}
  logFirSeq = apply(gMat,2,myfunc)
  
  fir = a/n * rowSums(sapply(
    1:n, function(i){logFirSeq[i] * eMat[,i] %*% t(zMat[,i]) %*% oldLambda}))
  sec = a * v * t(deltaMat) %*% (deltaMat %*% oldMu - oldKesi + 1/v * oldNu)
  newMu = oldMu + fir - sec
  return(newMu)
}


iterateBeta.function = function(k, oldBeta, newMu, oldLambda, yita){
  # print(paste("k:",k))
  oldBeta_k = oldBeta[k]
  simSeq = 1 + t(oldLambda) %*% gFunc(newMu, oldBeta)
  wikmSeq = t(oldLambda) %*% gFuncFir(n+k)
  
  logFirSeq = sapply(simSeq,logFuncFir)
  logSecSeq = sapply(simSeq,logFuncSec)
  
  numerator = sum(logFirSeq * wikmSeq)
  denominator = sum(logSecSeq * wikmSeq^2)
  
  
  newBeta_k = oldBeta_k - numerator/denominator
  newBeta = oldBeta
  newBeta[k] = newBeta_k
  #对beta某个分量的更新
  newLambda = sapply(1:r, iterateLambda.function, oldLambda, newMu, newBeta, yita)
  
  return(list(newBeta,newLambda))
}


iterateLambda.function = function(s, oldLambda, newMu, oldBeta, yita){
  # print(paste("s:",s))
  oldLambda_s = oldLambda[s]
  timSeq = 1 + t(oldLambda) %*% gFunc(newMu,oldBeta)
  gisThetaSeq = gFunc(newMu,oldBeta)[s,]
  
  logFirSeq = sapply(timSeq,logFuncFir)
  logSecSeq = sapply(timSeq,logFuncSec)
  
  numerator = sum(logFirSeq * gisThetaSeq) - 
    n*penalFuncFir(abs(oldLambda_s),yita)*ifelse(oldLambda_s>0,1,-1)
  denominator = sum(logSecSeq * gisThetaSeq^2) - 
    n*penalFuncSec(abs(oldLambda_s),yita)
  newLambda_s = oldLambda_s - numerator/denominator
  
  # newLambda = oldLambda
  # newLambda[s] = newLambda_s
  return(newLambda_s)
}

iterateDelta.function = function(newMu,oldNu){
  
  newDelta = deltaMat %*% newMu + 1/v * oldNu
  return(newDelta)
}

iterateNu.function = function(oldNu,newMu,newKesi){
  newNu = oldNu + v * (deltaMat %*% newMu - newKesi)
  return(newNu)
}


iterateKesi.function=function(newDelta,yita){
  ST=function(t,lambda){
    return(sign(t)*ifelse(abs(t)-lambda>0,abs(t)-lambda,0))
  }
  
  update.function=function(delta,yita){
    if(abs(delta)<=(yita+yita/v)){
      return(ST(delta,yita/v))
    }else if(abs(delta)<=b*yita){
      return(ST(delta,b*yita/((b-1)*v))/(1-1/((b-1)*v)))
    }else{
      return(delta)
    }
  }
  
  newKesi = sapply(newDelta, update.function, yita)
  return(matrix(newKesi,ncol=1))
}