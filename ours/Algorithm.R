#--------------------------- Iteration -------------------------------
computeBIC.function = function(yitaVec, initParams){
  # R = computeR.function(initParams)
  R = 1
  j=0
  oldParams = initParams
  
  while(R >= epsilon & j<1000){
    #使用j的原因是防止迭代过多次仍然不能收敛，用j控制一下让能跑出来结果
    newParams = iterate.function(oldParams,yitaVec)
    R = computeR.function(oldParams,newParams)
    oldParams = newParams
    j=j+1
    print(j)
  }
  
  k = sum(abs(newParams[[1]])>epsilon) # kesi向量的非零个数
  BIC = 2 * lFunc(newParams, yitaVec) + log(n) * (k + p)
  
  return(list(BIC,newParams))
} 


computeR.function = function(oldParams,newParams){
  oldKesi = oldParams[[1]]; newKesi = newParams[[1]]
  R = max(abs(oldKesi-newKesi))
  return(R)
}


iterate.function = function(oldParams,yitaVec){
  oldKesi = oldParams[[1]]; oldBeta = oldParams[[2]]
  oldLambda = oldParams[[3]]
  
  eps.tol_1 <- 0.005 
  iter.max <- 10
  
  for(count.num in 1:iter.max){
    temp=optim.kesibeta(oldKesi,oldBeta,oldLambda,yitaVec)
    oldKesi=temp[[1]]
    oldBeta=temp[[2]]
    oldLambda=temp[[3]]
    if(count.num>2){
      oldKesi[abs(oldKesi)<eps.tol_1]=0
    }
  }
  newKesi=oldKesi
  newBeta = oldBeta
  newLambda = oldLambda
  
  newParams = list(newKesi, newBeta, newLambda)
  return(newParams)
}


iterateKesi.function = function(k, oldKesi, oldBeta, oldLambda, yita){
  w <- 1e-10
  oldKesi_k = oldKesi[k]
  sijSeq = 1 + t(oldLambda) %*% gFunc(oldKesi, oldBeta)
  uijkSeq = t(oldLambda) %*% gFuncFir(k)
  
  logFirSeq = sapply(sijSeq,logFuncFir)
  logSecSeq = sapply(sijSeq,logFuncSec)
  
  numerator = sum(logFirSeq * uijkSeq) + 
    m*penalFuncFir(oldKesi_k,yita) + m*sum(DMat[,k])
  denominator = sum(logSecSeq * uijkSeq^2) +m*penalFuncSec(oldKesi_k,yita)
  
  if (abs(denominator) <= w){ #denominator can not be zero
    denominator <- denominator + w
  }
  temp.val=numerator/denominator
  if (abs(temp.val) > 0.001)  temp.val <- sign(temp.val) * 0.001
  
  newKesi_k = oldKesi_k - temp.val

  #newLambda = sapply(1:r, iterateLambda.function,
  #                   oldLambda, newKesi, oldBeta, yitaVec[2]) 
  
  return(newKesi_k)
}


iterateBeta.function = function(k, oldKesi, oldBeta, oldLambda){
  # print(paste("k:",k))
  w <- 1e-10
  oldBeta_k = oldBeta[k]
  sijSeq = 1 + t(oldLambda) %*% gFunc(oldKesi, oldBeta)
  wijkSeq = t(oldLambda) %*% gFuncFir(m+k)
  
  logFirSeq = sapply(sijSeq,logFuncFir)
  logSecSeq = sapply(sijSeq,logFuncSec)
  
  numerator = sum(logFirSeq * wijkSeq)
  denominator = sum(logSecSeq * wijkSeq^2)
  
  if (abs(denominator) <= w){ #denominator can not be zero
    denominator <- denominator + w
  }
  temp.val=numerator/denominator
  if (abs(temp.val) > 0.001)  temp.val <- sign(temp.val) * 0.001
  newBeta_k = oldBeta_k - temp.val

  #newLambda = sapply(1:r, iterateLambda.function,
  #                   oldLambda, oldKesi, newBeta, yita)
  
  return(newBeta_k)
}

optim.kesibeta <- function(Kesi,Beta,Lambda,yitaVec)  #iter.max
{
  iter.max <- 30  #convergence????? 30 60 90
  eps.tol_2 <- 0.005 
  num_solve=(n-1)*(n-2)/2
  for (iter.num in c(1:iter.max))
  {
    Kesi = sapply(1:m, iterateKesi.function,Kesi,Beta,Lambda,yitaVec[1])
    Kesi=matrix(Kesi,ncol=1)
    
    KesiOrigin=Kesi[n:m] #迭代得到的
    KesiSolve=DMat[1:num_solve,1:n-1]%*%Kesi[1:n-1] #根据通解计算得到
    supp_new=c(which(abs(KesiSolve-KesiOrigin)>0.01))
    if(length(supp_new)>0) Kesi[supp_new+n-1]=KesiSolve[supp_new]
    
    Beta = sapply(1:p, iterateBeta.function,Kesi,Beta,Lambda)
    Beta=matrix(Beta,ncol=1)
    Lambda=optim.lambda(Kesi,Beta,Lambda,yitaVec[2])
    Lambda=matrix(Lambda,ncol=1)
    #Lambda[abs(Lambda)<=eps.tol_2]=0
  }
  newParams = list(Kesi,Beta,Lambda)
  return(newParams)
}


iterateLambda.function = function(s, oldKesi, oldBeta, oldLambda, yita){
  # print(paste("s:",s))
  w <- 1e-10
  oldLambda_s = oldLambda[s]
  tijSeq = 1 + t(oldLambda) %*% gFunc(oldKesi, oldBeta)
  gisThetaSeq = gFunc(oldKesi, oldBeta)[s,]
  
  logFirSeq = sapply(tijSeq,logFuncFir)
  logSecSeq = sapply(tijSeq,logFuncSec)
  
  if(s<=m){
    numerator = sum(logFirSeq * gisThetaSeq)-m*penalFuncFir(oldLambda_s,yita)
    denominator = sum(logSecSeq * gisThetaSeq^2)-m*penalFuncSec(oldLambda_s,yita)
  }
  else{
    numerator = sum(logFirSeq * gisThetaSeq)
    denominator = sum(logSecSeq * gisThetaSeq^2)
  }
  
  if (abs(denominator) <= w){ #denominator can not be zero
    denominator <- denominator + w
  }
  temp.val=numerator/denominator
  if (abs(temp.val) > 0.001)  temp.val <- sign(temp.val) * 0.001
  newLambda_s = oldLambda_s + temp.val

  return(newLambda_s)
}

optim.lambda <- function(Kesi,Beta,Lambda,yita)  #iter.max
{
  iter.max <- 30  #convergence????? 30 60 90
  for (iter.num in c(1:iter.max))
  {
    Lambda = sapply(1:r,iterateLambda.function,Kesi,Beta,Lambda,yita)
    Lambda=matrix(Lambda,ncol=1)
  }
  return(Lambda)
}


iterateLambda2.function = function(s, oldKesi, oldBeta, oldLambda, yita){
  # print(paste("s:",s))
  w <- 1e-10
  oldLambda_s = oldLambda[s]
  tijSeq = 1 + t(oldLambda) %*% gFunc(oldKesi, oldBeta)
  gisThetaSeq = gFunc(oldKesi, oldBeta)[s,]
  
  logFirSeq = sapply(tijSeq,logFuncFir)
  logSecSeq = sapply(tijSeq,logFuncSec)
  
  if(s<=m){
    numerator = sum(logFirSeq * gisThetaSeq)-m*penalFuncFir(oldLambda_s,yita)
    denominator = sum(logSecSeq * gisThetaSeq^2)-m*penalFuncSec(oldLambda_s,yita)
  }
  else{
    numerator = sum(logFirSeq * gisThetaSeq)
    denominator = sum(logSecSeq * gisThetaSeq^2)
  }
  
  if (abs(denominator) <= w){ #denominator can not be zero
    denominator <- denominator + w
  }
  temp.val=numerator/denominator
  if (abs(temp.val) > 0.001)  temp.val <- sign(temp.val) * 0.001
  newLambda_s = oldLambda_s + temp.val
  
  return(newLambda_s)
}
