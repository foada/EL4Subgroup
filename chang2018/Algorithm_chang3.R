#--------------------------- Iteration -------------------------------
computeBIC.function = function(yitaVec, initParams){
  # R = computeR.function(initParams)
  R = 1
  j=0
  oldParams = initParams
  
  while(R >= epsilon & j<500){
    #使用j的原因是防止迭代过多次仍然不能收敛，用j控制一下让能跑出来结果
    newParams = iterate.function(oldParams,yitaVec)
    R = computeR.function(oldParams,newParams)
    oldParams = newParams
    j=j+1
    print(j)
  }
  
  k = sum(abs(newParams[[1]])>epsilon) # beta向量的非零个数
  BIC = 2 * lFunc(newParams, yitaVec) + log(n) * (k + p)
  
  return(list(BIC,newParams))
} 


computeR.function = function(oldParams,newParams){
  oldBeta = oldParams[[1]]; newBeta = newParams[[1]]
  R = max(abs(oldBeta-newBeta))
  return(R)
}


iterate.function = function(oldParams,yitaVec){
  oldBeta = oldParams[[1]]
  oldLambda = oldParams[[2]]
  
  eps.tol <- 0.005 
  iter.max <- 30
  # 迭代beta及lambda
  for(count.num in 1:iter.max){
    temp=optim.beta(oldBeta,oldLambda,yitaVec)
    #update beta_k
    oldBeta=temp[[1]]
    #the update of lambda while updating beta_k
    oldLambda=temp[[2]]
    if(count.num>5){
      oldBeta[abs(oldBeta)<eps.tol]=0
    }
  }
  newBeta = oldBeta
  newLambda = oldLambda
  
  newParams = list(newBeta, newLambda)
  return(newParams)
}


iterateBeta.function = function(k, oldBeta, oldLambda, yitaVec){
  w <- 1e-10
  oldBeta_k = oldBeta[k]
  sijSeq = 1 + t(oldLambda) %*% gFunc(oldBeta)
  wijkSeq = t(oldLambda) %*% gFuncFir(k)
  
  logFirSeq = sapply(sijSeq,logFuncFir)
  logSecSeq = sapply(sijSeq,logFuncSec)
  
  numerator = sum(logFirSeq * wijkSeq) + n*penalFuncFir(oldBeta_k,yitaVec[1])
  denominator = sum(logSecSeq * wijkSeq^2)+n*penalFuncSec(oldBeta_k,yitaVec[1])
  
  #if (abs(denominator) <= w){ #denominator can not be zero
  #  denominator <- denominator + w
  #}
  #temp.val=numerator/denominator
  temp.val=numerator
  
  if (abs(temp.val) > 0.01)  temp.val <- sign(temp.val) * 0.01
  newBeta_k = oldBeta_k - temp.val
  
  #if(abs(newBeta_k)<1e-03){
  #  newBeta_k=0
  #}
  #newBeta = oldBeta
  #newBeta[k] = newBeta_k

  #newLambda = sapply(1:r, iterateLambda.function,
  #                   oldLambda, newBeta, yitaVec[2])
  #newLambda[abs(newLambda)<1e-03]=0
  
  return(newBeta_k)
}

optim.beta <- function(Beta,Lambda,yitaVec)  #iter.max
{
  iter.max <- 30  #convergence????? 30 60 90
  for (iter.num in c(1:iter.max))
  {
    Beta = sapply(1:r, iterateBeta.function,Beta,Lambda,yitaVec[1])
    Lambda=optim.lambda(Beta,Lambda,yitaVec[2])
    Lambda[abs(Lambda)<=1e-03]=0
  }
  newParams = list(Beta,Lambda)
  return(newParams)
}

iterateLambda.function = function(s,oldBeta,oldLambda,yita){
  # print(paste("s:",s))
  w <- 1e-10
  oldLambda_s = oldLambda[s]
  tijSeq = 1 + t(oldLambda) %*% gFunc(oldBeta)
  gisThetaSeq = gFunc(oldBeta)[s,]
  
  logFirSeq = sapply(tijSeq,logFuncFir)
  logSecSeq = sapply(tijSeq,logFuncSec)
  
  numerator = sum(logFirSeq * gisThetaSeq) - n*penalFuncFir(oldLambda_s,yita)
  denominator = sum(logSecSeq * gisThetaSeq^2)-n*penalFuncSec(oldLambda_s,yita)
  
  #if (abs(denominator) <= w){ #denominator can not be zero
  #  denominator <- denominator + w
  #}
  #temp.val=numerator/denominator
  temp.val=numerator
  
  if (abs(temp.val) > 0.01)  temp.val <- sign(temp.val) * 0.01
  newLambda_s = oldLambda_s + temp.val
  
  #if(abs(newLambda_s)<1e-03){
  #  newLambda_s=0
  #}
  
  return(newLambda_s)
}

optim.lambda <- function(Beta,Lambda,yita)  #iter.max
{
  iter.max <- 30  #convergence????? 30 60 90
  for (iter.num in c(1:iter.max))
  {
    Lambda = sapply(1:r,iterateLambda.function,Beta,Lambda,yita)
  }
  return(Lambda)
}