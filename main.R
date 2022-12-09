library(MASS)
setwd("E:/360MoveData/Users/lx/Documents/My study/研究/基于经验似然的子群分析/代码")
source("init_1127.R");source("func_1127.R"); source("Algorithm_1127.R")
#load('initLambda.RData')
initParams = initParams.function()
#initParams[[3]]=matrix(initLambda,ncol=1)
initParams[[3]]=initParams[[3]]/10000000

temp = apply(yitaMat, 1, computeBIC.function, initParams)
BICSeq = temp[1,]; paramSeq = temp[2,]
params = paramSeq[which.min(BICSeq)]


yitaVec = yitaMat[1,]
R = 1 ; l = 0
oldParams = initParams
kesiList = oldParams[[1]]
betaList=oldParams[[2]]
lambdaList=oldParams[[3]]

while(R >= epsilon){
  newParams = iterate.function(oldParams,yitaVec)
  R = computeR.function(oldParams,newParams)
  oldParams = newParams
  l = l + 1; print(l)
  kesiList = cbind(kesiList,oldParams[[1]])
  betaList = cbind(betaList,oldParams[[2]])
  lambdaList = cbind(lambdaList,oldParams[[3]])
}
