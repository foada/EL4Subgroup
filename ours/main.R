#setwd("E:/360MoveData/Users/lx/Documents/My study/研究/基于经验似然的子群分析/代码")
source("init_1204_2.R");source("func_1204_2.R"); source("Algorithm_1204_2.R")
initParams = initParams.function()

#temp = apply(yitaMat, 1, computeBIC.function, initParams)
#BICSeq = temp[1,]; paramSeq = temp[2,]
#params = paramSeq[which.min(BICSeq)]


yitaVec = yitaMat[1,]
#R = 1 ; l = 0
oldParams = initParams
#kesiList = oldParams[[1]]
#betaList=oldParams[[2]]
#lambdaList=oldParams[[3]]

while(R >= epsilon & l<100){
  newParams = iterate.function(oldParams,yitaVec)
  R = computeR.function(oldParams,newParams)
  oldParams = newParams
  l = l + 1; print(R)
  kesiList = cbind(kesiList,oldParams[[1]])
  betaList = cbind(betaList,oldParams[[2]])
  lambdaList = cbind(lambdaList,oldParams[[3]])
}

start_t=proc.time()
newParams = iterate.function(oldParams,yitaVec)
end_t=proc.time()
print(end_t-start_t)
