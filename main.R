library(MASS); library(data.table)
source("init.R");source("func.R"); source("Algorithm.R")
initParams = initParams.function()

temp = apply(yitaMat, 1, computeBIC.function, initParams)
BICSeq = temp[1,]; paramSeq = temp[2,]
params = paramSeq[which.min(BICSeq)]


yitaVec = yitaMat[1,]
R = 1 ; k = 0
oldParams = initParams
kesiList = oldParams[[1]]

while(R >= epsilon){
  k = k + 1; print(k)
  newParams = iterate.function(oldParams,yitaVec)
  R = computeR.function(oldParams,newParams)
  oldParams = newParams
  kesiList = cbind(kesiList,oldParams[[1]])
}