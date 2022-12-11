setwd("E:/MyFiles/Programming/R/EL4Subgroup/ours/")
library(snowfall)
sfInit(parallel=TRUE, cpus=4); sfLibrary(snowfall)
sfSource("init.R");sfSource("func.R"); sfSource("Algorithm.R")
# source("init_1207.R");source("func_1207.R"); source("Algorithm_1207.R")
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

#while(R >= epsilon & l<100){
#  newParams = iterate.function(oldParams,yitaVec)
#  R = computeR.function(oldParams,newParams)
#  oldParams = newParams
#  l = l + 1; print(R)
#  kesiList = cbind(kesiList,oldParams[[1]])
#  betaList = cbind(betaList,oldParams[[2]])
#  lambdaList = cbind(lambdaList,oldParams[[3]])
#}

start_t=proc.time()
newParams = iterate.function(oldParams,yitaVec)
end_t=proc.time()
print(end_t-start_t)

#sfStop()

kesi1=newParams[[1]]
beta1=newParams[[2]]
lambda1=newParams[[3]]
kesiList=newParams[[4]]
betaList=newParams[[5]]
lambdaList=newParams[[6]]
r1=DMat%*%kesi1
max(abs(r1))