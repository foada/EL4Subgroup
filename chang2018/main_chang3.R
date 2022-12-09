source("init_chang3.R");source("func_chang2.R"); source("Algorithm_chang3.R")
initParams = initParams.function()

yitaVec = yitaMat[1,]
R = 1 ; l = 0
oldParams = initParams
betaList=oldParams[[1]]
lambdaList=oldParams[[2]]

while(R >= epsilon & l<300){#先尝试迭代100次
  newParams = iterate.function(oldParams,yitaVec)
  R = computeR.function(oldParams,newParams)
  oldParams = newParams
  l = l + 1; 
  #print(R,digits = 5)
  betaList = cbind(betaList,oldParams[[1]])
  lambdaList = cbind(lambdaList,oldParams[[2]])
}

while(R >= epsilon & l<300){#先尝试迭代100次
  newParams = iterate.function(oldParams,yitaVec)
  R = computeR.function(oldParams,newParams)
  oldParams = newParams
  l = l + 1; 
  #print(R,digits = 5)
  betaList = cbind(betaList,oldParams[[1]])
  lambdaList = cbind(lambdaList,oldParams[[2]])
}

start.t <- proc.time()
newParams = iterate.function(oldParams,yitaVec)
end.t <- proc.time()
op.time <- end.t - start.t
