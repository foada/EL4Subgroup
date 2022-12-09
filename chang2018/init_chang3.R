#---------------------------Init Function-------------------------------
library(MASS)
initXY.function = function(){
  # 模拟自变量X, 得到数值矩阵zMat
  set.seed(123)
  sd = 1 ; pho = 0.5
  V = (sd ^ 2) * sapply(1:p, function(x){pho ^ abs(x - c(1:p))})
  zMat = t( mvrnorm(n = n, mu = rep(0, p), V) )
  
  # 生成响应值yVec
  beta = rep(0,times=p)
  beta[1]=3; beta[2]=1.5; beta[5]=2
  #beta = runif(p, min = 0.5, max = 1) #均匀分布生成回归系数beta(真实值)
  #mu = rep(c(-2,2),each=n/2) #生成均值向量(真实值)
  err = rnorm(n, 0, 1) #生成随机误差项
  yVec = t(zMat) %*% beta + err
  
  return(list(zMat,yVec,beta,err))
}

initParams.function = function(){
  # 初始化 μ 及 β
  initBeta = ginv( zMat %*% t(zMat) ) %*% zMat %*% yVec
  #initMu = yVec - t(zMat) %*% initBeta +rnorm(n, 0, 0.5)
  #初始化kesi和theta
  #initKesi = deltaMat %*% initMu
  #initTheta = rbind(initKesi,initBeta)
  
  # 初始化λ
  initLambda = matrix(0, nrow=r, ncol=1)
  #optim.lambda <- function(Lambda,Beta,yita)  #iter.max
  #{
  #  iter.max <- 30  #convergence????? 30 60 90
  #  for (iter.num in c(1:iter.max))
  #  {
  #    Lambda = sapply(1:r, iterateLambda.function,
  #                       Lambda, Beta, yita)
  #  }
  #  return(Lambda)
  #}
  initLambda=optim.lambda(initBeta,initLambda,yitaMat[1,2])
  initLambda=matrix(initLambda,ncol=1)
  
  initParams = list(initBeta,initLambda)
  return(initParams)
}

#--------------------------- Initialization -------------------------------
## global variables
p = 50 ; n = 25 ; r = p ; 
l = 0 ; epsilon = 1e-04
b = 3.7; c = 1/n
#penal tuning para
#eta_1 for kesi, eta_2 for lambda
eta_1.max <- 1.5/sqrt(n) # 30 10 5 3 1.5
eta_1.min <- 0.3/sqrt(n) # 15 5 1 0.5 0.3
d1 <- 5 
eta_1.seq = exp(seq(log(eta_1.min),log(eta_1.max),len=d1)) 

eta_2.max <- 20/sqrt(n) #10 20
eta_2.min <- 10/sqrt(n)#1 5
d2 <-  5 # 5 15
eta_2.seq =  exp(seq(log(eta_2.min),log(eta_2.max),len=d2))  
yitaMat = matrix(c(rep(eta_1.seq,each=d2),rep(eta_2.seq,d1)),ncol=2)

temp = initXY.function()
zMat = temp[[1]] # p * n ; # xi = t(x[, i])
yVec = temp[[2]] # n * 1 ; # yi = y[i]
beta=temp[[3]] #beta真实值
err=temp[[4]]

#eMat = diag(n) # n维单位阵
#zMat = rbind(eMat, zMat) # r * n
## init Variables
# 格式
# initMu = matrix(0, nrow=n, ncol=1)
# initBeta = matrix(0, nrow=p, ncol=1)
# initLambda = matrix(0, nrow=r, ncol=1)
# initKesi = matrix(0, nrow=s, ncol=1);
# initNu = matrix(0, nrow=s, ncol=1) ; #ij = (i-1)*n - i*(i-1)/2 + (j-i)
# initTheta = rbind(initMu, initBeta) # matrix(r, 1) ; r = n+p
