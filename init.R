#---------------------------Init Function-------------------------------
library(MASS)
initXY.function = function(){
  # 模拟自变量X, 得到数值矩阵xMat
  set.seed(123)
  sd = 1 ; pho = 0.5
  V = (sd ^ 2) * sapply(1:p, function(x){pho ^ abs(x - c(1:p))})
  xMat = t( mvrnorm(n = n, mu = rep(0, p), V) )
  
  # 生成响应值yVec
  beta = runif(p, min = 0.5, max = 1) #均匀分布生成回归系数beta(真实值)
  mu = rep(c(-2,2),each=n/2) #生成均值向量(真实值)
  err = rnorm(n, 0, 0.5) #生成随机误差项
  yVec = mu + t(xMat) %*% beta + err
  
  return(list(xMat,yVec,beta,mu,err))
}

initParams.function = function(){
  # 初始化 μ 及 β
  initBeta = solve( xMat %*% t(xMat) ) %*% xMat %*% (yVec-mu)
  initMu = yVec - t(xMat) %*% initBeta +rnorm(n, 0, 0.5)
  #初始化kesi和theta
  initKesi = deltaMat %*% initMu
  initTheta = rbind(initKesi,initBeta)
  
  # 初始化λ
  initLambda = matrix(0, nrow=r, ncol=1)
  mFir=function(lambda){
    g = gFunc(initKesi,initBeta)
    logstarFir=sapply(1 + t(lambda) %*% g,logFuncFir)
    lambdaFir=-rowSums(sapply(1:m,function(i){logstarFir[i]*g[,i]}))
    return(matrix(lambdaFir,ncol=1))
  }
  mSec=function(lambda){
    g = gFunc(initKesi,initBeta)
    logstarSec=sapply(1 + t(lambda) %*% g,logFuncSec)
    M=array(0,dim=c(r,r))
    for(i in 1:m){
      M=M+logstarSec[i]*g[,i]%*%t(g[,i])
    }
    gc() # 清除冗余内存
    return(-M)
  }
  #Newton method
  newton_metod=function(oldLambda){
    r=1
    #k=0
    while (r>epsilon) {
      newLambda=oldLambda-ginv(mSec(oldLambda))%*%mFir(oldLambda)
      r= sqrt(sum(mFir(newLambda)^2))
      oldLambda=newLambda
      #k=k+1
    }
    return(oldLambda)
  }
  initLambda=newton_metod(initLambda)
  #initLambda=jieguo[[1]]
  #f=jieguo[[2]]
  
  initParams = list(initKesi,initBeta,initLambda)
  return(initParams)
}

#--------------------------- Initialization -------------------------------
## global variables
p = 10 ; n = 50 ; m = n*(n-1)/2 ; r = p + m ; q = r ;
l = 0 ; epsilon = 1e-04
b = 3.7; c = 1/n
D1 = seq(1, 3, 0.5); D2 = seq(1, 3, 0.5)
d1 = length(D1) ; d2 = length(D2)
yitaMat = matrix(c(rep(D1,each=d2),rep(D2,d1)),ncol=2)

temp = initXY.function()
xMat = temp[[1]] # p * n ; # xi = t(x[, i])
yVec = temp[[2]] # n * 1 ; # yi = y[i]
beta=temp[[3]] #beta真实值
mu=temp[[4]]  #mu真实值
err=temp[[5]]

deltaMat = matrix(0, nrow=m, ncol=n) # m * n
for(i in seq(1, n-1, 1)){
  for(j in seq(i+1, n, 1)){
    ij=(i-1)*n - i*(i-1)/2 + (j-i)
    deltaMat[ij, i]=1
    deltaMat[ij, j]=-1
  }
}

DMat = data.table::data.table()
for(i in seq(1, n-2, 1)){
  DMat_i = matrix(0, nrow=(n-i)*(n-i-1)/2, ncol=m)
  for(j in seq(i+1, n-1, 1)){
    ij = (i-1)*n - i*(i-1)/2 + (j-i)
    for(k in seq(j+1, n, 1)){
      ik = (i-1)*n - i*(i-1)/2 + (k-i)
      jk = (j-1)*n - j*(j-1)/2 + (k-j)
      ijk = (j-1-i)*(n-i) - (j-i)*(j-1-i)/2 + (k-j)
      DMat_i[ijk,c(ij,ik,jk)] = c(-1, 1, -1)
    }
  }
  DMat = rbind(DMat,DMat_i)
}
DMat = as.matrix(DMat) 

eMat = diag(m) # m维单位阵
zMat = rbind(eMat, t(deltaMat %*% t(xMat))) # r * m
## init Variables
# 格式
# initMu = matrix(0, nrow=n, ncol=1)
# initBeta = matrix(0, nrow=p, ncol=1)
# initLambda = matrix(0, nrow=r, ncol=1)
# initKesi = matrix(0, nrow=s, ncol=1);
# initNu = matrix(0, nrow=s, ncol=1) ; #ij = (i-1)*n - i*(i-1)/2 + (j-i)
# initTheta = rbind(initMu, initBeta) # matrix(r, 1) ; r = n+p
