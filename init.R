#---------------------------Init Function-------------------------------
initXY.function = function(){
  # 模拟自变量X, 得到数值矩阵xMat
  set.seed(123)
  sd = 1 ; pho = 0.5
  V = (sd ^ 2) * sapply(1:p, function(x){pho ^ abs(x - c(1:p))})
  xMat = t( mvrnorm(n = n, mu = rep(0, 10), V) )
  
  # 生成响应值yVec
  beta = runif(p, min = 0.5, max = 1) #均匀分布生成回归系数beta(真实值)
  mu = -2+4*rbinom(n, 1, 0.5) #生成均值向量(真实值)
  err = rnorm(n, 0, 0.5) #生成随机误差项
  yVec = mu + t(xMat) %*% beta + err
  
  return(list(xMat,yVec,beta,mu,err))
}

initParams.function = function(){
  # 初始化 μ 及 β
  initBeta = solve( xMat %*% t(xMat) ) %*% xMat %*% (yVec-mu)
  initMu = yVec - t(xMat) %*% initBeta
  
  # 初始化λ
  initLambda = matrix(0, nrow=r, ncol=1)
  mFir=function(lambda){
    g = sapply(1:length(err), 
               function(i){zMat[,i] * err[i]})
    logstarFir=sapply(1 + t(lambda) %*% g,logFuncFir)
    lambdaFir=-rowSums(sapply(1:n,function(i){logstarFir[i]*g[,i]}))
    return(matrix(lambdaFir,ncol=1))
  }
  mSec=function(lambda){
    g = sapply(1:length(err), 
               function(i){zMat[,i] * err[i]})
    logstarSec=sapply(1 + t(lambda) %*% g,logFuncSec)
    m=array(0,dim=c(r,r,n))
    for(i in 1:n){
      m[,,i]=logstarSec[i]*g[,i]%*%t(g[,i])
    }
    return(-apply(m,1:2,sum))
  }
  #Newton method
  newton_metod=function(oldLambda){
    r=1
    k=0
    while (r>epsilon) {
      newLambda=oldLambda-ginv(mSec(oldLambda))%*%mFir(oldLambda)
      r= sqrt(sum(mFir(newLambda)^2))
      oldLambda=newLambda
      k=k+1
    }
    return(list(oldLambda,k))
  }
  jieguo=newton_metod(initLambda)
  initLambda=jieguo[[1]]
  k=jieguo[[2]]
  
  initKesi=deltaMat %*% initMu
  initNu = matrix(0, nrow=s, ncol=1)
  
  
  initParams = list(initMu, initBeta, initLambda, initKesi, initNu)
  return(initParams)
}

#--------------------------- Initialization -------------------------------
## global variables
p = 10 ; n = 100 ; r = p + n ; q = r ; s = n*(n-1)/2
m = 0 ; epsilon = 1e-04
v = 1 ; a = 0.01 ; b = 3.7; c = 1/n
D1 = seq(1, 3, 0.5); D2 = seq(1, 3, 0.5)
d1 = length(D1) ; d2 = length(D2)
yitaMat = matrix(c(rep(D1,each=d2),rep(D2,d1)),ncol=2)

temp = initXY.function()
xMat = temp[[1]] # p * n ; # xi = t(x[, i])
yVec = temp[[2]] # n * 1 ; # yi = y[i]
beta=temp[[3]] #beta真实值
mu=temp[[4]]  #mu真实值
err=temp[[5]]
eMat = diag(n) # n维单位阵
zMat = rbind(eMat, xMat) # r * n
deltaMat = matrix(0, nrow=s, ncol=n) # s * n
for(i in seq(1, n-1, 1)){
  for(j in seq(i+1, n, 1)){
    ij=(i-1)*n - i*(i-1)/2 + (j-i)
    deltaMat[ij, i]=1
    deltaMat[ij, j]=-1
  }
}

## init Variables
# 格式
# initMu = matrix(0, nrow=n, ncol=1)
# initBeta = matrix(0, nrow=p, ncol=1)
# initLambda = matrix(0, nrow=r, ncol=1)
# initKesi = matrix(0, nrow=s, ncol=1);
# initNu = matrix(0, nrow=s, ncol=1) ; #ij = (i-1)*n - i*(i-1)/2 + (j-i)
# initTheta = rbind(initMu, initBeta) # matrix(r, 1) ; r = n+p
