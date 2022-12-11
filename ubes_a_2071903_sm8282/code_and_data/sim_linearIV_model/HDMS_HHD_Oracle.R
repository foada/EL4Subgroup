###################################################
##        oracle estimates in Table 2            ##
###################################################
library(MASS)
library(parallel)

##########data#######################################################

data_gen <- function(n, dw, dz, beta, dw2) #data6,7,8
{
  library(MASS)
  #dz <- 3
  betax<-beta[1]
  betaz<-beta[2:(dz+1)]
  
  z1<-rep(1,n) #intercept
  muz<-rep(0,dz-1)
  Sigmaz<-diag(dz-1)
  z2<-mvrnorm(n,muz,Sigmaz)
  z<-cbind(z1,z2)
  
  mue<-c(0,0)
  Sigma<-matrix(c(1,0.5,0.5,1),2,2) #控制内生性 0.5 0.6
  error<-mvrnorm(n,mue,Sigma)
  epsilon<-error[,1]
  mu<-error[,2]
  
  #dw2<-round(n*0.2) # 0.1 0.2
  #dw2<-round(log(n)) #0.1n 0.2n n^{1/3} n^{1/2} log(n)
  dw1<-dw-dw2
  
  #dw2<-3
  #dw1<-dw-dw2
  
  # dw1 valid iv
  muiv<-rep(0,dw1)
  Sigmaiv<-diag(dw1)
  iv<-mvrnorm(n,muiv,Sigmaiv) 
  #gammaiv<-c(0.8,rep(0.5,dw1-1))
  gammaiv<-c(0.8,seq(0.4,0.1,length.out = dw1-1))
  #gammaz<-rep(0.1,dz)
  gammaz<-rep(0.8,dz)
  x<-iv%*%gammaiv+z%*%gammaz+mu
  
  y<-x%*%betax+z%*%betaz+epsilon
  
  # dw2 invalid iv
  # delta.min <- 0.5
  # delta.max <- 0.7 
  # delta <- seq(delta.min,delta.max,len=dw2)
  # muiiv <- rep(0,dw2)
  # Sigmaiiv <- diag(dw2)
  # nu <- mvrnorm(n,muiiv,Sigmaiiv)
  # iiv <- matrix(rep(epsilon,dw2),nrow = n)%*%diag(delta) + nu
  
  #w<-cbind(iv,iiv)
  w <- iv
  
  return(list(y=y,x=x,z=z,w=w))
}

######init.B##########################################

init.B.2SLS.1 <- function(y, x, z, w)
{
  zw<-cbind(z,w)
  res<-lm(x~zw[,-1])
  x.hat<-res$fitted.values
  init.beta<-lm(y~cbind(x.hat, z[,-1]))$coef
  init.beta <- init.beta[c(2,1,3,4)]
  
  return(init.beta)
}

init.B.2SLS <- function(y, x, z, w)
{
  n <- length(y)
  w <- w[, 1:(n*0.1)]
  zw<-cbind(z,w)
  res<-lm(x~zw[,-1])
  x.hat<-res$fitted.values
  init.beta<-lm(y~cbind(x.hat, z[,-1]))$coef
  init.beta <- init.beta[c(2,1,3,4)]
  
  return(init.beta)
}

########rep.fun#########################

rep.fun <- function(index, para) #n init.B beta
{
  
  library(MASS)
  
  dz<-para$dz
  dw<-para$dw
  n <- para$n  # sample size
  beta <- para$beta
  dw2 <- para$dw2
  
  data <- data_gen(n, dw, dz, beta, dw2)
  y <- data$y
  x <- data$x
  z <- data$z
  w <- data$w
  
  init.beta <- init.B.2SLS(y,x,z,w)
  
  return(init.beta)
}

############demo###############################

start.t <- proc.time()

para <- list(n = n, dw = dw, dz = dz, beta = beta, dw2=dw2)

rep.seq <- c(1:rep.num)

#core.num<- detectCores(logical = FALSE)
core.num <- 4

cl <- makeCluster(core.num)
clusterSetRNGStream(cl, iseed = NULL)

clusterExport(cl, c('init.B.2SLS', 'data_gen'))
HDMS.esti <- parSapplyLB(cl, rep.seq, rep.fun, para = para, simplify = TRUE)

stopCluster(cl)

write.csv(HDMS.esti, file = paste("2SLS_esti_","n=", n, "dw=", dw, "s=", dw2, ".csv", sep= ""))

betaxz.mat <- t(HDMS.esti[-2, ])

RMSE <- sqrt(colMeans((betaxz.mat-beta[1])^2))
bias <- colMeans(betaxz.mat-beta[1])
std <- sqrt(RMSE^2 - bias^2)

RMSE2 <- sqrt(sum(RMSE[-1]^2))

error <- t(rbind(RMSE, bias, std, RMSE2))

print(error)

write.csv(error, file = paste("2SLS_error_","n=",n,"dw=",dw, "s=", dw2, ".csv", sep= ""))

######time#######
end.t <- proc.time()
op.time <- end.t - start.t
print(op.time)







