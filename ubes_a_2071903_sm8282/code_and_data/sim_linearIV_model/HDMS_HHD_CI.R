##############################################
### use PPEL to compute confidence regions ###
##############################################
library(MASS)
library(parallel)

#setwd("/home/jiazhang/HDMS_ZJ")  #
#setwd("C:/Users/jean/Desktop/res0419/data4")
source("HDMS_HHD_FUN.R") # HHD or HD

######init.B##########################################

init.B.2SLS <- function(y, x, z, w)
{
  zw1<-cbind(z,w[,1])
  res<-lm(x~zw1[,-1])
  x.hat<-res$fitted.values
  init.beta<-lm(y~cbind(x.hat,z[,-1]))$coef
  init.beta <- init.beta[c(2,1,3,4)]
  
  residual<-y-cbind(x,z)%*%init.beta
  residual<-as.vector(residual)
  residual<-diag(residual)
  
  #index <- sample(1:n, round(n*0.8), replace = FALSE)
  #init.xi<-colMeans((residual%*%w[,-1])[index,])
  
  init.xi<-colMeans(residual%*%w[,-1])
  init.B<-c(as.vector(init.beta),init.xi)
  return(init.B)
}

##########data#######################################################

data_gen <- function(n, dw, dz, beta, dw2) #data6,7,8
{
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
  
  dw1<-dw-dw2
  
  muiv<-rep(0,dw1)
  Sigmaiv<-diag(dw1)
  iv<-mvrnorm(n,muiv,Sigmaiv) 
  gammaiv<-c(0.8,seq(0.4,0.1,length.out = dw1-1))
  gammaz<-rep(0.8,dz)
  x<-iv%*%gammaiv+z%*%gammaz+mu
  
  y<-x%*%betax+z%*%betaz+epsilon
  
  # dw2 invalid iv
  delta.min <- 0.5
  delta.max <- 0.7 
  delta <- seq(delta.min,delta.max,len=dw2)
  muiiv <- rep(0,dw2)
  Sigmaiiv <- diag(dw2)
  nu <- mvrnorm(n,muiiv,Sigmaiiv)
  iiv <- matrix(rep(epsilon,dw2),nrow = n)%*%diag(delta) + nu
  # nu3<-rnorm(n)
  # iiv3<-epsilon*2+nu3
  # nu4<-rnorm(n)
  # iiv4<-epsilon*4+nu4
  
  #w<-cbind(iv,iiv3,iiv4)
  w<-cbind(iv,iiv)
  
  return(list(y=y,x=x,z=z,w=w))
}

########rep.fun#########################

# 2SLS & normal 
rep.fun <- function(index, para) #n init.B beta
{
  
  library(MASS)
  source("HDMS_HHD_FUN.R") # HHD or HD
  
  dz<-para$dz
  dw<-para$dw
  n <- para$n  # sample size
  beta <- para$beta
  dw2 <- para$dw2
  
  seed.num <- index*12345
  repeat{
    set.seed(seed.num)
    data <- data_gen(n, dw, dz, beta, dw2)
    y <- data$y
    x <- data$x
    z <- data$z
    w <- data$w
    
    seed.num <- seed.num + 10
    
    init.B <- init.B.2SLS(y,x,z,w)
    if (abs(init.B[1])<2) break
  }
  
  alpha <- c(0.1,0.05,0.01)
  normal.quantile.alpha <-qnorm(1-alpha/2)
  
  #2sls confidence interval
  library(sem)
  data <- cbind(y,x,z[,-1],w[,1])
  res <- summary(tsls(data[,1] ~ data[,2] + data[,3:4], ~ data[,5] + data[, 3:4], data=data))
  beta.2sls <- res$coefficients[2,1]
  beta.2sls.astd <- res$coefficients[2,2]
  tsls.ci.low <- beta.2sls - beta.2sls.astd*normal.quantile.alpha
  tsls.ci.up <- beta.2sls + beta.2sls.astd*normal.quantile.alpha
  tsls.ci <- cbind(tsls.ci.low, tsls.ci.up)
  #normal.ci <- tsls.ci
  
  #导入PPEL
  file <- paste("2SLS_HDMS_esti_n=",n,"dw=",dw, "s=", dw2, ".csv",sep="")
  HDMS.esti<-read.csv(file,header=TRUE)
  i <- index
  esti.mat <- HDMS.esti[,(5*(i-1)+2):(5*(i-1)+6)]
  esti.mat <- as.matrix(esti.mat)
  B.star <- esti.mat[, 1] 
  
  #normal confidence interval
  beta.debias <- esti.mat[1,2]
  beta.debias.astd <- esti.mat[1,5]
  normal.ci.low <- beta.debias - beta.debias.astd*normal.quantile.alpha
  normal.ci.up <- beta.debias + beta.debias.astd*normal.quantile.alpha
  normal.ci <- cbind(normal.ci.low, normal.ci.up)
  
  #chi confidence interval
  r <- dz + dw
  p <- length(B.star)
  zw<-cbind(z,w)
  grad.vect <- matrix(0.0, nrow = r, ncol = p)
  for(i in c(1:p))
  {
    index.B <- i
    if(index.B==1)
    {
      grad.mat<-diag(as.vector(x*(-1)))%*%zw
    }
    if(index.B>1 & index.B<(dz+2))
    {
      grad.mat<-diag(as.vector(z[,(index.B-1)]*(-1)))%*%zw  
    }
    if(index.B>(dz+1))
    {
      partial.g<-rep(0,r)
      partial.g[index.B]<--1
      partial.g<-rep(partial.g,n)
      partial.g<-matrix(partial.g,ncol=r,byrow=TRUE)
      grad.mat<-partial.g
    }
    
    grad.vect[,i] <-  colMeans(grad.mat)
  }
  
  xi1 <- c(1,rep(0,r-1))
  ttau <- 0.08*sqrt(log(p)/n) #0.5 0.1
  #ttau <- 0.2*(log(p)/n)
  
  library(lpSolve)
  objective.in <- rep(1,2*r)
  G <- t(grad.vect)
  const.mat <- rbind(cbind(G,-G),cbind(-G,G),-diag(2*r))
  const.dir <- "<="
  const.rhs <- c(rep(1,r)*ttau+xi1,rep(1,r)*ttau-xi1,rep(0,2*r))
  res <- lp(direction = "min",objective.in = objective.in, const.mat = const.mat,
            const.dir = const.dir, const.rhs = const.rhs)
  res.dir <- res$solution
  a.1.n <- res.dir[1:r] - res.dir[-(1:r)]
  
  para <- list(a.1.n=a.1.n, B.star=B.star, y=y, x=x, z=z, w=w)
  
  ell.star <- function(beta, alpha, para)
  {
    a.1.n <- para$a.1.n
    B.star <- para$B.star
    y <- para$y
    x <- para$x
    z <- para$z
    w <- para$w
    
    B <-c(beta,B.star[-1])  
    g.ee <- auxi.fun(B,y,x,z,w)$g.ee
    f.ee <- g.ee %*% a.1.n
    
    min.f.ee <- min(f.ee)
    max.f.ee <- max(f.ee)
    interval.low <- (-1+1/n)/max.f.ee
    interval.up <- (-1+1/n)/min.f.ee
    
    f.ee.lambda <- function(lambda,f.ee)
    {
      lambda <- matrix(lambda, ncol = 1) ##  (dz+dw) by 1 matrix
      f.ee <- as.matrix(f.ee)
      lambda.f <- 1.0 + f.ee %*% lambda ##  1+ lambda^t g(x,y,B)
      #log.f <- log.star(lambda.f)  #log
      log.f <- log(lambda.f)
      obj.fun <--1* sum(log.f)
      return(obj.fun)
    }
    
    res <-optimize(f=f.ee.lambda,interval=c(interval.low,interval.up),
                   f.ee=f.ee)
    esti.Chi <- -2*res$objective
    
    lh.ratio <- exp(-esti.Chi*0.5)
    Chi.quantile.alpha <- qchisq(1-alpha,1)
    obj.val <- lh.ratio - exp(-0.5* Chi.quantile.alpha)
    
    return(obj.val) 
  }
  
  res <- optimize(ell.star,c(0,1),alpha[1],para,maximum = TRUE)
  beta.max <- res$maximum
  
  B <-c(beta.max,B.star[-1])  
  g.ee <- auxi.fun(B,y,x,z,w)$g.ee
  f.ee <- g.ee %*% a.1.n
  top <- t(f.ee) %*% f.ee / n
  bottom <- n * (t(a.1.n) %*% grad.vect[,1])^2
  beta.max.astd <- sqrt(top/bottom)[1,1]
  
  chi.ci.low <- beta.max - beta.max.astd*normal.quantile.alpha
  chi.ci.up <- beta.max + beta.max.astd*normal.quantile.alpha
  chi.ci <- cbind(chi.ci.low, chi.ci.up)
  
  return(cbind(chi.ci,normal.ci,tsls.ci))
}

############demo###############################

start.t <- proc.time()

para <- list(n=n, dw = dw, dz = dz, beta = beta, dw2=dw2)

rep.seq <- c(1:rep.num)

# core.num<- detectCores(logical = FALSE)
core.num<-25

cl <- makeCluster(core.num)
clusterSetRNGStream(cl, iseed = NULL)

clusterExport(cl, c('init.B.2SLS', 'data_gen', "n"))
Chi.esti <- parSapplyLB(cl, rep.seq, rep.fun, para = para, simplify = "array")

stopCluster(cl)

write.csv(Chi.esti, file = paste("Chi_ci_","n=", n, "dw=", dw, "s=", dw2, ".csv",sep= ""))



