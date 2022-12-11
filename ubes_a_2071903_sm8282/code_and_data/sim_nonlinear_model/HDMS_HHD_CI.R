##################################################
### use PPEL to calculate confidence intervals ###
##################################################
library(MASS)
library(parallel)

#setwd("/home/jiazhang/HDMS_ZJ")  #
#setwd("C:/Users/jean/Desktop/res0419/data4")
source("HDMS_HHD_FUN.R") # HHD or HD

######init.B##########################################

init.B <- function(Y)
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

data_gen <- function(n, r, s, beta) ## return a n * (r+1) matrix
{
  beta1 <- beta[1]
  beta2 <- beta[2]
  
  alpha <- rnorm(n, 1, 1)  # with mean zero
  Alpha <- matrix(alpha, n, r+1)-1
  Error1 <- matrix(rnorm(n*(r+1), 0, 1/sqrt(2)), n, r+1)
  Error <- 1 / sqrt(2) * Alpha + Error1
  lambda.vec <- lam.fun(beta1, r)
  beta2_star <- beta2 + 0.3    # change: 0.5 0.3 0.1
  Beta2 <- cbind(matrix(beta2, n, r+1-s), matrix(beta2_star, n, s))
  Y <- alpha %*% t(lambda.vec) + Beta2 + Error
  
  return(Y)
}

########rep.fun#########################

# 2SLS & normal 
rep.fun <- function(index, para) #n init.B beta
{
  
  library(MASS)
  source("HDMS_HHD_FUN.R") # HHD or HD
  
  n <- para$n  # sample size
  r <- para$r
  r1 <- para$r1
  s <- para$s
  r2 <- r - r1
  p <- 2 + r2
  
  beta <- para$beta
  
  seed.num <- index * 12345
  Y <- data_gen(n, r, s, beta)
  
  # init.B <- c(beta, rep(0, r2-s), (lam.fun(beta[1], r)[(r+2-s):(r+1)] - 1)*0.5)  # 真值
  
  alpha <- c(0.1,0.05,0.01)
  normal.quantile.alpha <-qnorm(1-alpha/2)
  
  #2sls confidence interval
  # library(sem)
  # data <- cbind(y,x,z[,-1],w[,1])
  # res <- summary(tsls(data[,1] ~ data[,2] + data[,3:4], ~ data[,5] + data[, 3:4], data=data))
  # beta.2sls <- res$coefficients[2,1]
  # beta.2sls.astd <- res$coefficients[2,2]
  # tsls.ci.low <- beta.2sls - beta.2sls.astd*normal.quantile.alpha
  # tsls.ci.up <- beta.2sls + beta.2sls.astd*normal.quantile.alpha
  # tsls.ci <- cbind(tsls.ci.low, tsls.ci.up)
  #normal.ci <- tsls.ci
  
  #导入PPEL
  file <- paste("HDMS_esti_n=",n,"r=",r, "s=", s, ".csv",sep="")
  HDMS.esti<-read.csv(file,header=TRUE)
  i <- index
  esti.mat <- HDMS.esti[,(5*(i-1)+2):(5*(i-1)+6)]
  esti.mat <- as.matrix(esti.mat)
  B.star <- esti.mat[, 1] 
  beta1 <- B.star[1]
  beta2 <- B.star[2]
  
  #normal confidence interval
  beta.debias <- esti.mat[1,2]
  beta.debias.astd <- esti.mat[1,5]
  normal.ci.low <- beta.debias - beta.debias.astd*normal.quantile.alpha
  normal.ci.up <- beta.debias + beta.debias.astd*normal.quantile.alpha
  normal.ci <- cbind(normal.ci.low, normal.ci.up)
  
  # chi confidence interval
  # 计算投影方向
  grad.vect <- matrix(0, nrow = r, ncol = p)
  for(i in c(1:p))
  {
    index.B <- i
    if(index.B == 1)
    {
      grad.mat <- (beta2 - Y[, 1]) %*% t(der.lam.fun(beta1, r))
    }
    if(index.B == 2)
    {
      grad.mat <- matrix(lam.fun(beta1, r)[-1] - 1, n, r, byrow = T)
    }
    if(index.B > 2 )
    {
      partial.g<-rep(0,r)
      partial.g[(index.B - 2 + r1)] <- -1
      partial.g<-rep(partial.g,n)
      partial.g<-matrix(partial.g, ncol=r, byrow=TRUE)
      grad.mat<-partial.g
    }
    
    grad.vect[,i] <-  colMeans(grad.mat)
  }
  
  xi1 <- c(1, rep(0,p-1))
  ttau <- 0.08*sqrt(log(p)/n) # tuing: 0.5 0.1
  
  library(lpSolve)
  objective.in <- rep(1,2*r)
  G <- t(grad.vect)
  const.mat <- rbind(cbind(G,-G),cbind(-G,G),-diag(2*r))
  const.dir <- "<="
  const.rhs <- c(rep(1,p)*ttau+xi1,rep(1,p)*ttau-xi1,rep(0,2*r))
  res <- lp(direction = "min",objective.in = objective.in, const.mat = const.mat,
            const.dir = const.dir, const.rhs = const.rhs, scale = 196)
  res.dir <- res$solution
  a.1.n <- res.dir[1:r] - res.dir[-(1:r)]
  
  # projected-EL估计
  para <- list(a.1.n=a.1.n, B.star=B.star, Y = Y)
  ell.star <- function(beta, para)
  {
    a.1.n <- para$a.1.n
    B.star <- para$B.star
    Y <- para$Y
    
    B <-c(beta,B.star[-1])  
    g.ee <- auxi.fun(B,Y)
    f.ee <- g.ee %*% a.1.n
    
    min.f.ee <- min(f.ee)
    max.f.ee <- max(f.ee)
    interval.low <- (-1+1/n)/max.f.ee
    interval.up <- (-1+1/n)/min.f.ee
    
    f.ee.lambda <- function(lambda,f.ee)
    {
      lambda <- matrix(lambda, ncol = 1) 
      f.ee <- as.matrix(f.ee)
      lambda.f <- 1.0 + f.ee %*% lambda 
      #log.f <- log.star(lambda.f)  #log
      log.f <- log(lambda.f)
      obj.fun <--1* sum(log.f)
      return(obj.fun)
    }
    
    res <-optimize(f=f.ee.lambda,interval=c(interval.low,interval.up),
                   f.ee=f.ee)
    esti.Chi <- -2*res$objective
    
    lh.ratio <- exp(-esti.Chi*0.5)
    #Chi.quantile.alpha <- qchisq(1-alpha,1)
    #obj.val <- lh.ratio - exp(-0.5* Chi.quantile.alpha)
    obj.val <- lh.ratio
    
    return(obj.val) 
  }
  
  res <- optimize(ell.star,c(0,1),para,maximum = TRUE)
  beta.max <- res$maximum
  
  B <-c(beta.max,B.star[-1])  
  g.ee <- auxi.fun(B,Y)
  f.ee <- g.ee %*% a.1.n
  top <- t(f.ee) %*% f.ee / n
  bottom <- n * (t(a.1.n) %*% grad.vect[,1])^2
  beta.max.astd <- sqrt(top/bottom)[1,1]
  
  chi.ci.low <- beta.max - beta.max.astd*normal.quantile.alpha
  chi.ci.up <- beta.max + beta.max.astd*normal.quantile.alpha
  chi.ci <- cbind(chi.ci.low, chi.ci.up)
  
  tsls.ci <- normal.ci  # 可删除
  
  return(cbind(chi.ci,normal.ci,tsls.ci))
}

############demo###############################

start.t <- proc.time()

para <- list(n = n, r = r, r1 = r1, beta = beta, s = s)

rep.seq <- c(1:rep.num)

# core.num<- detectCores(logical = FALSE)
core.num<-50

cl <- makeCluster(core.num)
clusterSetRNGStream(cl, iseed = NULL)

clusterExport(cl, c('init.B', 'data_gen', "n"))
Chi.esti <- parSapplyLB(cl, rep.seq, rep.fun, para = para, simplify = "array")

stopCluster(cl)

write.csv(Chi.esti, file = paste("Chi_ci_", "n=", n, "r=", r, "s=", s, ".csv", sep= ""))



