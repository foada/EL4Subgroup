############################################################
###    PEL estimation in bootstrap samples               ###                                    
############################################################

library(MASS)
library(parallel)

#setwd("C:/Users/jean/Desktop/optim_code")
source("HDMS_HHD_FUN.R") #

#rep.tuning
rep.tuning <- function(index, para) #BIC
{
  library(MASS)
  source("HDMS_HHD_FUN.R")
  
  eps.tol <- 0.005  #????
  
  y <- para$y
  x <- para$x
  z <- para$z
  w <- para$w
  init.B <- para$init.B
  para.matrix <- para$para.matrix
  
  n <- nrow(w)
  dz<-ncol(z)
  dw<-ncol(w)
  r<-dz+dw
  r1<-dz+1
  
  tau <- para.matrix[index, 1]
  nu <- para.matrix[index, 2]
  
  esti <- main.iter.tun(tau, nu, init.B, y, x, z, w)
  
  B <- esti$B
  xi <- B[(dz+2):r]
  lambda <- esti$lambda
  g.ee <- esti$g.ee
  
  obj.val <- ee.lambda(lambda, g.ee, xi, nu, tau)
  #B.dof <- sum(abs(B[(dz+2):r]) >= eps.tol)   
  B.dof <- sum(abs(B[(dz+2):r]) >= eps.tol) + r1
  
  
  #obj.val <- obj.val + B.dof * (log(n) + 2*log(dw-1)) ## EBIC
  obj.val <- obj.val + B.dof * (log(n) + 2*log(1)) #BIC
  
  # write.table(c(index, tau, nu, obj.val), file = paste("tuning_para", index, sep= "_"))
  
  return(obj.val)
}

#init.B ×¢Òâ½Ø¾à
init.B.2SLS <- function(y, x, z, w)
{
  zw1<-cbind(z,w[,1])
  res<-lm(x~zw1[,-1])
  x.hat<-res$fitted.values
  init.beta<-lm(y~cbind(x.hat,z[,-1]))$coef
  init.beta <- init.beta[c(2,1,3)]
  
  residual<-y-cbind(x,z)%*%init.beta
  residual<-as.vector(residual)
  residual<-diag(residual)
  
  #index <- sample(1:n, round(n*0.8), replace = FALSE)
  #init.xi<-colMeans((residual%*%w[,-1])[index,])
  
  init.xi<-colMeans(residual%*%w[,-1])
  init.B<-c(as.vector(init.beta),init.xi)
  return(init.B)
}

#tun.para
tun.para <- function(y, x, z, w) #tuning para
{
  
  init.B <- init.B.2SLS(y,x,z,w)
  
  n <- nrow(w)
  dw <- ncol(w)
  dz <- ncol(z)
  r<-dz+dw
  
  eps.tol <- 1e-10
  
  tau.max <- 1.5/sqrt(n) # 30 10 5 3 1.5
  tau.min <- 0.3/sqrt(n) # 15 5 1 0.5 0.3
  ntau <- 19 #19
  tau.seq = exp(seq(log(tau.min),log(tau.max),len=ntau)) ## candidate tau, tuning parameter for B
  tau.seq <-  matrix(tau.seq, ncol = 1)
  
  nu.max <- 20/sqrt(n) #10 20
  nu.min <- 10/sqrt(n)#1 5
  nnu <-  5 # 5 15
  nu.seq =  exp(seq(log(nu.min),log(nu.max),len=nnu))   ## candidate nu
  nu.seq <- matrix(nu.seq, ncol = 1)
  
  tau.matrix <- tau.seq %*% matrix(1, nrow = 1, ncol = nnu)
  tau.vect <- matrix(tau.matrix, ncol = 1)
  nu.matrix <- matrix(1, nrow = ntau, ncol = 1) %*% t(nu.seq)
  nu.vect <- matrix(nu.matrix, ncol = 1)
  
  para.matrix <- cbind(tau.vect, nu.vect)
  
  rep.num <- nrow(para.matrix)
  rep.seq <- c(1:rep.num)
  
  para <- list(y = y, x = x, z=z, w=w, init.B = init.B, para.matrix = para.matrix)
  
  #core.num <- detectCores(logical = FALSE)
  #core.num <- min(rep.num, core.num)
  core.num<-30
  
  cl <- makeCluster(core.num)
  
  clusterSetRNGStream(cl, iseed = NULL)
  
  clusterExport(cl, c('optim.lambda', 'optim.B2', 'main.iter.tun', "n"))
  
  obj.vect <- parSapplyLB(cl, rep.seq, rep.tuning, para = para, simplify = TRUE)
  
  stopCluster(cl)
  
  order.obj <- order(obj.vect)
  min.ind <- order.obj[1]
  opt.tau <- para.matrix[min.ind, 1]
  opt.nu <- para.matrix[min.ind, 2]
  
  opt.para <- matrix(c(opt.tau, opt.nu), ncol = 1)
  
  return(opt.para)
}

#data <- read.csv("C:/Users/jean/Desktop/Real Data/data3.csv", header = TRUE)[,-c(1)]
data <- read.csv("data3.csv", header = TRUE)[,-c(1)]
data <- as.matrix(data)
Data <- na.omit(data)
N<-nrow(Data)

n <- ceiling(0.8*N) # size of subsample: N-1 0.9N 0.8N 

#subsample.num <- 500
subsample.num <- 500 # N 500
beta.mat <- matrix(0, nrow = subsample.num, ncol = 2)

for(i in 1 : subsample.num)
{
  subsample <- sample(c(1 : N), n, replace = FALSE)
  #subsample <- c(1 : N)[-i]  # for n=N-1
  data <- Data[subsample, ]
  
  y <- data[,1]
  x <- data[,3]
  z <- cbind(rep(1, n), data[, 4])
  w <- data[,5:ncol(data)]
  y <- matrix(y,ncol=1)
  x <- matrix(x,ncol=1)
  
  pen.para <- tun.para(y,x,z,w)
  opt.tau <- pen.para[1]
  opt.nu <- pen.para[2]
  
  init.B <- init.B.2SLS(y,x,z,w)
  ebic_fit2 <- main.iter(opt.tau, opt.nu, init.B, y, x, z, w)
  B <- ebic_fit2$B

  beta.mat[i, ] <- B[c(1, 2)]
  print(B[c(1, 2)])
}

write.csv(beta.mat,file = paste("beta_mat_", n, ".csv", sep=""))
print(summary(beta.mat[, 1]))
print(summary(beta.mat[, 2]))
