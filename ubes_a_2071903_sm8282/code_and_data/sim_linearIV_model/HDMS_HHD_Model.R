####################################################################
#######   calculate FP FN RMSE BIAS STD for Tables 1 and 2   #######
#######         change L191 for different signals            #######
####################################################################
library(MASS)
library(parallel)

#setwd("/home/jiazhang/HDMS_ZJ")  #
#setwd("C:/Users/jean/Desktop/optim_code")
source("HDMS_HHD_FUN.R") #

#######rep.tuning###################################################

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
  
  n <- nrow(x)
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
  
  #obj.val <- ee.lambda(lambda, g.ee, xi, nu, tau)
  obj.val <- 2 * ee.fun(lambda, g.ee)
  #B.dof <- sum(abs(B[(dz+2):r]) >= eps.tol)   
  B.dof <- sum(abs(B[(dz+2):r]) >= eps.tol) + r1
  
  
  #obj.val <- obj.val + B.dof * (log(n) + 2*log(dw-1)) ## EBIC
  obj.val <- obj.val + B.dof * (log(n) + 2*log(1)) #BIC
  
  # write.table(c(index, tau, nu, obj.val), file = paste("tuning_para", index, sep= "_"))
  
  return(obj.val)
}

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

#########tun.para##########################################

tun.para <- function(n, dw, dz, beta, dw2) #tuning para
{
  seed.num <- 1
  repeat{
    set.seed(seed.num * 12345)
    data <- data_gen(n, dw, dz, beta, dw2)
    y <- data$y
    x <- data$x
    z <- data$z
    w <- data$w
    
    seed.num <- seed.num + 1
    
    init.B <- init.B.2SLS(y,x,z,w)
    #if (abs(init.B[1]-0.5) < 0.5) break # 1
    if (abs(init.B[1]) < 2) break
  }
  
  r<-dz+dw
  
  eps.tol <- 1e-10
  
  tau.max <- 1.5/sqrt(n) # 30 10 5 3 1.5
  tau.min <- 0.3/sqrt(n) # 15 5 1 0.5 0.3
  ntau <- 20 #19
  tau.seq = exp(seq(log(tau.min),log(tau.max),len=ntau)) ## candidate tau, tuning parameter for B
  tau.seq <-  matrix(tau.seq, ncol = 1)
  
  nu.max <- 20/sqrt(n) #10 20
  nu.min <- 10/sqrt(n)#1 5
  nnu <-  10 # 5 15
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
  core.num<-25
  
  cl <- makeCluster(core.num)
  
  clusterSetRNGStream(cl, iseed = NULL)
  
  clusterExport(cl, c('optim.lambda', 'optim.B2', 'main.iter', "n"))
  
  obj.vect <- parSapplyLB(cl, rep.seq, rep.tuning, para = para, simplify = TRUE)
  
  stopCluster(cl)
  
  order.obj <- order(obj.vect)
  min.ind <- order.obj[1]
  opt.tau <- para.matrix[min.ind, 1]
  opt.nu <- para.matrix[min.ind, 2]
  
  opt.para <- matrix(c(opt.tau, opt.nu), ncol = 1)
  
  return(opt.para)
}

##########data#######################################################

data_gen <- function(n, dw, dz, beta, dw2) #data6,7,8
{
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
  delta.min <- 0.5  # change for different signals: 0.3-0.5, 0.7-0.9
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

rep.fun <- function(index, para) #n init.B beta
{
  
  library(MASS)
  source("HDMS_HHD_FUN.R")
  
  eps.tol <- 1e-10
  
  dz<-para$dz
  dw<-para$dw
  pen.para <- para$pen.para
  
  opt.tau <- pen.para[1]
  opt.nu <- pen.para[2]
  
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
    #if (abs(init.B[1]-0.5)<1) break
    if (abs(init.B[1])<2) break
  }
  
  r<-dz+dw
  
  ebic_fit2 <- main.iter(opt.tau, opt.nu, init.B, y, x, z, w, dw2)
  PEL2.ebic <- ebic_fit2$B
  biacor_PEL2.ebic <- ebic_fit2$biacor_B
  abias <- ebic_fit2$bias
  astd <- ebic_fit2$B.astd
  
  # if(index<11)
  # {
  #   xi.iter <- ebic_fit2$B.iter[,(r-5):r]
  #   write.csv(xi.iter, file = paste("xi_iter_", index, ".csv", sep = ""))
  # }
  
  esti.B <- cbind(PEL2.ebic, biacor_PEL2.ebic, init.B, abias, astd)
  
  if (index%%85==0) write.csv(esti.B, file = paste(n, "_", dw, "_", dw2, "-", index, ".csv", sep=""))
  
  return(esti.B)
}

############demo###############################

start.t <- proc.time()

pen.para <- tun.para(n, dw, dz, beta, dw2)

#write.csv(pen.para, file =paste("penpara_","n=",n,"dw=",dw,".csv",sep=""))

para <- list(n=n, dw = dw, dz = dz, pen.para = pen.para, beta = beta, dw2=dw2)

rep.seq <- c(1:rep.num)

#core.num<- detectCores(logical = FALSE)
core.num<-25

cl <- makeCluster(core.num)
clusterSetRNGStream(cl, iseed = NULL)

clusterExport(cl, c('init.B.2SLS', 'data_gen', 'optim.lambda', 'optim.B2','main.iter', "n"))
HDMS.esti <- parSapplyLB(cl, rep.seq, rep.fun, para = para, simplify = "array")

stopCluster(cl)

write.csv(HDMS.esti, file = paste("2SLS_HDMS_esti_","n=", n, "dw=", dw, "s=", dw2, ".csv",sep= ""))

r<-dz+dw

#dw2<-3
#dw2<-round(log(n))

#########FP FN#################
dw1<-dw-dw2
supp.B0<-c(rep(0,dw1-1),rep(1,dw2))
eps.tol <- 0.005 #0.005

dist.mat <- matrix(0, nrow = rep.num, ncol = 2) 
FP.mat <- matrix(0, nrow = rep.num, ncol = 2)  ## Falsely Positive Number
FN.mat <- matrix(0, nrow = rep.num, ncol = 2)  ## Falsely Negative

for (i in 1:rep.num){
  esti.mat <- HDMS.esti[ , , i]
  esti.mat <- as.matrix(esti.mat)
  
  PEL2.ebic <- esti.mat[, 1]
  PEL2.ebic.beta<-PEL2.ebic[1:(dz+1)]
  
  dist.ebic<- norm(PEL2.ebic.beta - beta, type = "F")
  dist.mat[i,1] <- dist.ebic
  
  act.ebic <- abs(PEL2.ebic[(dz+2):r])> eps.tol
  FP.mat[i,1] <- mean(act.ebic[1:(dw1-1)])
  FN.mat[i,1] <- -1 * mean(act.ebic[-(1:(dw1-1))]-supp.B0[-(1:(dw1-1))])
  
  biacor_PEL2.ebic <- esti.mat[, 2]
  biacor_PEL2.ebic.beta<-biacor_PEL2.ebic[1:(dz+1)]
  
  biacor_dist.ebic <- norm(biacor_PEL2.ebic.beta - beta, type = "F")
  dist.mat[i,2] <- biacor_dist.ebic
  
  biacor_act.ebic <- abs(biacor_PEL2.ebic[(dz+2):r])> eps.tol
  FP.mat[i,2] <- mean(biacor_act.ebic[1:(dw1-1)])
  FN.mat[i,2] <- -1 * mean(biacor_act.ebic[-(1:(dw1-1))]-supp.B0[-(1:(dw1-1))])
}

dist.mean <- colMeans(dist.mat)
dist.std <- sqrt(diag(cov(dist.mat)))

FP.mean <- colMeans(FP.mat)
FP.std <- sqrt(diag(cov(FP.mat)))

FN.mean <- colMeans(FN.mat)
FN.std <- sqrt(diag(cov(FN.mat)))

result.esti <- rbind(dist.mean, dist.std, FP.mean, FP.std, FN.mean, FN.std)

write.csv(result.esti, file = paste("n=",n,"dw=",dw, "s=", dw2, ".csv", sep= ""))

print(list(pen.para=pen.para, result.esti=result.esti))


#betax.mat <- matrix(0, nrow = rep.num, ncol = 2)
betax.mat <- matrix(0, nrow = rep.num, ncol = 3)
betax.std <- matrix(0, nrow = rep.num, ncol = 1)

for (i in 1:rep.num){
  esti.mat <- HDMS.esti[ , , i]
  #esti.mat <- HDMS.esti[,((2*i):(2*i+1))]
  esti.mat <- as.matrix(esti.mat)
  
  PEL2.ebic <- esti.mat[, 1]
  PEL2.ebic.betax <- PEL2.ebic[1]
  betax.mat[i,1] <- PEL2.ebic.betax
  
  #act.ebic <- abs(PEL2.ebic[(dz+2):r])> eps.tol
  
  biacor_PEL2.ebic <- esti.mat[, 2]
  biacor_PEL2.ebic.betax<-biacor_PEL2.ebic[1]
  betax.mat[i,2] <- biacor_PEL2.ebic.betax
  
  # biacor_act.ebic <- abs(biacor_PEL2.ebic[(dz+2):r])> eps.tol
  
  init.B <- esti.mat[, 3]
  init.B.betax<-init.B[1]
  betax.mat[i,3] <- init.B.betax
  
  B.astd <- esti.mat[,5]
  betax.std[i,1] <- B.astd[1]
}

RMSE <- sqrt(colMeans((betax.mat-beta[1])^2))
bias <- colMeans(betax.mat-beta[1])
std <- sqrt(RMSE^2 - bias^2)

error <- t(rbind(RMSE, bias, std))

print(error)

write.csv(error, file = paste("error_","n=",n,"dw=",dw, "s=", dw2, ".csv", sep= ""))

######time#######
end.t <- proc.time()
op.time <- end.t - start.t
print(op.time)







