###################################################
###     COMPUTE FP FN RMSE BIAS STD FOR PEL     ###
###################################################
library(MASS)
library(parallel)

#setwd("/home/jiazhang/HDMS_ZJ")  #
#setwd("C:/Users/jean/Desktop/optim_code")
source("HDMS_HHD_FUN.R") #

#######rep.tuning###################################################

rep.tuning <- function(index, para) # BIC
{
  library(MASS)
  source("HDMS_HHD_FUN.R")
  
  eps.tol <- 0.005  #????
  
  Y <- para$Y
  init.B <- para$init.B
  para.matrix <- para$para.matrix
  
  n <- nrow(Y)
  r <- ncol(Y) - 1
  p <- length(init.B)
  r2 <- p - 2
  r1 <- r - r2
  
  tau <- para.matrix[index, 1]
  nu <- para.matrix[index, 2]
  
  esti <- main.iter.tun(tau, nu, init.B, Y)
  
  B <- esti$B
  xi <- B[(1+2):p]
  lambda <- esti$lambda
  g.ee <- esti$g.ee
  
  obj.val <- 2 * ee.fun(lambda, g.ee)
  #obj.val <- ee.lambda(lambda, g.ee, xi, nu, tau)
  #B.dof <- sum(abs(B[(dz+2):r]) >= eps.tol)   
  B.dof <- sum(abs(B[(1+2):p]) >= eps.tol) + 2
  
  #obj.val <- obj.val + B.dof * (log(n) + 2*log(dw-1)) ## EBIC
  obj.val <- obj.val + B.dof * (log(n) + 2*log(1)) #BIC
  
  # write.table(c(index, tau, nu, obj.val), file = paste("tuning_para", index, sep= "_"))
  
  return(obj.val)
}

#########tun.para##########################################

tun.para <- function(n, r, r1, beta, s) #tuning para
{
  r2 <- r - r1
  p <- r2 + 2
  
  set.seed(12345)
  Y <- data_gen(n, r, s, beta)
  
  # init.B <- c(beta, rep(0, r2-s), (lam.fun(beta[1], r)[(r+2-s):(r+1)] - 1)*0.5)  # 真值作为初值
  init.B <- c(beta, rep(0, r2-s), (lam.fun(beta[1], r)[(r+2-s):(r+1)] - 1)*0.2) + rnorm(p, 0, 0.1)
  
  eps.tol <- 1e-10
  
  tau.max <- 1.5/sqrt(n)  # 30 10 5 3   1.5
  tau.min <- 0.3/sqrt(n)  # 15 5  1 0.5 0.3
  ntau <- 20 #19
  tau.seq <- exp(seq(log(tau.min),log(tau.max),len=ntau)) ## candidate tau, tuning parameter for B
  tau.seq <-  matrix(tau.seq, ncol = 1)
  
  nu.max <- 20/sqrt(n)  #10 20
  nu.min <- 10/sqrt(n)  #1  5
  nnu <-  10 # 5 15
  nu.seq <-  exp(seq(log(nu.min),log(nu.max),len=nnu))   ## candidate nu
  nu.seq <- matrix(nu.seq, ncol = 1)
  
  tau.matrix <- tau.seq %*% matrix(1, nrow = 1, ncol = nnu)
  tau.vect <- matrix(tau.matrix, ncol = 1)
  nu.matrix <- matrix(1, nrow = ntau, ncol = 1) %*% t(nu.seq)
  nu.vect <- matrix(nu.matrix, ncol = 1)
  
  para.matrix <- cbind(tau.vect, nu.vect)
  
  rep.num <- nrow(para.matrix)
  rep.seq <- c(1:rep.num)
  
  para <- list(Y = Y, init.B = init.B, para.matrix = para.matrix)
  
  #core.num <- detectCores(logical = FALSE)
  #core.num <- min(rep.num, core.num)
  core.num<-25
  
  cl <- makeCluster(core.num)
  
  clusterSetRNGStream(cl, iseed = NULL)
  
  clusterExport(cl, c('optim.lambda', 'optim.B2', 'main.iter.tun', "n", "r1"))
  
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

data_gen <- function(n, r, s, beta) ## return a n * (r+1) matrix
{
  beta1 <- beta[1]
  beta2 <- beta[2]
  
  alpha <- rnorm(n,1,1)  # with mean zero
  Alpha <- matrix(alpha, n, r+1) - 1
  Error1 <- matrix(rnorm(n*(r+1), 0, 1/sqrt(2)), n, r+1)
  Error <- 1 / sqrt(2) * Alpha + Error1
  lambda.vec <- lam.fun(beta1, r)
  # beta2_star <- beta2 + rnorm(1, 0.5, 0.3)    # change: 0.5
  beta2_star <- beta2 + 0.3
  Beta2 <- cbind(matrix(beta2, n, r+1-s), matrix(beta2_star, n, s))
  Y <- alpha %*% t(lambda.vec) + Beta2 + Error
  
  return(Y)
}

########rep.fun#########################

rep.fun <- function(index, para) #n init.B beta
{
  
  library(MASS)
  source("HDMS_HHD_FUN.R")
  
  eps.tol <- 1e-10
  
  r <- para$r
  r1 <- para$r1
  r2 <- r - r1
  p <- 2 + r2
  
  pen.para <- para$pen.para
  opt.tau <- pen.para[1]
  opt.nu <- pen.para[2]
  
  n <- para$n  # sample size
  beta <- para$beta
  s <- para$s
  
  seed.num <- index * 12345
  Y <- data_gen(n, r, s, beta)
  
  # init.B <- c(beta, rep(0, r2-s), (lam.fun(beta[1], r)[(r+2-s):(r+1)] - 1)*0.5)
  init.B <- c(beta, rep(0, r2-s), (lam.fun(beta[1], r)[(r+2-s):(r+1)] - 1)*0.3) + rnorm(p, 0, 0.1)
  
  ebic_fit2 <- main.iter(opt.tau, opt.nu, init.B, Y, s)
  PEL2.ebic <- ebic_fit2$B
  biacor_PEL2.ebic <- ebic_fit2$biacor_B
  abias <- ebic_fit2$bias
  astd <- ebic_fit2$B.astd
  
  if(index < 11)
  {
    B.iter <- ebic_fit2$B.iter
    write.csv(B.iter, file = paste("B_iter_", index, ".csv", sep = ""))
  }
  
  esti.B <- cbind(PEL2.ebic, biacor_PEL2.ebic, init.B, abias, astd)
  
  if (index %% 85 == 0) write.csv(esti.B, file = paste(n, "_", r, "_", s, "_", index, ".csv", sep=""))
  
  return(esti.B)
}

############ demo ###############################

start.t <- proc.time()

pen.para <- tun.para(n, r, r1, beta, s)

write.csv(pen.para, file =paste("penpara_", "n=", n, "r", r, "s=", s, ".csv", sep=""))

para <- list(n = n, r = r, r1 = r1, pen.para = pen.para, beta = beta, s = s)

rep.seq <- c(1:rep.num)

#core.num<- detectCores(logical = FALSE)
core.num<-25

cl <- makeCluster(core.num)
clusterSetRNGStream(cl, iseed = NULL)

clusterExport(cl, c('data_gen', 'optim.lambda', 'optim.B2','main.iter', "n"))
HDMS.esti <- parSapplyLB(cl, rep.seq, rep.fun, para = para, simplify = "array")

stopCluster(cl)

write.csv(HDMS.esti, file = paste("HDMS_esti_","n=", n, "r=", r, "s=", s, ".csv",sep= ""))

######### FP FN ##############################
r2 <- r - r1
supp.B0 <- c(rep(0, r2-s), rep(1, s))
eps.tol <- 0.005   # 0.005
p <- r2 + 2

dist.mat <- matrix(0, nrow = rep.num, ncol = 2) 
FP.mat <- matrix(0, nrow = rep.num, ncol = 2)  ## Falsely Positive Number
FN.mat <- matrix(0, nrow = rep.num, ncol = 2)  ## Falsely Negative

for (i in 1:rep.num)
  {
  esti.mat <- HDMS.esti[ , , i]
  esti.mat <- as.matrix(esti.mat)
  
  PEL2.ebic <- esti.mat[, 1]
  PEL2.ebic.beta <- PEL2.ebic[1:2]
  
  dist.ebic<- norm(matrix(PEL2.ebic.beta - beta, 2, 1), type = "F")
  dist.mat[i,1] <- dist.ebic
  
  act.ebic <- abs(PEL2.ebic[3:p]) > eps.tol
  FP.mat[i,1] <- mean(act.ebic[1:(r2-s)])
  FN.mat[i,1] <- -1 * mean(act.ebic[-(1:(r2-s))] - supp.B0[-(1:(r2-s))])
  
  biacor_PEL2.ebic <- esti.mat[, 2]
  biacor_PEL2.ebic.beta <- biacor_PEL2.ebic[1:2]
  
  biacor_dist.ebic <- norm(matrix(biacor_PEL2.ebic.beta - beta, 2, 1), type = "F")
  dist.mat[i,2] <- biacor_dist.ebic
  
  biacor_act.ebic <- abs(biacor_PEL2.ebic[3:p])> eps.tol
  FP.mat[i,2] <- mean(biacor_act.ebic[1:(r2-s)])
  FN.mat[i,2] <- -1 * mean(biacor_act.ebic[-(1:(r2-s))]-supp.B0[-(1:(r2-s))])
}

dist.mean <- colMeans(dist.mat)
dist.std <- sqrt(diag(cov(dist.mat)))

FP.mean <- colMeans(FP.mat)
FP.std <- sqrt(diag(cov(FP.mat)))

FN.mean <- colMeans(FN.mat)
FN.std <- sqrt(diag(cov(FN.mat)))

result.esti <- rbind(dist.mean, dist.std, FP.mean, FP.std, FN.mean, FN.std)

write.csv(result.esti, file = paste("n=",n,"r=",r, "s=", s, ".csv", sep= ""))

print(list(pen.para=pen.para, result.esti=result.esti))

########### RMSE #############################################

betax.mat <- matrix(0, nrow = rep.num, ncol = 3)
betax.std <- matrix(0, nrow = rep.num, ncol = 1)

for (i in 1:rep.num){
  esti.mat <- HDMS.esti[ , , i]
  esti.mat <- as.matrix(esti.mat)
  
  PEL2.ebic <- esti.mat[, 1]
  PEL2.ebic.betax <- PEL2.ebic[1]
  betax.mat[i,1] <- PEL2.ebic.betax
  
  biacor_PEL2.ebic <- esti.mat[, 2]
  biacor_PEL2.ebic.betax <- biacor_PEL2.ebic[1]
  betax.mat[i,2] <- biacor_PEL2.ebic.betax
  
  init.B <- esti.mat[, 3]
  init.B.betax <- init.B[1]
  betax.mat[i,3] <- init.B.betax
  
  B.astd <- esti.mat[, 5]
  betax.std[i,1] <- B.astd[1]
}

RMSE <- sqrt(colMeans((betax.mat - beta[1]) ^ 2))
bias <- colMeans(betax.mat - beta[1])
std <- sqrt(RMSE^2 - bias^2)

error <- t(rbind(RMSE, bias, std))

print(error)

write.csv(error, file = paste("error_","n=",n,"r=",r, "s=", s, ".csv", sep= ""))

######time#######
end.t <- proc.time()
op.time <- end.t - start.t
print(op.time)







