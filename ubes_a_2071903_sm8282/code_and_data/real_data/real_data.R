#############################################################
#  real dara: point estimation and confidence interval      #
#############################################################

library(MASS)
library(parallel)

#setwd("C:/Users/jean/Desktop/optim_code")
source("HDMS_HHD_FUN.R") #

#data <- read.csv("C:/Users/jean/Desktop/Real Data/data5_1.csv", header = TRUE)[,-c(1,3)]
data <- read.csv("data3.csv", header = TRUE)[,-c(1)]
data <- as.matrix(data)
data <- na.omit(data)
n<-nrow(data)
y <- data[,1]
x <- data[,3]
z <- cbind(rep(1,n),data[,4])
#z <- rep(1,n)
w <- data[,5:ncol(data)]
y <- matrix(y,ncol=1)
x <- matrix(x,ncol=1)
#z <- matrix(z,ncol=1)

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
  
  tau.max <- 1.0/sqrt(n) # 30 10 5 3 1.5
  tau.min <- 0.3/sqrt(n) # 15 5 1 0.5 0.3
  ntau <- 20 #19
  tau.seq = exp(seq(log(tau.min),log(tau.max),len=ntau)) ## candidate tau, tuning parameter for B
  tau.seq <-  matrix(tau.seq, ncol = 1)
  
  nu.max <- 20/sqrt(n) #10 20
  nu.min <- 10/sqrt(n)#1 5
  nnu <- 10 # 5 15
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
  core.num<-50
  
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

#rep.fun

pen.para <- tun.para(y,x,z,w)
opt.tau <- pen.para[1]
opt.nu <- pen.para[2]

dw <- ncol(w)
dz <- ncol(z)

eps.tol <- 1e-10

init.B <- init.B.2SLS(y,x,z,w)

ebic_fit2 <- main.iter(opt.tau, opt.nu, init.B, y, x, z, w)
PEL2.ebic <- ebic_fit2$B
biacor_PEL2.ebic <- ebic_fit2$biacor_B
abias <- ebic_fit2$bias
astd <- ebic_fit2$B.astd
esti.B <- cbind(PEL2.ebic, biacor_PEL2.ebic, init.B, abias, astd)

#inference

alpha <- c(0.1,0.05,0.01)
normal.quantile.alpha <-qnorm(1-alpha/2)

#2sls confidence interval
library(sem)
data <- cbind(y,x,z[,-1],w[,1])
res <- summary(tsls(data[,1] ~ data[,2] + data[,3], ~ data[,4] + data[, 3], data=data))
beta.2sls <- res$coefficients[2,1]
beta.2sls.astd <- res$coefficients[2,2]
tsls.ci.low <- beta.2sls - beta.2sls.astd*normal.quantile.alpha
tsls.ci.up <- beta.2sls + beta.2sls.astd*normal.quantile.alpha
tsls.ci <- cbind(tsls.ci.low, tsls.ci.up)

#print(res)

esti.mat <- esti.B
# esti.mat <- read.csv("esti.B.csv", header = TRUE)[,-1]
# esti.mat <- as.matrix(esti.mat)

#normal confidence interval
beta.debias <- esti.mat[1,2]
beta.debias.astd <- esti.mat[1,5]
normal.ci.low <- beta.debias - beta.debias.astd*normal.quantile.alpha
normal.ci.up <- beta.debias + beta.debias.astd*normal.quantile.alpha
normal.ci <- cbind(normal.ci.low, normal.ci.up)

#chi confidence interval
B.star <- esti.mat[, 1] 
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
ttau <- 0.16*sqrt(log(p)/n) # change the tuning parameter
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

# a <- seq(-3,3,len=100)
# b <- rep(0,100)
# for (i in 1:100) b[i] <- ell.star(a[i],alpha = 0.01, para)
# plot(a,b,type="l")

res <- optimize(ell.star,c(0,2),alpha[1],para,maximum = TRUE)
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

ci<-cbind(chi.ci,normal.ci,tsls.ci)

print(pen.para)
print(esti.B)
print(beta.max)
print(c(beta.max.astd,beta.debias.astd,beta.2sls.astd))
print(ci)