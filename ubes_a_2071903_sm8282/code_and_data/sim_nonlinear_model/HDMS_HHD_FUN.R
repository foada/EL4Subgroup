#############################################################
###  useful functions for HDMS_HHD_Model and HDMS_HHD_CI ###
#############################################################


############pseudo log##########################

log.star <- function(x, thresh = 1/n){ ## approximated logarithm function to overcome the bounded support problem
  x <- as.matrix(x)
  temp.x <- x
  temp.x[x < thresh] <- thresh
  case.1 <- log(temp.x)  ## for the case  x >= thresh
  case.2 <- log(thresh) - 1.5 + 2*x/thresh - 0.5 * (x/thresh)^2  ## for the case  x < thresh
  log.x <- (x >= thresh) * case.1 + (x < thresh) * case.2
  
  return(log.x)
}


der.log.star <- function(x, thresh = 1/n){
  x <- as.matrix(x)
  temp.x <- x
  temp.x[x < thresh] <- thresh
  deriv <- (x >= thresh) * (1.0/temp.x) + (x < thresh) * (2.0/thresh - x/thresh^2)
  
  return(deriv)
}


der2.log.star <- function(x, thresh = 1/n){
  x <- as.matrix(x)
  temp.x <- x
  temp.x[x < thresh] <- thresh
  deriv2 <- (x >= thresh) * (-1.0/temp.x^2) + (x < thresh) * (- 1/thresh^2)
  
  return(deriv2)
}

###########scad#################################

scad.penalty <- function(beta, lambda){
  ## where lambda is tuning parameter
  a <- 3.7 ## a = 3.7 suggested by Fan and Li(2001)
  beta <- as.matrix(beta, ncol=1)
  abs.beta <- abs(beta)
  temp.1 <- lambda * abs.beta  ## for the case abs.beta <= lambda
  temp.2 <- (2*a*lambda*abs.beta - abs.beta^2 - lambda^2)/(2*(a-1.0)) ## for the case abs.beta <= a*lambda
  temp.3 <- lambda^2 * (a+1.0)/2
  
  scad.val <- temp.1*(abs.beta <= lambda) + temp.2*(abs.beta > lambda)*(abs.beta <= a * lambda) + temp.3*(abs.beta > a * lambda)
  
  return(scad.val)
}

#对绝对值函数求导
scad.grad <- function(beta, lambda){
  ## where lambda is tuning parameter
  a <- 3.7 ## a = 3.7 suggested by Fan and Li(2001)
  beta <- as.matrix(beta, ncol=1)
  abs.beta <- abs(beta)
  
  #temp.1 <- lambda  ## for the case abs.beta <= lambda
  temp.1 <- lambda*sign(beta) #改
  #temp.2 <- (a*lambda - abs.beta)/(a-1.0) ## for the case abs.beta <= a*lambda
  temp.2 <- max( a*lambda - abs.beta,0)/(a-1.0)*sign(beta) #改
  
  grad.scad <- temp.1*(abs.beta <= lambda) + temp.2*(abs.beta > lambda)*(abs.beta <= a * lambda)
  
  return(grad.scad)
}

#二次近似及其倒数
scad.quad <- function(beta,lambda)
{
  beta <- as.matrix(beta)
  beta <- matrix(beta, ncol = 1)
  d <- length(beta)
  thresh <- lambda/4 #?
  
  init.beta <- beta
  index.1 <- (abs(init.beta) < 1e-4)
  index.2 <- (abs(init.beta) < thresh)
  index.3 <- as.logical((1-index.1)*index.2)
  init.beta[index.1] <- thresh #?
  init.beta[index.3] <- sign(beta[index.3])*thresh #?
  
  temp.1 <- scad.penalty(init.beta, lambda)
  temp.2 <- scad.grad(abs(init.beta), lambda)
  temp.3 <- (beta^2-init.beta^2) / 2
  
  scad.val <- temp.1 + temp.2 * temp.3 / abs(init.beta)
  scad.val <- sum(scad.val)
  
  scad.der <- temp.2 * beta / abs(init.beta)
  
  temp.vect<- temp.2 / abs(init.beta)
  scad.hessian <- diag(c(temp.vect), nrow = d, ncol = d)
  
  return(list(scad.val = scad.val, scad.der = scad.der, scad.hessian = scad.hessian))
}

##########lambda.fun##################################

lam.fun <- function(beta1, r) # return a (r+1)-dimensional vector
{
  r.vec <- c(0, 1:r)
  lambda.vec <- 2 / ( 1 + exp(- beta1 * r.vec ))  # sigmoid function
  return(lambda.vec)
}

der.lam.fun <- function(beta1, r) # return a (r)-dimensional vector
{
  r.vec <- 1:r
  der.lam.vec <- 2 * r.vec * exp(- beta1 * r.vec) / (1 + exp(- beta1 * r.vec))^2 
  return(der.lam.vec) 
}

### Extended moment function #################################

auxi.fun<-function(B, Y)  ## return a n * r matrix
{
 r <- ncol(Y) - 1
 p <- length(B)
 r2 <- p - 2
 r1 <- r - r2
 
 beta1 <- B[1]
 beta2 <- B[2]
 
 xi <- B[3:p]
 xizero <- rep(0, r1)
 xinew <- c(xizero, xi)
 xi.mat <- matrix(xinew, n, r, byrow=TRUE)
 
 lam.vec <- lam.fun(beta1, r)[-1]
 g.ee.2 <- matrix((lam.vec - 1) * beta2, n, r, byrow = TRUE)
 g.ee <- Y[, -1] + g.ee.2 - Y[, 1] %*% t(lam.vec) - xi.mat
 
 return(g.ee)
}

#### LD PPEL objective function #####################################

ee.lambda <- function(lambda, g.ee, xi, nu, tau)
{
  ## nu is tuning parameter for lambda, tau is tuning parameter for B
  lambda <- matrix(lambda, ncol = 1)  ##  r by 1 matrix
  g.ee <- as.matrix(g.ee)
  
  n <- nrow(g.ee)
  r <- ncol(g.ee)
  r2 <- length(xi)
  r1 <- r - r2
  
  lambda.g <- 1.0 + g.ee %*% lambda  ##  1+ lambda^t g(x,y,B)
  log.g <- log.star(lambda.g)
  obj.fun.1 <- sum(log.g)
  
  location <- c(rep(0,r1), rep(1,r2))
  lambda2 <- lambda[which(location == 1)]
  scad.lambda <- scad.quad(lambda2, nu)$scad.val
  obj.fun.2 <- n *  scad.lambda
  
  scad.xi <- scad.quad(xi, tau)$scad.val
  obj.fun.3 <- n * scad.xi
  
  obj.val <- obj.fun.1 - obj.fun.2 + obj.fun.3
  
  return(obj.val)
}

#### EL objective function ##################################3

ee.fun <- function(lambda, g.ee)
{ 
  lambda <- matrix(lambda, ncol = 1) 
  g.ee <- as.matrix(g.ee)
  
  lambda.g <- 1.0 + g.ee %*% lambda 
  log.g <- log.star(lambda.g)
  obj.val <- sum(log.g)
  
  return(obj.val)
}

##########grad.lambda#################

grad.lambda <- function(index.lambda, lambda, g.ee, nu, r1) 
{ 
  ## grad vector respect to lambda
  lambda <- matrix(lambda, ncol=1)
  g.ee <- as.matrix(g.ee)
  r <- ncol(g.ee)
  
  epsilon <- 1e-10
  
  scad.lambda <- scad.quad(beta = lambda[index.lambda], lambda = nu)
  scad.der <- scad.lambda$scad.der
  scad.hess <- scad.lambda$scad.hessian
  
  g.lambda <- 1.0 + g.ee %*% lambda
  
  sub.g <- g.ee[ , index.lambda]
  
  if(index.lambda < (r1+1))
  {
    g.grad <- sum(sub.g * der.log.star(g.lambda))
    g.hessian <-  sum(sub.g^2 * der2.log.star(g.lambda))  # not used
  }
  else
  {
    g.grad <- sum(sub.g * der.log.star(g.lambda)) - n * scad.der
    g.hessian <- sum(sub.g^2 * der2.log.star(g.lambda)) - n * scad.hess  # not used
  }
  
  if (abs(g.hessian) <= epsilon){
    g.hessian <- g.hessian + epsilon
  }
  
  return(list(g.grad = g.grad, g.hessian = g.hessian))
}

#######grad.B#############################################

grad.B <- function(B, supp.B, lambda, Y, tau)
{
  lambda <- matrix(lambda, ncol = 1)
  r <- ncol(Y) - 1
  p <- length(B)
  r2 <- p - 2
  r1 <- r - r2
  
  g.ee <- auxi.fun(B, Y)
  
  B <- matrix(B, ncol = 1)
  beta1 <- B[1, 1]
  beta2 <- B[2, 1]
  
  epsilon <- 1e-10
  
  len.B <- length(supp.B)
  grad.vect <- matrix(0.0, nrow = len.B)
  grad2.vect <-  grad.vect
  
  for(i in c(1:len.B))
  {
    index.B <- supp.B[i]
    lambda.g <- 1.0 + g.ee %*% lambda
    if(index.B == 1)
    {
      lambda.partial.g <- (beta2 - Y[, 1]) %*% t(der.lam.fun(beta1, r)) %*% lambda
    }
    if(index.B == 2)
    {
      lambda.partial.g <- matrix(lam.fun(beta1, r)[-1] - 1, n, r, byrow = T) %*% lambda
    }
    if(index.B > 2)
    {
      partial.g <- rep(0, r)
      partial.g[(index.B - 2 + r1)] <- -1
      partial.g <- rep(partial.g, n)
      partial.g <- matrix(partial.g, ncol=r, byrow=TRUE)
      lambda.partial.g <- partial.g %*% lambda  
    }
    grad.val <- sum(der.log.star(lambda.g) * lambda.partial.g)
    hess.val<-sum(der2.log.star(lambda.g)*(lambda.partial.g^2))  # revised later
    
    if(index.B > 2)
    {
      scad.B <- scad.quad(B[index.B], tau)
      scad.der <- scad.B$scad.der
      scad.hessian <- scad.B$scad.hessian
      grad.val <- grad.val + n * scad.der
      hess.val <- hess.val+n* scad.hessian
    }
    
    if (abs(hess.val) <= epsilon)
    {
      hess.val <- hess.val + epsilon
    }
    
    grad.vect[i] <- grad.val
    grad2.vect[i] <- hess.val
  }
  
  return(list(grad.vect = grad.vect, grad2.vect = grad2.vect))
}

#########bias corr##################################################

# for simulation
biacorr <- function(B, auxi, para, supp.B, s)
{
  lambda <- para$lambda
  lambda <- matrix(lambda, ncol = 1)
  supp.lambda <- para$supp.lambda
  
  beta1 <- B[1]
  beta2 <- B[2]
  
  Y <- para$Y
  Y <- as.matrix(Y)
  
  g.mat <- auxi
  g.mat <- as.matrix(g.mat)
  
  n <- nrow(Y)
  r <- ncol(Y) - 1
  p <- length(B)
  r2 <- p - 2
  r1 <- r - r2
  
  epsilon <- 0.05
  
  supp.B.2 <- supp.B
  len.B.2 <- length(supp.B.2)
  len.lambda <- length(supp.lambda)
  
  grad.vect <- matrix(0.0, nrow = len.lambda, ncol = len.B.2)
  for(i in c(1:len.B.2))
  {
    index.B <- supp.B.2[i]
    if(index.B == 1)
    {
      grad.mat <- (beta2 - Y[, 1]) %*% t(der.lam.fun(beta1, r)[supp.lambda])
    }
    if(index.B == 2)
    {
      grad.mat <- matrix(lam.fun(beta1, r)[-1] - 1, n, r, byrow = T)[ , supp.lambda]
    }
    if(index.B > 2)
    {
      partial.g<-rep(0,r)
      partial.g[(index.B - 2 + r1)] <- -1
      partial.g<-rep(partial.g,n)
      partial.g<-matrix(partial.g, ncol=r, byrow=TRUE)
      grad.mat<-partial.g[ , supp.lambda]
    }
    grad.vect[,i] <-  colMeans(grad.mat)
  }
  
  temp.vect <- g.mat %*% lambda
  #temp.vect[which(abs(temp.vect) >= 1e-6)] <- 1e-6 * sign(temp.vect[which(abs(temp.vect) >= 1e-6)])
  denom.vect <- 1.0 +  temp.vect
  
  temp.mat <- t(g.mat[, supp.lambda]) %*% g.mat[, supp.lambda]/n
  eig <- eigen(temp.mat)$values
  if (min(eig) <= max(eig) * epsilon) temp.mat <- temp.mat + epsilon * diag(len.lambda)
  Omega <- solve(temp.mat)
  
  Phi <- t(grad.vect) %*% Omega %*% grad.vect
  dim.phi <- nrow(Phi)
  svd.val <- svd(Phi)$d  #???
  if (min(abs(svd.val)) <= epsilon*max(abs(svd.val)))
  {
    Phi <- Phi + epsilon * diag(dim.phi)
  }
  
  eta <- colMeans(solve(diag(c(denom.vect), n)) %*% g.mat[, supp.lambda])
  eta[1:r1] <- 0
  
  Psi <- solve(Phi) %*% t(grad.vect) %*% Omega %*% eta
  
  ### 渐近方差
  B0 <- c(B[1:2], rep(0, p-2))
  g.mat <- auxi.fun(B0, Y)
  g.mat <- as.matrix(g.mat)
  
  A <- c((r1+1):(r-s))  # 真值
  Ip <- union(1:r1, intersect(supp.lambda, A))
  
  len.Ip <- length(Ip)
  grad.vect.Ip <- matrix(0, nrow = len.Ip, ncol = 2)
  for(i in c(1:2))
  {
    index.B <- i
    if(index.B == 1)
    {
      grad.mat <- (beta2 - Y[, 1]) %*% t(der.lam.fun(beta1, r)[Ip])
    }
    if(index.B == 2)
    {
      grad.mat <- matrix(lam.fun(beta1, r)[-1] - 1, n, r, byrow = T)[ , Ip]  
    }
    grad.vect.Ip[,i] <-  colMeans(grad.mat)
  }
  
  V.Ip <- t(g.mat[, Ip]) %*% g.mat[, Ip]/n
  eig <- eigen(V.Ip)$values
  if (min(eig) <= max(eig) * epsilon) V.Ip <- V.Ip + epsilon * diag(len.Ip)
  V.Ip.inv <- solve(V.Ip)
  
  J.Ip <- t(grad.vect.Ip) %*% V.Ip.inv %*% grad.vect.Ip
  eig <- eigen(J.Ip)$values
  if (min(eig) <= max(eig) * epsilon) J.Ip <- J.Ip + epsilon * diag(2)
  
  res <- eigen(J.Ip)
  eigen.value <- res$values
  eigen.vector <- res$vectors
  sqrt.J.Ip <- eigen.vector %*% diag(sqrt(eigen.value)) %*% solve(eigen.vector)
  b <- c(1, 0)
  alpha <- solve(sqrt.J.Ip, b)
  alpha.norm <- sqrt(sum(alpha^2))
  betax.astd <- alpha.norm / sqrt(n)
  
  return(list(Psi=Psi, betax.astd=betax.astd))
}

##########optim.lambda#####################################

#梯度上升，(不)考虑hessian
optim.lambda <- function(lambda, supp.lambda, g.ee, nu, r1)  #iter.max
{
  iter.max <- 30  #convergence????? 30 60 90
  
  for (iter.num in c(1:iter.max))
  {
    for (index in supp.lambda)
    {
      driv.lambda <- grad.lambda(index, lambda, g.ee, nu, r1) 
      driv.1 <- driv.lambda$g.grad
      driv.2 <- driv.lambda$g.hessian
      
      temp.val <- driv.1 
      #temp.val <- driv.1 / driv.2
      if (abs(temp.val) > 0.001)  temp.val <- sign(temp.val) * 0.001
      
      #lambda[index] <- lambda[index] - temp.val
      lambda[index] <- lambda[index] + temp.val
    }
  }
  
  return(lambda)
}

##########optim.B#######################################

optim.B2 <- function(B, supp.B, lambda, tau, nu, Y) #iter.max
{
  lambda <- matrix(lambda, ncol = 1)
  r <- ncol(Y) - 1
  p <- length(B)
  r2 <- p - 2
  r1 <- r - r2
  
  iter.max <- 20 #20 40 60
  
  for (iter.num in c(1:iter.max))
    {
    driv.B <- grad.B(B, supp.B, lambda, Y, tau)
    driv.1 <- driv.B$grad.vect
    driv.2 <- driv.B$grad2.vect
    
    temp.vect <- driv.1
    # temp.vect <- driv.1/driv.2
    temp.vect[abs(temp.vect) > 0.001] <- sign(temp.vect[abs(temp.vect) > 0.001])*0.001
    
    B[supp.B] <- B[supp.B] - temp.vect
    #B[supp.B] <- B[supp.B] + temp.vect
    
    g.ee <- auxi.fun(B, Y)
    supp.lambda <- 1:length(lambda)
    max.lambda <- optim.lambda(lambda, supp.lambda, g.ee, nu, r1) ### Coordinate descent for lambda
    lambda <- max.lambda
  }
  
  return(list(B = B, lambda = lambda))
}

#############main.iter####################################################

main.iter.tun <- function(tau, nu, init.B, Y) #iter.max
{
  n <- nrow(Y)
  r <- ncol(Y) - 1
  p <- length(init.B)
  r2 <- p - 2
  r1 <- r - r1
  
  eps.tol <- 0.005 
  
  init.lambda <- matrix(0.0, nrow = r, ncol = 1) #?
  lambda <- init.lambda
  
  iter.max <- 30  #30 60 90
  
  B <- init.B
  supp.B <- 1:p
  supp.lambda <- 1:r
  g.ee <- auxi.fun(B, Y)
  
  max.lambda <- optim.lambda(lambda, supp.lambda, g.ee, nu, r1) ### Coordinate descent for lambda
  lambda <- max.lambda
  
  for(count.num in 1:iter.max)
  {
    
    res <- optim.B2(B, supp.B, lambda, tau, nu, Y) ### Optimization function for B by direct derivate
    
    min.B <- res$B
    lambda <- res$lambda
    
    if(count.num > 5) # penalize part of parameter >5 20 60
    {
      temp.supp <- c(which(abs(min.B[(2+1):p]) > eps.tol))
      if (length(temp.supp) > 0)
      {
        xi.supp <- c(which(abs(min.B[(2+1):p]) > eps.tol)) + 2
        supp.B <- c(1:2, xi.supp)
      }
      else supp.B <- c(1:2)
    }
    
    B < -min.B
    B[-supp.B] <- 0
  }
  
  g.ee <- auxi.fun(B, Y)
  
  return(list(B =  B, lambda = lambda, g.ee = g.ee))
}

# for simulation
main.iter <- function(tau, nu, init.B, Y, s) #iter.max
{
  n <- nrow(Y)
  r <- ncol(Y) - 1
  p <- length(init.B)
  r2 <- p - 2
  r1 <- r - r2
  
  eps.tol <- 0.005 #????? #0.05 0.10 0.15 0.005 0.01
  
  init.lambda <- matrix(0, nrow = r, ncol = 1) #?
  lambda <- init.lambda
  
  iter.max <- 30  #30 60 90
  
  B <- init.B
  supp.B <- 1:p
  g.ee <- auxi.fun(B, Y)
  supp.lambda <- 1:r
  
  #记录迭代路径
  iter.num <- iter.max + 1
  B.iter <- matrix(0, iter.num, p)
  B.iter[1,] <- init.B
  
  max.lambda <- optim.lambda(lambda, supp.lambda, g.ee, nu, r1) ### Coordinate descent for lambda
  lambda <- max.lambda
  
  for(count.num in 1:iter.max)
  {
    res <- optim.B2(B, supp.B, lambda, tau, nu, Y) ### Optimization function for B by direct derivate
    
    min.B <- res$B
    lambda <- res$lambda
    
    if(count.num > 5) # penalize part of parameter >5 20 60
    {
      temp.supp <- c(which(abs(min.B[(2+1):p]) > eps.tol))
      if (length(temp.supp) > 0)
      {
        xi.supp <- c(which(abs(min.B[(2+1):p]) > eps.tol)) + 2
        supp.B <- c(1:2, xi.supp)
      }
      else supp.B <- c(1:2)
    }
    
    B <- min.B
    B[-supp.B] <- 0
    
    #记录迭代路径
    B.iter[(count.num + 1), ] <- B
  }
  
  auxi <- auxi.fun(B, Y)
  order.g <- order(abs(lambda), decreasing = TRUE)
  len.B <- length(supp.B)
  lambda.num <- max(round(n/log(n)),len.B+1)
  g.supp <- order.g[1:lambda.num]
  supp.lambda <- union(c(1:r1), g.supp)
  
  
  para <- list(Y = Y, lambda = lambda, supp.lambda = supp.lambda, tau = tau, nu = nu)
  
  bias <- rep(0, p)
  bias.supp.B <- biacorr(B, auxi, para, supp.B, s)
  bias[supp.B] <- bias.supp.B$Psi
  
  B.astd <- rep(0, p)
  B.astd[1] <- bias.supp.B$betax.astd
  
  biacor_B <- B - bias
  
  return(list(B =  B, biacor_B = biacor_B, B.iter = B.iter, bias = bias, B.astd= B.astd))
}

