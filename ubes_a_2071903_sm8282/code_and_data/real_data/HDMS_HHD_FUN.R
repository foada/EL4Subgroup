#############pseudo log##########################

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

scad.grad.2 <- function(beta, lambda){
  ## where lambda is tuning parameter
  a <- 3.7 ## a = 3.7 suggested by Fan and Li(2001)
  beta <- as.matrix(beta, ncol=1)
  abs.beta <- abs(beta)
  
  #temp.1 <- lambda  ## for the case abs.beta <= lambda
  temp.1 <- lambda*sign(beta) #改
  #temp.2 <- (a*lambda - abs.beta)/(a-1.0) ## for the case abs.beta <= a*lambda
  #temp.2 <- max( a*lambda - abs.beta,0)/(a-1.0)*sign(beta) #改
  temp.2 <-  (a*lambda - abs.beta)/(a-1.0)*sign(beta)
  
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

### Extended moment function #################################

auxi.fun<-function(B, y, x, z, w) 
{
  dz<-ncol(z)
  dw<-ncol(w)
  r<-dz+dw
  B <- matrix(B, ncol = 1)
  betax<-B[1,]
  betaz<-B[2:(dz+1),]
  xi<-B[(dz+2):r,]
  
  betaz<-matrix(betaz,ncol=1)
  xizero<-rep(0,dz+1)
  xinew<-c(xizero,xi)
  xinew<-rep(xinew,n)
  xinew<-matrix(xinew,nrow=n,byrow=TRUE)
  
  residual<-y-x*betax-z%*%betaz
  residual<-as.vector(residual)
  residual<-diag(residual)
  
  zw<-cbind(z,w)
  g.ee<-residual%*%zw-xinew
  return(list(g.ee=g.ee))
}

#### LD PPEL objective function #####################################

ee.lambda <- function(lambda, g.ee, xi, nu, tau)
{ ## nu is tuning parameter for lambda, tau is tuning parameter for B
  lambda <- matrix(lambda, ncol = 1) ##  (dz+dw) by 1 matrix
  g.ee <- as.matrix(g.ee)
  
  n <- nrow(g.ee)
  r<-ncol(g.ee)
  r2<-length(xi)
  r1<-r-r2
  
  lambda.g <- 1.0 + g.ee %*% lambda ##  1+ lambda^t g(x,y,B)
  log.g <- log.star(lambda.g)
  obj.fun.1 <- sum(log.g)
  
  location<-c(rep(0,r1),rep(1,r2))
  lambda2<-lambda[which(location==1)]
  scad.lambda <- scad.quad(lambda2, nu)$scad.val
  obj.fun.2 <- n *  scad.lambda
  
  scad.xi <- scad.quad(xi, tau)$scad.val
  obj.fun.3 <- n * scad.xi
  
  obj.val <- obj.fun.1 - obj.fun.2 + obj.fun.3
  
  return(obj.val)
}

##########grad.lambda#################

#HD, with penalty
grad.lambda <- function(index.lambda, lambda, g.ee, nu) #add dz or xi
{ ## grad vector respect to lambda
  lambda <- matrix(lambda, ncol=1)
  
  g.ee <- as.matrix(g.ee)
  r <- ncol(g.ee)
  #r2<-length(xi)
  dz<-3 #??????????
  r1<-dz+1
  #n <- nrow(g.ee)
  
  epsilon <- 1e-10
  
  scad.lambda <- scad.quad(beta=lambda[index.lambda], lambda=nu)
  scad.der <- scad.lambda$scad.der
  scad.hess <- scad.lambda$scad.hessian
  
  g.lambda <- 1.0 + g.ee %*% lambda
  
  sub.g <- g.ee[ , index.lambda]
  
  if(index.lambda<(r1+1))
  {
    #g.grad <- sum(sub.g / g.lambda) # ?why not use log.star
    g.grad <- sum(sub.g * der.log.star(g.lambda))
    
    #g.hessian <- (-1.0) * sum(sub.g^2 / g.lambda^2) 
    g.hessian <-  sum(sub.g^2 * der2.log.star(g.lambda)) 
  }
  else
  {
    #g.grad <- sum(sub.g / g.lambda)- n * scad.der
    g.grad <- sum(sub.g * der.log.star(g.lambda))- n * scad.der
    #g.grad <- sum(sub.g * der.log.star(g.lambda))
    
    #g.hessian <- (-1.0) * sum(sub.g^2 / g.lambda^2) - n * scad.hess
    g.hessian <- sum(sub.g^2 * der2.log.star(g.lambda)) - n * scad.hess
    #g.hessian <- sum(sub.g^2 * der2.log.star(g.lambda))
  }
  
  if (abs(g.hessian) <= epsilon){
    g.hessian <- g.hessian + epsilon
  }
  
  return(list(g.grad = g.grad, g.hessian = g.hessian))
}

#### EL objective function ##################################3

ee.fun <- function(lambda, g.ee)
{ 
  lambda <- matrix(lambda, ncol = 1) 
  g.ee <- as.matrix(g.ee)
  
  lambda.g <- 1.0 + g.ee %*% lambda ##  1+ lambda^t g(x,y,B)
  log.g <- log.star(lambda.g)
  obj.val <- sum(log.g)
  
  return(obj.val)
}


#######grad.B#############################################

grad.B <- function(B, supp.B, lambda, y, x, z, w, tau)
{
  lambda <- matrix(lambda, ncol = 1)
  dz<-ncol(z)
  dw<-ncol(w)
  r<-dz+dw
  r1<-dz+1
  
  auxi <- auxi.fun(B,y,x,z,w)
  g.ee <- auxi$g.ee
  
  B<-matrix(B, ncol = 1)
  
  epsilon <- 1e-10
  
  len.B <- length(supp.B)
  
  grad.vect <- matrix(0.0, nrow = len.B)
  grad2.vect <-  grad.vect
  
  for(i in c(1:len.B))
  {
    index.B <- supp.B[i]
    lambda.g <- 1.0 + g.ee %*% lambda
    if(index.B==1)
    {
      lambda.partial.g<-diag(as.vector(x*(-1)))%*%cbind(z,w)%*%lambda
    }
    if(index.B>1 & index.B<(dz+2))
    {
      lambda.partial.g<-diag(as.vector(z[,(index.B-1)]*(-1)))%*%cbind(z,w)%*%lambda  
    }
    if(index.B>(dz+1))
    {
      partial.g<-rep(0,r)
      partial.g[index.B]<--1
      partial.g<-rep(partial.g,n)
      partial.g<-matrix(partial.g,ncol=r,byrow=TRUE)
      lambda.partial.g<-partial.g%*%lambda  
    }
    grad.val<-sum(der.log.star(lambda.g)*lambda.partial.g)
    hess.val<-sum(der2.log.star(lambda.g)*(lambda.partial.g^2))
    
    if(index.B>r1)
    {
      scad.B <- scad.quad(B[index.B], tau)
      scad.der <- scad.B$scad.der
      scad.hessian <- scad.B$scad.hessian
      grad.val <- grad.val + n * scad.der
      hess.val<-hess.val+n* scad.hessian
    }
    
    ########################################
    
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

#for real data
biacorr <- function(B, auxi, para, supp.B)
{
  lambda <- para$lambda
  lambda <- matrix(lambda, ncol = 1)
  supp.lambda <- para$supp.lambda
  
  x <- para$x
  x <- as.matrix(x)
  
  y <- para$y
  y <- as.matrix(y)
  
  z <- para$z
  z <- as.matrix(z)
  
  w <- para$w
  w <- as.matrix(w)
  
  zw<-cbind(z,w)
  
  g.mat <- auxi$g.ee
  g.mat <- as.matrix(g.mat)
  
  n <- nrow(x)
  dz<-ncol(z)
  r1<-dz+1
  dw<-ncol(w)
  r<-dz+dw
  
  epsilon <- 0.05
  
  #supp.B.2 <- intersect(supp.B,supp.lambda)
  supp.B.2 <- supp.B
  
  len.B.2 <- length(supp.B.2)
  len.lambda <- length(supp.lambda)
  
  grad.vect <- matrix(0.0, nrow = len.lambda, ncol = len.B.2)
  
  for(i in c(1:len.B.2))
  {
    index.B <- supp.B.2[i]
    if(index.B==1)
    {
      grad.mat<-diag(as.vector(x*(-1)))%*%zw[ ,supp.lambda]
    }
    if(index.B>1 & index.B<(dz+2))
    {
      grad.mat<-diag(as.vector(z[,(index.B-1)]*(-1)))%*%zw[ ,supp.lambda]  
    }
    if(index.B>(dz+1))
    {
      partial.g<-rep(0,r)
      partial.g[index.B]<--1
      partial.g<-rep(partial.g,n)
      partial.g<-matrix(partial.g,ncol=r,byrow=TRUE)
      grad.mat<-partial.g[ ,supp.lambda]
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
  
  #acov <- ginv(Phi)
  
  Psi <- solve(Phi) %*% t(grad.vect) %*% Omega %*% eta
  
  #渐近方差
  B0 <- c(B[1:r1],rep(0,r-r1))
  auxi <- auxi.fun(B0,y,x,z,w)
  g.mat <- auxi$g.ee
  g.mat <- as.matrix(g.mat)
  
  #A <- c(5:(r-dw2))
  #Ip <- union(1:r1,intersect(supp.lambda, A))
  A.hat <- setdiff(c(1:r),supp.B.2)
  Ip <- union(1:r1,intersect(supp.lambda, A.hat))
  
  
  len.Ip <- length(Ip)
  grad.vect.Ip <- matrix(0.0, nrow = len.Ip, ncol = r1)
  for(i in c(1:r1))
  {
    index.B <- i
    if(index.B==1)
    {
      grad.mat<-diag(as.vector(x*(-1)))%*%zw[ ,Ip]
    }
    if(index.B>1 & index.B<(dz+2))
    {
      grad.mat<-diag(as.vector(z[,(index.B-1)]*(-1)))%*%zw[ ,Ip]  
    }
    
    grad.vect.Ip[,i] <-  colMeans(grad.mat)
  }
  V.Ip <- t(g.mat[, Ip]) %*% g.mat[, Ip]/n
  
  eig <- eigen(V.Ip)$values
  if (min(eig) <= max(eig) * epsilon) V.Ip <- V.Ip + epsilon * diag(len.Ip)
  
  V.Ip.inv <- solve(V.Ip)
  J.Ip <- t(grad.vect.Ip) %*% V.Ip.inv %*% grad.vect.Ip
  
  eig <- eigen(J.Ip)$values
  if (min(eig) <= max(eig) * epsilon) J.Ip <- J.Ip + epsilon * diag(r1)
  
  res <- eigen(J.Ip)
  eigen.value <- res$values
  eigen.vector <- res$vectors
  sqrt.J.Ip <- eigen.vector%*%diag(sqrt(eigen.value))%*%solve(eigen.vector)
  b <-c(1,rep(0,(r1-1)))
  alpha <- solve(sqrt.J.Ip, b)
  alpha.norm <- sqrt(sum(alpha^2))
  betax.astd <- alpha.norm/sqrt(n)
  
  return(list(Psi=Psi, betax.astd=betax.astd))
}

# for simulation
biacorr.1 <- function(B, auxi, para, supp.B, dw2)
{
  lambda <- para$lambda
  lambda <- matrix(lambda, ncol = 1)
  supp.lambda <- para$supp.lambda
  
  x <- para$x
  x <- as.matrix(x)
  
  y <- para$y
  y <- as.matrix(y)
  
  z <- para$z
  z <- as.matrix(z)
  
  w <- para$w
  w <- as.matrix(w)
  
  zw<-cbind(z,w)
  
  g.mat <- auxi$g.ee
  g.mat <- as.matrix(g.mat)
  
  n <- nrow(x)
  dz<-ncol(z)
  r1<-dz+1
  dw<-ncol(w)
  r<-dz+dw
  
  epsilon <- 0.05
  
  #supp.B.2 <- intersect(supp.B,supp.lambda)
  supp.B.2 <- supp.B
  #dw2 <- 3 #log(n)
  #supp.B.2 <- c(1:r1, (r-dw2+1):r)
  
  len.B.2 <- length(supp.B.2)
  len.lambda <- length(supp.lambda)
  
  grad.vect <- matrix(0.0, nrow = len.lambda, ncol = len.B.2)
  
  for(i in c(1:len.B.2))
  {
    index.B <- supp.B.2[i]
    if(index.B==1)
    {
      grad.mat<-diag(as.vector(x*(-1)))%*%zw[ ,supp.lambda]
    }
    if(index.B>1 & index.B<(dz+2))
    {
      grad.mat<-diag(as.vector(z[,(index.B-1)]*(-1)))%*%zw[ ,supp.lambda]  
    }
    if(index.B>(dz+1))
    {
      partial.g<-rep(0,r)
      partial.g[index.B]<--1
      partial.g<-rep(partial.g,n)
      partial.g<-matrix(partial.g,ncol=r,byrow=TRUE)
      grad.mat<-partial.g[ ,supp.lambda]
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
  
  #acov <- ginv(Phi)
  
  Psi <- solve(Phi) %*% t(grad.vect) %*% Omega %*% eta
  
  #渐近方差
  B0 <- c(B[1:r1],rep(0,r-r1))
  auxi <- auxi.fun(B0,y,x,z,w)
  g.mat <- auxi$g.ee
  g.mat <- as.matrix(g.mat)
  
  #dw2<-round(log(n))    # dw2
  #dw1<-dw-dw2
  A <- c(5:(r-dw2))
  Ip <- union(1:r1,intersect(supp.lambda, A))
  
  len.Ip <- length(Ip)
  grad.vect.Ip <- matrix(0.0, nrow = len.Ip, ncol = r1)
  for(i in c(1:r1))
  {
    index.B <- i
    if(index.B==1)
    {
      grad.mat<-diag(as.vector(x*(-1)))%*%zw[ ,Ip]
    }
    if(index.B>1 & index.B<(dz+2))
    {
      grad.mat<-diag(as.vector(z[,(index.B-1)]*(-1)))%*%zw[ ,Ip]  
    }
    
    grad.vect.Ip[,i] <-  colMeans(grad.mat)
  }
  V.Ip <- t(g.mat[, Ip]) %*% g.mat[, Ip]/n
  
  eig <- eigen(V.Ip)$values
  if (min(eig) <= max(eig) * epsilon) V.Ip <- V.Ip + epsilon * diag(len.Ip)
  
  V.Ip.inv <- solve(V.Ip)
  J.Ip <- t(grad.vect.Ip) %*% V.Ip.inv %*% grad.vect.Ip
  
  eig <- eigen(J.Ip)$values
  if (min(eig) <= max(eig) * epsilon) J.Ip <- J.Ip + epsilon * diag(r1)
  
  res <- eigen(J.Ip)
  eigen.value <- res$values
  eigen.vector <- res$vectors
  sqrt.J.Ip <- eigen.vector%*%diag(sqrt(eigen.value))%*%solve(eigen.vector)
  b <-c(1,rep(0,(r1-1)))
  alpha <- solve(sqrt.J.Ip, b)
  alpha.norm <- sqrt(sum(alpha^2))
  betax.astd <- alpha.norm/sqrt(n)
  
  return(list(Psi=Psi, betax.astd=betax.astd))
}

#真值
biacorr.2 <- function(B, auxi, para, supp.B)
{
  lambda <- para$lambda
  lambda <- matrix(lambda, ncol = 1)
  supp.lambda <- para$supp.lambda
  
  x <- para$x
  x <- as.matrix(x)
  
  y <- para$y
  y <- as.matrix(y)
  
  z <- para$z
  z <- as.matrix(z)
  
  w <- para$w
  w <- as.matrix(w)
  
  zw<-cbind(z,w)
  
  n <- nrow(x)
  dz<-ncol(z)
  r1<-dz+1
  dw<-ncol(w)
  r<-dz+dw
  
  #带入真值计算g.ee
  dw2<-round(log(n))    # dw2
  dw1<-dw-dw2
  
  delta.min <- 0.5
  delta.max <- 0.8 
  delta <- seq(delta.min,delta.max,len=dw2)
  xi0 <- delta
  #B0 <- c(B[1:r1],rep(0,r-r1-dw2),xi0)
  B0 <- c(rep(0.5,r1),rep(0,r-r1-dw2),xi0)
  auxi <- auxi.fun(B0,y,x,z,w)
  
  
  g.mat <- auxi$g.ee
  g.mat <- as.matrix(g.mat)
  
  # n <- nrow(x)
  # dz<-ncol(z)
  # r1<-dz+1
  # dw<-ncol(w)
  # r<-dz+dw
  
  epsilon <- 0.05
  
  supp.B0 <- c(1:r1,(r-dw2+1):r)
  supp.B.2 <- intersect(supp.B0,supp.lambda)
  #supp.B.2 <- intersect(supp.B,supp.lambda)
  #supp.B.2 <- supp.B
  #dw2 <- 3 #log(n)
  #supp.B.2 <- c(1:r1, (r-dw2+1):r)
  
  len.B.2 <- length(supp.B.2)
  len.lambda <- length(supp.lambda)
  
  grad.vect <- matrix(0.0, nrow = len.lambda, ncol = len.B.2)
  
  for(i in c(1:len.B.2))
  {
    index.B <- supp.B.2[i]
    if(index.B==1)
    {
      grad.mat<-diag(as.vector(x*(-1)))%*%zw[ ,supp.lambda]
    }
    if(index.B>1 & index.B<(dz+2))
    {
      grad.mat<-diag(as.vector(z[,(index.B-1)]*(-1)))%*%zw[ ,supp.lambda]  
    }
    if(index.B>(dz+1))
    {
      partial.g<-rep(0,r)
      partial.g[index.B]<--1
      partial.g<-rep(partial.g,n)
      partial.g<-matrix(partial.g,ncol=r,byrow=TRUE)
      grad.mat<-partial.g[ ,supp.lambda]
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
  
  #acov <- ginv(Phi)
  
  Psi <- solve(Phi) %*% t(grad.vect) %*% Omega %*% eta
  
  #渐近方差
  # B0 <- c(B[1:r1],rep(0,r-r1))
  # auxi <- auxi.fun(B0,y,x,z,w)
  # g.mat <- auxi$g.ee
  # g.mat <- as.matrix(g.mat)
  # 
  # dw2<-round(log(n))    # dw2
  #dw1<-dw-dw2
  
  A <- c(5:(r-dw2))
  Ip <- union(1:r1,intersect(supp.lambda, A))
  
  len.Ip <- length(Ip)
  
  grad.vect.Ip <- matrix(0.0, nrow = len.Ip, ncol = r1)
  for(i in c(1:r1))
  {
    index.B <- i
    if(index.B==1)
    {
      grad.mat<-diag(as.vector(x*(-1)))%*%zw[ ,Ip]
    }
    if(index.B>1 & index.B<(dz+2))
    {
      grad.mat<-diag(as.vector(z[,(index.B-1)]*(-1)))%*%zw[ ,Ip]  
    }
    
    grad.vect.Ip[,i] <-  colMeans(grad.mat)
  }
  
  V.Ip <- t(g.mat[, Ip]) %*% g.mat[, Ip]/n
  
  eig <- eigen(V.Ip)$values
  if (min(eig) <= max(eig) * epsilon) V.Ip <- V.Ip + epsilon * diag(len.Ip)
  
  V.Ip.inv <- solve(V.Ip)
  J.Ip <- t(grad.vect.Ip) %*% V.Ip.inv %*% grad.vect.Ip
  
  eig <- eigen(J.Ip)$values
  if (min(eig) <= max(eig) * epsilon) J.Ip <- J.Ip + epsilon * diag(r1)
  
  res <- eigen(J.Ip)
  eigen.value <- res$values
  eigen.vector <- res$vectors
  sqrt.J.Ip <- eigen.vector%*%diag(sqrt(eigen.value))%*%solve(eigen.vector)
  b <-c(1,rep(0,(r1-1)))
  alpha <- solve(sqrt.J.Ip, b)
  alpha.norm <- sqrt(sum(alpha^2))
  betax.astd <- alpha.norm/sqrt(n)
  
  return(list(Psi=Psi, betax.astd=betax.astd))
}

##########optim.lambda#####################################

#梯度上升，(不)考虑hessian
optim.lambda <- function(lambda, supp.lambda, g.ee, nu)  #iter.max
{
  #lambda.num <- length(supp.lambda)
  iter.max <- 30  #convergence????? 30 60 90
  #n <- nrow(g.ee)
  
  for (iter.num in c(1:iter.max))
  {
    for (index in supp.lambda)
    {
      driv.lambda <- grad.lambda(index, lambda, g.ee, nu) 
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

optim.lambda.3 <- function(lambda, supp.lambda, g.ee, nu)  #iter.max
{
  #lambda.num <- length(supp.lambda)
  iter.max <- 30  #convergence????? 30 60 90
  #n <- nrow(g.ee)
  
  for (iter.num in c(1:iter.max))
  {
    for (index in supp.lambda)
    {
      driv.lambda <- grad.lambda(index, lambda, g.ee, nu) 
      driv.1 <- driv.lambda$g.grad
      driv.2 <- driv.lambda$g.hessian
      
      temp.val <- driv.1 
      #temp.val <- driv.1 / driv.2
      if (abs(temp.val) > 0.001)  temp.val <- sign(temp.val) * 0.001
      
      #lambda[index] <- lambda[index] - temp.val
      lambda[index] <- lambda[index] + temp.val
    }
  }
  
  eps.tol <- 0.005
  supp.zero <- setdiff(which(lambda<eps.tol),c(1:r1))  #r1 needed
  lambda[supp.zero]<-0
  
  return(lambda)
}


##########optim.B#######################################

optim.B2 <- function(B, supp.B, lambda, tau, nu, y ,x, z, w) #iter.max
{
  lambda <- matrix(lambda, ncol = 1)
  dz<-ncol(z)
  r1<-dz+1
  
  #auxi <- auxi.fun(B,y,x,z,w)
  #g.ee <- auxi$g.ee
  #len.lambda <- ncol(g.ee)
  
  # lambda.num <- min(round(n/log(n)/2), len.lambda)
  #lambda.num <- len.lambda
  # lambda.num <- min(round(n/log(n)/1), len.lambda)
  
  #g.bar <- colSums(g.ee)
  
  #order.g <- order(abs(g.bar), decreasing = TRUE)
  #supp.lambda <- order.g[1:lambda.num]
  
  # supp.lambda <- union(supp.lambda,c((len.lambda-2):len.lambda))
  
  #sub.g <- g.ee[, supp.lambda]
  #sub.lambda <- matrix(lambda[supp.lambda], ncol=1)
  # lambda[-supp.lambda]<-0
  
  iter.max <- 20 #20 40 60
  
  # B.iter.2<-matrix(0,nrow=iter.max,ncol=r)
  
  for (iter.num in c(1:iter.max)){
    driv.B <- grad.B(B, supp.B, lambda, y, x, z, w, tau)
    driv.1 <- driv.B$grad.vect
    driv.2 <- driv.B$grad2.vect
    
    temp.vect <- driv.1
    # temp.vect <- driv.1/driv.2
    temp.vect[abs(temp.vect) > 0.001] <- sign(temp.vect[abs(temp.vect) > 0.001])*0.001
    
    B[supp.B] <- B[supp.B] - temp.vect
    #B[supp.B] <- B[supp.B] + temp.vect
    
    # B.iter.2[iter.num,]<-B
    
    auxi <- auxi.fun(B,y,x,z,w)
    g.ee <- auxi$g.ee
    #g.bar <- colSums(g.ee)
    
    #order.g <- order(abs(g.bar), decreasing = TRUE)
    #supp.lambda <- order.g[1:lambda.num]
    supp.lambda <- 1:length(lambda)
    
    # supp.lambda <- union(supp.lambda,c((len.lambda-2):len.lambda))
    
    #sub.g <- g.ee[, supp.lambda]
    #sub.lambda <- matrix(lambda[supp.lambda], ncol=1)
    #lambda[-supp.lambda]<-0 #??
    
    max.lambda <- optim.lambda(lambda, supp.lambda, g.ee, nu) ### Coordinate descent for lambda
    lambda <- max.lambda
  }
  
  return(list(B = B, lambda = lambda))
}

#############main.iter#####################################################
main.iter.tun <- function(tau, nu, init.B, y, x, z, w) #iter.max
{
  x <- as.matrix(x)
  y <- as.matrix(y)
  z <- as.matrix(z)
  w <- as.matrix(w)
  
  n <- nrow(x)
  dz<-ncol(z)
  r1<-dz+1
  dw<-ncol(w)
  r2<-dw-1
  r<-r1+r2
  
  eps.tol <- 0.005 
  
  #len.lambda <- r
  
  init.lambda <- matrix(0.0, nrow = r, ncol = 1) #?
  lambda <- init.lambda
  
  #lambda.num <- len.lambda
  
  iter.max <- 30  #30 60 90
  
  B <- init.B
  
  supp.B <- 1:r
  
  supp.lambda <- 1:r
  
  auxi <- auxi.fun(B,y,x,z,w)
  g.ee <- auxi$g.ee
  
  max.lambda <- optim.lambda(lambda, supp.lambda, g.ee, nu) ### Coordinate descent for lambda
  lambda <- max.lambda
  
  for(count.num in 1:iter.max)
  {
    
    res <- optim.B2(B, supp.B, lambda, tau, nu, y ,x, z, w)### Optimization function for B by direct derivate
    
    min.B <- res$B
    lambda <- res$lambda
    
    if(count.num > 5) # penalize part of parameter >5 20 60
    {
      temp.supp <- c(which(abs(min.B[(r1+1):r]) > eps.tol))
      if (length(temp.supp) > 0)
      {
        xi.supp <- c(which(abs(min.B[(r1+1):r]) > eps.tol))+r1
        supp.B <- c(1:r1,xi.supp)
      }
      else supp.B <- c(1:r1)
    }
    
    B<-min.B
    B[-supp.B]<-0
  }
  
  auxi <- auxi.fun(B,y,x,z,w)
  g.ee <- auxi$g.ee
  
  return(list(B =  B, lambda = lambda, g.ee = g.ee))
}

# for simulation
main.iter.1 <- function(tau, nu, init.B, y, x, z, w, dw2) #iter.max
{
  x <- as.matrix(x)
  y <- as.matrix(y)
  z <- as.matrix(z)
  w <- as.matrix(w)
  dw2 <- dw2
  
  n <- nrow(x)
  dz<-ncol(z)
  r1<-dz+1
  dw<-ncol(w)
  r2<-dw-1
  r<-r1+r2
  
  eps.tol <- 0.005 #????? #0.05 0.10 0.15 0.005 0.01
  
  init.lambda <- matrix(0.0, nrow = r, ncol = 1) #?
  lambda <- init.lambda
  
  iter.max <- 30  #30 60 90
  
  B <- init.B
  
  supp.B <- 1:r
  
  auxi <- auxi.fun(B,y,x,z,w)
  g.ee <- auxi$g.ee
  
  supp.lambda <- 1:r
  
  #记录迭代路径
  iter.num <- iter.max + 1
  B.iter<-matrix(0, iter.num, r)
  B.iter[1,] <- init.B
  
  max.lambda <- optim.lambda(lambda, supp.lambda, g.ee, nu) ### Coordinate descent for lambda
  lambda <- max.lambda
  
  for(count.num in 1:iter.max)
  {
    res <- optim.B2(B, supp.B, lambda, tau, nu, y ,x, z, w)### Optimization function for B by direct derivate
    
    min.B <- res$B
    lambda <- res$lambda
    
    # if(count.num > 5) # penalze all parameter
    #   {
    #   temp.supp <- c(which(abs(min.B) > eps.tol))
    #   if (length(temp.supp) > 0)
    #     {
    #     supp.B <- c(which(abs(min.B) > eps.tol))  ##
    #     }
    #   }
    
    if(count.num > 5) # penalize part of parameter >5 20 60
    {
      temp.supp <- c(which(abs(min.B[(r1+1):r]) > eps.tol))
      if (length(temp.supp) > 0)
      {
        xi.supp <- c(which(abs(min.B[(r1+1):r]) > eps.tol))+r1
        supp.B <- c(1:r1,xi.supp)
      }
      else supp.B <- c(1:r1)
    }
    
    B<-min.B
    B[-supp.B]<-0
    
    #记录迭代路径
    B.iter[(count.num + 1), ] <- B
  }
  
  auxi <- auxi.fun(B,y,x,z,w)
  #g.ee <- auxi$g.ee
  
  order.g <- order(abs(lambda), decreasing = TRUE)
  len.B <- length(supp.B)
  lambda.num <- max(round(n/log(n)),len.B+1)
  supp.lambda <- union(1:r1,order.g[1:lambda.num])
  
  
  para <- list(y = y, x = x, z=z, w=w, lambda = lambda, supp.lambda = supp.lambda, tau = tau, nu = nu)
  
  bias <- rep(0,r)
  bias.supp.B <- biacorr(B, auxi, para, supp.B, dw2)
  bias[supp.B] <- bias.supp.B$Psi
  
  B.astd <- rep(0,r)
  B.astd[1] <- bias.supp.B$betax.astd
  
  biacor_B <- B - bias
  
  return(list(B =  B, biacor_B = biacor_B, B.iter = B.iter, bias = bias, B.astd= B.astd))
}

#for real data
main.iter <- function(tau, nu, init.B, y, x, z, w) #iter.max
{
  x <- as.matrix(x)
  y <- as.matrix(y)
  z <- as.matrix(z)
  w <- as.matrix(w)
  #dw2 <- dw2
  
  n <- nrow(x)
  dz<-ncol(z)
  r1<-dz+1
  dw<-ncol(w)
  r2<-dw-1
  r<-r1+r2
  
  eps.tol <- 0.005 #????? #0.05 0.10 0.15 0.005 0.01
  
  init.lambda <- matrix(0.0, nrow = r, ncol = 1) #?
  lambda <- init.lambda
  
  iter.max <- 30  #30 60 90
  
  B <- init.B
  
  supp.B <- 1:r
  
  auxi <- auxi.fun(B,y,x,z,w)
  g.ee <- auxi$g.ee
  
  supp.lambda <- 1:r
  
  #记录迭代路径
  iter.num <- iter.max + 1
  B.iter<-matrix(0, iter.num, r)
  B.iter[1,] <- init.B
  
  max.lambda <- optim.lambda(lambda, supp.lambda, g.ee, nu) ### Coordinate descent for lambda
  lambda <- max.lambda
  
  for(count.num in 1:iter.max)
  {
    res <- optim.B2(B, supp.B, lambda, tau, nu, y ,x, z, w)### Optimization function for B by direct derivate
    
    min.B <- res$B
    lambda <- res$lambda
    
    # if(count.num > 5) # penalze all parameter
    #   {
    #   temp.supp <- c(which(abs(min.B) > eps.tol))
    #   if (length(temp.supp) > 0)
    #     {
    #     supp.B <- c(which(abs(min.B) > eps.tol))  ##
    #     }
    #   }
    
    if(count.num > 5) # penalize part of parameter >5 20 60
    {
      temp.supp <- c(which(abs(min.B[(r1+1):r]) > eps.tol))
      if (length(temp.supp) > 0)
      {
        xi.supp <- c(which(abs(min.B[(r1+1):r]) > eps.tol))+r1
        supp.B <- c(1:r1,xi.supp)
      }
      else supp.B <- c(1:r1)
    }
    
    B<-min.B
    B[-supp.B]<-0
    
    #记录迭代路径
    B.iter[(count.num + 1), ] <- B
  }
  
  auxi <- auxi.fun(B,y,x,z,w)
  #g.ee <- auxi$g.ee
  
  order.g <- order(abs(lambda), decreasing = TRUE)
  len.B <- length(supp.B)
  lambda.num <- min(r, max(round(n/log(n)),len.B+1))
  supp.lambda <- union(1:r1,order.g[1:lambda.num])
  
  
  para <- list(y = y, x = x, z=z, w=w, lambda = lambda, supp.lambda = supp.lambda, tau = tau, nu = nu)
  
  bias <- rep(0,r)
  #bias.supp.B <- biacorr(B, auxi, para, supp.B, dw2)
  bias.supp.B <- biacorr(B, auxi, para, supp.B)
  bias[supp.B] <- bias.supp.B$Psi
  
  B.astd <- rep(0,r)
  B.astd[1] <- bias.supp.B$betax.astd
  
  biacor_B <- B - bias
  
  return(list(B =  B, biacor_B = biacor_B, B.iter = B.iter, bias = bias, B.astd= B.astd))
}

main.iter.3 <- function(tau, nu, init.B, y, x, z, w) #iter.max
{
  x <- as.matrix(x)
  y <- as.matrix(y)
  z <- as.matrix(z)
  w <- as.matrix(w)
  #dw2 <- dw2
  
  n <- nrow(x)
  dz<-ncol(z)
  r1<-dz+1
  dw<-ncol(w)
  r2<-dw-1
  r<-r1+r2
  
  eps.tol <- 0.005 #????? #0.05 0.10 0.15 0.005 0.01
  
  init.lambda <- matrix(0.0, nrow = r, ncol = 1) #?
  lambda <- init.lambda
  
  iter.max <- 30  #30 60 90
  
  B <- init.B
  
  supp.B <- 1:r
  
  auxi <- auxi.fun(B,y,x,z,w)
  g.ee <- auxi$g.ee
  
  supp.lambda <- 1:r
  
  #记录迭代路径
  iter.num <- iter.max + 1
  B.iter<-matrix(0, iter.num, r)
  B.iter[1,] <- init.B
  
  max.lambda <- optim.lambda(lambda, supp.lambda, g.ee, nu) ### Coordinate descent for lambda
  lambda <- max.lambda
  
  for(count.num in 1:iter.max)
  {
    res <- optim.B2(B, supp.B, lambda, tau, nu, y ,x, z, w)### Optimization function for B by direct derivate
    
    min.B <- res$B
    lambda <- res$lambda
    
    # if(count.num > 5) # penalze all parameter
    #   {
    #   temp.supp <- c(which(abs(min.B) > eps.tol))
    #   if (length(temp.supp) > 0)
    #     {
    #     supp.B <- c(which(abs(min.B) > eps.tol))  ##
    #     }
    #   }
    
    if(count.num > 5) # penalize part of parameter >5 20 60
    {
      temp.supp <- c(which(abs(min.B[(r1+1):r]) > eps.tol))
      if (length(temp.supp) > 0)
      {
        xi.supp <- c(which(abs(min.B[(r1+1):r]) > eps.tol))+r1
        supp.B <- c(1:r1,xi.supp)
      }
      else supp.B <- c(1:r1)
    }
    
    B<-min.B
    B[-supp.B]<-0
    
    #记录迭代路径
    B.iter[(count.num + 1), ] <- B
  }
  
  auxi <- auxi.fun(B,y,x,z,w)
  #g.ee <- auxi$g.ee
  
  order.g <- order(abs(lambda), decreasing = TRUE)
  len.B <- length(supp.B)
  lambda.num <- max(round(n/log(n)),len.B+1)
  supp.lambda <- union(1:r1,order.g[1:lambda.num])
  
  
  return(list(B =  B, B.iter = B.iter))
}
