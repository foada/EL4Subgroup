############################################################
######              RMSE of beta_z in Table S1          #####
############################################################
rep.num<-500 
beta <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol=1)

setwd("C:/Users/Administrator/Desktop/HDMS/Rcode and sim res/sim_res_linear_bic/0.7_0.9")

n <- 200
dw <- 240
dw2 <- 17

file <- paste("2SLS_HDMS_esti_n=", n, "dw=", dw, "s=", dw2, ".csv", sep = "")

HDMS.esti<-read.csv(file, header=TRUE) #¸Ä

beta.mat <- matrix(0, nrow = rep.num, ncol = 6)

for (i in 1:rep.num){
  esti.mat <- HDMS.esti[,(5*(i-1)+2):(5*(i-1)+6)]
  esti.mat <- as.matrix(esti.mat)
  
  # for beta1
  PEL2.ebic <- esti.mat[, 1]
  PEL2.ebic.beta1 <- PEL2.ebic[3]
  beta.mat[i,1] <- PEL2.ebic.beta1
  
  biacor_PEL2.ebic <- esti.mat[, 2]
  biacor_PEL2.ebic.beta1 <- biacor_PEL2.ebic[3]
  beta.mat[i,2] <- biacor_PEL2.ebic.beta1
  
  init.B <- esti.mat[, 3]
  init.B.beta1 <- init.B[3]
  beta.mat[i,3] <- init.B.beta1
  
  # for beta2
  PEL2.ebic.beta2 <- PEL2.ebic[4]
  beta.mat[i,4] <- PEL2.ebic.beta2
  
  biacor_PEL2.ebic.beta2 <- biacor_PEL2.ebic[4]
  beta.mat[i,5] <- biacor_PEL2.ebic.beta2
  
  init.B.beta2 <- init.B[4]
  beta.mat[i,6] <- init.B.beta2
}

RMSE <- sqrt(colMeans((beta.mat-beta[1])^2))
bias <- colMeans(beta.mat-beta[1])
std <- sqrt(RMSE^2 - bias^2)

error <- t(rbind(RMSE, bias, std))
print(error)

MSE <- RMSE^2
RMSE2 <- rep(0, 3)
RMSE2[1] <- sqrt(MSE[1] + MSE[4])
RMSE2[2] <- sqrt(MSE[2] + MSE[5])
RMSE2[3] <- sqrt(MSE[3] + MSE[6])
print(RMSE2)

write.csv(RMSE2, file = paste("ci2_n=", n, "dw=", dw, "s=", dw2, ".csv", sep = ""))
