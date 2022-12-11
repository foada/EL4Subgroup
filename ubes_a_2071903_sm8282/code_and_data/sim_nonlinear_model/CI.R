
###########################################################################
###  use results from HDMS_HHD_CI to calculate coverage and ci length   ###
##########################################################################

rep.num<-500 #500
beta <- matrix(c(0.5, -2), ncol=1)

setwd("C:/Users/Administrator/Desktop/HDMS/Rcode and sim res/sim_res_nonlinear_bic/0.2")
n <- 200
r <- 240
s <- 17 

file <- paste("Chi_ci_n=", n, "r=", r, "s=", s, ".csv", sep ="")

HDMS.esti<-read.csv(file, header=TRUE) #¸Ä

chi.coverage <- matrix(0, nrow = rep.num, ncol = 3)
chi.ci.len <- matrix(0, nrow = rep.num, ncol = 3)

normal.coverage <- matrix(0, nrow = rep.num, ncol = 3)
normal.ci.len <- matrix(0, nrow = rep.num, ncol = 3)

tsls.coverage <- matrix(0, nrow = rep.num, ncol = 3)
tsls.ci.len <- matrix(0, nrow = rep.num, ncol = 3)

for (i in 1:rep.num){
  ci.mat <- HDMS.esti[,(6*(i-1)+2):(6*(i-1)+7)]
  ci.mat <- as.matrix(ci.mat)
  
  chi.ci <- ci.mat[,1:2]
  normal.ci <- ci.mat[,3:4]
  tsls.ci <- ci.mat[,5:6]
  
  for(j in 1:3)
  {
    if (chi.ci[j,1]<=beta[1]&chi.ci[j,2]>=beta[1]) chi.coverage[i,j]<-1 
    if (normal.ci[j,1]<=beta[1]&normal.ci[j,2]>=beta[1]) normal.coverage[i,j]<-1
    if (tsls.ci[j,1]<=beta[1]&tsls.ci[j,2]>=beta[1]) tsls.coverage[i,j]<-1
    
    
    chi.ci.len[i,j] <- chi.ci[j,2] - chi.ci[j,1]
    normal.ci.len[i,j] <- normal.ci[j,2] - normal.ci[j,1]
    tsls.ci.len[i,j] <- tsls.ci[j,2] - tsls.ci[j,1]
  }
}


res <- matrix(0,6,3)
res[1,] <- colMeans((chi.coverage))
res[2,] <- colMeans((normal.coverage))
res[3,] <- colMeans((tsls.coverage))

res[4,] <- colMeans((chi.ci.len))
res[5,] <- colMeans((normal.ci.len))
res[6,] <- colMeans((tsls.ci.len))

print(res)

write.csv(res, file = paste("ci_n=", n, "r=", r, "s=", s, ".csv", sep = ""))

# par(mfcol=c(1,3))
# index <- which(abs(tsls.ci.len[,3])<1)
# len.90 <- cbind(chi.ci.len[,1], normal.ci.len[,1], tsls.ci.len[,1])[index,]
# len.95 <- cbind(chi.ci.len[,2], normal.ci.len[,2], tsls.ci.len[,2])[index,]
# len.99 <- cbind(chi.ci.len[,3], normal.ci.len[,3], tsls.ci.len[,3])[index,]
# 
# boxplot(len.90,names = c("Projected-PPEL", "PPEL", "2SLS"),
#         col = c("red","green","blue"))
# boxplot(len.95,names = c("Projected-PPEL", "PPEL", "2SLS"),
#         col = c("red","green","blue"))
# boxplot(len.99,names = c("Projected-PPEL", "PPEL", "2SLS"),
#         col = c("red","green","blue"))
# 
# par(mfcol=c(1,1))
# index <- which(abs(tsls.ci.len[,3])<1)
# boxplot(cbind(chi.ci.len,normal.ci.len,tsls.ci.len)[index,],
#         names = c("", "Projected-EL", "","", "PEL", "",
#                   "", "2SLS", ""),
#         col = c("red","red","red","green","green","green","blue","blue","blue"))

par(mfcol=c(1,1))
# index <- which(abs(tsls.ci.len[,3])<1)
boxplot(cbind(chi.ci.len,normal.ci.len),
        names = c("", "PPEL", "","", "PEL", ""),
        col = c("white","white","white","white","white","white"))
