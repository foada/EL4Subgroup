####################################################
###     run for different parameter schemes      ###
####################################################

##### case 100-120 s=3 ##################

n<-100 #100 200 500 1000 50
dw<-50 #50 100 200

dw2<- 6 #log(n)

rep.num <- 500

dz<-3
beta <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol=1) #data_3

#运行主程序
source("HDMS_HHD_CI.R")

rm(list=ls())

#### case 100-120 s=3 ##################

n<-100 #100 200 500 1000 50
dw<-120 #50 100 200

dw2<- 6 #log(n)

rep.num <- 500

dz<-3
beta <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol=1) #data_3

#运行主程序
source("HDMS_HHD_CI.R")

rm(list=ls())


##### case 100-120 s=log(n) ##################

n<-100 #100 200 500 1000 50
dw<-120 #50 100 200

dw2<- 8 #log(n)

rep.num <- 500

dz<-3
beta <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol=1) #data_3

#运行主程序
source("HDMS_HHD_CI.R")

rm(list=ls())


#####case 100-120 s=2*round(n^(1/5))##################

n<-100 #100 200 500 1000 50
dw<-120 #50 100 200

dw2<- 13 #log(n)

rep.num <- 500

dz<-3
beta <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol=1) #data_3

#运行主程序
source("HDMS_HHD_CI.R")

rm(list=ls())

##### case 100-120 s=3 ##################

n<-200 #100 200 500 1000 50
dw<-100 #50 100 200

dw2<- 6 #log(n)

rep.num <- 500

dz<-3
beta <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol=1) #data_3

#运行主程序
source("HDMS_HHD_CI.R")

rm(list=ls())

##### case 100-120 s=3 ##################

n<-200 #100 200 500 1000 50
dw<-240 #50 100 200

dw2<- 6 #log(n)

rep.num <- 500

dz<-3
beta <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol=1) #data_3

#运行主程序
source("HDMS_HHD_CI.R")

rm(list=ls())


##### case 100-120 s=log(n) ##################

n<-200 #100 200 500 1000 50
dw<-240 #50 100 200

dw2<- 12 #log(n)

rep.num <- 500

dz<-3
beta <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol=1) #data_3

#运行主程序
source("HDMS_HHD_CI.R")

rm(list=ls())


#####case 200-20 s=2*round(n^(1/5))##################

n<-200 #100 200 500 1000 50
dw<-240 #50 100 200

dw2<- 17 #log(n)

rep.num <- 500

dz<-3
beta <- matrix(c(0.5, 0.5, 0.5, 0.5), ncol=1) #data_3

#运行主程序
source("HDMS_HHD_CI.R")

rm(list=ls())
