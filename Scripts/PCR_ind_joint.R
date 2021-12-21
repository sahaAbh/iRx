#setwd("./Work/IRX")
library("corrplot")
library(R.matlab)

find_K_M1<-function(Yc,Yp, VarE=0.8){
  lim <- min(dim(Yc)[1], dim(Yp)[1]) -10
  sv <- svd(Yc)
  sv1 <- svd(Yp)
  varExp <- sv$d^2/sum(sv$d^2)
  varExp1 <- sv1$d^2/sum(sv1$d^2)
  K_0 <- min(which(cumsum(varExp) >= VarE))
  K_1 <- min(which(cumsum(varExp1) >= VarE))
  K_max <- max(K_0,K_1)
  if(K_max > lim)K_max=lim 
  K_max
}

pcr_pred<- function(rdata, nc, np){

rt<- load(rdata)  

Tlarge <- cbind(C,P)
Tscale <- scale( t(Tlarge), center=T,scale=T)

Yc <-  Tscale[1:nc,]
Yp <-  Tscale[(nc+1):(nc+np),]

Yc_pca <- prcomp(Yc) 
Yp_pca <- prcomp(Yp)   
K_max <- find_K_M1(Yc,Yp, VarE=0.8)

rt_m1 <- data.frame(D=D, Yc_pca$x[,1:K_max])
pcr_model <- lm(D ~ ., data = rt_m1)
data_x_pred <- data.frame(Yp_pca$x[,1:K_max])

D_est1<- predict(pcr_model,data_x_pred)
D_est1_std <- (D_est1 - mean(D_est1))/sd(D_est1) 
return(list(PCR_1=D_est1, PCR_1_std = D_est1_std))
}
# correlation between loadings of PCsc

pcr_pred2<- function(rdata, nc, np){

  rt<- load(rdata)  
  
  Tlarge <- cbind(C,P)
  Tscale <- scale( t(Tlarge), center=T,scale=T)
  
  Yc <-  Tscale[1:nc,]
  Yp <-  Tscale[(nc+1):(nc+np),]
  
Y<- rbind(Yc,Yp)
Y_pca <- prcomp(Y)
sv <- svd(Y)
varExp <- sv$d^2/sum(sv$d^2)
K_max <- min(which(cumsum(varExp) >= 0.8))

rt_m2 <- data.frame(D=D, Y_pca$x[1:nc,1:K_max])
pcr_model <- lm(D ~ ., data = rt_m2)
data_x_pred <- data.frame(Y_pca$x[(nc+1):(nc+np),1:K_max])

D_est2<- predict(pcr_model,data_x_pred)
D_est2_std <- (D_est2 - mean(D_est2))/sd(D_est2) 
return(list(PCR_2=D_est2,PCR_2_std = D_est2_std))
}
