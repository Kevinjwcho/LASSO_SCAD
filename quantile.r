# install.packages("glmnet")
# install.packages("MASS")
# install.packages("Rglpk")

rm(list=ls(all=TRUE))

library(glmnet)
library(MASS)
library(Rglpk)

Pprime<-function(theta,lambda)
{
  a=3.7
  y=(theta <= lambda)*lambda + (a*lambda - theta)/(a-1) *((theta>lambda)&(theta<a*lambda)) + 0*(theta>a*lambda)
  return(y)
}

denzero<-function(vec)
{
  return(sum(vec!=0))
}

quant_loss <- function(x, tau){0.5*abs(x) + (tau - 0.5)*x}

lcomb <- function(n, x) {
  # prod((n-x+1):n) / factorial(x)
  lgamma(n+1)-lgamma(x+1)-lgamma(n-x+1)
}


new_C <- function(n, gam, p, denz){
  gam*lcomb(p, denz)/n
}


#ptm<-proc.time()
#set.seed(1010)
#n=1000;p=100
#nzc=trunc(p/10)
#x=matrix(rnorm(n*p),n,p)
#beta=rnorm(nzc)
#fx= x[,seq(nzc)] %*% beta
#eps=rnorm(n)*5
#y=drop(fx+eps)
#px=exp(fx)
#px=px/(1+px)
#ly=rbinom(n=length(px),prob=px,size=1)
#set.seed(1011)
#cvob1=cv.glmnet(x,y)
#plot(cvob1)
#coef(cvob1)
#predict(cvob1,newx=x[1:5,], s="lambda.min")
#title("Gaussian Family",line=2.5)
#proc.time()-ptm



# p <- 100
p_list <- c(100, 200, 400)

n <- 200

rho <- 0.5

tau<- 0.5

IterN <- 200




quantile_hbic <- function(IterN, n, p, mu, vmatirx, tmpB, cfe){
  # selection result
  
  res1<-c()  #log(p)
  res2<-c()  #new_C(n, 0, p, denz)
  res3<-c()  #new_C(n, 0.5, p, denz)
  res4<-c()  #new_C(n, 1, p, denz)
  
  # Correct and Incorrect
  NC_log<-c()
  NC_newC_gam_0<-c()
  NC_newC_gam_0.5<-c()
  NC_newC_gam_1<-c()
  
  NI_log<-c()
  NI_newC_gam_0<-c()
  NI_newC_gam_0.5<-c()
  NI_newC_gam_1<-c()
  
  for (iter in 1:IterN)
  {
    if(iter%%10 == 0){cat("Iteration is now ", iter, " of ", IterN, "\n")}
    X <-mvrnorm(n,mu,vmatrix)
    e <- rt(n,2) # t error
    
    Y <- X%*%tmpB+e
    
    cvob1=cv.glmnet(X,Y, intercept=FALSE)
    # plot(cvob1)
    beta_i<-as.numeric(coef(cvob1))
    
    l_min<-l_width<-0.005
    l_max<-0.35
    lam_grid<-seq(l_min, l_max, l_width)
    
    
    beta_sc<-rep(NA,p+1)
    for (i in 1:length(lam_grid)){
      lam_vec<-n*Pprime(abs(beta_i[-1]),lam_grid[i]);
      # ind<-(1:p)[lam_vec==0]
      
      obj <- c(rep(tau,n), rep((1-tau),n), c(0,lam_vec), c(0,lam_vec))
      mat <- cbind(diag(n), -diag(n), cbind(1,X), -cbind(1,X))
      dir <- rep("==", n)
      rhs <- Y
      max <- FALSE
      
      opt <- Rglpk_solve_LP(obj,mat,dir,rhs,max=max)$solution
      betaq <- opt[(2*n+1):(2*n+p+1)]-opt[(2*n+p+2):(2*n+2*p+2)]
      
      #print(fit1)
      beta_sc<-rbind(beta_sc, betaq)
    }
    beta_sc<-beta_sc[-(1),]
    
    
    BICh1_i<-c() # log(p)
    BICh2_i<-c() # new_C(n, 0, p, denz)
    BICh3_i<-c() # new_C(n, 0.5, p, denz)
    BICh4_i<-c() # new_C(n, 1, p, denz)
    
    L<-dim(beta_sc)[1]
    for (i in 1:L){
      beta_si<-beta_sc[i,]    
      # BICh1_i<-c(BICh1_i,log(sum(quant_loss((Y - X %*% beta_si[-1]),tau))/n) + (log(n)+log(p))*denzero(beta_si[-1])/n)
      # BICh2_i<-c(BICh2_i,log(sum(quant_loss((Y - X %*% beta_si[-1]),tau))/n) + log(p)*log(n)*denzero(beta_si[-1])/n)
      BICh1_i<-c(BICh1_i,2*log(sum(quant_loss((Y - X %*% beta_si[-1]), tau))) + log(p)*log(n)*denzero(beta_si[-1])/n)
      BICh2_i<-c(BICh2_i,2*log(sum(quant_loss((Y - X %*% beta_si[-1]), tau))) + 2*new_C(n, gam=0, p, denzero(beta_si[-1]))+log(n)*denzero(beta_si[-1])/n)
      BICh3_i<-c(BICh3_i,2*log(sum(quant_loss((Y - X %*% beta_si[-1]), tau))) + 2*new_C(n, gam=0.5, p, denzero(beta_si[-1]))+log(n)*denzero(beta_si[-1])/n)
      BICh4_i<-c(BICh4_i,2*log(sum(quant_loss((Y - X %*% beta_si[-1]), tau))) + 2*new_C(n, gam=1, p, denzero(beta_si[-1]))+log(n)*denzero(beta_si[-1])/n)
    }
    
    sel1<-which(beta_sc[which.min(BICh1_i),-1]!=0)
    sel2<-which(beta_sc[which.min(BICh2_i),-1]!=0)
    sel3<-which(beta_sc[which.min(BICh3_i),-1]!=0)
    sel4<-which(beta_sc[which.min(BICh4_i),-1]!=0)
    
    tindex<-which(tB!=0)  
    
    temp1 <-rep(0,p)
    temp1[sel1]<-1
    if (sum(temp1[tindex])<length(tindex)) { res1<-c(res1,"U") }# 2 stands for underfit
    if ((sum(temp1[tindex])==length(tindex)) && (sum(temp1[-tindex])==0)) { res1<-c(res1,"C") }# 1 stands for correct fit
    if ((sum(temp1[tindex])==length(tindex)) && (sum(temp1[-tindex])>0))  { res1<-c(res1,"O") }# 3 stands for over fit
    NC_log<-c(NC_log,sum(temp1[tindex]))
    NI_log<-c(NI_log,sum(temp1[-tindex]))
    
    
    temp2 <-rep(0,p)
    temp2[sel2]<-1
    if (sum(temp2[tindex])<length(tindex)) { res2<-c(res2,"U") }# 2 stands for underfit
    if ((sum(temp2[tindex])==length(tindex)) && (sum(temp2[-tindex])==0)) { res2<-c(res2,"C") }# 1 stands for correct fit
    if ((sum(temp2[tindex])==length(tindex)) && (sum(temp2[-tindex])>0))  { res2<-c(res2,"O") }# 3 stands for over fit
    NC_newC_gam_0<-c(NC_newC_gam_0,sum(temp2[tindex]))
    NI_newC_gam_0<-c(NI_newC_gam_0,sum(temp2[-tindex]))
    
    temp3 <-rep(0,p)
    temp3[sel3]<-1
    if (sum(temp3[tindex])<length(tindex)) { res3<-c(res3,"U") }# 2 stands for underfit
    if ((sum(temp3[tindex])==length(tindex)) && (sum(temp3[-tindex])==0)) { res3<-c(res3,"C") }# 1 stands for correct fit
    if ((sum(temp3[tindex])==length(tindex)) && (sum(temp3[-tindex])>0))  { res3<-c(res3,"O") }# 3 stands for over fit
    NC_newC_gam_0.5<-c(NC_newC_gam_0.5,sum(temp3[tindex]))
    NI_newC_gam_0.5<-c(NI_newC_gam_0.5,sum(temp3[-tindex]))
    
    temp4 <-rep(0,p)
    temp4[sel4]<-1
    if (sum(temp4[tindex])<length(tindex)) { res4<-c(res4,"U") }# 2 stands for underfit
    if ((sum(temp4[tindex])==length(tindex)) && (sum(temp4[-tindex])==0)) { res4<-c(res4,"C") }# 1 stands for correct fit
    if ((sum(temp4[tindex])==length(tindex)) && (sum(temp4[-tindex])>0))  { res4<-c(res4,"O") }# 3 stands for over fit
    NC_newC_gam_1<-c(NC_newC_gam_1,sum(temp4[tindex]))
    NI_newC_gam_1<-c(NI_newC_gam_1,sum(temp4[-tindex]))
    
    
    cat("Sel1=",sel1,"\n")
    cat("Sel2=",sel2,"\n")
    cat("Sel3=",sel3,"\n")
    cat("Sel4=",sel4,"\n")
    cat("------------------------------------", "\n")
    
  }
  
  C <- rbind(length(which(res1 == "C")), length(which(res2 == "C")), length(which(res3 == "C")), length(which(res4 == "C")))/IterN
  
  O <- rbind(length(which(res1 == "O")), length(which(res2 == "O")), length(which(res3 == "O")), length(which(res4 == "O")))/IterN
  
  U <- rbind(length(which(res1 == "U")), length(which(res2 == "U")), length(which(res3 == "U")), length(which(res4 == "U")))/IterN
  
  NC <- rbind(mean(NC_log), mean(NC_newC_gam_0), mean(NC_newC_gam_0.5), mean(NC_newC_gam_1))
  
  NI <- rbind(mean(NI_log), mean(NI_newC_gam_0), mean(NI_newC_gam_0.5), mean(NI_newC_gam_1))
  
  result <- data.frame(Complete = C, Overfit = O, Underfit = U, Nc= NC, NI=NI)
  return(result)
}


ptm <- Sys.time()
set.seed(1)
for(i in 1:length(p_list)){
  
  p <- p_list[i]
  cat("p =", p, "\n")
  vmatrix <- rho^abs(outer(1:(p),1:(p),"-"))
  mu<-rep(0,p)
  cfe<-2
  tmpB <- rep(0,p)
  rind<-c(1,2,5)
  tmpB[rind]<-c(3,1.5,2)
  tB<-tmpB # the true coefficient vector
  
  result <- quantile_hbic(IterN, n, p, mu, vmatirx, tmpB, cfe)
  w <- paste0("quantile_result_", p)
  assign(w, result)
  write.csv(get(w), paste0("C:/Users/Cho/Google 드라이브/Paper working/이은령교수님 cowork/LASSO_SCAD/results/results_0916/n_200",w,".csv"))
  
}
Sys.time()-ptm

quantile_result_100
quantile_result_200
quantile_result_400

