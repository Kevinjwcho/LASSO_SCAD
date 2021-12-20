# install.packages("glmnet")
# install.packages("MASS")

rm(list=ls(all=TRUE))

library(glmnet)
library(MASS)

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

C_list <- c(3)

n <- 200

rho <- 0.5

IterN<-200




logistic_hbic <- function(IterN, n, p, mu, vmatirx, tmpB, cfe){
  # selection result
  
  res1<-c()  #log(p)
  res2<-c()  #new_C(n, 0, p, denz)
  res3<-c()  #new_C(n, 0.5, p, denz)
  res4<-c()  #new_C(n, 1, p, denz)
  res5<-c()  #new_C(n, 1, p, denz)
  
  # Correct and Incorrect
  NC_log<-c()
  NC_newC_gam_0<-c()
  NC_newC_gam_0.5<-c()
  NC_newC_gam_1<-c()
  NC_newC_gam_log<-c()
  
  NI_log<-c()
  NI_newC_gam_0<-c()
  NI_newC_gam_0.5<-c()
  NI_newC_gam_1<-c()
  NI_newC_gam_log<-c()
  
  NCC_log <- numeric(length(which(tmpB!=0)))
  NCC_newC_gam_0 <- numeric(length(which(tmpB!=0)))
  NCC_newC_gam_0.5 <- numeric(length(which(tmpB!=0)))
  NCC_newC_gam_1 <- numeric(length(which(tmpB!=0)))
  NCC_newC_gam_log <- numeric(length(which(tmpB!=0)))
  
  
  
  for (iter in 1:IterN)
  {
    if(iter%%10 == 0){cat("Iteration is now ", iter, " of ", IterN, "\n")}
    # X_old <-mvrnorm(n,mu,vmatrix)
    # X <- pnorm(X_old)-0.5
    X <-mvrnorm(n,mu,vmatrix)
    fx<-X%*%tmpB
    px<-exp(fx)/(1+exp(fx))
    
    Y=rbinom(n=length(px),prob=px,size=1)
    
    cvob1=cv.glmnet(X,Y, family="binomial", intercept=FALSE)
    #cvob1=cv.glmnet(X,Y, family="binomial")
    # plot(cvob1)
    beta_i<-as.numeric(coef(cvob1))
    
    # l_min<-0.01
    # l_width<-10^(-3)
    # l_max<-0.065
    l_min<-l_width<-0.005
    l_max<-0.35
    lam_grid<-seq(l_min, l_max, l_width)
    
    
    beta_sc<-rep(NA,p+1)
    for (i in 1:length(lam_grid)){
      lam_vec<-Pprime(abs(beta_i[-1]),lam_grid[i]);
      ind<-(1:p)[lam_vec==0]
      
      temp<-X%*%diag(lam_grid[i]/lam_vec)
      temp[,ind]<-X[,ind]
      fit1=glmnet(temp,Y,family="binomial", intercept=FALSE)
      #print(fit1)
      beta_scad<-as.numeric(coef(fit1,s=lam_grid[i])) # extract coefficients at a single value of lambda
      beta_sc<-rbind(beta_sc, beta_scad)
    }
    beta_sc<-beta_sc[-1,]
    
    
    BICh1_i<-c() # log(p)
    BICh2_i<-c() # new_C(n, 0, p, denz)
    BICh3_i<-c() # new_C(n, 0.5, p, denz)
    BICh4_i<-c() # new_C(n, 1, p, denz)
    BICh5_i<-c() # new_C(n, 1, p, denz)
    
    L<-dim(beta_sc)[1]
    for (i in 1:L){
      beta_si<-beta_sc[i,]
      hfx<- c(X %*% beta_si[-1])
      # BICh1_i<-c(BICh1_i,2*sum(-Y*hfx + log(1+exp(hfx)))/n + (log(n)+log(p))*log(n)*denzero(beta_si[-1])/n)
      BICh1_i<-c(BICh1_i,2*sum(-Y*hfx + log(1+exp(hfx)))/n + log(p)*log(n)*denzero(beta_si[-1])/n)
      BICh2_i<-c(BICh2_i,2*sum(-Y*hfx + log(1+exp(hfx)))/n + 2*new_C(n, gam=0, p, denzero(beta_si[-1]))+log(n)*denzero(beta_si[-1])/n)
      BICh3_i<-c(BICh3_i,2*sum(-Y*hfx + log(1+exp(hfx)))/n + 2*new_C(n, gam=0.5, p, denzero(beta_si[-1]))+log(n)*denzero(beta_si[-1])/n)
      BICh4_i<-c(BICh4_i,2*sum(-Y*hfx + log(1+exp(hfx)))/n + 2*new_C(n, gam=1, p, denzero(beta_si[-1]))+log(n)*denzero(beta_si[-1])/n)
      BICh5_i<-c(BICh5_i,2*sum(-Y*hfx + log(1+exp(hfx)))/n + 2*new_C(n, gam=log(p), p, denzero(beta_si[-1]))+log(n)*denzero(beta_si[-1])/n)
    }
    
    sel1<-which(beta_sc[which.min(BICh1_i),-1]!=0)
    sel2<-which(beta_sc[which.min(BICh2_i),-1]!=0)
    sel3<-which(beta_sc[which.min(BICh3_i),-1]!=0)
    sel4<-which(beta_sc[which.min(BICh4_i),-1]!=0)
    sel5<-which(beta_sc[which.min(BICh5_i),-1]!=0)
    
    tindex<-which(tB!=0)  
    
    temp1 <-rep(0,p)
    temp1[sel1]<-1
    if (sum(temp1[tindex])<length(tindex)) { res1<-c(res1,"U") }# 2 stands for underfit
    if ((sum(temp1[tindex])==length(tindex)) && (sum(temp1[-tindex])==0)) { res1<-c(res1,"C") }# 1 stands for correct fit
    if ((sum(temp1[tindex])==length(tindex)) && (sum(temp1[-tindex])>0))  { res1<-c(res1,"O") }# 3 stands for over fit
    NC_log<-c(NC_log,sum(temp1[tindex]))
    NI_log<-c(NI_log,sum(temp1[-tindex]))
    NCC_log <- NCC_log + temp1[tindex]
    
    
    temp2 <-rep(0,p)
    temp2[sel2]<-1
    if (sum(temp2[tindex])<length(tindex)) { res2<-c(res2,"U") }# 2 stands for underfit
    if ((sum(temp2[tindex])==length(tindex)) && (sum(temp2[-tindex])==0)) { res2<-c(res2,"C") }# 1 stands for correct fit
    if ((sum(temp2[tindex])==length(tindex)) && (sum(temp2[-tindex])>0))  { res2<-c(res2,"O") }# 3 stands for over fit
    NC_newC_gam_0<-c(NC_newC_gam_0,sum(temp2[tindex]))
    NI_newC_gam_0<-c(NI_newC_gam_0,sum(temp2[-tindex]))
    NCC_newC_gam_0 <- NCC_newC_gam_0 + temp2[tindex]
    
    temp3 <-rep(0,p)
    temp3[sel3]<-1
    if (sum(temp3[tindex])<length(tindex)) { res3<-c(res3,"U") }# 2 stands for underfit
    if ((sum(temp3[tindex])==length(tindex)) && (sum(temp3[-tindex])==0)) { res3<-c(res3,"C") }# 1 stands for correct fit
    if ((sum(temp3[tindex])==length(tindex)) && (sum(temp3[-tindex])>0))  { res3<-c(res3,"O") }# 3 stands for over fit
    NC_newC_gam_0.5<-c(NC_newC_gam_0.5,sum(temp3[tindex]))
    NI_newC_gam_0.5<-c(NI_newC_gam_0.5,sum(temp3[-tindex]))
    NCC_newC_gam_0.5 <- NCC_newC_gam_0.5 + temp3[tindex]
    
    temp4 <-rep(0,p)
    temp4[sel4]<-1
    if (sum(temp4[tindex])<length(tindex)) { res4<-c(res4,"U") }# 2 stands for underfit
    if ((sum(temp4[tindex])==length(tindex)) && (sum(temp4[-tindex])==0)) { res4<-c(res4,"C") }# 1 stands for correct fit
    if ((sum(temp4[tindex])==length(tindex)) && (sum(temp4[-tindex])>0))  { res4<-c(res4,"O") }# 3 stands for over fit
    NC_newC_gam_1<-c(NC_newC_gam_1,sum(temp4[tindex]))
    NI_newC_gam_1<-c(NI_newC_gam_1,sum(temp4[-tindex]))
    NCC_newC_gam_1 <- NCC_newC_gam_1 + temp4[tindex]
    
    
    temp5 <-rep(0,p)
    temp5[sel5]<-1
    if (sum(temp5[tindex])<length(tindex)) { res5<-c(res5,"U") }# 2 stands for underfit
    if ((sum(temp5[tindex])==length(tindex)) && (sum(temp5[-tindex])==0)) { res5<-c(res5,"C") }# 1 stands for correct fit
    if ((sum(temp5[tindex])==length(tindex)) && (sum(temp5[-tindex])>0))  { res5<-c(res5,"O") }# 3 stands for over fit
    NC_newC_gam_log<-c(NC_newC_gam_log,sum(temp5[tindex]))
    NI_newC_gam_log<-c(NI_newC_gam_log,sum(temp5[-tindex]))
    NCC_newC_gam_log <- NCC_newC_gam_log + temp5[tindex]
    
    
    cat("Sel1=",sel1,"\n")
    cat("Sel2=",sel2,"\n")
    cat("Sel3=",sel3,"\n")
    cat("Sel4=",sel4,"\n")
    cat("Sel5=",sel5,"\n")
    cat("------------------------------------", "\n")
    
    
  }
  
  C <- rbind(length(which(res1 == "C")), length(which(res2 == "C")), length(which(res3 == "C")), length(which(res4 == "C")), length(which(res5 == "C")))/IterN
  
  O <- rbind(length(which(res1 == "O")), length(which(res2 == "O")), length(which(res3 == "O")), length(which(res4 == "O")), length(which(res5 == "O")))/IterN
  
  U <- rbind(length(which(res1 == "U")), length(which(res2 == "U")), length(which(res3 == "U")), length(which(res4 == "U")), length(which(res5 == "U")))/IterN
  
  NC <- rbind(mean(NC_log), mean(NC_newC_gam_0), mean(NC_newC_gam_0.5), mean(NC_newC_gam_1), mean(NC_newC_gam_log))
  
  NI <- rbind(mean(NI_log), mean(NI_newC_gam_0), mean(NI_newC_gam_0.5), mean(NI_newC_gam_1), mean(NI_newC_gam_log))
  
  NCC <- rbind(NCC_log, NCC_newC_gam_0, NCC_newC_gam_0.5, NCC_newC_gam_1 , NCC_newC_gam_log)/IterN
  
  result <- data.frame(Complete = C, Overfit = O, Underfit = U, Nc= NC, NI=NI)
  
  result2 <- as.data.frame(NCC)
  colnames(result2) <- apply(as.data.frame(tindex), 1, function(x) paste0("X_", x))
  
  result <- cbind(result, result2)
  return(result)
}


ptm <- Sys.time()
set.seed(1)
for(j in 1:length(C_list)){
  C_beta <- C_list[j]
  cat("C_beta =",C_beta,"\n")
for(i in 1:length(p_list)){
  p <- p_list[i]
  cat("p =",p,"\n")
  vmatrix <- rho^abs(outer(1:(p),1:(p),"-"))
  mu<-rep(0,p)
  cfe<-2
  tmpB <- rep(0,p)
  rind<-c(1:5)
  tmpB[rind]<-c(C_beta/(1:5))
  tB<-tmpB # the true coefficient vector
  
  result <- logistic_hbic(IterN, n, p, mu, vmatirx, tmpB, cfe)
  w <- paste0("logistic_result_ver2_", "p", p, "_C", C_beta)
  assign(w, result)
  write.csv(get(w), paste0("C:/Users/Cho/Google 드라이브/Paper working/이은령교수님 cowork/LASSO_SCAD/results/results_0916/n_200_",w,".csv"))
  
}
}
Sys.time()-ptm


n <- 400
ptm <- Sys.time()
set.seed(1)
for(j in 1:length(C_list)){
  C_beta <- C_list[j]
  cat("C_beta =",C_beta,"\n")
  for(i in 1:length(p_list)){
    p <- p_list[i]
    cat("p =",p,"\n")
    vmatrix <- rho^abs(outer(1:(p),1:(p),"-"))
    mu<-rep(0,p)
    cfe<-2
    tmpB <- rep(0,p)
    rind<-c(1:5)
    tmpB[rind]<-c(C_beta/(1:5))
    tB<-tmpB # the true coefficient vector
    
    result <- logistic_hbic(IterN, n, p, mu, vmatirx, tmpB, cfe)
    w <- paste0("logistic_result_ver2_", "p", p, "_C", C_beta)
    assign(w, result)
    write.csv(get(w), paste0("C:/Users/Cho/Google 드라이브/Paper working/이은령교수님 cowork/LASSO_SCAD/results/results_0916/n_400_",w,".csv"))
    
  }
}
Sys.time()-ptm

# 
# logistic_result_ver2_p100_C1
# logistic_result_ver2_p200_C1
# logistic_result_ver2_p400_C1
# 
# logistic_result_ver2_p100_C3
# logistic_result_ver2_p200_C3
# logistic_result_ver2_p400_C3

