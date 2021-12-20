library(MASS)

p <- 10

n <- 1000

rho <- 0.5



vmatrix <- rho^abs(outer(1:(p),1:(p),"-"))
mu<-rep(0,p)

tmpB <- rep(0,p)
rind<-c(1,2,5)
tmpB[rind]<-c(3,1.5,2)
tB<-tmpB # the true coefficient vector

X <-mvrnorm(n,mu,vmatrix)
fx<-X%*%tmpB
px<-exp(fx)/(1+exp(fx))

Y=rbinom(n=length(px),prob=px,size=1)

pdf("C:/Users/Cho/Google 드라이브/Paper working/이은령교수님 cowork/LASSO_SCAD/Scatter_and_histogram.pdf")

for(i in 1:10){
plot(X[,i], Y, main=paste0("X_", i, " scatter plot"))
hist(X[,i], main=paste0("Histogram for X_", i))
}

graphics.off()