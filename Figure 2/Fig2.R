rm(list=ls())
#setwd("/gpfs/ysm/project/zhang_heping/mc3526/Ball_Impurity/")
library(Rcpp)
sourceCpp("BallDistanceVector.cpp")
sourceCpp("BallImpurity.cpp")
pvec=c(0.1, 0.74, 1.47)

n=500
rep=100
sd=0.1

d=0.4
set.seed(202303)
results=NULL
for(p in pvec){
  pheno_list = list()
  Z=runif(n,0,1)*p+pi*rbinom(n,size=1,prob=0.5)
  Yx=rep(0,n)
  Yy=rep(0,n)
  for(i in 1:n){
    eps=rnorm(2,0,sd)
    #theta=seq(0,2*pi,by=2*pi/(y_length-1))
    x=cos(Z[i])+eps[1]
    y=sin(Z[i])+eps[2]
    Yx[i]=x
    Yy[i]=y
    pheno_list[[i]] = matrix(cbind(x,y),ncol = 2)
  }
  pdf(paste0("cut_points",p,".pdf"),width =5,height=5)
  plot(Yx,Yy,pch=20,xlab="",ylab="",xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),cex.lab=2.2, cex.axis=2.2)
  dev.off()

  alphavec=seq(0,pi,by=0.01)
  PG_vec=rep(0,length(alphavec))
  D=BallDistanceVector(pheno_list)#dist(pheno_list)
  for(k in 1:length(alphavec)){
    alpha=alphavec[k]
    set1=which(Yy-tan(alpha)*Yx>=0)
    set2=setdiff(1:n,set1)
    m1=length(set1)
    m2=length(set2)
    if(m1==0 | m2==0){
      PG_vec[k]=0
    }else{
    BI_full=BallImpurity(n,D,1:n,d)
    BI_left=BallImpurity(n,D,set1,d)
    BI_right=BallImpurity(n,D,set2,d)
    PG_vec[k]=BI_full-m1/n*BI_left-m2/n*BI_right
  }}
  results=rbind(results,PG_vec)
}
for(j in 1:3){
  pdf(paste0("new_PG_demo_md3_p=",pvec[j],".pdf"),width=5,height=5)
  plot(alphavec,results[j,],xlab=expression(alpha),ylab="",yaxt="n",pch=20,cex.lab=2.2,cex.axis=2)
  lines(alphavec,results[j,], col="red")
  dev.off()
}
 