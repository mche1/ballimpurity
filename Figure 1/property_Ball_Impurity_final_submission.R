##Property of Ball Impurity
##Menglu Che 20231119
rm(list=ls())
library(shapes)
library(Rcpp)
library(Ball)    # for Ball Correlation
library(energy)  # for Distance Correlation

#setwd("/Users/mengluche/Documents/Research/Ball_impurity/")
sourceCpp("BallDistanceVector.cpp")
sourceCpp("BallImpurity.cpp")
repeat.time=100

#Gaussion
n=500
p=2
d=0.5
e_seq=c(20:50)/50
BI_seq_mean=rep(0,length(e_seq))
Dcov_seq_mean=rep(0,length(e_seq))
Bcov_seq_mean=rep(0,length(e_seq))
Sum_var_seq_mean=rep(0,length(e_seq))

for(r in 1:repeat.time)
{
  BI_seq=rep(NA,length(e_seq))
  Dcov_seq=rep(NA,length(e_seq))
  Bcov_seq=rep(NA,length(e_seq))
  Sum_var_seq=rep(NA,length(e_seq))
  for(i in 1:length(e_seq))
  {
    e=e_seq[i] #noise level sou
    Y=matrix(rnorm(n*p,mean=0, sd=e),nrow=n)
    #plot(Y)
    distance_matrix=BallDistanceVector(Y)
    BI_seq[i]=BallImpurity(n,distance_matrix,1:n,d)
    Dcov_seq[i]=dcov(Y,Y)
    Bcov_seq[i]=bcov(Y,Y)
    Sum_var_seq[i]=var(Y[,1])+var(Y[,2])
  }
  BI_seq_mean=BI_seq_mean+BI_seq
  Dcov_seq_mean=Dcov_seq_mean+Dcov_seq
  Bcov_seq_mean=Bcov_seq_mean+Bcov_seq
  Sum_var_seq_mean=Sum_var_seq_mean+Sum_var_seq
}

result=cbind(BI_seq_mean,Bcov_seq_mean, Dcov_seq_mean,Sum_var_seq_mean)
result=as.data.frame(scale(result))
min.y=min(result)
max.y=max(result)
pdf("Simulation3_nov24.pdf",width=7,height=7)
plot(e_seq,result$BI_seq_mean,pch=1,type="b",ylim=c(min.y,max.y),xlab="noise",ylab="",yaxt="n",cex.lab=2.5,cex.axis=2.5)
lines(e_seq,result$Bcov_seq_mean,pch=2,type="b",col="blue")
lines(e_seq,result$Dcov_seq_mean,pch=3,type="b",col="red")
lines(e_seq,result$Sum_var_seq_mean,pch=4,type="b",col="green")
legend("topleft", legend=c("BI","BCov", "DCov","SumVar"),
       col=c("black","blue", "red","green"), lty=1,pch=1:4, cex=2,
       box.lty=0)
dev.off()

##Mixtured model 1 Gaussion
n=500
p=2
Y=matrix(NA,nrow=n,ncol=p)
e=0.1
d=0.2
Pi_seq=seq(0,1,0.05)
n_seq=Pi_seq*n
BI_seq_mean=rep(0,length(Pi_seq))
Bcov_seq_mean=rep(0,length(Pi_seq))
Dcov_seq_mean=rep(0,length(Pi_seq))
Sum_var_seq_mean=rep(0,length(Pi_seq))

for(r in 1:repeat.time)
{
  BI_seq=rep(NA,length(Pi_seq))
  Bcov_seq=rep(NA,length(Pi_seq))
  Dcov_seq=rep(NA,length(Pi_seq))
  Sum_var_seq=rep(NA,length(Pi_seq))
  #### sample 2n Gaussian distributed points, and truncate 
  Y1=matrix(rep(0,n*4),nrow=2*n)
  Y2=matrix(rep(0,n*4),nrow=2*n)
  mu1=c(2,2)
  mu2=c(2,-2)
  for(j in 1:2){
    Y1[,j]=rnorm(2*n,mean=mu1[j],sd=e)
    Y2[,j]=rnorm(2*n,mean=mu2[j],sd=e)
  }
  d1=apply(Y1-rep(mu1,each=nrow(Y1)),1,norm,type = "2")
  d2=apply(Y2-rep(mu2,each=nrow(Y2)),1,norm,type = "2")
  id1=which(d1<=0.5)
  id2=which(d2<=0.5)
  for(i in 1:length(Pi_seq))
  {
    n1=n*Pi_seq[i]
    n2=n-n1
    Y.t=rbind(Y1[sample(id1,size=n1),],Y2[sample(id2,size=n2),])
    Y=list()
    for(k in 1:nrow(Y.t)){Y[[k]]<-Y.t[k,]}
    distance_matrix=BallDistanceVector(Y)
    BI_seq[i]=BallImpurity(n,distance_matrix,1:n,d)
    Dcov_seq[i]=dcov(Y.t,Y.t)
    Bcov_seq[i]=bcov(Y.t,Y.t)
    Sum_var_seq[i]=var(Y.t[,1])+var(Y.t[,2])
    #plot(Y)
  }
  print(BI_seq)
  BI_seq_mean=BI_seq_mean+BI_seq
  Bcov_seq_mean=Bcov_seq_mean+Bcov_seq
  Dcov_seq_mean=Dcov_seq_mean+Dcov_seq
  Sum_var_seq_mean=Sum_var_seq_mean+Sum_var_seq
}
result=cbind(BI_seq_mean,Bcov_seq_mean,Dcov_seq_mean,Sum_var_seq_mean)
result=as.data.frame(scale(result))
min.y=min(result)
max.y=max(result)
pdf("Simulation4_nov25.pdf",width=7,height=7)
plot(Pi_seq,result$BI_seq_mean,pch=1,type="b",ylim=c(min.y,max.y),xlab=expression(pi),ylab="",yaxt="n",cex.lab=2.5, cex.axis=2.5)
lines(Pi_seq,result$Bcov_seq_mean,pch=2,type="b",col="blue")
lines(Pi_seq,result$Dcov_seq_mean,pch=3,type="b",col="red")
lines(Pi_seq,result$Sum_var_seq_mean,pch=4,type="b",col="green")
legend("bottom", legend=c("BI","Bcov", "Dcov","SumVar"),
       col=c("black","blue","red","green"), lty=1,pch=1:4, cex=1.8,
       box.lty=0)
dev.off()

##Mixtured model 3 
 n=500
 p=2
 Y=matrix(NA,nrow=n,ncol=p)
 e=0.1
 d=0.2
 Pi_seq=seq(0,pi/2,length.out = 20)
 
 BI_seq_mean=rep(0,length(Pi_seq))
 
 Bcov_seq_mean=rep(0,length(Pi_seq))
 Dcov_seq_mean=rep(0,length(Pi_seq))
 Sum_var_seq_mean=rep(0,length(Pi_seq))
 
 for(r in 1:repeat.time)
 {
   BI_seq=rep(NA,length(Pi_seq))
   Bcov_seq=rep(NA,length(Pi_seq))
   Dcov_seq=rep(NA,length(Pi_seq))
   Sum_var_seq=rep(NA,length(Pi_seq))
   for(i in 1:length(Pi_seq))
   {
     Pi=runif(n,min=0,max=Pi_seq[i])+rbinom(n,1,0.5)*pi
     Y[,1]=cos(Pi)+rnorm(n,mean=0,sd=e)
     Y[,2]=sin(Pi)+rnorm(n,mean=0,sd=e)
     plot(Y,xlim=c(-1,1),ylim=c(-1,1))
     distance_matrix=BallDistanceVector(Y)
     BI_seq[i]=BallImpurity(n,distance_matrix,1:n,d)
     Dcov_seq[i]=dcov(Y,Y)
     Bcov_seq[i]=bcov(Y,Y)
     Sum_var_seq[i]=var(Y[,1])+var(Y[,2])
   }
   BI_seq_mean=BI_seq_mean+BI_seq
   Bcov_seq_mean=Bcov_seq_mean+Bcov_seq
   Dcov_seq_mean=Dcov_seq_mean+Dcov_seq
   Sum_var_seq_mean=Sum_var_seq_mean+Sum_var_seq
 }
 
 
 result=cbind(BI_seq_mean,Bcov_seq_mean,Dcov_seq_mean,Sum_var_seq_mean)
 result=as.data.frame(scale(result))
 min.y=min(result)
 max.y=max(result)
 pdf("Simulation5_nov24.pdf",width=7,height=7)
 plot(Pi_seq,result$BI_seq_mean,pch=1,type="b",ylim=c(min.y,max.y+0.5),xlab=expression(psi),ylab="",yaxt="n",cex.lab=2.5,cex.axis=2.5)
 lines(Pi_seq,result$Bcov_seq_mean,pch=2,type="b",col="blue")  
 lines(Pi_seq,result$Dcov_seq_mean,pch=3,type="b",col="red")
 lines(Pi_seq,result$Sum_var_seq_mean,pch=4,type="b",col="green")
 legend("topleft", legend=c("BI","BCov", "DCov","SumVar"),
        col=c("black","blue", "red","green"), lty=1,pch=1:4, cex=1.8,
        box.lty=0)
 dev.off()

