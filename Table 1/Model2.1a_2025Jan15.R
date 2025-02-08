#Simulation of Ball Impurity vector version
#baased on Ting Li 20200722
rm(list=ls())
library(Rcpp)
library(energy) #for Distance Correlation
#library(SBISIS) #for Ball Correlation
library(Ball)
library(VGAM)
setwd("/Users/mengluche/Documents/Research/Ball_impurity/simulation/Euclid")
sourceCpp("BallDistanceVector.cpp")
sourceCpp("BallImpurity.cpp")
# calculate the frequency of the true signals being selected 
#Sel.rate<-function(n,d=10,true.v,rank.mtx) {
  # Input
  # n        :  the sample size
  # d        :  coeficient of cutoffs
  # true.v   :  the true variables index
  # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
  #             each column corresponds the ranked index in one replication.
  # Output
  # rate     :  the proportions that every single active predictor is selected 
  #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
  rank.mtx.sel<-rank.mtx[1:d,]
  r<-min(dim(rank.mtx)[2],length(rank.mtx))
  p0<-length(true.v)
  R<-matrix(0,p0,r)
  rate<-c()
  for (i in 1:p0) {
    for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i]))==0) }
    rate[i]<-mean(R[i,])
  }
  return(rate)
}
Sel.rate<-function(n,d=10,true.v,rank.mtx) {
  # Input
  # n        :  the sample size
  # d        :  coeficient of cutoffs
  # true.v   :  the true variables index
  # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
  #             each column corresponds the ranked index in one replication.
  # Output
  # rate     :  the proportions that every single active predictor is selected 
  #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
  rank.mtx.sel<-rank.mtx[1:d,]
  r<-min(dim(rank.mtx)[2],length(rank.mtx)) #repeat times
  p0<-nrow(true.v)
  R<-matrix(0,p0,r)
  rate<-rep(NA,p0)
  for (i in 1:p0) {
    for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i,j]))==0) }
    rate[i]<-mean(R[i,])
  }
  return(rate)
}
# calculate the frequency all true signals being selected.
Sel.rate.all<-function(n,d=10,true.v,rank.mtx) {
  # Input
  # n        :  the sample size
  # d        :  coeficient of cutoffs
  # true.v   :  the true variables index
  # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
  #             each column corresponds the ranked index in one replication.
  # Output
  # rate     :  the proportions that every single active predictor is selected 
  #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
  rank.mtx.sel<-rank.mtx[1:d,]
  r<-min(dim(rank.mtx)[2],length(rank.mtx))
  p0<-nrow(true.v)
  R<-matrix(0,p0,r)
  rate<-rep(NA,p0)
  for (i in 1:p0) {
    for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i,j]))==0) }
  }
  rate<-mean(apply(R,2,min))
  return(rate)
}

n=200
p=100
e_X=0.05 #noise level on X
e_Y=0.05 #noise level on Y
b_X=4  ##coefficient generating Z
b_Y=3  ## coefficient generating Y
repeat_times=100 #repeat times

IG_result_weighted_sum=matrix(NA,nrow=p,ncol=repeat_times)
BCor_result=matrix(NA,nrow=p,ncol=repeat_times)
Dcor_result=matrix(NA,nrow=p,ncol=repeat_times)
true_set_rec=matrix(NA,nrow=4,ncol=repeat_times)

for (r in 1:repeat_times)
{
  #Generate X
  set.seed(20250124+r)
  X=matrix(rbinom(n*p,1,0.5),nrow=n, ncol=p)
  true_set=sample(1:100,size=4)
  #add noise on X
  Z=apply(X,c(1,2),function(x) b_X*x+rnorm(1,mean=0,sd=e_X))
  Z=scale(Z)
  #Generate Y
  Y=list()
  for(i in 1:n)
  {
    Y1=b_Y*Z[i,true_set[1]]*Z[i,true_set[2]]*Z[i,true_set[3]]+b_Y*Z[i,true_set[4]]+rnorm(1,mean=0,sd=e_Y)
    
    Y[[i]]=c(Y1)
  }
  #true_set=c(1,5,10,15)
  distance_matrix=BallDistanceVector(Y)

  Y_matrix=matrix(unlist(Y),nrow=n, ncol=length(Y[[1]]), byrow=T)
  Y_matrix=scale(Y_matrix)
  for(i in 1:n)
  {
    Y[[i]]=Y_matrix[i,]
  }
  
  d=quantile(distance_matrix,0.5,na.rm=T)
  BI0=BallImpurity(n,distance_matrix,1:n,d)
  #  IG_max_abs=rep(NA,p) #Information Gain
  IG_weighted_sum=rep(NA,p) #Information Gain
  BCor=rep(NA,p) #Ball Correlation    
  Dcor=rep(NA,p) #Distance Correlation
  for(j in 1:p)
  {
    tmp1=which(X[,j]==1)
    tmp0=which(X[,j]==0)
    #Information Gain on Ball Impurity: BI_0(we take 1 here for convenience)-1/n(BI_x1*n1+BI_x2)
    IG_weighted_sum[j]=BI0-(BallImpurity(n,distance_matrix,tmp1,d)*length(tmp1)/n+BallImpurity(n,distance_matrix,tmp0,d)*length(tmp0)/n)

    BCor[j]= bcor(X[,j],Y_matrix)
    Dcor[j]= dcor(X[,j],Y_matrix)
  }
  IG_result_weighted_sum[,r]=order(IG_weighted_sum,decreasing=T)
  BCor_result[,r]=order(BCor,decreasing=T)
  Dcor_result[,r]=order(Dcor,decreasing=T)
  true_set_rec[,r]=true_set
}

result=list()

select_number=10 #2*floor(n/log(n))
result<-c(result,
          list("Indicator P_m of BI"=c(Sel.rate(n,select_number,true_set_rec,IG_result_weighted_sum),"Indicator P_a of BI"=Sel.rate.all(n,select_number,true_set_rec,IG_result_weighted_sum))),
          list("Indicator P_m of BCor"=c(Sel.rate(n,select_number,true_set_rec,BCor_result),"Indicator P_a of BCor"=Sel.rate.all(n,select_number,true_set_rec,BCor_result))),
          list("Indicator P_m of Dcor"=c(Sel.rate(n,select_number,true_set_rec,Dcor_result),"Indicator P_a of Dcor"=Sel.rate.all(n,select_number,true_set_rec,Dcor_result))))
result
save.image(file="Model2.1a_Jan2025_v2.Rdata")
#IG describes the difference of 3 distributions!!!
