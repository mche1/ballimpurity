rm(list=ls())
source("/utility_functions/simu_rank_matrix.R")
library(snow)
library(snowfall)
library(tictoc)
library(movMF)
library(MASS)
rep=100
n=200          #### sample size
p_X <- 500    #### number of columns in X
n_sig=3
p_Y<-5

startingseed=202411
znoise_level=0.2
ynoise_level=0.3
B_True=matrix(c(0,1,1,1,0,0,2,2,0,0,0,3,0,0,0),nrow=3,byrow=T)
mu0=c(1,0,0,0,0)#
outlier_ratio=0.01 #to be varied
rank_mat=list(nh_BI4=matrix(rep(0,rep*n_sig),nrow=rep),nh_BCor=matrix(rep(0,rep*n_sig),nrow=rep),nh_DCor=matrix(rep(0,rep*n_sig),nrow=rep),
              max_BI4=matrix(rep(0,rep*n_sig),nrow=rep),max_BCor=matrix(rep(0,rep*n_sig),nrow=rep),max_DCor=matrix(rep(0,rep*n_sig),nrow=rep))

tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
sfLibrary("Ball")
rank_mat$nh_BI4=sfSapply(seq(1:100),sim_BI4)
sfStop()
toc()

tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
sfLibrary(Ball)
rank_mat$nh_BCor=sfSapply(seq(1:100),sim_BCor)
sfStop()
toc()

tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
sfLibrary(Ball)
rank_mat$nh_DCor=sfSapply(seq(1:100),sim_DCor)
sfStop()
toc()



tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
sfLibrary(Ball)
rank_mat$max_BI4=sfSapply(seq(1:100),sim_maxBI4)
sfStop()
toc()

tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
rank_mat$max_BI8=sfSapply(seq(1:100),sim_maxBI8)
sfStop()
toc()

tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
rank_mat$max_BCor=sfSapply(seq(1:100),sim_maxBCor)
sfStop()
toc()

tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
rank_mat$max_DCor=sfSapply(seq(1:100),sim_maxDCor)
sfStop()
toc()

save.image("md31_outlier_1.RData")

thres=21
rowSums(rank_mat$nh_BI4<thres)
sum(colSums(rank_mat$nh_BI4<thres)>2)
rowSums(rank_mat$nh_BCor<thres)
sum(colSums(rank_mat$nh_BCor<thres)>2)
rowSums(rank_mat$nh_DCor<thres)
sum(colSums(rank_mat$nh_DCor<thres)>2)




rowSums(rank_mat$max_BI4<thres)
sum(colSums(rank_mat$max_BI4<thres)>2)
rowSums(rank_mat$max_BCor<thres)
sum(colSums(rank_mat$max_BCor<thres)>2)
rowSums(rank_mat$max_DCor<thres)
sum(colSums(rank_mat$max_DCor<thres)>2)



sim_BI4=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
  Y=matrix(NA,nrow=n,ncol=p_Y)
  out_id=sample(1:n,size=floor(n*outlier_ratio))
  disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
  if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
    disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
  }else{
    disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
  disturb=disturb_sample[disturb_filtered_id,]
  for(i in 1:dim(X)[1])
  {
    mu=Z[i,,drop=FALSE]%*%B_True
    if(i %in% out_id){
      mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
    }
    Y[i,]=rmovMF(1,theta=mu)
  }
  #### Add noise
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
  res_r=simu_rank_matrix(Y,X,sig,method="BI",dweight=4)$rank
  return(res_r)
}
sim_BCor=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
  # B_True=matrix(c(0,1,0,1,0,0,1,1,0,0)*1,nrow=2,byrow=T)
  #mu0=c(1,0,0,0,0)#*signal_level
  Y=matrix(NA,nrow=n,ncol=p_Y)
  out_id=sample(1:n,size=floor(n*outlier_ratio))
  disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
  if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
    disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
  }else{
    disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
  disturb=disturb_sample[disturb_filtered_id,]
  for(i in 1:dim(X)[1])
  {
    mu=Z[i,,drop=FALSE]%*%B_True
    if(i %in% out_id){
      mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
    }
    Y[i,]=rmovMF(1,theta=mu)
  }
  #### Add noise
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
  res_r=simu_rank_matrix(Y,X,sig,method="BCor")$rank
  return(res_r)
}
sim_DCor=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
  # B_True=matrix(c(0,1,0,1,0,0,1,1,0,0)*1,nrow=2,byrow=T)
  #mu0=c(1,0,0,0,0)#*signal_level
  Y=matrix(NA,nrow=n,ncol=p_Y)
  out_id=sample(1:n,size=floor(n*outlier_ratio))
  disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
  if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
    disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
  }else{
    disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
  disturb=disturb_sample[disturb_filtered_id,]
  for(i in 1:dim(X)[1])
  {
    mu=Z[i,,drop=FALSE]%*%B_True
    if(i %in% out_id){
      mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
    }
    Y[i,]=rmovMF(1,theta=mu)
  }
  #### Add noise
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
  res_r=simu_rank_matrix(Y,X,sig,method="DCor")$rank
  return(res_r)
}

sim_maxBI4=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
  # B_True=matrix(c(0,1,0,1,0,0,1,1,0,0)*1,nrow=2,byrow=T)
  #mu0=c(1,0,0,0,0)#*signal_level
  Y=matrix(NA,nrow=n,ncol=p_Y)
  out_id=sample(1:n,size=floor(n*outlier_ratio))
  disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
  if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
    disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
  }else{
    disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
  disturb=disturb_sample[disturb_filtered_id,]
  for(i in 1:dim(X)[1])
  {
    mu=Z[i,,drop=FALSE]%*%B_True
    if(i %in% out_id){
      mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
    }
    Y[i,]=rmovMF(1,theta=mu)
  }
  #### Add noise
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
  res_r=simu_rank_matrix_euclid(Y,X,sig,method="BI",dweight=4,dist_method = "maximum")$rank
  return(res_r)
}
sim_maxBCor=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
  # B_True=matrix(c(0,1,0,1,0,0,1,1,0,0)*1,nrow=2,byrow=T)
  #mu0=c(1,0,0,0,0)#*signal_level
  Y=matrix(NA,nrow=n,ncol=p_Y)
  out_id=sample(1:n,size=floor(n*outlier_ratio))
  disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
  if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
    disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
  }else{
    disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
  disturb=disturb_sample[disturb_filtered_id,]
  for(i in 1:dim(X)[1])
  {
    mu=Z[i,,drop=FALSE]%*%B_True
    if(i %in% out_id){
      mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
    }
    Y[i,]=rmovMF(1,theta=mu)
  }
  #### Add noise
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
  res_r=simu_rank_matrix_euclid(Y,X,sig,method="BCor",dist_method = "maximum")$rank
  return(res_r)
}
sim_maxDCor=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
  # B_True=matrix(c(0,1,0,1,0,0,1,1,0,0)*1,nrow=2,byrow=T)
  #mu0=c(1,0,0,0,0)#*signal_level
  Y=matrix(NA,nrow=n,ncol=p_Y)
  out_id=sample(1:n,size=floor(n*outlier_ratio))
  disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
  if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
    disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
  }else{
    disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
  disturb=disturb_sample[disturb_filtered_id,]
  for(i in 1:dim(X)[1])
  {
    mu=Z[i,,drop=FALSE]%*%B_True
    if(i %in% out_id){
      mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
    }
    Y[i,]=rmovMF(1,theta=mu)
  }
  #### Add noise
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
  res_r=simu_rank_matrix_euclid(Y,X,sig,method="DCor",dist_method = "maximum")$rank
  return(res_r)
}

