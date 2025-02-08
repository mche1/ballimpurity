rm(list=ls())
source("/utility_functions/simu_rank_matrix.R")
library(snow)
library(Ball)
library(snowfall)
library(tictoc)
library(movMF)
library(MASS)
rep=100
n=200          #### sample size
p_X <- 500    #### number of columns in X
n_sig=2
p_Y<-5

outlier_ratio=0.1
startingseed=202479
B_True=matrix(c(0,2/3,-1/3,1,0,0,-1/3,-1/3,4/3,0)*1.5,nrow=2,byrow=T)

rank_mat=list(nh_BI4=matrix(rep(0,rep*n_sig),nrow=rep),nh_BCor=matrix(rep(0,rep*n_sig),nrow=rep),nh_DCor=matrix(rep(0,rep*n_sig),nrow=rep),
              max_BI4=matrix(rep(0,rep*n_sig),nrow=rep),max_BCor=matrix(rep(0,rep*n_sig),nrow=rep),max_DCor=matrix(rep(0,rep*n_sig),nrow=rep))

tic()
sfInit(cpus=4,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
sfLibrary(Ball)
rank_mat$nh_BI4[1:30,]=sfSapply(seq(1:30),sim_BI4)
sfStop()
toc()
sum(rank_mat$nh_BI4[1:100,1]<11)
sum(rank_mat$nh_BI4[1:100,2]<11)
sum(rowSums(rank_mat$nh_BI4[1:100,]<11)>1)

tic()
sfInit(cpus=4,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
sfLibrary(Ball)
rank_mat$nh_BCor[1:100,]=sfSapply(seq(1:100),sim_BCor)
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
sfInit(cpus=4,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
rank_mat$eu_BI4[1:100,]=sfSapply(seq(1:100),sim_euBI4)
sfStop()
toc()
sum(rank_mat$eu_BI4[1:100,1]<11)
sum(rank_mat$eu_BI4[1:100,2]<11)
sum(rowSums(rank_mat$eu_BI4[1:100,]<11)>1)

tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
rank_mat$eu_DCor[1:100,]=sfSapply(seq(1:100),sim_euDCor)
sfStop()
toc()
sum(rank_mat$eu_DCor[1:100,1]<11)
sum(rank_mat$eu_DCor[1:100,2]<11)
sum(rowSums(rank_mat$eu_DCor[1:100,]<11)>1)

tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
rank_mat$eu_BCor[1:100,]=sfSapply(seq(1:100),sim_euBCor)
sfStop()
toc()
sum(rank_mat$eu_BCor[,1]<11)
sum(rank_mat$eu_BCor[,2]<11)
sum(rowSums(rank_mat$eu_BCor[],]<11)>1)


tic()
sfInit(cpus=4,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
#rank_mat$nh_BCor=sfSapply(seq(1:100),sim_BCor)
rank_mat$max_BI4[1:100,]=sfSapply(seq(1:100),sim_maxBI4)
sfStop()
toc()

sum(rank_mat$max_BI4[1:100,1]<11)
sum(rank_mat$max_BI4[,2]<11)
sum(rowSums(rank_mat$max_BI4[,]<11)>1)


tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
rank_mat$max_BCor[1:100,]=sfSapply(seq(1:100),sim_maxBCor)
sfStop()
toc()

tic()
sfInit(cpus=3,parallel=TRUE) 
sfExportAll()
sfLibrary(MASS)
sfLibrary(movMF)
rank_mat$max_DCor[1:100,]=sfSapply(seq(1:100),sim_maxDCor)
sfStop()
toc()

save.image("md32_outlier_10.RData")

sim_BI4=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
  mu0=c(1,0,0,0,0)#*signal_level
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
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
  D=nhdist(Y)
  res_r=simu_rank_matrix(Y,X,sig,D,method="BI",dweight=4)$rank
  return(res_r)
}
sim_BCor=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
  mu0=c(1,0,0,0,0)#*signal_level
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
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
  D=nhdist(Y)
  res_r=simu_rank_matrix(Y,X,sig,D,method="BCor")$rank
  return(res_r)
}
sim_DCor=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
  mu0=c(1,0,0,0,0)#*signal_level
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
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
  D=nhdist(Y)
  res_r=simu_rank_matrix(Y,X,sig,D,method="DCor")$rank
  return(res_r)
}
sim_euBI4=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
  mu0=c(1,0,0,0,0)#*signal_level
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
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
  D=as.matrix(dist(Y,method = "euclidean"))
  res_r=simu_rank_matrix(Y,X,sig,D,method="BI",dweight=4)$rank
  return(res_r)
}
sim_euBCor=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
  mu0=c(1,0,0,0,0)#*signal_level
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
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
  D=as.matrix(dist(Y,method = "euclidean"))
  res_r=simu_rank_matrix(Y,X,sig,D,method="BCor")$rank
  return(res_r)
}
sim_euDCor=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
  mu0=c(1,0,0,0,0)#*signal_level
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
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
  D=as.matrix(dist(Y,method = "euclidean"))
  res_r=simu_rank_matrix(Y,X,sig,D,method="DCor")$rank
  return(res_r)
}

sim_maxBI4=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
  mu0=c(1,0,0,0,0)#*signal_level
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
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
  D=as.matrix(dist(Y,method = "maximum"))
  res_r=simu_rank_matrix(Y,X,sig,D,method="BI",dweight=4)$rank
  return(res_r)
}
sim_maxBCor=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
  mu0=c(1,0,0,0,0)#*signal_level
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
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
  D=as.matrix(dist(Y,method = "maximum"))
  res_r=simu_rank_matrix(Y,X,sig,D,method="BCor")$rank
  return(res_r)
}
sim_maxDCor=function(r){
  seed=startingseed+r
  set.seed(seed)
  X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
  sig = sample(1:p_X,size=n_sig,replace=F)
  #print(sig)
  Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
  mu0=c(1,0,0,0,0)#*signal_level
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
  Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
  D=as.matrix(dist(Y,method = "maximum"))
  res_r=simu_rank_matrix(Y,X,sig,D,method="DCor")$rank
  return(res_r)
}
