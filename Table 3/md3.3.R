rm(list=ls())
library(Ball)
library(BallImpurityFunc)
library(movMF)
source("/utility_functions/biRegressionTree_NE.R")
source("/utility_functions/simu_rank_matrix.R")
source("/utility_functions/PruneTree.R")
rep=100
n=200          #### sample size
p_X <- 200    #### number of columns in X
n_sig=3
p_Y<-5
#noise_level <- 0.2
signal_level1 <- 3
signal_level2 <- 3
noise_level1 <- 0.2
noise_level2 <- 0.3

outlier_ratio_vector=c(0.10,0.05,0.01)
results=NULL#
for(outlier_ratio in outlier_ratio_vector){
  rank_mat=NULL
  for(r in 1:rep){
    seed=202279+r
    set.seed(seed)
    X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
    sig = sample(1:p_X,size=n_sig,replace=F)
    print(sig)
    Z=X[,sig]*signal_level1+matrix(rnorm(n*n_sig,mean=0,sd=noise_level1),nrow=n)
    transformed_cov=cbind(Z[,1]*Z[,2],Z[,3])
    q_Z=ncol(transformed_cov)
    B_True=matrix(c(0,0,1,2,2,0,0,1,0,0),nrow=2,byrow = F)#c(0,rep(0,p_Y-q_Z-1),rep(1,q_Z))*signal_level2
    mu0=c(1,0,0,0,0)#*signal_level
    B_True_vec <- as.vector(B_True*signal_level2)
    mu0_vec<-as.vector(mu0)
    ##Pheno
    library(movMF)
    Y=matrix(NA,nrow=dim(transformed_cov)[1],ncol=p_Y)
    out_id=sample(1:n,size=floor(n*outlier_ratio))
    library(MASS)
    disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
    if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
      disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))}
    else{disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
    disturb=disturb_sample[disturb_filtered_id,]
    for(i in 1:dim(transformed_cov)[1])
    {
      #mu=colSums(t(transformed_cov[i,,drop = FALSE]) %*% t(as.matrix(B_True_vec)))+mu0_vec
      mu=transformed_cov[i,,drop=FALSE]%*%B_True
      if(i %in% out_id){
        mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
      }
      Y[i,]=rmovMF(1,theta=mu)}
    dist_Y=nhdist(Y)
    biTree<- buildTree(list(X=X,Y=Y), leafFunc=regLeaf, method='BallImpurity',D=dist_Y, ops = c(0, 4), maxDepth=4, leafDepth=1, lastSplitFeat=NULL, lossChoice='loss2')
    newTree=PruneTreeComplete(biTree,0.15,depth=4)
    bisig=unique(getnodes(newTree))
    if(length(bisig)>10){bisig=bisig[1:10]}
    rank_matrix_BI=c(as.integer(sig[1]%in%bisig),as.integer(sig[2]%in%bisig),as.integer(sig[3]%in%bisig))
    rank_mat=rbind(rank_mat,c(rank_matrix_BI))
  }
  ratesi=apply((rank_mat<=10)&(rank_mat>0),2,sum)
  ratesallbi=sum(apply(rank_mat[,(1:n_sig)]>0,1,sum)>=3)
  results=rbind(results,c(ratesi,ratesallbi))#,ratesallbi4))
}
save.image("dir_md3.3.RData")
plotdata<-t(results)
cols<-NULL
for(i in 1:length(outlier_ratio_vector)){
  cols<-c(cols,paste0("outlier",outlier_ratio_vector[i]*100,"%"))
} 
colnames(plotdata)<-cols
colnames(plotdata)<-c("S1","S2","S3","all")
barplot(plotdata,col=c("lightblue","mistyrose","cornsilk"),main=paste0("Rates of all the active predictors being selected in ", rep," reps"),
        legend=rownames(plotdata),args.legend = list(x="topright",inset=0,box.lty=0),beside=TRUE,ylim=c(0,rep+10))
dev.off()