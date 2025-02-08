rm(list=ls())
parent.path="/gpfs/gibbs/project/zhang_heping/mc3526/Ball_Imp/"
library(Ball)
library(BallImpurityFunc)
library(Rcpp)
sourceCpp(paste0(parent.path,"code/BallImpurity.cpp"))
sourceCpp(paste0(parent.path,"code/BallDistanceVector.cpp"))
source(paste0(parent.path,"simulation/Euclid/simu_rank_matrix_euclid.R"))
source(paste0(parent.path,"simulation/Euclid/biRegressionTree_E.R"))
source(paste0(parent.path,"code/PruneTree.R"))
source(paste0(parent.path,"code/splitDataset.R"))
rep=100
n=200          #### sample size
34724-
#noise_level <- 0.2
signal_level1 <- 1
signal_level2 <- 0.5
noise_level1 <- 0.2
noise_level2 <- 0.25
beta0=c(0.5,1)*signal_level2
results=NULL#matrix(rep(0,(n_sig+1)*length(outlier_ratio_vector)),nrow=length(outlier_ratio_vector))
#for(k in 1:length(alpha_vector)){
rank_mat=NULL
  for(i in 1:rep){
    seed=202279+i
    set.seed(seed)
    X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
    #### the original predictors are binary with p=0.5
    sig = sample(1:p_X,size=n_sig,replace=F)
    #### randomly select n_sig of them to be active predictors
    #print(sig)
    Z=X[,sig]*signal_level1+matrix(rnorm(n*n_sig,mean=0,sd=noise_level1),nrow=n)
    ### transform X into Z by magnifying by signal_level1 and adding perturbation 
    ## with sd=noise_level1
    transformed_cov=cbind(Z[,1]*Z[,2],Z[,3])
    
    Y=transformed_cov%*%beta0
    rank_matrix_BCor=simu_rank_matrix_euclid(Y,X,sig,method="BCor")$rank
    bcorsig=sig[which(rank_matrix_BCor<=10)]
    biTree<- buildTree(list(X=X,Y=Y), leafFunc=regLeaf, method='BallImpurity', ops=c(0, minLeafSize=20), maxDepth=5, lossChoice="loss2")
    newTree=PruneTreeComplete(biTree,0.10,depth=8)
    bisig=unique(getnodes(newTree))
    if(length(bisig)>10){bisig=bisig[1:10]}
    rank_matrix_BT=c(as.integer(sig[1]%in%btsig),as.integer(sig[2]%in%btsig),as.integer(sig[3]%in%btsig))
    rank_mat=rbind(rank_mat,c(rank_matrix_BT))
    #dcorsig=sig[which(rank_matrix_DCor<=10)]
    #rank_mat=rbind(rank_mat,c(length(intersect(sig,dcorsig)), length(intersect(sig,bcorsig)), length(intersect(sig, btsig))))#rank_mat=rbind(rank_mat,c(rank_matrix_BCor,rank_matrix_DCor))
}
  ratesi=apply((rank_mat<=10)&(rank_mat>0),2,sum)
  
  ratesallbt=sum(apply((rank_mat[,(1:n_sig)]>0)&(rank_mat[,(1:n_sig)]<=10),1,sum)>=3)
  #ratesalldcor=sum(apply(rank_mat[,(3*n_sig+1):(4*n_sig)]<=10,1,sum)>=3)
  
  results=rbind(results,c(ratesi,ratesallbt))
#}
save.image("euclid_md1_tree_20230617.RData")
result=matrix(results,byrow=F,nrow=3)
rownames(result)<-c("dcor","bcor","bi","btree")
colnames(result)<-c("S1","S2","S3","all")
pdf("euclid_md1_tree_20230613.pdf")
barplot(result,col=c("lightblue","mistyrose","cornsilk","lightgreen"),cex.main=2.4, cex.axis=2, cex.name=2,legend=rownames(result),args.legend = list(x="topleft",inset=0,box.lty=0,cex=2),beside=TRUE,ylim=c(65,115),xpd=FALSE)
#barplot(result,col=c("lightblue","mistyrose","cornsilk"),legend=rownames(result),args.legend = list(x="topright",inset=0,box.lty=0),beside=TRUE,ylim=c(0,115))
dev.off()
#for(j in 1:n_sig){
#   pdf(paste0("may5_euclid_md1.pdf"),width=4,height=4)
#   plotdata=matrix(results[,(c(0:2)*n_sig+j)],nrow=4,byrow=T) 
#   
#   cols<-NULL
#   for(i in 1:length(alpha_vector)){
#     cols<-c(cols,paste0("alpha=",alpha_vector[i]))
#   }
#   colnames(plotdata)<-cols
#   plotdata=plotdata#[c(1,3,4),c(2:4)]
#   rownames(plotdata)<-c("BCor","BI,d=q/16","BI,d=q/8","DCor")
#   barplot(plotdata,col=c("lightblue","mistyrose","cornsilk","lavender"),main=paste0("Rates of the ",j,"-th active predictor is selected in 100 reps"),legend=rownames(plotdata),args.legend = list(x="topright",inset=0,box.lty=0),beside=TRUE,ylim=c(0,115))
#   dev.off()
# }
# pdf(paste0("sept30_alpha_covariate_all.pdf"))
# plotdata=matrix(results[,(n_sig*4+c(1:4))],nrow=4,byrow=T) 
# cols<-NULL
# for(i in 1:length(alpha_vector)){
#   cols<-c(cols,paste0("alpha=",alpha_vector[i]))
# }
# colnames(plotdata)<-cols
# plotdata=plotdata#[c(1,3,4),c(2:4)]
# rownames(plotdata)<-c("BCor","BI,d=q/16","BI,d=q/8","DCor")
# barplot(plotdata,col=c("lightblue","mistyrose","cornsilk","lavender"),main="Rates of all active predictors are selected in 100 reps",legend=rownames(plotdata),args.legend = list(x="topright",inset=0,box.lty=0),beside=TRUE,ylim=c(0,100))
# dev.off()
#.libPaths( c( .libPaths(), "/vast/palmer/apps/avx2/software/R/4.2.0-foss-2020b/lib64/R/library") )
#libPaths("/vast/palmer/apps/avx2/software/R/4.2.0-foss-2020b/lib64/R/library")
#install.packages("/gpfs/ysm/project/zhang_heping/mc3526/Ball_Impurity/BallImpurityFunc_0.1.0.tar.gz",repos=NULL,type="source",lib="/vast/palmer/apps/avx2/software/R/4.2.0-foss-2020b/lib64/R/library")
#source("/gpfs/ysm/project/zhang_heping/mc3526/Ball_Impurity/simulation/biRegressionFunc_test.R")

#biTree<- buildTree(list(X=X,Y=Y), leafFunc=regLeaf, method='BallImpurity', ops=c(0, minLeafSize=20), maxDepth=5, lossChoice="loss2")
