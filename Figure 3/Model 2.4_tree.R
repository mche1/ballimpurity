rm(list=ls())
parent.path="/Users/mengluche/Documents/Research/Ball_impurity/jasa_code_submission/"
library(Ball)
library(Rcpp)
sourceCpp(paste0(parent.path,"utility_functions/BallImpurity.cpp"))
sourceCpp(paste0(parent.path,"utility_functions/BallDistanceVector.cpp"))
#source(paste0(parent.path,"simulation/Euclid/simu_rank_matrix_euclid.R"))
source(paste0(parent.path,"utility_functions/biRegressionTree_Euclidean.R"))
source(paste0(parent.path,"utility_functions/PruneTree.R"))
source(paste0(parent.path,"utility_functions/splitDataset.R"))
rep=1

n=200          #### sample size
p_X <- 3    #### number of columns in X
n_sig=3
#noise_level <- 0.2
signal_level1 <- 1
signal_level2 <- 3
e_Y=2
beta0=c(1,1)*signal_level2
results=NULL#matrix(rep(0,(n_sig+1)*length(outlier_ratio_vector)),nrow=length(outlier_ratio_vector))
#for(k in 1:length(alpha_vector)){
rank_mat=NULL
  for(i in 1:rep){
    seed=202279+i
    set.seed(seed)
    X=mvrnorm(n=200,mu=c(5,10,15),Sigma = matrix(c(3,1,-1,1,4,-0.5,1,-0.5,2),ncol=3))#matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
    #### the original predictors are binary with p=0.5
    #sig = sample(1:p_X,size=n_sig,replace=F)
    #### randomly select n_sig of them to be active predictors
    #print(sig)
    Z=cbind(X[,1]>=5,X[,2]>=11,X[,3]>=14)
    transformed_cov=cbind(Z[,1]*Z[,2],Z[,3])
    Y=transformed_cov%*%beta0+rnorm(1,mean=0,sd=e_Y)
    Y=scale(Y)
    D=BallDistanceVector(Y)
    biTree<- buildTree(list(X=X,Y=Y), leafFunc=regLeaf, method='BallImpurity',D, ops=c(0, minLeafSize=10), maxDepth=4, lossChoice="loss2")
}
biTree<-PruneTreeComplete(biTree,1.5,4)
save.image("model2.4_20250120.RData")

