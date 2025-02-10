PurityGainFunc=function(Y,X,D,d,subsample=FALSE,subrep=10){
  ### Y is a list or vector of responses
  ### X is a vector of covariates
  ### D is the distance matrix of all the Y values
  ### d is the parameter in BI
  ### optional: subsample method can be called, with default number of rep=10
  ### default is to use the full sample
 
  setwd("/Users/mengluche/Documents/Research/Ball_impurity")
  source("splitDataset.R")
  source("code/code_1.0/BImp_new.R")
  #load("/gpfs/ysm/project/zhang_heping/mc3526/Ball_Impurity/realdata/full_dist_mat_dis.RData")
  #library(Rcpp)
  #sourceCpp("/Users/mengluche/Documents/Research/Ball_impurity/BallDistanceVector.cpp")
  #sourceCpp("/Users/mengluche/Documents/Research/Ball_impurity/BallImpurity.cpp")
  N=length(X)
  all_var=unique(X)
  if(length(all_var)==1){  ##only one variant present
    break
    return(0)
  }
  n.cut=length(all_var)-1
  dataset=list(X=as.matrix(X),Y=Y)
  n.cut=length(all_var)-1 ##possible cuts
  cutpoints=sort(all_var,decreasing=TRUE)
  print(cutpoints)
  if(subsample==FALSE){
    BI_full=BImp(N,D,1:N,1:N,d)
    PG.vec=rep(0,n.cut)###prepare the PG vector, each entry corresponding to a cut
    for(cut in 1:n.cut){
      split=splitDataset(dataset,ops1=1,ops2=cutpoints[cut])
      m1=length(split$leftInd)
      m2=N-m1
      BI_left=BImp(N,D,1:N,split$leftInd,d)
      BI_right=BImp(N,D,1:N,split$rightInd,d)#BallImpurity(n,distance_matrix,split$rightInd,d)
      PG_k=BI_full-m1/n*BI_left-m2/n*BI_right
      PG.vec[cut]<-PG_k}
  }
 return(PG.vec)
}