simu_rank_matrix=function(pheno_list,X,sig,D,method="BI",dweight=4){
  ### for response Y, computes the BI/BCor/DCor of all columns of X with Y
  ### returns the rank of active predictors which is the sig_idx columns of X
  # pheno_list: list of response Y
  # X: matrix of predictors X
  # sig: the true signals
  # D: distance matrix of Y. It can be a matrix or an nhdist object.
  # method: one of "BI", "BCor" or "DCor"
  # dweight: the reciprocal of the coefficient to be used of the parameter d. 
  #          d is set to median(D)/dweight
  p_X=ncol(X)
  n=nrow(X)
  n_sig=length(sig)
  dist_pheno=as.matrix(D)
  library(Ball)
  Cor=rep(0,p_X)
  METHODS<-c("BI","BCor","DCor")
  methodIndex<-pmatch(method,METHODS)
  if(is.na(method))
    stop("invalid correlation method")
  if(method==-1)
    stop("ambiguous correlation method")
  method<-METHODS[methodIndex]
  if(method %in% c("BI")){
    d=quantile(dist_pheno,0.5)/dweight
    BI_full=BImp_bin(n,dist_pheno,1:n,1:n,d=d)
    for(i in 1:p_X){
      Cor[i]=max(PurityGainFunc(pheno_list,X[,i],dist_pheno,BI_full,d=d))
    } 
  }
  else if (method %in% c("BCor")){
    for(i in 1:p_X){
      Cor[i]=bcor(dist_pheno,dist(X[,i]),distance=TRUE)
    }
  }
  else if (method %in% c("DCor")){
    library(energy)
    for(i in 1:p_X){
      Cor[i]=dcor(as.dist(dist_pheno),dist(X[,i]))
      }
  }
  o=order(Cor,decreasing = T)
  rank.Cor=rep(0,n_sig)
  for(i in 1:n_sig){
      rank.Cor[i]=which(o==sig[i])
  }
  return(list(rank=rank.Cor,order=o))
}
PurityGainFunc=function(Y,X,D,BI_full,d){
  ### Y is a list or vector of responses
  ### X is a vector of covariates
  ### D is the distance matrix of all the Y values
  ### BI_full is the BI of the root set
  ### d is the parameter in BI
  ### optional: subsample method can be called, with default number of rep=10
  ### default is to use the full sample
  setwd("/Users/mengluche/Documents/Research/Ball_impurity")
  source("splitDataset.R")
  N=length(X)
  all_var=unique(X)
  if(length(all_var)==1){  ##only one variant present
    break
    return(0)
  }
  n.cut=length(all_var)-1
  dataset=list(X=as.matrix(X),Y=Y)
  cutpoints=sort(all_var,decreasing=TRUE)
  PG=NULL
  for(cut in 1:n.cut){
      split=splitDataset(dataset,ops1=1,ops2=cutpoints[cut])
      m1=length(split$leftInd)
      m2=N-m1
      BI_left=BImp_bin(N,D,1:N,split$leftInd,d)
      BI_right=BImp_bin(N,D,1:N,split$rightInd,d)#BallImpurity(n,distance_matrix,split$rightInd,d)
      PG_k=BI_full-m1/n*BI_left-m2/n*BI_right
      PG<-c(PG,PG_k)}
  return(PG)
}
BImp_bin=function(n,D,setid,splitid,d){ 
  ######BImp_bin computes the Ball Impurity of sample of n out of totally N individuals, 
  ##### with a set of split ids which are subset of these n individuals. 
  ##### It is based on binary search and reduces the time complexity of the original ball impurity function.
  ##### distance matrix is given as D, of dimension n by N; and the sample indices are 
  ####  given in setid; split indices are given in splitid (note, split ids are the positions
  ###   they are in the sample. E.g. the sample is (10,20,...,19000)-th of the original
  ##    set and the split contains its first individual, then the first element in splitid is 1, not 10.
  ##### the parameter d of BImp needs to be specified as a positive constant.
  if(n!=length(setid)){print("dimensions do not match!")
    #  break
    return(0)
    break
  }
  bin_search = function(v, t, eps) {
    lo <- 1; hi <- length(v)
    while (lo <= hi) {
      mid <- as.integer(round((lo + hi) / 2)) # always even!
      if (abs(v[mid] - t) <= eps) {
        return(mid)
      } else if (v[mid] < t) { # C style would be OK
        lo <- mid + 1
      } else {
        hi <- mid - 1
      }
    }
    return(lo)
  }
  bin_search2 = function(v, t,eps) {
    lo <- 1; hi <- length(v)
    if(v[lo]>t){return(0)}
    else if(v[hi]<t){return(hi)}
    else{
      while (lo <= hi) {
        midl= as.integer(floor((lo + hi) / 2))
        midr= as.integer(ceiling((lo+hi)/2))
        if(v[midl]<=t && v[midr]>t){return(midl)}
        #else if(abs(v[midr]-t)<=eps){return(midr)}
        else if(v[midr]<t){lo<-midr}
        else {hi<-midl}
      }
      return(lo)
    }
  }
  BI=0
  N=ncol(D)
  m=length(splitid)
  Dsub=D[,setid]
  for(i in 1:n){ ###do the sum for each individual i
    x_i=Dsub[i,] ### the distance vector of (x1,xi) i=1,...,n 
    x_i_split=x_i[splitid]
    x.s=sort(x_i_split,index.return=T)
    
    x.is=x.s$x
    x.ii=x.s$ix
    
    for(j in 1:N){ ###find how many k satisfy D(i,k)<D(i,j)+d
      rhs=D[i,j]+d
      count=bin_search2(x.is,rhs,1e-9)
      BI=BI+(count/m)*(1-count/m)
    }  
  }
  BI=BI/(2*N^2)
  return(BI)
}
