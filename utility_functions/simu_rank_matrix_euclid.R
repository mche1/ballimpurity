simu_rank_matrix_euclid=function(pheno_list,X,sig,D,method="BI",dweight=4){
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
    source("/utility_functions/PurityGainFunc.R")
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