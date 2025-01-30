SnpIR_full=function(pheno,snp.info,D,BI_fullsample){
  # used to compute the impurity reduction of the phenotype variable with respect to the value of a given snp.
  # pheno: list of pheotype objects
  # snp.info the 0,1,2 value of the SNP
  # D: distance matrix of the phenotype objects
  # BI_fullsample: pre-computed BI of the full sample. 
  source("../utility_functions/splitDataset.R")
  source("../utility_functions/BImp_subsample.R")
  library(Rcpp)
  sourceCpp("../utility_functions/BallDistanceVector.cpp")
  sourceCpp("../utility_functions/Ball_Impurity/BallImpurity.cpp")
    N=length(snp.info)
    all_var=unique(snp.info)
    print(all_var)
    if(length(all_var)==1){  ##only one variant present
      break
      return(0)
    }
    n.cut=length(all_var)-1
    PG.full=rep(0,n.cut)###prepare the PG.full vector, each entry is a different cut (max 2 cuts for a snp with values 0,1,2)
    dataset=list(X=as.matrix(snp.info),Y=pheno)
    d=3
    print("begin BImp calculation:")
    print(Sys.time())
      BI_fullsample=BallImpurity(N,D,1:N,d)
      print(paste0("completed BImp calculation, BImp full sample equals ",BI_fullsample,"."))
      print(Sys.time())
      n.cut=length(all_var)-1 ##possible cuts

      cutpoints=sort(all_var,decreasing=TRUE)
      for(cut in 1:n.cut){
        split=splitDataset(dataset,ops1=1,ops2=cutpoints[cut])
        m1=length(split$leftInd)
        m2=N-m1
        print(c(m1,m2))
        #compute the BI on left and right leaf sets
        BI_left=BallImpurity(N,D,split$leftInd,d) 
        BI_right=BallImpurity(N,D,split$rightInd,d)
        PG_k=BI_fullsample-m1/n*BI_left-m2/n*BI_right
        PG.full[cut]<-PG_k
      }
    print(PG.full)
    return(PG.full)
}