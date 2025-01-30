SNP_IR_subsample=function(pheno,snp.info,subsize,rep,D){
  ### pheno is the list of pheno image of all subjects (list of N matrices)
  ### snp.info is the snp matrix (N by 1 snpmatrix)
  ### we take a subset of 1/subsize of the original set (appox)
  ### rep: the number of replicates
  ### finally the mean is taken for all replicates.
  ### The maximum of impurity reduction between both possible cuts are returned. 
  source("../utility_functions/splitDataset.R")
    print(table(snp.info))
    n=length(snp.info)
    min.snp=min(as.numeric(table(snp.info))) ###the min frequent snp type
   
    print(paste0("min freq snp:",min.snp))      
    pheno.info=pheno
    if(min.snp<2*subsize){ ### one of the variants is rare, ignore.
      all_var=as.vector(na.omit(unique(snp.info)))
      rare_variant=all_var[which.min(tabulate(match(snp.info,all_var)))]
      incl_id=which(snp.info!=rare_variant) ## exclude those with this rare variant.
    
      excl_id=sort(setdiff(1:n,incl_id))
      snp.info=snp.info[incl_id]
      #D=full_dist_mat_dis[incl_id,incl_id] 
      if(length(excl_id)>0){
        for(q in 1:length(excl_id)){
          pheno.info[[(excl_id[q]-q+1)]]<-NULL
        }
      }
    }
    
    N=length(snp.info)
    distance_matrix=D
    all_var=unique(snp.info)
    print(paste0("all variants:", all_var))
    if(length(all_var)==1){  ##only one variant present
      break
      return(0)
    }
    n.cut=length(all_var)-1
    PG.sub=matrix(rep(0,n.cut*rep),nrow=rep)  ##prepare the PG.sub matrix, each row is one replicate
    for(replicate in 1:rep){###replicate the subsample computation for 10 different subsamples
                          ### and take the average PurityGain
      sub.id=NULL         ##start to select the subsample
      #print(all_var)
      set.seed(replicate+1000)
      
      sub.id=sample(1:N,size=ceiling(N/subsize),replace=F)
      print(sub.id)
      if(length(unique(snp.info[sub.id]))<length(all_var)){###if not all variants are selected, re-sample until it includes all variants
        resamp=0
        while(length(unique(snp.info[sub.id]))<length(all_var)){
          set.seed(replicate+1000+resamp)
          sub.id=sample(1:N,size=ceiling(N/subsize),replace=F)
          resamp=resamp+1
        }
      }
      #print(sub.id)
      n.sub=length(sub.id) 
      D=distance_matrix[sub.id,sub.id] ######get the n by N distance matrix
      snp.sub=snp.info[sub.id]
      pheno.sub=pheno.info[sub.id]
      dataset=list(X=as.matrix(snp.sub),Y=pheno.sub)
      d=3
      print("begin BImp calculation:")
      print(Sys.time())
      #print(n.sub)
      BI_fullsubsample=BImp_subsample(n.sub,D,1:n.sub,1:n.sub,d)
      print(paste0("completed BImp calculation, BImp full subsample equals ",BI_fullsubsample,"."))
      print(Sys.time())
      n.cut=length(all_var)-1 ##possible cuts
      
      cutpoints=sort(all_var,decreasing=TRUE)
      for(cut in 1:n.cut){
        split=splitDataset(dataset,ops1=1,ops2=cutpoints[cut])
        m1=length(split$leftInd)
        m2=n.sub-m1
        print(c(m1,m2))
        print("starts BI_left")
        BI_left=BImp_subsample(n.sub,D,1:n.sub,split$leftInd,d)#BallImpurity(n,distance_matrix,split$leftInd,d)
        print("BI_left:")
        print(BI_left)
        print("starts BI_right")
        BI_right=BImp_subsample(n.sub,D,1:n.sub,split$rightInd,d)#BallImpurity(n,distance_matrix,split$rightInd,d)
        print("BI_right:")
        print(BI_right)
        PG_k=BI_fullsubsample-m1/n.sub*BI_left-m2/n.sub*BI_right
        print(PG_k)
        PG.sub[replicate,cut]<-PG_k
        
        #if(PG_k>PG){PG<-PG_k} ##we record the PG of the best split.
      }
    }
    print(PG.sub)
    PG=apply(PG.sub,2,mean)
    print(PG)
    return(max(PG))
}