PruneTreeLeaf<-function(tree,thres){
    if(length(tree)==4){### leaf node reached
      print("No branches!")
      #flag<<-0
      #break
      return(tree)
    }else if((length(tree$left)==4)&(length(tree$right)==4)){####the last branch
      flag=as.numeric(tree$splitstat1>thres)
      if(flag==1){### the last branch needs to be merged
        tree$splitInd<-NULL;tree$splitVal<-NULL;tree$splitImpurity<-NULL;tree$splitstat1<-NULL;
        tree$leafValue=(tree$right$leafValue*tree$right$nsample+tree$left$leafValue*tree$left$nsample)/tree$nsample
        tree$left<-NULL;tree$right<-NULL}
      #break #
      return(tree)
    }else{ ### not the last branch, keep pruning. 
      tree$left=PruneTreeLeaf(tree$left,thres); tree$right=PruneTreeLeaf(tree$right,thres)
      return(tree)}
}
PruneTreeComplete<-function(tree,thres,depth){
  for(d in 1:depth){
    tree=PruneTreeLeaf(tree,thres)
  }
  return(tree)
}
