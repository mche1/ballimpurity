lossFunc2 <- function(N,D,ind,d) {
  loss=BallImpurityFunc::BallImpurity(N,D,ind,d)
  return(loss)
}
PruningStat1<-function(D,ind1){ ###S_\tau
  SA=sum(D^2)-sum(diag(D^2))
  D1=D[ind1,ind1]
  S1=sum(D1^2)-sum(diag(D1^2))
  ind2=setdiff(seq(1:nrow(D)),ind1)
  D2=D[ind2,ind2]
  S2=sum(D2^2)-sum(diag(D2^2))
  df1=choose(nrow(D),2)
  numerator=SA/df1
  #denom=S1/choose(length(ind1),2)+S2/choose(nrow(D)-length(ind1),2)
  df2=(choose(length(ind1),2)+choose(nrow(D)-length(ind1),2))
  denom=(S1+S2)/df2
  fval=numerator/denom
  stat=1-pf(fval,df1,df2)
  return(stat)
}
regLeaf <- function(target, method = 'BallImpurity') {
  # 返回叶子节点的值
  # input: 被划分为同一区域的y值
  # output: 同一区域y值的均值，即为叶子节点的均值
  if (method == 'BallImpurity') {
    if (sum(class(target) == 'matrix') >= 1) { return(apply(target, 2, median)) }
    if (class(target) == 'numeric') {
      if (length(unique(target)) <= 2) {return(as.numeric(names(table(target))[table(target) == max(table(target))])[1])}
      else { return(median(target)) }
    }
    # 可以用mshape()
    
    if (class(target) == 'array') { 
      n <- dim(target)[3]
      if (n == 1) {return(target)}
      temp <- lapply(1:n, function(i, Y_=target){Y_[,,i]})
      temp <- list(data = temp, size = dim(target[,,1]), name = 'landmark')
      class(temp) <- 'riemdata'
      temp <- km_median(temp) # 'intrinsic', 'extrinsic'
      return(temp$median)
      
      # return(apply(simplify2array(target), 1:2, median))
    }
    if (class(target) == 'list') {
      if (length(target) == 1) {return(target[[1]])}
      temp <- list(data = target, size = dim(target[[1]]), name = 'landmark')
      class(temp) <- 'riemdata'
      temp <- km_median(temp) # 'intrinsic', 'extrinsic'
      return(temp$median)
      
      # n <- length(target); k <- nrow(target[[1]]); m <- ncol(target[[1]])
      # tmp <- array(as.numeric(unlist(target)), dim=c(k, m, n))
      # return(apply(simplify2array(tmp), 1:2, median))
    }
    
  } else {
    if (sum(class(target) == 'matrix') >= 1) {return(colMeans(target))}
    if (class(target) == 'numeric') {
      if (length(unique(target)) <= 2) {return(as.numeric(names(table(target))[table(target) == max(table(target))])[1])}
      else { return(mean(target)) }
    }
    # 可以用mshape()
    if (class(target) == 'array') { return(apply(simplify2array(target), 1:2, mean)) }
    if (class(target) == 'list') {
      n <- length(target); k <- nrow(target[[1]]); m <- ncol(target[[1]])
      tmp <- array(as.numeric(unlist(target)), dim=c(k, m, n))
      return(apply(simplify2array(tmp), 1:2, mean))
    }
  }
}
bestSplit <- function(dataset, D, leafFunc=regLeaf, method='BallImpurity', ops=c(0, 4), lastSplitFeat=NULL, lossChoice="loss2") {
  # CART 寻找最优切分特征与分裂值
  # input: dataset, leafFunc - 计算叶节点的方式，method - 损失函数计算方式，ops - 控制分裂停止的参数,  lastSplitFeat - 存放上一次分裂的节点
  # output: 返回切分的特征与分裂值；如果没有进一步分裂则返回预测值
  tolLoss <- ops[1]; tolN <- ops[2]
  X <- dataset$X; Y <- dataset$Y
  if (sum(class(Y) == 'matrix') >= 1) {tmpLen <- dim(unique(Y, MARGIN = 1))[1]
  }else if ((class(Y) == 'numeric') || (class(Y) == 'list')) {tmpLen <- length(unique(Y))
  }else if (class(Y) == 'array') {tmpLen <- dim(unique(Y, MARGIN = 3))[3]}
  if (tmpLen == 1) {return(leafFunc(Y, method = method))}
  if (lossChoice=='loss1') {
    errFunc = lossFunc1
  } else {
    errFunc = lossFunc2
    #D <- BallDistanceVector(Y) 
    n <- dim(D)[1]; d <- median(D)#(range(sort(D)[(n+1):n^2])[2]-range(sort(D)[(n+1):n^2])[1])/2
    #if ((class(Y) == 'array') || (class(Y) == 'list')) { D <- shapeDist(Y) } else { D <- BallDistanceVector(Y) }
  }
  if (sum(class(X) == 'matrix') >= 1) {
    n <- dim(X)[1]; p <- dim(X)[2] # X为欧氏数据
  } else if (class(X) == 'array') {
    p <- dim(X)[length(dim(X))]; n <- dim(X)[3]; D_X <- list() # X为形状数据
    for (i in 1:p) { D_X[[i]] <- shapeDist(X[,,,i]) }
  }
  
  ## 计算整棵树的情况
  if (lossChoice=='loss1') {YHat <- leafFunc(Y, method = method); loss <- errFunc(YHat, Y, method)} # loss1
  if (lossChoice=='loss2') {YHat <- leafFunc(Y, method = method); loss <-  BallImpurity(n,D,1:n,d)} # loss2
  chosen_feature <- c(1:p)
  # ## 寻找每个特征下最优分裂值
  # if (!rf) chosen_feature <- c(1:p)
  # if (rf == TRUE) chosen_feature <- sample(1:p,ceiling(sqrt(p)))
  # # if (rf == TRUE) chosen_feature <- filt_feature
  # # if((!is.null(lastSplitFeat))&&(length(chosen_feature)>2)&&(!rf)) chosen_feature <- chosen_feature[-lastSplitFeat]
  # 
  ##############
  ## 欧氏数据，求解使得损失函数最小的分裂特征及对应分裂值
  if (sum(class(X) == 'matrix') >= 1) {
    bestFeat <- chosen_feature[1]; bestVal <- X[1,bestFeat]; bestLoss <- 1e5
    for (i in 1:length(chosen_feature)) {
      feat <- chosen_feature[i]; x <- X[, feat]
      #print(paste0("now computing feature ",feat))
      #print(class(x))
      xSort <- sort(x); xSort <- unique(xSort)
      if (length(xSort) > 1) {
        cutVal <- (xSort[1:(length(xSort)-1)] + xSort[2:length(xSort)] )/2
        for (j in 1:length(cutVal)) {
          #     print(paste0("the ",j,"-th cut value"))
          tmpSubset <- splitDataset(dataset, feat, cutVal[j])
          leftSet <- tmpSubset[[1]]; rightSet <- tmpSubset[[2]]
          leftInd <- tmpSubset[[3]]; rightInd <- tmpSubset[[4]]
          leftLen <- length(leftInd); rightLen <- length(rightInd)
          leftY <- leftSet$Y; rightY <- rightSet$Y
          if ((leftLen >= tolN) && (rightLen >= tolN)) {
            if (lossChoice=='loss1') {
              leftYHat <- leafFunc(leftY, method = method); rightYHat <- leafFunc(rightY, method = method)
              nowLoss <- leftLen/n*errFunc(leftYHat, leftY, method) + rightLen/n*errFunc(rightYHat, rightY, method)
            } else if (lossChoice=='loss2') {
              D_left <- D[leftInd, leftInd]; D_right <- D[rightInd, rightInd]
              nowLoss <- leftLen/n*BallImpurity(n,D,leftInd,d) + rightLen/n*BallImpurity(n,D,rightInd,d)
              #  print(paste0("new loss: ",nowLoss))
            }
            canSplit <- TRUE
          } else {canSplit <- FALSE}
          
          if ((canSplit) && (nowLoss < bestLoss)) {
            bestLoss <- nowLoss; bestFeat <- feat; bestVal <- cutVal[j]
          }
        }
      }
      #print(paste0("best feature is ",bestFeat))
    }
    
    
    tmpSubset <- splitDataset(dataset, bestFeat, bestVal)
    leftSet <- tmpSubset[[1]]; rightSet <- tmpSubset[[2]]
    leftInd <- tmpSubset[[3]]; rightInd <- tmpSubset[[4]]
    leftLen <- length(leftInd); rightLen <- length(rightInd)
    leftY <- leftSet$Y; rightY <- rightSet$Y
    if ((loss - bestLoss) < tolLoss) {print("loss not improving"); print(paste0("leftLen",leftLen));return(YHat)}
    
    if ((leftLen < tolN) || (rightLen < tolN)) {print("min leaf size reached");return(YHat)}
    splitstat1  <- PruningStat1(D,leftInd)
    #splitstat2 <- PruningStat2(D,leftInd,rightInd,rep=100,d)
    #print(splitstat1)
    #print(splitstat2)
    #if(splitstat1>0.05){print("pruning stat reached");return(YHat)}
    return(c(bestFeat, bestVal, bestLoss,splitstat1))
  }
  
  ##############
  ## 形状数据，循环，求解使得损失函数最小的分裂特征及对应分裂中心点
  if (length(class(X))==1&&class(X)=="array") {
    bestFeat <- chosen_feature[1]; bestLoss <- 1e5; bestD_X <- D_X[[bestFeat]]
    #print(bestFeat)
    candidateVars <- which(bestD_X == max(bestD_X), arr.ind = TRUE)
    candidateVars <- getUniqueInd(candidateVars); bestCand <- candidateVars[1,]
    bestVal <- list(ind = bestCand, c1 = X[,,,bestFeat][,,bestCand[1]], c2 = X[,,,bestFeat][,,bestCand[2]]) 
    nCand <- nrow(candidateVars)
    
    for (i in 1:length(chosen_feature)) {
      feat <- chosen_feature[i]; tmpD_X <- D_X[[feat]]
      tempDx <- unique(as.numeric(tmpD_X)); tempDx <- tempDx[-which(tempDx == 0)]
      if (length(unique(X[,,,feat])) > 1){
        for (td in tempDx) {
          candidateVars <- which(tmpD_X == td, arr.ind = TRUE)
          candidateVars <- getUniqueInd(candidateVars)
          nCand <- nrow(candidateVars)
          for (j in 1:nCand) {
            candVar <- candidateVars[j,]
            tmpSubset <- splitDataset(dataset, candVar, tmpD_X)
            leftSet <- tmpSubset[[1]]; rightSet <- tmpSubset[[2]]
            leftInd <- tmpSubset[[3]]; rightInd <- tmpSubset[[4]]
            leftLen <- length(leftInd); rightLen <- length(rightInd)
            leftY <- leftSet$Y; rightY <- rightSet$Y
            if ((leftLen >= tolN) && (rightLen >= tolN)) {
              if (lossChoice=='loss1') {
                leftYHat <- leafFunc(leftY, method = method); rightYHat <- leafFunc(rightY, method = method)
                nowLoss <- leftLen/n*errFunc(leftYHat, leftY, method) + rightLen/n*errFunc(rightYHat, rightY, method)
              } else if (lossChoice=='loss2') {
                D_left <- D[leftInd, leftInd]; D_right <- D[rightInd, rightInd]
                nowLoss <- leftLen/n*BallImpurity(n,D,leftInd,d) + rightLen/n*BallImpurity(n,d,rightInd,d)
              }
              canSplit <- TRUE
            } else {canSplit <- FALSE}
            
            if ((canSplit) && (nowLoss < bestLoss)) {
              bestLoss <- nowLoss; bestFeat <- feat; bestD_X <- tmpD_X
              bestVal <- list(ind = candVar, c1 = X[,,,bestFeat][,,candVar[1]], c2 = X[,,,bestFeat][,,candVar[2]]) 
            }
          }
        }
      }
      
    }
    
    if ((loss - bestLoss) < tolLoss) {print("loss not improving");print(paste0("leftLen",leftLen));return(YHat)}
    
    tmpSubset <- splitDataset(dataset, bestVal$ind, bestD_X)
    leftSet <- tmpSubset[[1]]; rightSet <- tmpSubset[[2]]
    leftInd <- tmpSubset[[3]]; rightInd <- tmpSubset[[4]]
    leftLen <- length(leftInd); rightLen <- length(rightInd)
    leftY <- leftSet$Y; rightY <- rightSet$Y
    splitstat1  <- PruningStat1(D,leftInd)
    #    splitstat2 <- PruningStat2(D,leftInd,rightInd,rep=100,d)
    print(paste0("splitstat1=",splitstat1))
    if ((leftLen < tolN) || (rightLen < tolN)) {print("leaf size criterion met");return(YHat)}
    return(list(bestFeat, bestVal, bestLoss, splitstat1))
  }
  
}
buildTree <- function(dataset, leafFunc=regLeaf, method='BallImpurity',D, ops = c(0, 4), maxDepth=4, leafDepth=1, lastSplitFeat=NULL, lossChoice='loss2') {
  # 树的构建函数
  # input: 
  #   dataset, leafFunc - 计算叶节点的方式，errFunc - 损失函数，ops - 控制分裂停止的参数(ops[1]控制loss，ops[2]控制样本数量),
  #   since this is non-Euclidean data, need to supply distance matrix D.
  #   maxDepth - 控制树的深度, rf - 是否按随机森林的方式随机抽特征数， leafDepth - 当前节点的深度，lastSplitFeat - 上一次分裂的特征
  # output: 返回生成的树b
  
  # 选择最优划分条件，若满足停止条件，返回对应值
  nodeset<-NULL
  X <- dataset$X; Y <- dataset$Y
  splitResult <- bestSplit(dataset, D, leafFunc, method, ops, lastSplitFeat, lossChoice)
  #print(splitResult)
  if ((sum(class(X) == 'array')==1) && (sum(class(X)=='matrix')!=1)) {n <- dim(X)[3]} else {n <- dim(X)[1]}
  if (lossChoice=='loss1') {
    errFunc = lossFunc1
  } else {
    errFunc = lossFunc2
    #    if ((class(Y) == 'array') || (class(Y) == 'list')) { D <- shapeDist(Y) } else { 
    #D <- BallDistanceVector(Y);
    n <- dim(D)[1]; d <-median(D)# (range(sort(D)[(n+1):n^2])[2]-range(sort(D)[(n+1):n^2])[1])/2#}
  }
  
  if (lossChoice=='loss1') {yHat <- leafFunc(dataset$Y, method = method); loss <- errFunc(yHat, dataset$Y, method)}
  if (lossChoice=='loss2') {loss <-  BallImpurity(n,D,1:n,d)} 
  if (length(splitResult)>=3){if(length(splitResult[[2]])>1) {flag<-TRUE} else {flag<-FALSE}}
  if ((length(splitResult)==1) || (!((length(splitResult)>=3)&&flag) && (class(splitResult) != "numeric")) || (leafDepth >= maxDepth)) {
    yHat <- leafFunc(dataset$Y, method = method)
    resTree <- list()
    resTree$depth <- leafDepth
    resTree$nsample <- n
    resTree$leafValue <- yHat
    resTree$impur <- n * loss
    return(resTree)
  }
  
  if (length(class(X))==1&&class(X)=="array") {
    feat <- splitResult[[1]]; val <- splitResult[[2]]; impur <- splitResult[[3]]
  } else {
    feat <- splitResult[1]; val <- splitResult[2]; impur <- splitResult[3]; stat1<-splitResult[4]#;stat2<-splitResult[5]
  }
  resTree <- list()
  resTree$depth <- leafDepth
  resTree$nsample <- n
  resTree$splitInd <- feat
  resTree$splitVal <- val
  resTree$splitImpurity <- impur
  resTree$splitstat1<-stat1
  #resTree$splitstat2<-stat2
  # resTree$impur <- n * (loss - impur)
  resTree$impur <- n * loss
  #resTree$stat1 <- PruningStat1D,ind1)
  # 将数据集划分为左右子树，递归调用继续划分
  if (length(class(dataset$X))==1&&class(dataset$X)=="array") {
    tmpD_X <- shapeDist(dataset$X[,,,feat]); candVar <- val$ind
    subDataset <-splitDataset(dataset, candVar, tmpD_X)
  } else {subDataset <- splitDataset(dataset, feat, val)}
  lSet <- subDataset[[1]]; lD=D[subDataset$leftInd,subDataset$leftInd]; rSet <- subDataset[[2]]; rD=D[subDataset$rightInd,subDataset$rightInd]
  
  resTree$right <- buildTree(rSet,leafFunc, method,rD, ops, maxDepth, leafDepth+1, feat, lossChoice)
  resTree$left <- buildTree(lSet, leafFunc, method,lD, ops, maxDepth, leafDepth+1, feat, lossChoice)
  return(resTree)
}
getnodes<-function(tree){
  nodes<-NULL
  subtreelist<-list(tree)
  while(length(subtreelist[[1]])>0){
    nodes<-c(nodes,subtreelist[[1]]$splitInd)
    subtreelist<-c(subtreelist,list(subtreelist[[1]]$left,subtreelist[[1]]$right))
    subtreelist[[1]]<-NULL
    #return(c(tree$splitInd, getnodes(tree$left), getnodes(tree$right)))
  }
  return(nodes)
}