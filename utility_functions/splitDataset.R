splitDataset <- function(dataset, ops1, ops2) {
  # split dataset into left and right child nodes
  # input: 
  #   dataset contains covariates X and outcome Y
  #   When X is Euclidean，ops1 corresponds to the index of the feature to be used for split (X[, ops1]).
  #   ops2 is the feature value
  #  
  # output: left and right child nodes info
  leftSet <- list(); rightSet <- list()
  if (sum(class(dataset$X) == 'matrix') >= 1) {
    feat <- ops1; val <- ops2
    leftInd <- which(dataset$X[,feat] < val); rightInd <- which(dataset$X[,feat] >= val)
    leftSet$X <- dataset$X[leftInd,]; rightSet$X <- dataset$X[rightInd,]
  } else if (class(dataset$X) == 'array') {
    candidateVar <- ops1; D <- ops2
    n <- nrow(D); c1 <- candidateVar[1]; c2 <- candidateVar[2]; cInds <- c(1:n)[-c(c1,c2)]
    leftInd <- c(c1); rightInd <- c(c2)
    for (i in cInds) { if (D[c1, i] < D[c2, i]) {leftInd <- cbind(leftInd, i)} else {rightInd <- cbind(rightInd, i)} }
    leftSet$X <- dataset$X[,,leftInd,]; rightSet$X <- dataset$X[,,rightInd,]
    
    if (length(dim(leftInd)) == 1) {leftSet$X <- array(leftSet$X, c(dim(dataset$X)[1], dim(dataset$X)[2], 1, dim(dataset$X)[4]))}
    if (length(dim(rightInd)) == 1) {rightSet$X <- array(rightSet$X, c(dim(dataset$X)[1], dim(dataset$X)[2], 1, dim(dataset$X)[4]))}
    
    if (length(dim(leftSet$X)) < 4) {leftSet$X <- array(leftSet$X, c(dim(dataset$X)[1], dim(dataset$X)[2], dim(leftSet$X)[3], dim(dataset$X)[4]))}
    if (length(dim(rightSet$X)) < 4) {rightSet$X <- array(rightSet$X, c(dim(dataset$X)[1], dim(dataset$X)[2], dim(rightSet$X)[3], dim(dataset$X)[4]))}
  }
  
  if (sum(class(dataset$Y) == 'matrix') >= 1) {
    leftSet$Y <-dataset$Y[leftInd,]; rightSet$Y <-dataset$Y[rightInd,]
  } else if (class(dataset$Y) == 'numeric') {
    leftSet$Y <-dataset$Y[leftInd]; rightSet$Y <-dataset$Y[rightInd]
  } else if (class(dataset$Y) == 'array') {
    leftSet$Y <-dataset$Y[,,leftInd]; rightSet$Y <-dataset$Y[,,rightInd]
    if (length(dim(leftInd)) == 1) {leftSet$Y <- array(leftSet$Y, c(dim(dataset$Y)[1], dim(dataset$Y)[2], 1))}
    if (length(dim(rightInd)) == 1) {rightSet$Y <- array(rightSet$Y, c(dim(dataset$Y)[1], dim(dataset$Y)[2], 1))}
  } else if (class(dataset$Y) == 'list') {
    leftSet$Y=dataset[[2]][leftInd]
    rightSet$Y=dataset[[2]][rightInd]
       
  } 
  return(list(leftSet=leftSet, rightSet=rightSet, leftInd=leftInd, rightInd=rightInd))
  
}