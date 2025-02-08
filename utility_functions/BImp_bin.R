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
