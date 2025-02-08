library(movMF)
library(circular)
library(rgl)

kappa_vector=c(50,20,10,5,2)
j=1 ### plot the case for kappa=50. Change the index to get plots for other kappa values
n_sig=2
n=200
p_X=100
p_Y=3
signal_level1=0.5;noise_level1=0.05
set.seed(123)
X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
#### the original predictors are binary with p=0.5
sig = sample(1:p_X,size=n_sig,replace=F)
#### randomly select n_sig of them to be active predictors
#print(sig)
Z=X[,sig]+matrix(rnorm(n*n_sig,mean=0,sd=noise_level1),nrow=n)
### transform X into Z by magnifying by signal_level1 and adding perturbation 
## with sd=noise_level1
transformed_cov=Z%*%matrix(c(2/3,1/3,-1/3,-1/3),nrow=2,byrow = T)
### transform Z into a more compact, 2-dim covariate vector transformed_cov
### 
q_Z=ncol(transformed_cov)
B_True=matrix(c(0,1,1,0,-1,1),nrow=2,byrow=T)  
#B_True=c(0,rep(0,p_Y-q_Z-1),rep(1,q_Z))
mu0=c(1,0,0)#signal_level
mu0_vec<-as.vector(mu0)
##Pheno
library(movMF)
Y=matrix(NA,nrow=dim(transformed_cov)[1],ncol=length(B_True))
for(i in 1:dim(transformed_cov)[1]){
  mu=c(0,transformed_cov[i,,drop = FALSE]) #%*% B_True
  mu=mu+mu0_vec#mu=colSums(t(transformed_cov[i,,drop = FALSE]) %*% t(as.matrix(B_True_vec)))+mu0_vec
  mu1=mu*kappa_vector[j]
  mu_final=mu/norm(mu,type="2")*kappa_vector[j]
  
  Y[i,]=rmovMF(1,mu_final)
}
X_vec=rbind(c(0,0),c(0,1),c(1,0),c(1,1))
open3d()
spheres3d(0,0,0,radius=1.0,lit=TRUE,color="lavender",front="lines")
observer3d(0,0,7)
for(k in 1:4){
  x=Y[which(X[,sig[1]]==X_vec[k,1]&X[,sig[2]]==X_vec[k,2]),]
  print(colMeans(x))
  points3d(x[,1],x[,2],x[,3],col=rainbow(4)[k],mode="point",size=6.5)
 #plot(circular(x[,1]),cex=0.8,bin=50,stack=TRUE,shrink=1.3,main=paste0("kappa=",kappa[i]))
}
par3d(windowRect = c(0, 0, 400, 400))

mat <- par3d("userMatrix")
# Rotate around the 3 axis to better show the points
mat <- rotate3d(mat, angle = pi/2, x = 1, y = 0, z = 0)
par3d(userMatrix = mat)
mat <- rotate3d(mat, angle = pi/2, x = 0, y = -1, z = 0)
par3d(userMatrix = mat)
mat <- rotate3d(mat, angle = pi/16, x = 0, y = 0, z = 1)
par3d(userMatrix = mat)

legend3d("topright", legend = paste('Signal=', c('(0,0)', '(0,1)','(1,0)','(1,1)')), pch = 16, col = rainbow(4), cex=1.4, inset=c(0.01))
highlevel()
rgl.snapshot(filename="nov11_kappa50_4groups.png",fmt='png')