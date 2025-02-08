rm(list = ls())
p_snp <- 0.1
error <- FALSE
iter=100
n <- 1000
print(paste0("n=",n))
num.snps.per.blk <- 30
num.blks <- 1
noise_level <- 0.5
set.seed(20220614)
library(snpStats)
##load part of chr19 data (0-3MB) of UKB
data.chr.blk <- read.plink("/../geno_simu")

# extract geno type and conduct QC -----------------------------------------------

geno <- data.chr.blk$genotypes
member <- data.chr.blk$fam$member

geno <- geno[sample(1:nrow(geno), n, replace = F),, drop = FALSE]
member <- as.numeric(rownames(geno))

# per-SNP quality control for call rate and HWE test 
snpsum <- col.summary(geno)
# remove SNPs with call rate less than 90% and MAF < 0.01
call.cut <- 0.9
minor <- 0.05
call.keep <- with(snpsum, ((!is.na(MAF)) & (MAF > minor)) & ((!is.na(Call.rate) & Call.rate>= call.cut)))
call.keep[is.na(call.keep)] <- FALSE
geno <- geno[, call.keep]
map.call.keep <- data.chr.blk$map[call.keep, ]
snpsum.call.keep <- snpsum[call.keep, ]

# remove SNPs with HWE test p-value less than 10-6
HWE.cut.off <- 10^-6
hwe.keep <- with(snpsum.call.keep, (!is.na(z.HWE) & (abs(z.HWE)<abs(qnorm(HWE.cut.off/2)))))
hwe.keep[is.na(hwe.keep)] <- FALSE
geno <- geno[, hwe.keep]
map <- map.call.keep[hwe.keep, ]
# snpsum <- snpsum.call.keep[hwe.keep, ]

rm(list = setdiff(ls(), c("geno", "map", "member", 
                          "n", "p_snp", "num.snps.per.blk", 
                          "noise_level", "num.blks", "error", "iter")))

# Randomly sample 500 SNPs for a single block
snp.sample <- sample(ncol(geno), num.snps.per.blk * num.blks)
print(snp.sample)
geno <- geno[, snp.sample]
map <- map[snp.sample, ]

# code checking
snpsum <- col.summary(geno)
geno <- t(as(geno, "numeric"))
coding.test <- snpsum[, 8] - snpsum[, 6]
geno[which(coding.test>0), ] <- 2 - geno[which(coding.test>0), ]
geno <- as.data.frame(t(geno))
rm(list = setdiff(ls(), c("geno", "map", "member", "snpsum", 
                          "num.snps.per.blk", "n", "p_snp", 
                          "noise_level", "num.blks", "error", "iter")))

# impute NAs with MAF
for(j in 1:ncol(geno)){
  geno_imp = sample(c(0,1,2),size = sum(is.na(geno[, j])), prob = unlist(c(snpsum[j,c(8,7,6)])), replace =T)
  geno[which(is.na(geno[,j])),j] = geno_imp
}

rm(list=setdiff(ls(), c("geno", 
                        "map", "member", "snpsum", "num.snps.per.blk", "num.blks",
                        "n", "p_snp", "noise_level", "error", "iter")))

#### simulation true B #####
###smile
B_True = as.matrix(read.csv("emoji_smile_adverse.csv", header = T))
B_True[is.na(B_True)] = 0
dimnames(B_True) = NULL
B_True_vec <- as.vector(B_True)

rm(list=setdiff(ls(), c("geno", "map", "member","snpsum", "num.snps.per.blk", "num.blks", "error",
                        "B_True","B_True_vec",
                        "n", "p_snp", "noise_level", "iter")))


############## Set True Signals ################
sig = matrix(sample(ncol(geno), size = round(ncol(geno) * p_snp), replace = F),
             ncol = round(ncol(geno) * p_snp), byrow = T)
B_sig = runif(nrow(sig), 1.5, 2.5)
#Transform interactive snps
transformed_geno=NULL
for(i in 1:dim(sig)[1])
{
  transformed_geno=cbind(transformed_geno,rowSums(as.matrix(geno[,sig[i,], drop = FALSE])) >  0 * B_sig[i])
}

rm(list=setdiff(ls(), c("geno","sig", "B_sig", "transformed_geno", "map", "member", "num.snps.per.blk", "num.blks", "error", "iter",
                      "snpsum", "n", "p_snp", "noise_level", "B_True","B_True_vec")))

##Pheno
Y=matrix(NA,nrow=dim(transformed_geno)[1],ncol=length(B_True_vec))
for(i in 1:dim(transformed_geno)[1])
{
  Y[i,]=colSums(t(transformed_geno[i,,drop = FALSE]) %*% t(as.matrix(B_True_vec)))
}
rm(B_True_vec)


#############################################################
################### Add NOISE ######################
#############################################################
library(Ball)
library(Matrix)
library(Rcpp)
source('/utility_functions/0-tarv-transform-simu.R')
sourceCpp("/utility_functions/BallDistanceVector.cpp")

true_sig=unlist(c(sig))

##Adding noise
ee <- matrix(rnorm(dim(Y)[1]*dim(Y)[2], mean=0, sd=noise_level),nrow=dim(Y)[1])
ee <- t(apply(ee, 1, function(ex){
  e_mat = matrix(ex, dim(B_True)[1], dim(B_True)[1])
  e_mat = as.matrix(forceSymmetric(e_mat))
  return(as.vector(e_mat))
}))
if(error){
  Y = matrix(0, nrow = dim(geno)[1], ncol = dim(B_True)[1]*dim(B_True)[1])
}
Y_noise=Y+ee
Y_noise[is.na(Y_noise)]=0

pheno_list = list()
for(i in 1:dim(Y_noise)[1]){
  pheno_list[[i]] = matrix(Y_noise[i, ],
                           nrow = sqrt(dim(Y_noise)[2]),
                           ncol = sqrt(dim(Y_noise)[2]))
}
rm(ee, Y, Y_noise)

dist_pheno = BallDistanceVector(pheno_list)

#########################Ball Impurity Computation#############
source("/utility_functions/SNP_IR_subsample.R")
source("/utility_functions/SNP_IR_full.R")
source("/utility_functions/BImp_subsample.R")
sourceCpp("/utility_functions/BallImpurity.cpp")
BI_fullsubsample=BallImpurity(n,dist_pheno,1:n,d=3)
BallImp_full=rep(0,num.snps.per.blk)
BallImp_subsample=rep(0,num.snps.per.blk)
for(i in 1:num.snps.per.blk){
BallImp_full[i]=max(SNP_IR_full(pheno_list,geno[,i],dist_pheno,BI_fullsubsample))
BallImp_subsample[i]=max(SNP_IR_subsample(pheno.dis,geno.dis[,k],10,10,full_dist_mat_dis))
}
o1=order(BallImp_full,decreasing=T)
o2=order(BallImp_subsample,decreasing=T)
print("For full data method, true signals rank ")
print(which(o1%in%true_sig))
print(paste0("out of ",length(o1), " signals."))

print("For subsample method, true signals rank ")
print(which(o2%in%true_sig))
print(paste0("out of ",length(o2), " signals."))


pdf("2024Dec_simu_1k_strongSNR.pdf",width=8,height=6)
par(mar=c(6,6,2.5,2.5))
plot(map$position,BallImp_full_scaled,xlab="SNP position",ylab="Scaled BI reduction, BCor, DCor",ylim=c(min(BallImp_full_scaled),max(c(Cor_scaled,BallImp_full_scaled,DCor_scaled))),pch=1,cex=1.5,cex.main=2,cex.lab=2)
points(map$position,Cor_scaled,pch=2,cex=1.5)
points(map$position,DCor_scaled,pch=0,cex=1.5)
points(map$position[true_sig],BallImp_full_scaled[true_sig],pch=19,cex=2,col="red")
points(map$position[true_sig],Cor_scaled[true_sig],pch=17,cex=2,col="blue")
points(map$position[true_sig],DCor_scaled[true_sig],pch=15,cex=2,col="green")
abline(h=quantile(BallImp_full_scaled,0.9),col = 'red',lwd=2, lty = 2)
abline(h=quantile(Cor_scaled,0.9),col="blue",lwd=2,lty=2)
abline(h=quantile(DCor_scaled,0.9),col="green",lwd=2,lty=2)
legend("topleft", legend = c("BI reduction","BCor","DCor","BI reduction of true signal","BCor of true signal","DCor of true signal","90% quantile of Impurity reductions"), 
       col=c('black','black','black','red','blue','green','steelblue'),pch=c(1,2,0,19,17,15,NA),lty = c(NA,NA,NA,NA,NA,NA,2),cex=1.5)
dev.off()

pdf("2024Dec_simu_1k_strongSNR_subsample.pdf",width=8,height=6)
par(mar=c(6,6,2.5,2.5))
plot(map$position,BallImp_subsample,xlab="SNP position",ylab="Impurity reduction",pch=1,cex=1,cex.main=2,cex.lab=2)
points(map$position[true_sig],BallImp_subsample[true_sig],pch=19,cex=2,col="red")
abline(h=BallImp_subsample[o2[10]],col = 'steelblue', lty = 2)
legend("topleft", legend = c("BI reduction","BI reduction of true signal","90% quantile of BI reductions"), col=c('black','red','steelblue'),pch=c(1,19,NA),lty = c(NA,NA,2),cex=1.5)
dev.off()
