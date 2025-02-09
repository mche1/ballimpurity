rm(list = ls()); gc(reset = TRUE)
source("../utlitiy_functions/SNP_IR_full.R")
source("../utlitiy_functions/SNP_IR_subsample.R")
library(snpStats)

chr <- 
blk <- 
Sys.time()

dis.data<-paste0("../chr",chr,"/dis_geno_transed_chr",chr,"_",blk,".RData")
load(dis.data)
ver.data<-paste0("../chr",chr,"/ver_geno_transed_chr",chr,"_",blk,".RData")
load(ver.data)
snpsum <- col.summary(geno)
# remove SNPs with call rate less than 90% and MAF < 0.01
call.cut <- 0.9
minor <- 0.01
call.keep <- with(snpsum, ((!is.na(MAF)) & (MAF > minor)) & ((!is.na(Call.rate) & Call.rate>= call.cut)))
call.keep[is.na(call.keep)] <- FALSE
geno.dis <- geno.dis[, call.keep]
geno.ver <- geno.ver[,call.keep]
map.call.keep <- map.use[call.keep, ]
snpsum.call.keep <- snpsum[call.keep, ]

# remove SNPs with HWE test p-value less than 10-6
HWE.cut.off <- 10^-6
hwe.keep <- with(snpsum.call.keep, (!is.na(z.HWE) & (abs(z.HWE)<abs(qnorm(HWE.cut.off/2)))))
hwe.keep[is.na(hwe.keep)] <- FALSE
geno.dis <- geno.dis[, hwe.keep]
geno.ver <- geno.ver[, hwe.keep]
map <- map.call.keep[hwe.keep, ]
snpsum <- col.summary(rbind(geno.dis,geno.ver))


geno.dis <- t(as(geno.dis, "numeric"))
coding.test <- snpsum[, 8] - snpsum[, 6]
geno.dis[which(coding.test>0), ] <- 2 - geno.dis[which(coding.test>0), ]
geno.dis <- as.data.frame(t(geno.dis))


geno.ver <- t(as(geno.ver, "numeric"))
coding.test <- snpsum[, 8] - snpsum[, 6]
geno.ver[which(coding.test>0), ] <- 2 - geno.ver[which(coding.test>0), ]
geno.ver <- as.data.frame(t(geno.ver))

# impute NAs with MAF
for(j in 1:ncol(geno.dis)){
  geno_imp = sample(c(0,1,2),size = sum(is.na(geno.dis[, j])), prob = unlist(c(snpsum[j,c(8,7,6)])), replace =T)
  geno.dis[which(is.na(geno.dis[,j])),j] = geno_imp
}
for(j in 1:ncol(geno.ver)){
  geno_imp = sample(c(0,1,2),size = sum(is.na(geno.ver[, j])), prob = unlist(c(snpsum[j,c(8,7,6)])), replace =T)
  geno.ver[which(is.na(geno.ver[,j])),j] = geno_imp
}
rm(list=c("base.dis","base.ver","call.cut","call.keep","coding.test","dis.data","geno_imp","hwe.keep","HWE.cut.off","j","map.call.keep","map.use","minor","snpsum","snpsum.call.keep"))
# snpsum <- snpsum.call.keep[hwe.keep, ]
n.snp=nrow(map)
print(paste0("After QC, there are ", n.snp, " snps in the ",blk,"-th block of Chromesome ",chr,"."))
###########initialize the BI Reduction results matrix
PG.mat=matrix(rep(0,2*n.snp),nrow=n.snp)
row.names(PG.mat)<-row.names(map)
colnames(PG.mat)<-c("discovery","verification")
###########compute the BI reduction on the discovery set
load("../dis_pheno.RData")
load("../full_dist_mat_dis.RData")
for(k in 1:n.snp){
PG.mat[k,1]<-SNP_IR_subsample(pheno.dis,geno.dis[,k],10,10,full_dist_mat_dis)
}
rm(list=c("geno.dis","full_dist_mat_dis","pheno.dis"))
##################on the verification set
load("../ver_pheno.RData")
load("../full_dist_mat_ver.RData")
for(k in 1:n.snp){
  PG.mat[k,2]<-SNP_IR_subsample(pheno.ver,geno.ver[,k],10,10,full_dist_mat_ver)
}
rm(list=c("geno.ver","pheno.ver","full_dist_mat_ver"))
save.image(file=paste0("PG_results_chr",chr,"_blk",blk,".RData"))
