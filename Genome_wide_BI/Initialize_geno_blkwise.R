rm(list = ls()); gc(reset = TRUE)
require(snpStats) 
library(igraph)
library(data.table)
#require(genio)
setwd("/gpfs/ysm/project/zhang_heping/mc3526/basket-50595")
maintable=fread("../ukb50595.csv", select = c("eid","31-0.0","34-0.0","52-0.0","53-0.0","53-2.0","22001-0.0","22006-0.0","22009-0.1",
                                           "22009-0.2","22009-0.3","22009-0.4","22009-0.5"))

bulk=read.table("../ukb50595.bulk")

load("kinship_3degree.Rdata")
kinship_3degree[,1]=as.character(kinship_3degree[,1])
kinship_3degree[,2]=as.character(kinship_3degree[,2])
Kinship_network=simplify(graph_from_edgelist(as.matrix(kinship_3degree),directed=F))

delete_relative <- function(data)
{
  cases_kinship_network=induced_subgraph(Kinship_network, v=V(Kinship_network)$name %in% data$eid)
  S_G=delete.vertices(cases_kinship_network,degree(cases_kinship_network)==0)
  #rm(cases_kinship_network)
  cases_kikedout_kinship_id=as.numeric(V(S_G)$name)
  cases_id=setdiff(data$eid,cases_kikedout_kinship_id)
  new_data=data[data$eid %in% cases_id,]
  return(new_data)
}
maintable=maintable[which(maintable$eid %in% bulk$V1),]
colnames(maintable)<-c("eid","sex","birth_year","birth_month","assessment_date","image_date","genetic_sex","ethnicity","pc1","pc2","pc3",
                       "pc4","pc5")#,"pc6","pc7","pc8","pc9","pc10")
main<-maintable[-which(is.na(maintable$genetic_sex)|maintable$sex!=maintable$genetic_sex),]
main=delete_relative(main)
main$age=year(main$image_date)-main$birth_year+(month(main$image_date)-main$birth_month)/12
#main<-maintable[-which(maintable$sex!=maintable$genetic_sex),]
#rm(maintable,kinship_3degree,Kinship_network)
#setwd("/gpfs/ysm/scratch60/zhang_heping/mc3526/Ball_Impurity/data/")
blk_info=read.csv("/gpfs/ysm/scratch60/zhang_heping/mc3526/Ball_Impurity/data/num_blk.csv",header=TRUE)
#geno.transed <<- NULL
basechar=main[,c(1,2,8,14,9:13)]
set.seed(202204)
G1=basechar[which(basechar$age<=65&basechar$sex==0),]
G2=basechar[which(basechar$age>65&basechar$sex==0),]
G3=basechar[which(basechar$age<=65&basechar$sex==1),]
G4=basechar[which(basechar$age>65&basechar$sex==1),]
d1=sample(G1$eid,size=ceiling(nrow(G1)/2),replace = FALSE)
d2=sample(G2$eid,size=ceiling(nrow(G2)/2),replace = FALSE)
d3=sample(G3$eid,size=ceiling(nrow(G3)/2),replace = FALSE)
d4=sample(G4$eid,size=ceiling(nrow(G4)/2),replace = FALSE)
v1=setdiff(G1$eid,d1)
v2=setdiff(G2$eid,d2)
v3=setdiff(G3$eid,d3)
v4=setdiff(G4$eid,d4)
discover=sort(c(d1,d2,d3,d4))
verify=sort(c(v1,v2,v3,v4))
base.dis=basechar[which(basechar$eid%in%discover),] ###baseline characteristics of the discovery set
base.ver=basechar[which(basechar$eid%in%verify),] ####baseline charac. of the verification set
for(i in 1:22){
  wd=paste0("../data/divided/chr",i,"/")
  setwd(wd)
  num_blk=blk_info$num_blk[i]
  for(j in 1:num_blk){
    bedfile=paste0("chr",i,"_",j,".bed")
    if(file.exists(bedfile)){
      geno.use<<-list()
      map.use <<- list()
      bimfile=paste0("chr",i,"_",j,".bim")
      famfile=paste0("chr",i,"_",j,".fam")
      data.chr<-read.plink(bedfile,bimfile,famfile,na.strings="0")
      geno0=data.chr$genotypes
      geno=geno0[which(rownames(geno0)%in%basechar$eid),]
      #geno.ver=geno0[which(rownames(geno0)%in%base.ver$eid),]
      eid_keep=rownames(geno) ###left 40123 subjects
  
      
      
      ###baseline characteristics matrix, with columns: eid, sex, ethnicity, age, pc1-5.
      ########Quality Control below 
      
      snpsum=col.summary(geno)
      ##remove snps with call rate less than 90%
      call.cut<-0.9
      call.keep <- with(snpsum, (!is.na(Call.rate) & Call.rate>= call.cut))
      call.keep[is.na(call.keep)] <- FALSE
      geno <- geno[, call.keep]
      map.call.keep <- data.chr$map[call.keep, ]
      snpsum.call.keep <- snpsum[call.keep, ]
      # remove SNPs with HWE test p-value less than 10-7
      HWE.cut.off <- 10^-7
      snpsum.control <- col.summary(geno)
      hwe.keep <- with(snpsum.control, (!is.na(z.HWE) & (abs(z.HWE)<abs(qnorm(HWE.cut.off/2)))))
      hwe.keep[is.na(hwe.keep)] <- FALSE
      geno.use <- geno[, hwe.keep]
      map.use <<- map.call.keep[hwe.keep, ]
      geno.dis=geno.use[which(rownames(geno.use)%in%base.dis$eid),]
      geno.ver=geno.use[which(rownames(geno.use)%in%base.ver$eid),]
      save(base.dis, base.ver, geno.dis,geno.ver,geno.use,map.use,file = paste0(wd, "geno_transed_chr",i,"_",j,".RData"))
     }
  }
}
