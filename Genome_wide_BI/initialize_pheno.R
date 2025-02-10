rm(list = ls()); gc(reset = TRUE)
require(snpStats) 
library(igraph)
library(data.table)
wd="" #working directory
maintable=fread("/../ukb50595.csv", select = c("eid","31-0.0","34-0.0","52-0.0","53-0.0","53-2.0","22001-0.0","22006-0.0","22009-0.1",
                                                                                              "22009-0.2","22009-0.3","22009-0.4","22009-0.5"))

bulk=read.table("/../ukb50595.bulk")

load("/gpfs/ysm/project/zhang_heping/mc3526/basket-50595/kinship_3degree.Rdata")
kinship_3degree[,1]=as.character(kinship_3degree[,1])
kinship_3degree[,2]=as.character(kinship_3degree[,2])
Kinship_network=simplify(graph_from_edgelist(as.matrix(kinship_3degree),directed=F))

delete_relative <- function(data)
{
  cases_kinship_network=induced_subgraph(Kinship_network, v=V(Kinship_network)$name %in% data$eid)
  S_G=delete.vertices(cases_kinship_network,degree(cases_kinship_network)==0)
  rm(cases_kinship_network)
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
rm(maintable,kinship_3degree,Kinship_network)
blk_info=read.csv("/../num_blk.csv",header=TRUE)
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
###########intialize the responses#############
#bulk.dis=bulk[bulk$V1%in%base.dis$eid,]
n.dis=nrow(base.dis)
#bulk.ver=bulk[bulk$V1%in%base.ver$eid,]
n.ver=nrow(base.ver)
pheno.dis=array(0,dim=c(55,55,n.dis))
pheno.ver=array(0,dim=c(55,55,n.ver))
#field="25753_2_0"
for(m in 1:n.dis){
  id=base.dis$eid[m]
  field=bulk$V2[which(bulk$V1==id)]
  field=field[1]
  filename=paste0("/gpfs/ysm/project/zhang_heping/mc3526/basket-50595/field25753/",id,"_",field,".txt")
  partial_cor=fread(filename)
  B <- diag(0,55)
  B[upper.tri(B,diag = F)]<- as.numeric(partial_cor)
  B <- B+t(B)
  pheno.dis[,,m]<-B
}

for(m in 1:n.ver){
  id=base.ver$eid[m]
  field=bulk$V2[which(bulk$V1==id)]
  field=field[1]
  filename=paste0("/gpfs/ysm/project/zhang_heping/mc3526/basket-50595/field25753/",id,"_",field,".txt")
  partial_cor=fread(filename)
  B <- diag(0,55)
  B[upper.tri(B,diag = F)]<- as.numeric(partial_cor)
  B <- B+t(B)
  pheno.ver[,,m]<-B
}
save(pheno.dis,file = paste0(wd, "dis_pheno.RData"))
save(pheno.ver,file = paste0(wd, "ver_pheno.RData"))

full_dist_mat_dis=BallDistanceVector(pheno.dis)
full_dist_mat_ver=BallDistanceVector(pheno.ver)
save(full_dist_mat_dis,file = paste0(wd, "full_dist_mat_dis.RData"))
save(full_dist_mat_ver,file = paste0(wd, "full_dist_mat_ver.RData"))