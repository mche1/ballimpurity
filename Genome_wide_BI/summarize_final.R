rm(list = ls()); gc(reset = TRUE)
parent.path<-"/../"
PG_final<-NULL

###chr1-22 #####
for(chr in 1:22){
  setwd(paste0(parent.path,"/chr",chr))
  for(j in 1:250){
    filename<-paste0("PG_results_chr",chr,"_blk",j,".RData")
    if(file.exists(filename)){
      load(filename)
      snps<-map
      n_j=nrow(PG.mat)
      PG_final<-rbind(PG_final,cbind(rep(chr,n_j),rep(j,n_j),snps$position,PG.mat))
    }
  }
}

n=floor(nrow(PG_final)/100)
top.dis=sort(PG_final[,4],decreasing=T)[1:n]
top.ver=sort(PG_final[,5],decreasing=T)[1:n]

topsnps=rownames(PG_final)[which(PG_final[,4]>=top.dis[n]&PG_final[,5]>=top.ver[n])]
top_snps=PG_final[which(PG_final[,4]>=top.dis[n]&PG_final[,5]>=top.ver[n]),]
top_snps=cbind(rownames(top_snps),top_snps)
#dis_s
colnames(top_snps)<-c("SNP","CHR","BP","position","BI_dis","BI_ver")
write.table(topsnps,file="top_snps_2perc.txt")


TS=as.data.frame(top_snps)
TS$CHR=as.numeric(TS$CHR)
TS$BP=as.numeric(TS$BP)
TS$position=as.numeric(TS$position)
TS$BI_dis=as.numeric(TS$BI_dis)
TS$BI_ver=as.numeric(TS$BI_ver)

colnames(PG_final)<-c("chr","blk","bp","discovery","verification")
PG_final_table=as.data.frame(PG_final)
data_cum <- PG_final_table %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(chr, bp_add)
########annotate the supervariants
superv=data.table(chr=c(1,3,7,8,9,9,10,12,15,16,18,19),blk=c(119,151,139,110,26,120,30,34,65,61,71,55),bp_start=rep(0,12),bp_end=rep(0,12))
for(i in 1:12){
  superv$bp_start[i]=data_cum$bp_add[superv$chr[i]]+(superv$blk[i]-1)*1e6
  superv$bp_end[i]=data_cum$bp_add[superv$chr[i]]+(superv$blk[i])*1e6
}

#####sample 10k random snps 
sam=sample(seq(from=1,to=nrow(PG_final)),size=10000)
sam=sort(sam)
chrsam=PG_final[sam,1]
rownames(chrsam)<-NULL
########get plot data 
plotdata=data.frame(chr=as.numeric(c(chrsam,TS$CHR)),bp=as.numeric(c(PG_final[sam,3],TS$position)),p = c(PG_final[sam,4],as.numeric(TS$BI_dis)),snp=c(rownames(PG_final)[sam],TS$SNP))

plotdata<-plotdata %>%
  group_by(chr) %>%
  select(chr,bp,p,snp)
plotdata <- plotdata %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)
top10.dis=sort(plotdata$p,decreasing=T)[10]
chrtops <- plotdata %>% 
  #group_by(chr) %>% 
  filter(p>=top10.dis)
genes.dis=c("Gene RAB3C","Gene MKLN1","Gene RAD54B", "","", "Gene MAML2","","Gene TTC14-DT", "LINC02994","LINC02994")
chrtops=cbind(chrtops,genes.dis)
colnames(chrtops)[7]<-"Genes"
axis_set <- plotdata %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))
library(karyoploteR)
library(tidyverse)
library(ggtext)

axisdf = plotdata %>% group_by(chr) %>% summarize(center=( max(bp_cum) + min(bp_cum) ) / 2 )
pdf("manhattan_dis_feb1.pdf",width=14,height=8)
manhplot <- ggplot(plotdata, aes(x=bp_cum, y=p)) +
  xlab("Chromesome") + ylab("Purity Gain")+
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.2) +
  geom_point(data=chrtops,colour="red",size=3)+
  ggrepel::geom_text_repel(data = chrtops, aes(label = chrtops$snp))+
  annotate("text", x =as.numeric(chrtops[1,6])-2e7, y = 1.1e-4, label = "Gene RAB3C", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[1,6]), y =as.numeric(chrtops[1,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  annotate("text", x =as.numeric(chrtops[2,6])-2e7, y = 1.1e-4, label = "Gene MKLN1", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[2,6]), y =as.numeric(chrtops[2,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  annotate("text", x =as.numeric(chrtops[3,6])-2e7, y = 1.1e-4, label = "Gene RAD54B", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[3,6]), y =as.numeric(chrtops[3,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  annotate("text", x =as.numeric(chrtops[6,6])-2e7, y = 1.1e-4, label = "Gene MAML2", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[6,6]), y =as.numeric(chrtops[6,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  annotate("text", x =as.numeric(chrtops[8,6])-2e7, y = 1.1e-4, label = "Gene TTC14-DT", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[8,6]), y =as.numeric(chrtops[8,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  annotate("text", x =as.numeric(chrtops[9,6])-2e7, y = 1.1e-4, label = "LINC02994", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[9,6]), y =as.numeric(chrtops[9,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  annotate("rect", xmin = superv$bp_start[1], xmax = superv$bp_end[1], ymin = 2e-5, ymax = 1e-4,
             alpha = .1,fill = "blue")+#}
  scale_color_manual(values = rep(c("black", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 1e-5) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
print(manhplot)
dev.off()



plotdata=data.frame(chr=as.numeric(c(chrsam,TS$CHR)),bp=as.numeric(c(PG_final[sam,3],TS$position)),p = c(PG_final[sam,5],as.numeric(TS$BI_ver)),snp=c(rownames(PG_final)[sam],TS$SNP))

plotdata<-plotdata %>%
  group_by(chr) %>%
  select(chr,bp,p,snp)
# data_cum <- plotdata %>% 
#   group_by(chr) %>% 
#   summarise(max_bp = max(bp)) %>% 
#   mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
#   select(chr, bp_add)
plotdata <- plotdata %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)

top10.ver=sort(plotdata$p,decreasing=T)[10]
chrtops <- plotdata %>% 
  #group_by(chr) %>% 
  filter(p>=top10.ver)

axis_set <- plotdata %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))
axisdf = plotdata %>% group_by(chr) %>% summarize(center=( max(bp_cum) + min(bp_cum) ) / 2 )
pdf("manhattan_ver_feb1.pdf",width=14,height=8)
manhplot <- ggplot(plotdata, aes(x=bp_cum, y=p)) +
  xlab("Chromesome") + ylab("Purity Gain")+
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.2) +
  geom_point(data=chrtops,colour="red",size=3)+
  ggrepel::geom_text_repel(data = chrtops, aes(label = chrtops$snp))+
  #annotate("text", x =as.numeric(chrtops[1,6])-2e7, y = 1.1e-4, label = "Gene RAB3C", angle=90)+
  #annotate("pointrange", x = as.numeric(chrtops[1,6]), y =as.numeric(chrtops[1,3]), ymin = 2e-5, ymax = 1.2e-4,
  #         colour = "green", size=0.2)+
  annotate("text", x =as.numeric(chrtops[2,6])-2e7, y = 1.1e-4, label = "Gene TRIM71", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[2,6]), y =as.numeric(chrtops[2,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size = 0.2)+
  annotate("text", x =as.numeric(chrtops[3,6])+2e7, y = 1.1e-4, label = "Gene CDCP1", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[3,6]), y =as.numeric(chrtops[3,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size = 0.2)+
  annotate("text", x =as.numeric(chrtops[4,6])-2e7, y = 1.1e-4, label = "LINC02380", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[4,6]), y =as.numeric(chrtops[4,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size = 0.2)+
  annotate("text", x =as.numeric(chrtops[5,6])-2e7, y = 1.1e-4, label = "Gene HIVEP1", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[5,6]), y =as.numeric(chrtops[5,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  annotate("text", x =as.numeric(chrtops[7,6])-2e7, y = 1.1e-4, label = "Gene TEK", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[7,6]), y =as.numeric(chrtops[7,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  annotate("text", x =as.numeric(chrtops[8,6])-2e7, y = 1.1e-4, label = "Gene GABBR2", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[8,6]), y =as.numeric(chrtops[8,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  annotate("text", x =as.numeric(chrtops[9,6])-2e7, y = 1.1e-4, label = "LOC105376637", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[9,6]), y =as.numeric(chrtops[9,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  annotate("text", x =as.numeric(chrtops[10,6])-2e7, y = 1.1e-4, label = "LOC107984706", angle=90)+
  annotate("pointrange", x = as.numeric(chrtops[10,6]), y =as.numeric(chrtops[10,3]), ymin = 2e-5, ymax = 1.2e-4,
           colour = "green", size=0.2)+
  scale_color_manual(values = rep(c("black", "skyblue"), 22 )) +
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 1e-5) ) +     # remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
print(manhplot)
dev.off()
