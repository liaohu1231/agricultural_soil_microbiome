data=read.csv("all.ADVA_agri.rpkm",header = 1,row.names = 1,check.names = F,sep='\t')
group=read.csv("metadata.txt",header = 1,stringsAsFactors = T,sep="\t",row.names = 1)

library(ggplot2)
library(ggsci)
relative_abundance=function(d){
  dta_sum=apply(d,2,function(x){x/sum(x)})
}
abu=relative_abundance(data)
abu=as.data.frame(abu)
library(ggbeeswarm)
##PCoA
library(dplyr)
library(vegan)
results5=t(abu)
results5=as.matrix(results5)
library(vegan)
distance <- as.matrix(vegdist(results5, method= "bray",na.rm = T))
colnames(distance)=rownames(results5)
rownames(distance)=rownames(results5)
distance2=matrix(NA,nrow(distance),ncol(distance))
for(i in 1:nrow(distance)){
  for(j in 1:ncol(distance)){
    if ( i != j){
      distance2[i,j]=distance[i,j]
    }
  }
}
rownames(distance2)=rownames(distance)
colnames(distance2)=colnames(distance)

write.table(x = distance2,"viral_distance.tsv",quote = F,row.names = T,sep='\t')
library(pheatmap)
distance2=1-distance2
pheatmap(distance2,cluster_rows = T,cluster_cols = T,
         clustering_method = "complete",
         na_col = "lightblue",border_color = F,cellwidth = 8,cellheight = 6,
         color = colorRampPalette(c("white","gray","#40CCA7"))(100))

pcoa_v <- cmdscale(distance, k = (nrow(results5) - 1), eig = TRUE)

mds=metaMDS(distance)
point_v=data.frame(mds$point)

#point_v <- data.frame(pcoa_v$point)
pcoa_eig_v <- (pcoa_v$eig)[1:2] / sum(pcoa_v$eig)

sample_site_v <- data.frame (point_v)[1:2]
sample_site_v$names <- rownames(sample_site_v)
names(sample_site_v)[1:2] <- c('PCoA1', 'PCoA2')
sample_site_v$group = group[rownames(sample_site_v),9]

library(vegan)
dist = vegdist(t(abu), method = 'bray',na.rm = T)

p_value=anosim(x=dist,grouping = sample_site_v$group,permutations = 9999)
summary(p_value)
p_value_1=adonis2(formula = dist~group,sample_site_v,permutations = 9999)
p_value_1

sample_site_v$depth=group[rownames(sample_site_v),8]
p_value=anosim(x= dist,grouping = sample_site_v$depth,permutations = 9999)
summary(p_value)
p_value_1=adonis2(formula = dist~depth,sample_site_v,permutations = 9999)
p_value_1

p <- ggplot(sample_site_v, aes(PCoA1, PCoA2)) +
  theme(panel.grid = element_line(color = 'gray', linetype = 2), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(aes(fill=group,shape=depth),size = 2.5, alpha = 0.8) + 
  scale_fill_d3()+
  scale_color_d3()+
  scale_shape_manual(values = c(21,22,23,24))+
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig_v[1], 2), '%'), 
       y = paste('PCoA axis2: ', round(100 * pcoa_eig_v[2], 2), '%'))+
  #stat_ellipse(data = sample_site_v,mapping = aes(PCoA1, PCoA2,group = depth,color=depth),
               #level = 0.9, show.legend = TRUE,inherit.aes = F)+
  stat_ellipse(data = sample_site_v,mapping = aes(PCoA1, PCoA2,group = site),level = 0.975, show.legend = TRUE,inherit.aes = F)+
  annotate('text', colour="#8766d")+
  guides(fill = guide_legend(override.aes = list(size = 5, shape = 21)))
p
write.table(sample_site_v,"PCoA_table.csv",sep="\t",quote = F)

##MDS
mds=metaMDS(distance)
point_v=data.frame(mds$point)
pcoa_eig_v <- (mds$eig)[1:2] / sum(mds$eig)
sample_site_v <- data.frame (point_v)[1:2]
sample_site_v$names <- rownames(sample_site_v)
names(sample_site_v)[1:2] <- c('PCoA1', 'PCoA2')
sample_site_v$group = group[rownames(sample_site_v),9]
sample_site_v$depth=group[rownames(sample_site_v),8]

p <- ggplot(sample_site_v, aes(PCoA1, PCoA2)) +
  theme(panel.grid = element_line(color = 'gray', linetype = 2), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(aes(fill=group),shape=21,size = 2.5, alpha = 0.8) + 
  scale_fill_d3()+
  scale_color_d3()+
  labs(x = paste('MDS1'), 
       y = paste('MDS2'))+
  #stat_ellipse(data = sample_site_v,mapping = aes(PCoA1, PCoA2,group = depth,color=depth),
  #level = 0.9, show.legend = TRUE,inherit.aes = F)+
  annotate('text', colour="#8766d")+
  guides(fill = guide_legend(override.aes = list(size = 5, shape = 21)))
p
write.table(sample_site_v,"PCoA_table.csv",sep="\t",quote = F)
library(ggpubr)
ggplot(sample_site_v,aes(group,PCoA2))+
  geom_boxplot(aes(fill=group),outliers = F)+
  geom_point(aes(fill=group),position = position_jitter(width = 0.2),
             alpha=0.8,shape=21)+
  scale_fill_d3()+
  theme_classic()+
  stat_compare_means(comparisons = list(c("paddy","crop")))+
  labs(y="MDS2")
  
ggsave("box_mds2.pdf",device = "pdf",height = 3.5,width = 3.3)


##alpha diversity
library(vegan)
library(dplyr)
shannon=diversity(t(abu[,1:24]))%>%as.data.frame()
colnames(shannon)="shannon_index"
shannon$group=group[rownames(shannon),9]
shannon$depth=group[rownames(shannon),8]
shannon_aver=rowsum(shannon$shannon_index,shannon$group)/3
shannon_aver=as.data.frame(shannon_aver)

library(ggplot2)
library(ggpubr)
ggplot(shannon,aes(x = group,y=shannon_index))+
  geom_boxplot(aes(fill=depth))+
  theme_classic()+
  scale_fill_d3()

write.table(shannon,"shannon_index.tsv",sep="\t",quote = F)

##mantel
mental_matrix=function(x,y){
  nd=matrix(data = NA,ncol(x),2)
  for (i in 1:ncol(x)){
    print(i)
    scale.env=scale(x[,i],center = T,scale = T)
    dist_env=dist(scale.env,method="euclidean")
    abund_env=mantel(xdis = y,ydis = dist_env,method = "spearman",permutations = 9999,na.rm = T)
    R=cbind(abund_env$statistic,abund_env$signif)
    print(R)
    nd[i,]=R
  }
  rownames(nd)=colnames(x)
  colnames(nd)=c("R","p")
  nd=nd[order(nd[,1],decreasing = T),]
  nd=as.data.frame(nd)
  nd$p_adjust=p.adjust(nd$p,method = "fdr")
  nd=as.data.frame(nd)
  nd$env=as.factor(rownames(nd))
  nd$env=factor(nd$env,levels = as.factor(nd$env))
  return(nd)
}
group$depth_value=rep(c(25,50,75,100),6)
group1=group[,-c(8:13,20:62)]
mantel_1=mental_matrix(group1,dist)
hosts2[,8:31]=abu[as.character(hosts2$viral_contigs),1:24]
host_linked=aggregate.data.frame(hosts2[,8:31],by = list(hosts2$class),FUN = sum)
host_linked$Group.1[1]="unclassified"
rownames(host_linked)=host_linked$Group.1
host_linked=host_linked[,-1]
##host and bacterial distance
microbial=read.table("../prokaryotic/ADVD_MAGs_abu.tsv",sep='\t',header = T)
library(reshape2)
microbial=microbial[,c(1:3)]
microbial=dcast(microbial,formula = Sample_file~Genome_file)
rownames(microbial)=microbial$Sample_file
microbial=microbial[,-1]
microbial[is.na(microbial)]=0

dist_abund=vegdist(microbial,method = "bray")%>%as.matrix()
dist_host_abund=vegdist(t(host_linked),method = "bray")%>%as.matrix()
dist_population = vegdist(t(abu), method = 'bray',na.rm = T)%>%as.matrix()
##as.numeric of matrix
mantal_viral_bac=mantel(dist_population,dist_abund,
                        method = "spearman",
                        permutations = 999,na.rm = T)
mantal_viral_bac$signif

mantal_viral_host_linked=mantel(dist_population,dist_host_abund,
                                method = "spearman",
                                permutations = 999,na.rm = T)
mantal_viral_host_linked$signif

mantal_viral_statistic=rbind(mantal_viral_bac$statistic,
                             mantal_viral_host_linked$statistic)%>%as.data.frame()

colnames(mantal_viral_statistic)="mantal_statistic_R"
mantal_group=c("Host-linked","Bacterial")
mantal_viral_statistic=cbind(mantal_viral_statistic,mantal_group)

###ggcor
distance=as.matrix(vegdist(microbial, method= "bray",na.rm = T))
pcoa_v <- cmdscale(distance, k = (nrow(microbial) - 1), eig = TRUE)
point_v <- data.frame(pcoa_v$point)
sample_site_v <- data.frame ({pcoa_v$point})[1]
colnames(sample_site_v)=c("bacterial_PCoA1")
spec=cbind(t(host_linked),microbial)
env_f=cbind(group1,sample_site_v)
library(ggcor)
mantel_1=mantel_test(spec = spec,env= env_f,
                     spec.select = list(
                       bacteria_abundance=55:511
                     ))
mantel_2=mantel_test(spec = spec,env= env_f,
                     spec.select = list(
                       host_linked_abundance=1:54
                     ))
mantel_3=mantel_test(spec = t(abu),env= env_f,
                     spec.select = list(
                       viral_communities=1:nrow(abu)
                     ))
mantel_1$p_adjust=p.adjust(mantel_1$p.value,method = "fdr")

mantel_2$p_adjust=p.adjust(mantel_2$p.value,method = "fdr")
mantel_3$p_adjust=p.adjust(mantel_3$p.value,method = "fdr")
mantel=rbind(mantel_1,mantel_2,mantel_3)

mantel=mantel%>%mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                                labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
                       pd = cut(p_adjust, breaks = c(-Inf, 0.01, 0.05, Inf),
                                labels = c("< 0.01", "< 0.05", "0.05 - 0.1")))

mantel=mantel[mantel$p.value<0.1,]
library(ggplot2)
library(ggraph)
col2 <- colorRampPalette(c("#77C034","white" ,"#C388FE"),alpha = TRUE)

p=quickcor(env_f,type="upper")+
  geom_circle2()+
  add_link(df = mantel,mapping = aes(colour = pd,size=rd))+
  scale_fill_gradient2(midpoint = 0, low = "#1C5999", 
                       mid = "white",high = "#720021", space = "Lab" )+
  scale_size_manual(values = c(0.2, 1.5, 3))+
  scale_colour_manual(values = c("#D95F02", "#1B9E77","red")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 4), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 5))
p
ggsave(filename = "quickcor2.pdf",width = 10,height = 8,device = 'pdf')
#env

library(ggplot2)
library(dplyr)
library(ggsci)
library(ggpubr)
plot_data=group %>%
  
  group_by(group, Depth) %>%
  
  summarize(
    Mean_Value = mean(DNA.concerntration),
    SE_Value = sd(DNA.concerntration) / sqrt(n())
  )

ggplot(plot_data,aes(Depth,Mean_Value,fill=group))+
  geom_bar(stat='identity',position="dodge", alpha=.9, 
           colour='black', width = .6)+
  geom_errorbar(aes(ymin=Mean_Value-SE_Value, 
                    ymax=Mean_Value+SE_Value),
                stat = 'identity', 
                position = position_dodge(0.6), 
                width=.1)+
  theme_bw()+
  scale_fill_d3()+
  geom_signif(y_position=13.5, xmin=c(2.85), 
              xmax=c(3.15),
              annotation=c("**"),
              tip_length=0.03, size=0.7, 
              textsize = 5,vjust = 0.05)+
  ylim(0,180)+
  labs(y="DNA concern-\ntration (ng/uL)")+
  theme(axis.text.x = element_text(angle = 40,hjust = 0.8))
plot_data=group %>%
  group_by(group, Depth) %>%
  summarize(
    Mean_Value = mean(TOC.mg.g),
    SE_Value = sd(TOC.mg.g) / sqrt(n())
  )

ggplot(plot_data,aes(Depth,Mean_Value,fill=group))+
  geom_bar(stat='identity',position="dodge", alpha=.9, 
           colour='black', width = .6)+
  geom_errorbar(aes(ymin=Mean_Value-SE_Value, 
                    ymax=Mean_Value+SE_Value),
                stat = 'identity', 
                position = position_dodge(0.6), 
                width=.1)+
  theme_bw()+
  scale_fill_d3()+
  geom_signif(y_position=13.5, xmin=c(2.85), 
              xmax=c(3.15),
              annotation=c("**"),
              tip_length=0.03, size=0.7, 
              textsize = 5,vjust = 0.05)+
  ylim(0,25.5)+
  labs(y="TOC (mg/g)")+
  theme(axis.text.x = element_text(angle = 40,hjust = 0.8))
plot_data=group %>%
  
  group_by(group, Depth) %>%
  
  summarize(
    Mean_Value = mean(C.mg.g),
    SE_Value = sd(C.mg.g) / sqrt(n())
  )

ggplot(plot_data,aes(Depth,Mean_Value,fill=group))+
  geom_bar(stat='identity',position="dodge", alpha=.9, 
           colour='black', width = .6)+
  geom_errorbar(aes(ymin=Mean_Value-SE_Value, 
                    ymax=Mean_Value+SE_Value),
                stat = 'identity', 
                position = position_dodge(0.6), width=.1)+
  theme_bw()+
  scale_fill_d3()+
  geom_signif(y_position=23.5, xmin=c(0.85), xmax=c(1.15),
              annotation=c("*"),
              tip_length=0.03, size=0.7, 
              textsize = 5,vjust = 0.05)+ylim(0,25.5)+
  labs(y="Carbon (mg/g)")+
  theme(axis.text.x = element_text(angle = 40,hjust = 0.8))
  
group[grep("25",group$depth_value),1]
t.test(group[group$depth_value == 25 & 
               group$group =="crop",2],
            group[group$depth_value == 25 & 
                    group$group =="paddy",2],
       alternative = "less")


plot_data=group %>%
  
  group_by(group, Depth) %>%
  
  summarize(
    Mean_Value = mean(N.mg.g),
    SE_Value = sd(N.mg.g) / sqrt(n())
  )

ggplot(plot_data,aes(Depth,Mean_Value,fill=group))+
  geom_bar(stat='identity',position="dodge", alpha=.9, 
           colour='black', width = .6)+
  geom_errorbar(aes(ymin=Mean_Value-SE_Value, 
                    ymax=Mean_Value+SE_Value),
                stat = 'identity', 
                position = position_dodge(0.6), width=.1)+
  theme_bw()+
  scale_fill_d3()+
  geom_signif(y_position=1.7, xmin=c(0.85), xmax=c(1.15),
              annotation=c("*"),
              tip_length=0.03, size=0.7, 
              textsize = 5,vjust = 0.05)+ylim(0,2)+
  geom_signif(y_position=0.7, xmin=c(2.85), xmax=c(3.15),
              annotation=c("*"),
              tip_length=0.03, size=0.7, 
              textsize = 5,vjust = 0.05)+ylim(0,2)+
  labs(y="Nitrogen (mg/g)")+
  theme(axis.text.x = element_text(angle = 40,hjust = 0.8))
library(epitools)
wilcox.test(group[group$depth == "subsoil" & 
               group$group =="crop",3],
       group[group$depth == "subsoil" & 
               group$group =="paddy",3],
       alternative = "greater")

plot_data=group %>%
  group_by(group, Depth) %>%
  summarize(
    Mean_Value = mean(S.mg.g),
    SE_Value = sd(S.mg.g) / sqrt(n())
  )

ggplot(plot_data,aes(Depth,Mean_Value,fill=group))+
  geom_bar(stat='identity',position="dodge", alpha=.9, 
           colour='black', width = .6)+
  geom_errorbar(aes(ymin=Mean_Value-SE_Value, 
                    ymax=Mean_Value+SE_Value),
                stat = 'identity', 
                position = position_dodge(0.6), width=.1)+
  theme_bw()+
  scale_fill_d3()+
  labs(y="Sulfur (mg/g)")+
  theme(axis.text.x = element_text(angle = 40,hjust = 0.8))
library(epitools)
t.test(group[group$depth_value == 50 & 
               group$group =="crop",4],
       group[group$depth_value == 50& 
               group$group =="paddy",4],
       alternative = "less")
plot_data=group %>%
  
  group_by(group, Depth) %>%
  
  summarize(
    Mean_Value = mean(C.mg.g),
    SE_Value = sd(C.mg.g) / sqrt(n())
  )

ggplot(plot_data,aes(Depth,Mean_Value,fill=group))+
  geom_bar(stat='identity',position="dodge", alpha=.9, 
           colour='black', width = .6)+
  geom_errorbar(aes(ymin=Mean_Value-SE_Value, 
                    ymax=Mean_Value+SE_Value),
                stat = 'identity', 
                position = position_dodge(0.6), width=.1)+
  theme_bw()+
  scale_fill_d3()+
  geom_signif(y_position=23.5, xmin=c(0.85), xmax=c(1.15),
              annotation=c("*"),
              tip_length=0.03, size=0.7, 
              textsize = 5,vjust = 0.05)+ylim(0,25.5)+
  labs(y="Carbon (mg/g)")+
  theme(axis.text.x = element_text(angle = 40,hjust = 0.8))

##host
hosts=read.table("../host/viruses_host_uniq.tsv",sep='\t',header = F)
colnames(hosts)=c("viral_contigs","MAGs_contigs","value","methods","MAGs")
hosts=as.data.frame(hosts)
gtdbtk=read.csv("../prokaryotic/gtdbtk.bac120.summary.tsv",sep='\t')
tmp=strsplit(gtdbtk$classification,split = ";")
tmp=do.call(rbind,tmp)
gtdbtk$phylum=tmp[,2]
gtdbtk$class=tmp[,3]
rownames(gtdbtk)=gtdbtk$user_genome
hosts$phylum=gtdbtk[hosts$MAGs,21]
hosts$phylum=gsub("p__","",hosts$phylum)
hosts$class=gtdbtk[hosts$MAGs,22]
hosts$class=gsub("c__","",hosts$class)
hosts2=hosts[!duplicated(hosts[,c(1,5)]),]
write.table(hosts2,"Table_S2_hosts2.tsv",sep='\t',quote = F,row.names = T)
nrow(hosts[!duplicated(hosts[,5]),])

host_number=as.data.frame(table(hosts2$viral_contigs))
host_number2=as.data.frame(table(hosts2$viral_contigs,hosts2$phylum))
host_number=host_number2[host_number2$Freq>0,]
host_number2=host_number2[host_number2$Freq>3,]
host_number3=matrix(data = NA,nrow = 0,ncol = ncol(host_number))
for(i in 1:nrow(host_number2)){
  tmp=host_number[grep(host_number2[i,1],host_number$Var1),]
  host_number3=rbind(host_number3,tmp)
}


host_number3=host_number3[order(host_number3$Freq,decreasing = T),]
host_number2=host_number2[order(host_number2$Freq,decreasing = T),]
host_number2$Var1=factor(as.factor(host_number2$Var1),
                         levels = host_number2$Var1)
host_number3$Var1=factor(as.factor(host_number3$Var1),
                             levels = host_number2$Var1)

library(ggplot2)
library(ggsci)
library(tidyplots)
ggplot(host_number3,aes(Var1,Freq,fill=Var2))+
  geom_bar(stat = 'identity',width = .7)+
  theme_classic()+
  scale_fill_brewer(palette = "Set3")+
  labs(y="number of host MAGs",
       x="viral contigs")+
  theme(axis.text.x = element_blank(),
        legend.key.size = unit(0.35, "cm"))
  
library(ggplot2)
phy_number=as.data.frame(table(hosts2$phylum))
phy_number=phy_number[order(phy_number$Freq,decreasing = T),]
phy_number=phy_number[phy_number$Freq>10,]
phy_number$Var1=factor(phy_number$Var1,levels = as.factor(phy_number$Var1))
ggplot(phy_number,aes(Var1,Freq))+
  geom_bar(stat = 'identity',fill="#EFD7E7")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  labs(y="number of vOTUs")

ggplot(phy_number[phy_number$Freq<100,],
       aes(Var1,Freq))+
  geom_bar(stat = 'identity',fill="#EFD7E7")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

pates=hosts[grep(pattern = "Pates",hosts$phylum),]
pates=pates[!duplicated(pates$viral_contigs),]
write.table(pates,file = "Pates_virues.tsv",sep='\t',
            row.names = F,quote = F)
pates[,7:30]=abu[pates$viral_contigs,]
pates=data.frame(pates=apply(pates[,7:30],2,sum))
pates=pates*100
pates$depth=group[rownames(pates),8]
pates$group=group[rownames(pates),9]
library(dplyr)
plot_data=pates %>%
  group_by(group, depth) %>%
  summarize(
    Mean_Value = mean(pates),
    SE_Value = sd(pates) / sqrt(n())
  )

ggplot(plot_data,aes(depth,Mean_Value,fill=group))+
  geom_bar(stat='identity',position="dodge", alpha=.9, 
           colour='black', width = .6)+
  geom_errorbar(aes(ymin=Mean_Value-SE_Value, 
                    ymax=Mean_Value+SE_Value),
                stat = 'identity', 
                position = position_dodge(0.6), width=.1)+
  theme_bw()+
  scale_fill_d3()+
  labs(y="Relative abundance (%) \nof Patescibacteria viruses")+
  theme(axis.text.x = element_text(angle = 40,
                                   hjust = 0.8))
t.test(pates[pates$depth == "75-100 cm" & 
                    pates$group =="crop",1],
            pates[pates$depth == "75-100 cm" & 
                    pates$group =="paddy",1],
            alternative = "less")

##Halo_viruses
Halo=hosts[grep(pattern = "Halo",hosts$phylum),]
Halo[,7:30]=abu[Halo$viral_contigs,]
Halo=Halo[!duplicated(Halo$viral_contigs),]
write.table(Halo,"ANME_2D_virues.tsv",sep='\t',quote = F,row.names = F)

rownames(Halo)=Halo$viral_contigs
Halo=Halo[,7:30]*100
library(reshape2)
Halo=melt(as.matrix(Halo))
Halo=Halo[Halo$value>0,]
Halo$depth=group[Halo$Var2,8]
Halo$group=group[Halo$Var2,9]
library(dplyr)

library(ggpubr)
ggplot(Halo,aes(depth,log10(value),fill=group))+
  geom_boxplot(outliers = F,alpha=0.6)+
  geom_point(shape=21,stat='identity',alpha=0.6,
             position = position_dodge(width = 1))+
  theme_bw()+
  scale_fill_d3()+
  labs(y="Log10(Relative abundance (%) \nof ANME-2d viruses)")+
  theme(axis.text.x = element_text(angle = 40,
                                   hjust = 0.8))+
  stat_compare_means(comparisons = list(c("0-25 cm",
                                          "50-75 cm"),
                                        c("50-75 cm",
                                          "75-100 cm")))+
  ylim(c(-3,0.5))

##
Halo=hosts[grep(pattern = "Methylomirabilia",
                 hosts$class),]
Halo[,8:31]=abu[Halo$viral_contigs,1:24]
Halo=Halo[!duplicated(Halo$viral_contigs),]
write.table(Halo,"Methylomirabilia_virues.tsv",sep='\t',quote = F,row.names = F)

rownames(Halo)=Halo$viral_contigs
Halo=Halo[,8:31]*100
library(reshape2)
Halo=data.frame(Methylomirabilia=apply(Halo,2,sum))
Halo$depth=group[rownames(Halo),8]
Halo$group=group[rownames(Halo),9]
library(dplyr)
library(ggsci)
library(ggpubr)
ggplot(Halo,aes(depth,Methylomirabilia,fill=group))+
  geom_boxplot(outliers = F,alpha=0.6)+
  geom_point(shape=21,stat='identity',alpha=0.6,
             position = position_dodge(width = 0.8))+
  theme_bw()+
  scale_fill_d3()+
  labs(y="Relative abundance (%) of 
       all Methylomirabilia viruses")+
  theme(axis.text.x = element_text(angle = 40,
                                   hjust = 0.8))+
  stat_compare_means(comparisons = list(c("0-25 cm",
                                          "50-75 cm")))
##aob_viruses

Halo=hosts[grep(pattern = "Nitros",
                 hosts$class),]
Halo[,8:31]=abu[Halo$viral_contigs,1:24]
Halo=Halo[!duplicated(Halo$viral_contigs),]
write.table(Halo,"AOB_virues.tsv",sep='\t',quote = F,row.names = F)

rownames(Halo)=Halo$viral_contigs
Halo=Halo[,8:31]*100
library(reshape2)
Halo=data.frame(AOB=apply(Halo,2,sum))
Halo$depth=group[rownames(Halo),8]
Halo$group=group[rownames(Halo),9]
library(dplyr)
library(ggsci)
library(ggpubr)
ggplot(Halo,aes(depth,Methylomirabilia,fill=group))+
  geom_boxplot(outliers = F,alpha=0.6)+
  geom_point(shape=21,stat='identity',alpha=0.6,
             position = position_dodge(width = 0.8))+
  theme_bw()+
  scale_fill_d3()+
  labs(y="Relative abundance (%) of 
       all Nitrospirae viruses")+
  theme(axis.text.x = element_text(angle = 40,
                                   hjust = 0.8))+
  stat_compare_means(comparisons = list(c("0-25 cm",
                                          "25-50 cm")))


###lifestyle
lifestyle=read.table("phatyp_prediction.tsv",sep = '\t',header = 1,row.names = 1)
lifestyle[lifestyle$PhaTYPScore==0,2]="unclassified"
integrase=read.table("../functions/ADVD_lysogenic.tsv",row.names = 1)
lifestyle$inte=integrase[rownames(lifestyle),1]
virues=as.data.frame(abu)
virues$lifestyle=lifestyle[rownames(virues),2]
table(lifestyle$TYPE)

abu_life=aggregate.data.frame(virues[,1:24],by = list(virues$lifestyle),sum)
rownames(abu_life)=abu_life$Group.1
library(reshape2)
abu_life=melt(abu_life)
library(ggplot2)
plot_data=abu_life

plot_data$depth_value=as.factor(group[plot_data$variable,66])
plot_data$landuses=as.factor(group[plot_data$variable,9])
plot_data=plot_data[grep("virulent",plot_data$Group.1),]
mean(plot_data[grep("1",plot_data$variable),3])
mean(plot_data[grep("4",plot_data$variable),3])/mean(plot_data[grep("1",plot_data$variable),3])
library(ggsci)
library(ggpubr)
p9=ggplot(plot_data,aes(depth_value,value),
)+
  geom_boxplot(aes(fill=landuses))+
  theme_classic()+
  geom_point(aes(fill=depth_value),shape=21,
             position = position_jitter(width =0.2))+
  #facet_wrap(~depth_value,ncol = 6)+
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 40,hjust = 1))+
  labs(y="Relative abundance \nof lytic phage (%)")+
  stat_compare_means(paired = T,
                     comparisons = list(c("25","50"),
                                        c("25","100"),
                                        c("50","100")))+
  theme(legend.position = "none")+
  labs(x="Depth (cm)")

p9
rownames(plot_data)=plot_data$variable
group$lytic_abundance=plot_data[rownames(group),3]
dependent_vars <- colnames(group)[1:66]
model_list <- list()
for (var in dependent_vars) {
  formula <- as.formula(paste(var, "~lytic_abundance"))
  model_list[[var]] <- lm(formula, data = group)
}

# 创建一个空数据框用于存储统计量
results_df <- data.frame(
  model = character(),
  r_squared = numeric(),
  adj_r_squared = numeric(),
  f_statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# 循环提取每个模型的摘要信息
for (model_name in names(model_list)) {
  model <- model_list[[model_name]]
  s <- summary(model)
  
  results_df <- rbind(results_df, data.frame(
    model = model_name,
    r_squared = s$r.squared,
    adj_r_squared = s$adj.r.squared,
    f_statistic = s$fstatistic[1],
    p_value = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE),
    stringsAsFactors = FALSE
  ))
}

ggplot(group,aes(depth_value,lytic_abundance))+
  geom_point(aes(fill=group),stat = 'identity',size=4,shape = 21)+
  stat_smooth(method = "lm")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()

p1=ggplot(group,aes(TOC.mg.g,lytic_abundance))+
  geom_point(aes(fill=depth),stat = 'identity',size=4,shape = 21)+
  stat_smooth(aes(fill=depth),method = "lm")+
  #stat_smooth(method = "loess")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()+
  #theme(legend.position = "none")+
  labs(x="TOC (mg/g)")+
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), 
                      sep = "~`, `~")),
    label.x = 5,  # 标签的X轴位置
    label.y = 0.37,# 标签的Y轴位置
    size=3.5
  )
p1
library(ggplot2)
library(dplyr)
library(purrr)
library(broom)
library(tidyverse)
#group$lytic_abundance=group$lytic_abundance/100
model_stats <- group %>%
  group_by(depth) %>% # 分组变量：cyl
  nest() %>%
  mutate(
    # 拟合与 `geom_smooth(method=“lm”)` 相同的模型
    model = map(data, ~ lm(TOC.mg.g ~ lytic_abundance, data = .x)),
    glance = map(model, glance)
  ) %>%
  unnest(glance) %>%
  select(depth, r.squared, p.value) %>% # 选择关键统计量
  # 创建用于geom_text的标签列（parse=TRUE实现数学格式）
  mutate(
    label = sprintf("italic(R)^2 == %.3f\nitalic(p) == %.3f", 
                    r.squared, p.value)
  )

# 2. 为每个标签计算或指定一个合适的坐标位置
# 方法A：基于每组数据范围动态计算（推荐，自动适应数据）
label_pos <- group %>%
  group_by(depth) %>%
  summarise(
    x_pos = max(TOC.mg.g),   # 例如，放在每组x最大值处
    y_pos = min(lytic_abundance)   # 例如，放在每组y最小值处
  ) %>%
  left_join(model_stats, by = "depth") # 合并统计量标签

p10=ggplot(group,aes(TOC.mg.g,lytic_abundance))+
  geom_point(aes(fill=depth),stat = 'identity',size=4,shape = 21)+
  stat_smooth(aes(fill=depth),method = "lm")+
  #stat_smooth(method = "loess")+
  scale_fill_brewer(palette = "Set2")+
  theme_classic()+
  #theme(legend.position = "none")+
  labs(x="TOC (mg/g)",y="Relative abundance \n virulent phage (%)")+
  geom_text(
    data = label_pos, # 指定标注数据源
    aes(x = x_pos, y = y_pos, 
        label = label, 
        color = factor(depth)), # 颜色映射须与主图一致
    parse = TRUE,             # 解析数学表达式
    show.legend = FALSE,      # 文本不重复出现在图例
    hjust = 1, vjust = 0,     # 调整文本对齐（右上角：hjust=1, vjust=0）
    size = 6
  )

p10
p2=ggplot(group,aes(DNA.concerntration,lytic_abundance))+
  geom_point(aes(fill=group),stat = 'identity',size=4,shape = 21)+
  stat_smooth(method = "lm")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="DNA (ng/ml)")+
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), 
                      sep = "~`,`~")),
    label.x = 7,  # 标签的X轴位置
    label.y = 0.37,# 标签的Y轴位置
    size=3.5
  )
p2

p3=ggplot(group,aes(N.mg.g,lytic_abundance))+
  geom_point(aes(fill=group),stat = 'identity',
             size=4,shape = 21)+
  stat_smooth(method = "lm")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="N (mg/g)")+
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), 
                      sep = "~`,`~")),
    label.x = 0.4,  # 标签的X轴位置
    label.y = 0.37,# 标签的Y轴位置
    size=3.5
  )
p3
p4=ggplot(group,aes(Sn,lytic_abundance))+
  geom_point(aes(fill=group),stat = 'identity',size=4,shape = 21)+
  stat_smooth(method = "lm")+
  stat_smooth(aes(fill=depth),method = "lm")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()+
  theme(legend.position = "none") +
  labs(x="Sn (ppm)")+
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), 
                      sep = "~`,`~")),
    label.x = 7,  # 标签的X轴位置
    label.y = 37,# 标签的Y轴位置
    size=3.5
  )
p4
p5=ggplot(group,aes(S042.,lytic_abundance))+
  geom_point(aes(fill=depth),stat = 'identity',size=4,shape = 21)+
  stat_smooth(method = "gam")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="Sulfate (mg/Kg)")+
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), 
                      sep = "~`,`~")),
    label.x = 100,  # 标签的X轴位置
    label.y = 37,# 标签的Y轴位置
    size=3.5
  )
p5
p11=ggplot(group,aes(	
  Al,lytic_abundance))+
  geom_point(aes(fill=depth),stat = 'identity',size=4,shape = 21)+
  stat_smooth(method = "lm")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="Cr (ppm)")+
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), 
                      sep = "~`,`~")),
    label.x = 0.1,  # 标签的X轴位置
    label.y = 37,# 标签的Y轴位置
    size=3.5
  )
p11

##多重线形相关模型
model=lm(formula = lytic_abundance~Sn+TOC.mg.g+	
     Ag.mg.kg.+Al+pH+Th+As+La+depth_value+Cu,data = group)
summary(model)
par(mfrow = c(2, 2))
plot(model)

library(ggplot2)
library(gridExtra)

# 创建数据框
diagnostic_data <- data.frame(
  fitted = fitted(model),
  residuals = residuals(model),
  standardized = rstandard(model),
  cooks = cooks.distance(model)
)

# 残差 vs 拟合值图
p1 <- ggplot(diagnostic_data, aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "残差 vs 拟合值", x = "拟合值", y = "残差")

# Q-Q图
p2 <- ggplot(diagnostic_data, aes(sample = standardized)) +
  stat_qq() + stat_qq_line(color = "red") +
  labs(title = "正态Q-Q图", x = "理论分位数", y = "样本分位数")+theme_bw()

# 尺度-位置图
p3 <- ggplot(diagnostic_data, aes(x = fitted, y = sqrt(abs(standardized)))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "尺度-位置图", x = "拟合值", y = "标准化残差的平方根")

# Cook's距离图
p4 <- ggplot(diagnostic_data, aes(x = seq_along(cooks), y = cooks)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Cook's distance", x = "观测序号", y = "Cook's distance")
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol = 2)

library(ggplot2)
ggplot(data.frame(actual = group$lytic_abundance, predicted = fitted(model)), 
       aes(x = actual, y = predicted)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", alpha = 0.5) +
  labs(title = "Actual vs predicted value",
       subtitle = paste("R² =", round(summary(model)$r.squared, 3)),
       x = "Actual abudnace", y = "Predicted abundance") +
  theme_minimal()



model=lm.fit(formula = lytic_abundance~Sn+TOC.mg.g+	
     Ag.mg.kg.+Al+pH+Th+As+La+depth_value+Cu,data = group)%>%summary()
loess(formula = lytic_abundance~Sn+TOC.mg.g+	
        Ag.mg.kg.+Al,data = group)%>%summary()
ggplot(group,aes(NO3.,lytic_abundance))+
  geom_point(aes(fill=group),stat = 'identity',size=4,shape = 21)+
  stat_smooth(method = "lm")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()+
  theme(legend.position = "none")

p6=ggplot(group,aes(pH,lytic_abundance))+
  geom_point(aes(fill=group),stat = 'identity',size=4,shape = 21)+
  stat_smooth(method = "lm")+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()+
  theme(legend.position = "none")+
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), 
                      sep = "~`,`~")),
    label.x = 6,  # 标签的X轴位置
    label.y = 0.37,# 标签的Y轴位置
    size=3.5
  )
p6

##proviruses
virues=as.data.frame(abu)
proviruses=read.table("proviruses_name.txt",row.names = 1)
virues$proviruses=proviruses[rownames(virues),1]
abu_pro=virues[grep("proviruses",virues$proviruses),]
abu_pro=data.frame(proviruses=apply(abu_pro[,1:24],2,sum))
abu_pro$depth_value=as.factor(group[rownames(abu_pro),66])
abu_pro$landuses=group[rownames(abu_pro),9]
cor.test(abu_pro$proviruses,as.numeric(abu_pro$depth_value))
cor.test(abu_pro$proviruses[1:12],as.numeric(abu_pro$depth_value)[1:12])
cor.test(abu_pro$proviruses[13:24],as.numeric(abu_pro$depth_value)[13:24])
lm(formula = abu_pro$proviruses~
     as.numeric(abu_pro$depth_value))%>%summary()
lm(abu_pro$proviruses[1:12]~as.numeric(abu_pro$depth_value)[1:12])%>%summary()
lm(abu_pro$proviruses[13:24]~as.numeric(abu_pro$depth_value)[13:24])%>%summary()

library(ggpubr)
p7=ggplot(abu_pro,aes(depth_value,proviruses))+
  geom_boxplot(aes(fill=landuses))+
  geom_point(aes(fill=landuses),shape=21,
             position = position_jitter(width =0.2))+
  scale_fill_brewer(palette = "Set1")+
  theme_classic()+
  labs(x="Depth (cm)",
       y="Relative abundance \nof proviruses (%)")+
  stat_compare_means(paired = T,
                     comparisons = list(c("25","50"),
                                        c("25","100"),
                                        c("50","100")))+
  theme(legend.position = "none")
p7
abu_pro$depth_value=group[rownames(abu_pro),66]

p8=ggplot(abu_pro,aes(depth_value,proviruses,
                   group=landuses))+
  geom_point(aes(fill=landuses),shape=21,size=4,
             position = position_jitter(width =0.2))+
  scale_fill_brewer(palette = "Set1")+
  scale_color_brewer(palette = "Set1")+
  geom_smooth(aes(color=landuses,fill=landuses),
              method = "lm",alpha=0.3)+
  theme_classic()+
  labs(x="Depth (cm)",
       y="Relative abundance \nof proviruses (%)")+
  theme(legend.position = "top")
p8
library(gridExtra)
grid.arrange(p7,p8,p9,p4,p1,p2,p3,p5,p6,ncol=3)
##host-linked
host=read.table("../host/viruses_host_uniq.tsv")
host$classification=gtdbtk[as.character(host$V5),2]
host$lifestyle=virues[as.character(host$V1),25]
library(stringr)
library(dplyr)
library(reshape2)
library(ggpubr)
library(ggplot2)
host$phylum<- str_split(host$classification, ";", simplify = TRUE)[,3]
host$phylum<- str_split(host$phylum, "__", simplify = TRUE)[,2]
arch=host[grep("Archaea",host$classification),]

Methylomirabilia=host[grep("Methylomirabilia",host$classification),]
Methylomirabilia=Methylomirabilia[!duplicated(Methylomirabilia$V1),]
Methylomirabilia[,9:33]=virues[as.character(Methylomirabilia$V1),1:24]
Methylomirabilia=aggregate.data.frame(Methylomirabilia[,9:32],list(Methylomirabilia$lifestyle),sum)
Methylomirabilia=melt(Methylomirabilia)
Methylomirabilia$depth_value=as.factor(group[as.character(Methylomirabilia$variable),66])
Methylomirabilia$landuses=group[as.character(Methylomirabilia$variable),9]
Methylomirabilia$Group.1=factor(Methylomirabilia$Group.1,levels = c("temperate","virulent","unclassified"))
ggplot(Methylomirabilia[grep("virulent",Methylomirabilia$Group.1),],aes(depth_value,value,fill=Group.1))+
  geom_boxplot(outliers = F,whisker.linewidth = 0.7,staplewidth = 0.5)+
  theme_classic()+
  #facet_grid(~landuses)+
  scale_fill_manual(values = c("#F9E2E6","#872C45","#6A4C9C"))+
  labs(y='Relative abundance (%)')+
  stat_compare_means(comparisons = list(c("25","75")))
ggplot(Methylomirabilia,aes(depth_value,value,fill=Group.1))+
  geom_boxplot(outliers = F,whisker.linewidth = 0.7,staplewidth = 0.5)+
  theme_classic()+
  facet_grid(~landuses)+
  scale_fill_manual(values = c("#F9E2E6","#872C45","#6A4C9C"))+
  labs(y='Relative abundance (%)')
ggsave("Methylomirabila_lifestyle.pdf",device = "pdf",width = 5,height = 2)

Methylomirabilia=host[grep("Nitro",host$classification),]
Methylomirabilia=Methylomirabilia[!duplicated(Methylomirabilia$V1),]
Methylomirabilia[,9:33]=virues[as.character(Methylomirabilia$V1),]
Methylomirabilia=aggregate.data.frame(Methylomirabilia[,9:32],list(Methylomirabilia$lifestyle),sum)
Methylomirabilia=melt(Methylomirabilia)
Methylomirabilia$depth_value=as.factor(group[as.character(Methylomirabilia$variable),66])
Methylomirabilia$landuses=group[as.character(Methylomirabilia$variable),9]
Methylomirabilia$Group.1=factor(Methylomirabilia$Group.1,levels = c("temperate","virulent","unclassified"))
ggplot(Methylomirabilia,aes(depth_value,value,fill=Group.1))+
  geom_boxplot(outliers = F,whisker.linewidth = 0.7,staplewidth = 0.5)+
  theme_classic()+
  #facet_grid(~landuses)+
  scale_fill_manual(values = c("#F9E2E6","#872C45","#6A4C9C"))+
  labs(y='Relative abundance (%)')
ggsave("Nitro_lifestyle.pdf",device = "pdf",width = 5,height = 2)


Methylomirabilia=host[grep("Halo",host$classification),]
Methylomirabilia=Methylomirabilia[!duplicated(Methylomirabilia$V1),]
Methylomirabilia[,9:33]=virues[as.character(Methylomirabilia$V1),]
Methylomirabilia=aggregate.data.frame(Methylomirabilia[,9:32],list(Methylomirabilia$lifestyle),sum)
Methylomirabilia=melt(Methylomirabilia)
Methylomirabilia$depth_value=as.factor(group[as.character(Methylomirabilia$variable),66])
Methylomirabilia$landuses=group[as.character(Methylomirabilia$variable),9]
Methylomirabilia$Group.1=factor(Methylomirabilia$Group.1,levels = c("temperate","virulent","unclassified"))
ggplot(Methylomirabilia,aes(depth_value,value,fill=Group.1))+
  geom_boxplot(outliers = F,whisker.linewidth = 0.3,staplewidth = 0.5)+
  theme_classic()+
  #facet_grid(~landuses)+
  scale_fill_manual(values = c("#F9E2E6","#872C45","#6A4C9C"))+
  labs(y='Relative abundance (%)')

ggsave("ANME_lifestyle.pdf",device = "pdf",width = 4,height = 2)
