quality=read.table("quality_summary_filter10k.tsv",sep='\t',header = T,row.names = 1)
quality[grep("Yes",quality$provirus),7]="proviruses"
quality$checkv_quality=factor(quality$checkv_quality,levels = rev(c("Complete","High-quality",
                                                                    "Medium-quality","proviruses",
                                                                    "Low-quality","Not-determined")))
library(ggplot2)
ggplot(quality,aes(checkv_quality,contig_length))+
  stat_boxplot(geom = "errorbar",
               width=0.3)+
  geom_boxplot(fill="#008B00",outlier.colour = "gray",
               outlier.shape = 21,outlier.fill = "white",outlier.size = 1)+
  theme_light()+ scale_y_log10()+ coord_flip()+labs(y="length (bp)",x="")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 60,size = 10,hjust = 1),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 10))+
  theme(panel.grid = element_line(color = 'white', linetype = 2, size = 0), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'))

quality=read.table("quality_summary.tsv",sep='\t',header = T,row.names = 1)
quality[grep("Yes",quality$provirus),7]="proviruses"
quality$checkv_quality=factor(quality$checkv_quality,levels = rev(c("Complete","High-quality",
                                                                    "Medium-quality","proviruses",
                                                                    "Low-quality","Not-determined")))
ggplot(quality,aes(checkv_quality,contig_length))+
  stat_boxplot(geom = "errorbar",
               width=0.3)+
  geom_boxplot(fill="#008B00",outlier.colour = "gray",
               outlier.shape = 21,outlier.fill = "white",outlier.size = 1)+
  theme_light()+ scale_y_log10()+ coord_flip()+labs(y="length (bp)",x="")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 60,size = 10,hjust = 1),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 10))+
  theme(panel.grid = element_line(color = 'white', linetype = 2, size = 0), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'))
ggsave("quality_summary.pdf",device = "pdf",width = 3,height = 3)

quality_table=as.data.frame(table(quality$checkv_quality))
quality_table$Var1=factor(quality_table$Var1,levels = rev(c("Complete","High-quality",
                                                                "Medium-quality","proviruses",
                                                                "Low-quality","Not-determined")))
ggplot(quality_table,aes(Var1,Freq))+
  geom_bar(fill="#008B00",stat = "identity")+coord_flip()+
  theme_light()+ scale_y_log10()+ coord_flip()+labs(y="Number",x="")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 60,size = 10,hjust = 1),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 10))+
  theme(panel.grid = element_line(color = 'white', linetype = 2, size = 0), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'))
ggsave("quality_summary_nember.pdf",device = "pdf",width = 2.5,height = 3)

quality_filter1=read.table("quality_summary_filter10k.tsv",sep='\t',header = 1)
quality_filter2=read.table("quality_summary_proviruses.tsv",sep = '\t',header = F)
quality_filter2=quality_filter2[,-c(11,14)]
colnames(quality_filter2)=colnames(quality_filter1)
quality_filter=rbind(quality_filter1,quality_filter2)
phatyp=read.table("phatyp_prediction.tsv",sep = '\t',header = 1)
quality_filter=quality_filter[!duplicated(quality_filter$contig_id),]
rownames(quality_filter)=quality_filter$contig_id
rownames(phatyp)=phatyp$Accession
quality_filter$lifestyle=phatyp[rownames(quality_filter),3]
write.table(quality_filter,"ASVDD_quality_filter.tsv",sep='\t',quote = F,row.names = F)

