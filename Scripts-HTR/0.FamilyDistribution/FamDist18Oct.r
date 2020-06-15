
library(ggplot2)
library(reshape2)
library(extrafont)
args=commandArgs(trailingOnly = T)

args[1]="dataCSV20181010.famknown2018Oct10"
args[2]="dataCSV20181010.famunid2018Oct10"
args[3]="dataCSV20181010.famhtr2018Oct10"
known=read.table(args[1],sep="\t")
unidentified=read.table(args[2],sep="\t")
htr=read.table(args[3],sep="\t")

outbase="RBPfamDist"
tstring=gsub(" .*$","",Sys.time(),perl=TRUE)
outname=paste(outbase,tstring,".pdf",sep="")

kn_htr=merge(known,htr,by="V1")
kn_htr_unid=merge(kn_htr,unidentified,by="V1")
colnames(kn_htr_unid)=c("V1","V2","V3","V4")
kn_htr_unid$V5=kn_htr_unid$V2 - kn_htr_unid$V3
colnames(kn_htr_unid)=c("name","known","HTR","unknown","previous")
ntab=data.frame(kn_htr_unid[,1],kn_htr_unid[,4],kn_htr_unid[,5],kn_htr_unid[,3])
colnames(ntab)=c("name","unknown","previous","HTR")
ntab=ntab[order(ntab$HTR+ntab$unknown+ntab$previous,decreasing=T),]
nntab=melt(ntab)

p=ggplot(nntab,aes(name,value,fill=variable))+geom_col(width = 0.618)+
  scale_x_discrete("RBP Family",limits=ntab$name)+ scale_fill_brewer(limits=c("HTR","previous","unknown"),labels=c("HTR","Previous","Unknown"),palette="Pastel1")+
  scale_y_continuous("RBP Counts")+ggtitle("RBPs protein class")+
  theme(panel.background=element_blank(),axis.line.x=element_line(size=1,color="grey",linetype="solid"),axis.line.y=element_line(size=1,color="grey",linetype="solid"),plot.title=element_text(size=rel(1.5),hjust = 0.5),axis.text.y=element_text(size=rel(1.2)),axis.text.x=element_text(angle=60,size=rel(1.2),vjust=0.5),axis.title=element_text(size=rel(1.5)),legend.text=element_text(size=rel(1.2)),legend.position=c(0.85,0.9),legend.title=element_blank(),legend.background=element_blank(),text = element_text(family = "Arial Narrow")) 

ggsave(p,filename = outname,height = 10,width = 10,units = "cm")