#for all stem-loop pfms, give out the overall and perbase Bits
library(ggplot2)
library(extrafont)
args <- commandArgs(trailingOnly = TRUE)

#font_import()
#loadfonts()
slmtable=args[1]
perbasebin=as.numeric(args[2])

outbase="SLMBits"
tstring=gsub(" .*$","",Sys.time(),perl=TRUE)
overalloutname=paste(outbase,"_overall_",tstring,".pdf",sep="")
perbaseoutname=paste(outbase,"_perbase_",perbasebin,tstring,".pdf",sep="")


slminfo=read.table(slmtable,sep="\t")

infostem=subset(slminfo,V3=="STEM")
infoloop=subset(slminfo,V3=="LOOP")

#minfostems=aggregate(infostem$V2,by=list(infostem$V1),FUN=mean)
#minfoloops=aggregate(infoloop$V2,by=list(infoloop$V1),FUN=mean)
#nminfostems=data.frame(infostems,"STEM")
#nminfoloops=data.frame(infoloops,"LOOP")
colnames(slminfo)=c("motif","Bits","Position")

#colnames(nminfostems)=c("motif","Bits","Position")
#colnames(nminfoloops)=c("motif","Bits","Position")

#allminfo=rbind(nminfostems,nminfoloops)
#pdf(perbaseoutname,width=6,height=6)
q=ggplot(slminfo,aes(Bits,fill=Position)) +geom_histogram(alpha=0.5, linetype=0,binwidth=perbasebin,stat="bin",position="identity")+scale_x_continuous(paste("Bits(bin=",perbasebin,")",sep=""))+scale_fill_brewer(palette="Pastel1")+scale_y_continuous("Counts")+ggtitle("Perbase information content")+
theme(text = element_text(family = "Arial Narrow"), panel.background=element_blank(),axis.line.x=element_line(size=1,color="grey",linetype="solid"),axis.line.y=element_line(size=1,color="grey",linetype="solid"),legend.position=c(0.5,0.95),legend.direction="horizontal",legend.title=element_blank(),legend.background=element_blank(),title=element_text(size=14),axis.text=element_text(size=12),axis.title=element_text(size=14),plot.title=element_text(hjust=0.5,size = 14),legend.text=element_text(size=12))
#dev.off()
ggsave(filename = perbaseoutname,width = 10,height =10,units = "cm",plot = q)


infoloop$V1=gsub(pattern = ".pfm.*",replacement = ".pfm",infoloop$V1)
infostem$V1=gsub(pattern = ".pfm.*",replacement = ".pfm",infostem$V1)
infostems=aggregate(infostem$V2,by=list(infostem$V1),FUN=sum)
infoloops=aggregate(infoloop$V2,by=list(infoloop$V1),FUN=sum)

ninfostems=data.frame(infostems,"STEM")
ninfoloops=data.frame(infoloops,"LOOP")
colnames(ninfoloops)=c("motif","Bits","Position")
colnames(ninfostems)=c("motif","Bits","Position")

snslminfo=rbind(ninfoloops,ninfostems)

overallbin=1
p=ggplot(snslminfo,aes(Bits,fill=Position)) +geom_histogram(binwidth=overallbin,stat="bin",position="identity",alpha=0.5,linetype=0)+scale_x_continuous(paste("Bits(bin=",overallbin,")",sep=""))+scale_fill_brewer(palette="Pastel1")+scale_y_continuous("Counts")+ggtitle("Overall information content")+
theme(panel.background=element_blank(),axis.line.x=element_line(size=1,color="grey",linetype="solid"),axis.line.y=element_line(size=1,color="grey",linetype="solid"),legend.position=c(0.5,0.95),legend.direction="horizontal",legend.title=element_blank(),legend.background=element_blank(),plot.title=element_text(size=14,hjust = 0.5),axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=12),text = element_text(family = "Arial Narrow"))

ggsave(filename = overalloutname,width = 10,height =10,units = "cm",plot = p)

#pdf(overalloutname,width=6,height=6)

#dev.off()
