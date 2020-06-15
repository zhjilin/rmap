library(ggplot2)
library(extrafont)
args=commandArgs(trailingOnly = TRUE)


outbase="STEMBits"
tstring=gsub(" .*$","",Sys.time(),perl=TRUE)

pwmtab= args[1]
ditab=args[2]

pwm=read.table(pwmtab,sep = "\t")
di=read.table(ditab,sep = "\t")
diname=di$V1
ndiname=gsub(".di.pfm","_short.pfm",diname)
#ndipre=cbind(di,ndiname)
di$V1=ndiname

overalloutname=paste(outbase,"_overall_",tstring,".pdf",sep="")
perbaseoutname=paste(outbase,"_perbase_",tstring,".pdf",sep="")
distem=subset(di,V3=="STEM")
pwmstem=subset(pwm,V3=="STEM")


#dimean=aggregate(distem$V2,by=list(distem$V1),FUN=mean)
#pwmmean=aggregate(pwmstem$V2,by=list(pwmstem$V1),FUN=mean)
#colnames(dimean)=c("motif","Di")
#colnames(pwmmean)=c("motif","Mono")
colnames(distem)=c("motif","Di","Tag")
colnames(pwmstem)=c("motif","Mono","Tag")
allstem=merge(distem,pwmstem,by="motif")
q=ggplot(allstem,aes(Mono,Di/2)) +geom_point(col="royalblue",alpha=0.5,size=2)+scale_fill_brewer(palette="Pastel1")+scale_y_continuous("Di")+scale_x_continuous("Mono")+ggtitle("Perbase IC of stem in Mono & Di matrix")+
  theme(panel.background=element_blank(),axis.line.x=element_line(size=1,color="grey",linetype="solid"),axis.line.y=element_line(size=1,color="grey",linetype="solid"),legend.position=c(0.5,0.95),legend.direction="horizontal",legend.title=element_blank(),legend.background=element_blank(),plot.title=element_text(size=14,hjust = 0.5),axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=12),text = element_text(family = "Arial Narrow"))

ggsave(filename = perbaseoutname,width = 10,height =10,units = "cm",plot = q,useDingbats=FALSE)

# the overall
dnn=distem$motif
mnn=pwmstem$motif
distem$motif=gsub("pfm.*","pfm",dnn)
pwmstem$motif=gsub("pfm.*","pfm",mnn)
disum=aggregate(distem$Di,by=list(distem$motif),FUN=sum)
pwmsum=aggregate(pwmstem$Mono,by=list(pwmstem$motif),FUN=sum)
colnames(disum)=c("motif","Di")
colnames(pwmsum)=c("motif","Mono")
allsum=merge(disum,pwmsum,by="motif")
p=ggplot(allsum,aes(Mono,Di/2)) +geom_point(col="blue",alpha=0.5,size=2)+scale_x_continuous(limits=c(0,8))+scale_y_continuous("Di",limits=c(0,8))+ggtitle("Overall IC of stem in Mono & Di matrix")+
  theme(panel.background=element_blank(),axis.line.x=element_line(size=1,color="grey",linetype="solid"),axis.line.y=element_line(size=1,color="grey",linetype="solid"),legend.position=c(0.5,0.95),legend.direction="horizontal",legend.title=element_blank(),legend.background=element_blank(),plot.title=element_text(size=14,hjust = 0.5),axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=12),text = element_text(family = "Arial Narrow"))
ggsave(filename = overalloutname,width = 10,height =10,units = "cm",plot = p,useDingbats=FALSE)

overallIncreaseBits=mean(allsum$Di -allsum$Mono)
write.table(overallIncreaseBits,file=paste("OverallIncreasedBits",tstring,".txt",sep = ""),col.names = F,row.names = F,quote = F)


