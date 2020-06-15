library(ggplot2)
library(scales)
library(extrafont)
args <- commandArgs(trailingOnly = TRUE)
filename=args[1]
outnamebase="PriSecPie"
tstring=gsub(" .*$","",Sys.time(),perl=TRUE)
outname=paste(outnamebase,tstring,"pdf",sep=".")
pie=read.table(filename,sep="\t")

#pdf(outname)
p=ggplot(pie,aes(x="",y=V2,fill=V1))+geom_bar(width=1,stat="identity")+coord_polar("y",start=pi/3.14)+scale_x_discrete("")+scale_y_continuous("")+scale_fill_manual(breaks=c("Primary","Secondary"), values=c("#FBB4AE","#B3CDE3"),labels=c("Primary motif only","With secondary motif"))+
geom_text(aes(y = c(0, cumsum(pie$V2)[-length(pie$V2)]),label = percent(pie$V2/sum(pie$V2))), size=8)+
theme(text = element_text(family = "Arial Narrow"), plot.title=element_text(hjust = 0.5,size=14),axis.text=element_blank(),axis.title=element_text(size=14),axis.line=element_blank(),panel.background=element_blank(),legend.text=element_text(size=12),legend.direction = "horizontal",legend.position = "top",legend.background=element_blank(),legend.title=element_blank())+ggtitle("RBPs with secondary motifs")
ggsave(p,filename = outname,height = 10,width = 10,units = "cm")
#dev.off()

