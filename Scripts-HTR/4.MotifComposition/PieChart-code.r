#script to generate the pie chart for Fig. 4G 
#update to barplot 2017-05-02

library(ggplot2)
library(extrafont)
library(RColorBrewer)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

rhttable=args[1]
#pubtable=args[2]

rht=read.table(rhttable,sep="\t")

RdBu=brewer.pal(n=11,"RdBu")
RdBu16=colorRampPalette(RdBu)(16)
RdBu4=brewer.pal(n=4,"RdBu")


tstring=gsub(" .*$","",Sys.time(),perl=TRUE)
mononame=paste("mono-",tstring,".pdf",sep="")
diname=paste("di-",tstring,".pdf",sep="")
rhttabn=data.frame(rht,"RHT")
names(rhttabn)=c("base","value","Type","Experiment")
monorht=rhttabn[grepl(pattern="mono",rhttabn$Type,perl=T),]
dirht=rhttabn[grepl(pattern="di",rhttabn$Type,perl=T),]
monorhtn=monorht
dirhtn=dirht



dirhtn$base=gsub(pattern="T","U",dirhtn$base,perl=T)
monorhtn$base=gsub(pattern="T","U",monorhtn$base,perl=T)
dp=ggplot(dirhtn,aes(x="",y=value,fill=base))+geom_bar(width = 1,stat = "identity")+coord_polar(theta="y",start=pi/3.14)+scale_fill_manual(values = rev(RdBu16))+scale_x_discrete("")+scale_y_continuous("")+
  geom_text(aes(x=1.28,label = percent(dirhtn$value/sum(dirhtn$value))), size=4,position = position_stack(vjust = 0.5))+
  theme(text = element_text(family = "Arial Narrow"), plot.title=element_text(hjust = 0.5,size=14),axis.text=element_blank(),axis.title=element_text(size=14),axis.line=element_blank(),panel.background=element_blank(),legend.text=element_text(size=12),legend.background=element_blank(),legend.title=element_blank())+ggtitle("di-nucleotides")

mp=ggplot(monorhtn,aes(x="",y=value,fill=base))+geom_bar(width = 1,stat = "identity")+coord_polar("y",start=pi/3.14)+scale_fill_manual(values = rev(RdBu4))+scale_x_discrete("")+scale_y_continuous("Percent")+
  geom_text(aes(x=1.28,label = percent(monorhtn$value/sum(monorhtn$value))), size=4,position = position_stack(vjust = 0.5))+
  theme(text = element_text(family = "Arial Narrow"), plot.title=element_text(hjust = 0.5,size=14),axis.text=element_blank(),axis.title=element_text(size=14),axis.line=element_blank(),panel.background=element_blank(),legend.text=element_text(size=12),legend.background=element_blank(),legend.title=element_blank())+ggtitle("mono-nucleotide")


ggsave(plot=mp,height=12,width=12, filename=mononame, useDingbats=FALSE,units="cm")
ggsave(plot=dp,height=12,width=12, filename=diname, useDingbats=FALSE,units="cm")
