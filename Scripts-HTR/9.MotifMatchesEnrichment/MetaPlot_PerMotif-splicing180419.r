#script to generate the position specific meta plot
#update 20180120 to support the permotif mode
#update 20180303 to support the selected region for smoothing
#update 20180419 to get the desired order of the RBP
library(ggplot2)
library(plyr)
library(extrafont)
library(pagenum)
library(gridExtra)

#change the span size from 0.1 to 0.05, which means the only 5 % of the points used for local smoothing 
per_motif=function(x,sv,tag){
#   x=alltab
#   sv=0.01
#   tag="spliceDONOR"
  labeltext=" "
  if(any(grepl("DONOR",tag))){
    labeltext="exon    intron"
    deduction=1000
  }else if(any(grepl("ACCEPTOR",tag))){
    labeltext="intron    exon"
    deduction=1002
  }else{
    deduction=1000
  }
  ggplot(x,aes(position-deduction,count))+geom_point(aes(colour=factor(strand),fill=factor(strand)),alpha=0.7)+geom_smooth(aes(linetype=strand),method="loess",span=sv,se=F,lwd=1.2,col="#252525")+
    scale_color_manual(breaks=c("antisense","sense"),values = c("#deebf7","#4292c6"))+
    scale_fill_manual(breaks=c("antisense","sense"),values = c("#deebf7","#4292c6"))+
    scale_linetype_manual(breaks=c("antisense","sense"),values = c("dotted","solid"),labels=c("antisense","sense"))+
    scale_y_continuous(limits=c(0,max(x$count)*1.05),"normalized count")+
    scale_x_continuous("position")+
    theme(plot.title = element_text(hjust = 0.5,size=16,family = "Arial"),axis.title = element_text(size=16,family="Arial"),legend.text = element_text(size=14,family = "Arial"),axis.text = element_text(size=14,family="Arial"),axis.line = element_line(color="black"),panel.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_line(color="grey"),legend.background = element_blank(),legend.key = element_blank(),legend.title = element_blank(),legend.direction = "horizontal",legend.position = "top",legend.box = "horizontal")+
    annotate("text",label=labeltext,x=0,y=0.8*max(x$count),family="Arial",size=14)
  
#old plot style

#   ggplot(x,aes(position-deduction,count))+geom_point(aes(shape=factor(strand)),alpha=0.75)+geom_smooth(aes(linetype=strand),method="loess",span=sv,se=F,lwd=1.2)+
#     scale_linetype_manual(breaks=c("antisense","sense"),values = c("dotted","solid"),labels=c("anti-sense","sense"))+
#     scale_y_continuous(limits=c(0,max(x$count)*1.05),"normalized count")+
#     scale_x_continuous("position")+
#     theme(plot.title = element_text(hjust = 0.5,size=16,family = "Arial"),axis.title = element_text(size=16,family="Arial"),legend.text = element_text(size=14,family = "Arial"),axis.text = element_text(size=14,family="Arial"),axis.line = element_line(color="black"),panel.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_line(color="grey"),legend.background = element_blank(),legend.key = element_blank(),legend.title = element_blank(),legend.direction = "horizontal",legend.position = "top",legend.box = "horizontal")+
#     annotate("text",label=labeltext,x=0,y=0.8*max(x$count),family="Arial",size=14)
}

args=commandArgs(trailingOnly = TRUE)
fsense=args[1]
fanti=args[2]
totalnum=args[3]
tag=args[4]
name_map=args[5]
spanval=as.numeric(args[6])
winsize=as.numeric(args[7])
order=args[8]

tstring=Sys.Date()

# fanti="DONOR-antisense-180120.freq"
# fsense="DONOR-sense-180120.freq"
# fsense="START-sense-180120.freq"
# fanti="START-antisense-180120.freq"
# totalnum="total_number180120.table"
# name_map="name_file.map.17Apr28.tsv"
# tag="STARTcodon"
# spanval=0.1
# winsize=25
# name_map="../../../../name_file.map.17Apr28.tsv"
# tag="DONOR"
# tag="TSS"


tabsense=read.table(fsense,sep="\t",col.names = c("motif","position","count","feature"),quote = "")
tabanti=read.table(fanti,sep="\t",col.names = c("motif","position","count","feature"),quote = "")
tnum=read.table(totalnum,sep="\t",col.names = c("motif","number"),quote = "")
proteinname=read.table(name_map,sep="\t",quote = "",col.names = c("protein","barcode","batch","seed","multi","cycle","anno","motif"))

orderfile=read.table(order,sep="\t",quote = "")

motif_list=as.character(orderfile$V1)

tabsense$strand="sense"
tabanti$strand="antisense"
rtabsense=merge(tabsense,tnum,by="motif")
rtabanti=merge(tabanti,tnum,by="motif")

# outname=paste(tag,"-",tstring,".pdf",sep="")

sectype=c("antisense","lincRNA","protein_coding","processed_transcript","nonsense_mediated_decay","tRNA")


nrtabsense=aggregate((rtabsense[,3]/rtabsense[,6]),list(rtabsense[,1],rtabsense[,2],rtabsense[,5]),sum)
nrtabanti=aggregate((rtabanti[,3]/rtabanti[,6]),list(rtabanti[,1],rtabanti[,2],rtabanti[,5]),sum)
colnames(nrtabsense)=c("motif","position","strand","count")
colnames(nrtabanti)=c("motif","position","strand","count")

# nrtabsense=aggregate((rtabsense[,3]/rtabsense[,6]),list(rtabsense[,1],rtabsense[,2],rtabsense[,4],rtabsense[,5]),sum)
# nrtabanti=aggregate((rtabanti[,3]/rtabanti[,6]),list(rtabanti[,1],rtabanti[,2],rtabanti[,4],rtabanti[,5]),sum)
# 
# colnames(nrtabsense)=c("motif","position","feature","strand","count")
# colnames(nrtabanti)=c("motif","position","feature","strand","count")


pdf(paste(tag,"-PerMotif-range",winsize,"-",tstring, "span-",spanval,".pdf",sep=""),height=12,width=10.5)
pnum=1
i=1
plot=list()
actualwin=winsize+1
for ( n in 1:length(motif_list)){
  outname=gsub("_short.pfm","", motif_list[n],perl = TRUE)
  namepro=paste(proteinname[grep(outname,proteinname$motif),1],proteinname[grep(outname,proteinname$motif),4],sep="-")
  hsense=subset(nrtabsense[nrtabsense$motif == motif_list[n],], ((position- 1000) >= -actualwin) & ((position- 1000) <= actualwin))
  hanti=subset(nrtabanti[nrtabanti$motif == motif_list[n],], ((position- 1000) >= -actualwin) & ((position- 1000) <= actualwin))
  
    
  # hsense=subset(nrtabsense[nrtabsense$motif == motif_list[n],], (position>=500 & position <=1502))
  # hanti=subset(nrtabanti[nrtabanti$motif == motif_list[n],], (position>=500 & position <=1502))
  

  alltab=rbind(hsense,hanti)
#  if(tag=="TSS" || tag=="TTS"){
#        alltab_sub=subset(alltab,feature =="protein_coding" )
#        alltab=alltab_sub
#  }
  
  pn=per_motif(alltab,spanval,tag)+ggtitle(namepro)
  plot[[i]]=pn

  if (i %% 6 == 0) { ## print 6 plots on a page
    print (do.call(grid.arrange,  plot))
    plot = list() # reset plot 
    i = 0 # reset index
    pagenum(num="", text=pnum,x=0.5,y=0.01,cex=1.5)
    pnum=pnum+1
  }
  i = i + 1
}

if ( length(plot) != 0 ) { 
  print (do.call(grid.arrange,  plot))
  pagenum(num="", text=pnum,x=0.5,y=0.01,cex=1.5)
}
dev.off()
