#upadte 20200419, shift the coordinate of motif matches on minus strand to the first base of the motif to generate the distance to feature histogram 
library(ggplot2)
library(plyr)
library(extrafont)
library(pagenum)
library(gridExtra)

#change the span size from 0.1 to 0.05, which means the only 5 % of the points used for local smoothing 
per_motif=function(x,tag, barlen){
  
  labeltext=" "
  if(any(grepl("DONOR",tag,ignore.case = T))){
    labeltext="exon    intron"
    deduction=1000
  }else if(any(grepl("ACCEPTOR",tag,ignore.case = T))){
    labeltext="intron    exon"
    deduction=1002
  }else{
    deduction=1000
  }
  
  ggplot(x,aes(position-deduction,count,fill=strand))+geom_bar(stat="identity",position="identity",alpha=0.4)+scale_color_brewer(aesthetics = "fill",palette = "Set1",direction = -1)+
    scale_x_continuous("Distance to feature(base)",limits=c(-50,50))+
	  scale_y_continuous("Count")+
    geom_segment(x=as.numeric(-barlen),xend=0,y=0.9*max(x$count),yend=0.9*max(x$count),col="#abdda4")+
    theme(plot.title = element_text(hjust = 0.5,size=16,family = "Arial",face = "bold",vjust = 0),axis.title = element_text(size=16,family="Arial"),legend.text = element_text(size=14,family = "Arial"),axis.text = element_text(size=14,family="Arial"),axis.line = element_line(color="black"),panel.background = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_line(color="grey"),legend.background = element_blank(),legend.key = element_blank(),legend.title = element_blank(),legend.direction = "horizontal",legend.position =c(0.5,0.98),legend.box = "horizontal")+
    annotate("text",label=labeltext,x=0,y=0.8*max(x$count),family="Arial",size=14)
}

args=commandArgs(trailingOnly = TRUE)
fsense=args[1]
fanti=args[2]
tag=args[3]
name_map=args[4]
winsize=as.numeric(args[5])
order=args[6]
motif_length=args[7]
#motif_length="/mnt/raid5/project/rnaselex/figures/GR-revision2020/distanceHistogram/motif_length.txt"

tstring=Sys.Date()

# fanti="DONOR-antisense-180120.freq"
# fsense="DONOR-sense-180120.freq"
# fsense="START-sense-180108.freq"
# # fanti="START-antisense-180108.freq"
# totalnum="total_number180120.table"
# name_map="../../../../name_file.map.17Apr28.tsv"
# tag="DONOR"
# tag="TSS"


tabsense=read.table(fsense,sep="\t",col.names = c("motif","position","count","feature"),quote = "")
tabanti=read.table(fanti,sep="\t",col.names = c("motif","position","count","feature"),quote = "")
#tnum=read.table(totalnum,sep="\t",col.names = c("motif","number"),quote = "")
proteinname=read.table(name_map,sep="\t",quote = "",col.names = c("protein","barcode","batch","seed","multi","cycle","anno","motif"))
motiflength=read.table(motif_length,sep="\t",col.names = c("motif","length"))

orderfile=read.table(order,sep="\t",quote = "")

motif_list=as.character(orderfile$V1)


tabsense$strand="sense"
tabanti$strand="antisense"


sectype=c("antisense","lincRNA","protein_coding","processed_transcript","nonsense_mediated_decay","tRNA")


pdf(paste(tag,"-distance-",tstring, ".pdf",sep=""),height=12,width=10)
pnum=1
i=1
plot=list()
actualwin=winsize+1
for ( n in 1:length(motif_list)){
  
  outname=gsub("_short.pfm","", motif_list[n],perl = TRUE)
  namepro=paste(proteinname[grep(outname,proteinname$motif),1],proteinname[grep(outname,proteinname$motif),4],sep="-")
  mlen=as.numeric(motiflength[grep(pattern =motif_list[n] , motiflength[,1]),2])
  hsense=subset(tabsense[tabsense$motif == motif_list[n],], ((position- 1000) >= -actualwin) & ((position- 1000) <= actualwin))
  hanti=subset(tabanti[tabanti$motif == motif_list[n],], ((position- 1000)-(mlen-1) >= -actualwin) & ((position- 1000)-(mlen-1) <= actualwin))
  
#   outname=gsub("_short.pfm","", motif_list[n],perl = TRUE)
#   namepro=paste(proteinname[grep(outname,proteinname$motif),1],proteinname[grep(outname,proteinname$motif),4],sep="-")
  
#   hsense=subset(tabsense[tabsense$motif == motif_list[n],], (position>=900 & position <=1102))
#   hanti=subset(tabanti[tabanti$motif == motif_list[n],], (position>=900 & position <=1102))
  

  alltab=rbind(hsense,hanti)

  
  pn=per_motif(alltab,tag,mlen)+ggtitle(namepro)
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
