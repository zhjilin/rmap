#script to generate the histogram of length for loop and stem of the motifs in Fig. 2E 

args <- commandArgs(trailingOnly = TRUE)
library(ggplot2)
library(extrafont)
library(reshape2)
filehtselex=args[1]
#filehtselex="motif_SL_length.tsv"
outprefix="SL-length_dist"
tstring=gsub(" .*$","",Sys.time(),perl=TRUE)
outname=paste(outprefix,tstring,".pdf",sep="")
rhtselex=read.table(filehtselex,sep="\t")
colnames(rhtselex)=c("motif","Stem","Loop")
nrhtselex=melt(rhtselex)
colnames(nrhtselex)=c("motif","tag","Length")


  ggplot(nrhtselex,aes(Length,fill=tag))+geom_histogram(alpha=0.5,stat="bin",binwidth=1,position="identity",linetype=0)+scale_fill_brewer(palette="Pastel1",direction = -1)+scale_y_continuous("Counts")+ggtitle("Length distribution of stem and loop")+
  theme(legend.position=c(0.8,0.9),legend.direction = "vertical",legend.title=element_blank(),panel.background=element_blank(),axis.line.y=element_line(size=1,colour="grey",linetype="solid"),axis.line.x=element_line(size=1,colour="grey",linetype="solid"),axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=12),legend.background=element_blank(),title=element_text(size=14),plot.title=element_text(hjust=0.5),text=element_text(family = "Arial Narrow"))
ggsave(filename = outname,plot = p,units = "cm",width = 10,height = 10)
#dev.off()
