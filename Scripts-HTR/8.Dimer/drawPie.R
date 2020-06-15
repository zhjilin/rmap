
#A script to generate pie chart to show the proportion of RBPs that can form tandem  motifs

library(ggplot2)
library(scales)
outnamebase="Pie-Dimer"
pie=read.table("pie.table",sep="\t")
tstring=gsub(" .*$","",Sys.time(),perl=TRUE)
outname=paste(outnamebase,tstring,"pdf",sep=".")
pdf(outname)
ggplot(pie,aes(x="",y=V2,fill=V1))+geom_bar(width=1,stat="identity")+coord_polar("y",start=pi/3.14)+scale_x_discrete("")+scale_y_continuous("")+scale_fill_brewer(palette="Pastel1")+
geom_text(aes(y = c(0, cumsum(pie$V2)[-length(pie$V2)]),label = percent(pie$V2/sum(pie$V2))), size=8)+
theme(panel.background=element_blank(),title=element_text(size=rel(1.8)),axis.text=element_blank(),axis.title=element_text(size=rel(1.5)),legend.text=element_text(size=rel(1.2)),legend.position=c(0.9,.9),legend.background=element_blank(),legend.title=element_blank())+ggtitle("RBPs with dimeric binding specificities")
dev.off()
