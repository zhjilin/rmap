library(png)
# library(extrafont)

args <- commandArgs(trailingOnly = TRUE)

# setwd("Downloads/Project/RNASELEX/Figures/HTR-RNAcompete")
tstring=gsub(" .*$","",Sys.time(),perl=TRUE)

outname=paste("htr-rnacompete.alnscore",tstring,".pdf",sep="")
#filename="Downloads/Project/RNASELEX/Figures/HTR-RNACom/ALL_newlist.dist.sort"

filename=args[1]
proname=args[2]
secname=args[3]
PSlist=args[4]

# # filename="ALL_newlist.dist.sort"
# PSlist="primary_secondary.list"
# # secname="ALL_listsec.dist"
# proname="protein.name.txt"
# filename="primary.score.tab.sort"
# secname="secondary.score.tab"

nsstat=read.table(filename,sep="\t",col.names = c("HTR","RNACompete","Score"))
pname=read.table(proname,sep = "\t")
nsecmotif=read.table(secname,sep="\t",col.names = c("HTR","RNACompete","Score"))
priseclist=read.table(PSlist,sep="\t")
colnames(pname)=c("HTR","protein")
colnames(priseclist)=c("primary","secondary")


# colnames(nsecmotif)=c("HTR","RNACompete","max","covariance")
# colnames(nsstat)=c("HTR","RNACompete","max","covariance")
nsstat$HTR=sub("pwm","pfm",nsstat$HTR)
nsecmotif$HTR=sub("pwm","pfm",nsecmotif$HTR)
pname$HTR=sub("pwm","pfm",pname$HTR)
snsstat=merge(nsstat,pname,by="HTR")
snsstat=snsstat[order(snsstat$Score,decreasing = T),]
Alist=paste("../dominating17Apr28/network_motifs/",snsstat$HTR,".png",sep="")
Blist=paste("converted_pfm/",snsstat$RNACompete,".png",sep="")
#Slist=paste("../dominating17Apr28/network_motifs/",priseclist$secondary,".png,")
pinix=2
piniy=28
sinix=14
siniy=28

pdf(outname,width=8,height=12)
plot(NA,xlim=c(0,30),ylim=c(0,28), xaxt="n",yaxt="n",bty="n")
text(pinix,piniy+0.5,"HTR-SELEX",cex=0.5)
text(sinix,siniy+0.5,"RNACompete",cex=0.5)

for(i in 1:length(Alist)){

  imgA=readPNG(Alist[i])
  imgB=readPNG(Blist[i])
  imgA_width=ncol(imgA)/100
  imgB_width=ncol(imgB)/100
  text(pinix-2,piniy-0.2,snsstat[i,4],cex=0.5)
  text(pinix+16,piniy-0.2,snsstat[i,3],cex=0.5)
  rasterImage(imgA,pinix,piniy-0.4,pinix+imgA_width,piniy)
  rasterImage(imgB,sinix,siniy-0.4,sinix+imgB_width,siniy)

  if( any(priseclist$primary == snsstat$HTR[i])){
    temS=paste("../dominating17Apr28/network_motifs/",priseclist[grep(snsstat$HTR[i],priseclist$primary,perl = T),2],".png",sep="")
    imgS=readPNG(temS)
    imgS_width=ncol(imgS)/100
    rasterImage(imgS,pinix +6,piniy-0.4,pinix+6+imgS_width,piniy)
    textS=nsecmotif[grep(priseclist[grep(snsstat$HTR[i],priseclist$primary,perl = T),2],nsecmotif$HTR),3]
    text(pinix+19,piniy-0.2,textS,cex=0.5)
  }
  piniy=piniy- 0.5
  siniy=siniy-0.5
}
dev.off()
