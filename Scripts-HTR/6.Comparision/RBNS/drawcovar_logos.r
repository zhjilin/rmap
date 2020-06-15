library(png)
# library(extrafont)

args <- commandArgs(trailingOnly = TRUE)

# setwd("Downloads/Project/RNASELEX/Figures/HTR-RNAcompete")
tstring=gsub(" .*$","",Sys.time(),perl=TRUE)

outname=paste("htr-rbns",tstring,".pdf",sep="")

htr=args[1]
htrsec=args[2]
rbns=args[3]
rbnssec=args[4]

setwd("/Users/zhangjl/Downloads/Project/RNASELEX/Figures/HTR_RBNS")

htr="RBNS_HTR_overlap_motiflogo/shrinkedLogo/htr_list"
htrsec="RBNS_HTR_overlap_motiflogo/shrinkedLogo/htr-prisec.list"
rbnssec="RBNS_HTR_overlap_motiflogo/shrinkedLogo/rbns-prisec.protein"


htrmotif=read.table(htr,sep="\t",col.names = c("protein","motif"))
rbnsmotif=data.frame(protein=htrmotif$protein,motif=paste(htrmotif$protein,".png",sep=""))

htrmotif_sec=read.table(htrsec,sep="\t",col.names = c("protein","primary","secondary"))
rbnsmotif_sec=read.table(rbnssec,sep="\t",col.names = c("protein"))


htrmotif=htrmotif[order(htrmotif$protein,decreasing = F),]
rbnsmotif=rbnsmotif[order(rbnsmotif$protein,decreasing = F),]
Alist=paste("../dominating17Apr28/network_motifs/",htrmotif$motif,sep="")
Blist=paste("RBNS_HTR_overlap_motiflogo/shrinkedLogo/primary/",rbnsmotif$motif,sep="")

pinix=2
piniy=28
sinix=14
siniy=28

pdf(outname,width=8,height=12)
plot(NA,xlim=c(0,30),ylim=c(0,28), xaxt="n",yaxt="n",bty="n")
text(pinix,piniy+0.5,"HTR-SELEX",cex=0.5)
text(sinix,siniy+0.5,"RBNS",cex=0.5)

for(i in 1:length(Alist)){

  imgA=readPNG(Alist[i])
  imgB=readPNG(Blist[i])
  imgA_width=ncol(imgA)/100
  imgB_width=ncol(imgB)/600
  text(pinix-2,piniy-0.2,htrmotif[i,1],cex=0.5)
  # text(pinix+16,piniy-0.2,snsstat[i,3],cex=0.5)
  rasterImage(imgA,pinix,piniy-0.4,pinix+imgA_width,piniy)
  rasterImage(imgB,sinix,siniy-0.4,sinix+imgB_width,siniy)

  if(any(as.character(htrmotif_sec$protein) == as.character(htrmotif$protein[i]))){
    temS=paste("../dominating17Apr28/network_motifs/",htrmotif_sec[grep(as.character(htrmotif$protein[i]),htrmotif_sec$protein,fixed=T),3],".png",sep="")
    imgS=readPNG(temS)
    imgS_width=ncol(imgS)/100
    rasterImage(imgS,pinix +6,piniy-0.4,pinix+6+imgS_width,piniy)
  #   textS=nsecmotif[grep(priseclist[grep(snsstat$HTR[i],priseclist$primary,perl = T),2],nsecmotif$HTR),3]
  #   text(pinix+19,piniy-0.2,textS,cex=0.5)
  }
  
  if(any(as.character(rbnsmotif_sec$protein) == as.character(htrmotif$protein[i]))){
    temR=paste("RBNS_HTR_overlap_motiflogo/shrinkedLogo/secondary/",as.character(htrmotif$protein[i]),"-2.png",sep="")
    imgR=readPNG(temR)
    imgR_width=ncol(imgR)/600
    rasterImage(imgR,sinix +6,siniy-0.4,sinix+6+imgR_width,siniy)
  }
  piniy=piniy- 0.5
  siniy=siniy-0.5
}
dev.off()
