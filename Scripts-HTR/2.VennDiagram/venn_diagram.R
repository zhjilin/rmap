#script to draw the Venn diagram
#update 20180627: add one more category-RNA bind N seq data

args=commandArgs(trailingOnly=TRUE)
library(Vennerable)

rhtset=args[1]
selexset=args[2]
competeset=args[3]
rnbsset=args[4]

# rhtset="RHT_SELEX_2016Dec13.table"
# selexset="SELEX_cook.table"
# competeset="RNA_Compete.table"
# rnbsset="RBNS.tsv"
outprefix="VennDiagram"
tstring=gsub(" .*$","",Sys.time(),perl=TRUE)
outname=paste(outprefix,tstring,".pdf",sep="")

rhtselex=read.table(rhtset,sep="\t",header=T,allowEscapes=F)
selex=read.table(selexset,sep="\t",header=T,allowEscapes=F)
rnacompete=read.table(competeset,sep="\t",header=T,allowEscapes=F)
rnbs=read.table(rnbsset,sep="\t",header=T,allowEscapes=F)

expcom=list("HTR-SELEX"=rhtselex$Approved.symbol,"RNAcompete"=rnacompete$Approved.symbol,"SELEX"=selex$Approved.symbol,"RNA Bind-N-Seq"=rnbs$Approved.symbol)
VSE=Venn(expcom)
V3=compute.Venn(VSE, doWeights = FALSE,type="ellipses")
pdf(outname)
grid.newpage()
plot(V3,show = list( SetLabels = TRUE,Faces = FALSE))
dev.off()
