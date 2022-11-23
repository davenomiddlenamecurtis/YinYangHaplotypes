#!/share/apps/R-4.0.3/bin/Rscript
library(dplyr)

AllChr=c(c(1:22),"X") 

AllResultsFile="/home/rejudcu/yinYang/AllYYResults.20221107.txt"
AncSummaryFile="/home/rejudcu/yinYang/Anc.chrAll.summ.20221101.txt"
SupTableFile="/home/rejudcu/yinYang/SupplementaryTable.20221101.txt"

Results=read.table(AllResultsFile,header=TRUE,stringsAsFactors=FALSE,sep="\t",quote="")

Results=Results[c("Coord","Chr","Start","End","Number","Length","Variant","VarName","CountAA","CountAB","CountBB","Fst","NonAfrFST","AFR","AMR","EAS","EUR","SAS","HetMean0","HetMean2")]
# AllCountAB etc. are to do with read depths
colnames(Results)[9:11]=c("CountRR","CountRA","CountAA")
colnames(Results)[19:20]=c("HetMeanRR","HetMeanAA")
Ancs=read.table(AncSummaryFile,header=TRUE,stringsAsFactors=FALSE,sep="\t",quote="")
Ancs$Coord=sprintf("%s:%s",Ancs$Chr,Ancs$Start)
Ancs=Ancs[c("Coord","BothAlt","BothRef","ChimpAlt","NeandAlt")]
NumAncs=rowSums(Ancs[,2:5])
Ancs$Derived=(NumAncs>=5 & (Ancs$BothAlt==NumAncs | Ancs$BothRef==NumAncs))*1
Ancs$Ancient=(NumAncs>=5 & (Ancs$ChimpAlt==NumAncs | Ancs$NeandAlt==NumAncs))*1
Ancs$Mixed=(NumAncs>=5 & (Ancs$ChimpAlt>=1 | Ancs$NeandAlt>=1))*1
sprintf("Number of human derived haplotypes = %d",sum(Ancs$Derived))
Results=merge(Results,Ancs,by="Coord",all.x=TRUE)
Results=select(Results,-c(Coord))
Results=arrange(Results,factor(Chr,levels=AllChr),Start)
Results[is.na(Results)]=0 # X chromosome 
write.table(Results,SupTableFile,sep="\t",quote=FALSE,row.names = FALSE)
