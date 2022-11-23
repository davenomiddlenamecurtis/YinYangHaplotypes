# run GO analysis on genes overlapping YY haplotypes

wd="C:/Users/dave/OneDrive/sharedseq/yinYang"

AllGenesFile="allGeneRegions.txt"
YYGenesFile="YYGenes.uniq.txt"
ResultsFile="GO.20220809.txt"
setwd(wd)

library(goseq)
AllGenes=read.table(AllGenesFile,header=FALSE,stringsAsFactors=FALSE) 
YYGenes=read.table(YYGenesFile,header=FALSE,stringsAsFactors=FALSE)[,1]
Lengths=AllGenes[,3]-AllGenes[,2]
genes=as.integer(AllGenes[,4]%in%YYGenes)
names(genes)=AllGenes[,4]


pwf=nullp(genes,bias.data=Lengths)
GO=goseq(pwf,"hg19","geneSymbol")
write.table(GO,ResultsFile,sep="\t",quote=FALSE,row.names = FALSE)

