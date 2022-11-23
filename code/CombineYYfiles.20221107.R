# combine all YY results files and separate out the ones which look like CNVs
library(dplyr)

SummaryFile="/cluster/project9/bipolargenomes/yinYang/allChrsSummary.20220727.txt"
HWEFile="/cluster/project9/bipolargenomes/yinYang/getYYHWE.20220801.txt"
HetFile="/cluster/project9/bipolargenomes/yinYang/HetResults.20220727.txt"
FstsFile="/cluster/project9/bipolargenomes/yinYang/allFsts.20220727.txt"
GeneFile="/cluster/project9/bipolargenomes/yinYang/mg.matchGenes.20220727.txt"
DepthFile="/cluster/project9/bipolargenomes/yinYang/results.allDepths.20220727.txt"
VEPFile="/cluster/project9/bipolargenomes/yinYang/VEP.annots.20220727.txt"
GWASHitsFile="/cluster/project9/bipolargenomes/yinYang/YY.GWAScat.20220727.tsv"
AllGenesFile="/home/rejudcu/yinYang/allGeneRegions.txt"
ResultsFile="/home/rejudcu/yinYang/AllYYResults.20221107.txt"
GOResultsFile="GO.20221107.txt"
# sink("CombineYYFiles.20221107.log")

AllChr=c(c(1:22),"X")

SummaryTable=read.table(SummaryFile,sep="",header=FALSE,stringsAsFactors=FALSE)
SummaryTable$Coord=sprintf("%s:%s",SummaryTable[,1],SummaryTable[,2])
SummaryTable=SummaryTable[,c(6,1:5)]
colnames(SummaryTable)=c("Coord","Chr","Start","End","Number","Length")
head(SummaryTable)

CommStr=sprintf("tail -n +52 %s > VEP.temp.annots.txt",VEPFile) # cannot use header line because starts with hash, FFS
print(CommStr)
system(CommStr)
VEPTable=read.table("VEP.temp.annots.txt",sep="\t",header=FALSE,stringsAsFactors=FALSE)
SNPsTable=VEPTable[,c(1,2,13)]
colnames(SNPsTable)=c("Variant","Coord","VarName")
head(SNPsTable)

HWETable=read.table(HWEFile,sep="",header=FALSE,stringsAsFactors=FALSE)
HWETable=HWETable[,c(1,2,9,11,8)]
HWETable$Coord=sprintf("%s:%s",HWETable[,1],HWETable[,2])
HWETable=HWETable[,c(6,1:5)]
HWETable=HWETable[-c(2:3)]
HWETable[,c(7,6,5)]=as.numeric(do.call(rbind,strsplit(HWETable[,4],"/")))
colnames(HWETable)=c("Coord","Het","Pval","GenosBBABAA","CountAA","CountAB","CountBB")
head(HWETable)

DepthTable=read.table(DepthFile,sep="",header=FALSE,stringsAsFactors=FALSE)
DepthTable$Coord=sprintf("%s:%s",DepthTable[,1],DepthTable[,2])
DepthTable=DepthTable[,c(24,1:23)]
colnames(DepthTable)=c("Coord","Chr","Start","AllCountAA","AllCountAB","AllCountBB",
	"MeanAA_A","MeanAA_B","SDAA_A","SDAA_B","AllBalAA","SDAllBalAA",
	"MeanAB_A","MeanAB_B","SDAB_A","SDAB_B","AllBalAB","SDAllBalAB",
	"MeanBB_A","MeanBB_B","SDBB_A","SDBB_B","AllBalBB","SDAllBalBB")
DepthTable=DepthTable[-c(2:3)]
head(DepthTable)

FstsTable=read.table(FstsFile,sep="",header=FALSE,stringsAsFactors=FALSE)
FstsTable$Coord=sprintf("%s:%s",FstsTable[,1],FstsTable[,2])
FstsTable=FstsTable[,c(10,1:9)]
colnames(FstsTable)=c("Coord","Chr","Pos","Fst","NonAfrFST","AFR","AMR","EAS","EUR","SAS")
FstsTable=FstsTable[,-c(2,3)]
head(FstsTable)

HetTable=read.table(HetFile,sep="",header=TRUE,stringsAsFactors=FALSE)
HetTable=HetTable[-c(2:4)]
head(HetTable)

GeneTable=read.table(GeneFile,sep="\t",header=FALSE,stringsAsFactors=FALSE)
GeneTable$Coord=sprintf("%s:%s",GeneTable[,1],GeneTable[,2])
GeneTable=GeneTable[,c(6,1:5)]
colnames(GeneTable)=c("Coord","Chr","Start","End","GeneNumber","Genes")
GeneTable=GeneTable[-c(2:4)]
head(GeneTable)

GWASHits=read.table(GWASHitsFile,sep="\t",header=FALSE,stringsAsFactors=FALSE,quote="")
GWASHits=GWASHits[,c(2,8,22)]
colnames(GWASHits)=c("PubMed","Pheno","SNP")
head(GWASHits)

AllYYResults=merge(SummaryTable,SNPsTable,by="Coord")
AllYYResults=merge(AllYYResults,HWETable,by="Coord")
# AllYYResults=merge(SummaryTable,HWETable,by="Coord")
AllYYResults=merge(AllYYResults,DepthTable,by="Coord")
AllYYResults=merge(AllYYResults,FstsTable,by="Coord")
AllYYResults=merge(AllYYResults,HetTable,by="Coord")
AllYYResults=merge(AllYYResults,GeneTable,by="Coord")
for (r in 1:nrow(AllYYResults)) {
	AllYYResults$Pheno[r]=paste(GWASHits$Pheno[GWASHits$SNP==AllYYResults$VarName[r]],collapse="; ")
	AllYYResults$PubMed[r]=paste(GWASHits$PubMed[GWASHits$SNP==AllYYResults$VarName[r]],collapse="; ")
}
head(AllYYResults)

AllYYResults=arrange(AllYYResults,factor(Chr,levels=AllChr),Start)

TotalTyped=AllYYResults$CountAA+AllYYResults$CountAB+AllYYResults$CountBB
AllYYResults=AllYYResults[TotalTyped>2000,]

IsRepeat=AllYYResults$Het>0.53
print(sprintf("Number of apparent repeats is %d",sum(IsRepeat)))
HetDepth=AllYYResults$MeanAB_A+AllYYResults$MeanAB_B
print(sprintf("Read depth of repeats: mean = %.1f, sd = %.2f",mean(HetDepth[IsRepeat]),sd(HetDepth[IsRepeat])))
print(sprintf("Read depth of other SNPs: mean = %.1f, sd = %.2f",mean(HetDepth[-IsRepeat]),sd(HetDepth[-IsRepeat])))
print(sprintf("Allele balance of repeats: mean = %.1f, sd = %.2f",mean(AllYYResults$AllBalAB[IsRepeat]),sd(AllYYResults$AllBalAB[IsRepeat])))
print(sprintf("Allele balance of other SNPs: mean = %.1f, sd = %.2f",mean(AllYYResults$AllBalAB[-IsRepeat]),sd(AllYYResults$AllBalAB[-IsRepeat])))
print(sprintf("Range of allele balance of repeats: min = %.2f, max= %.2f",min(AllYYResults$AllBalAB[IsRepeat]),max(AllYYResults$AllBalAB[IsRepeat])))

print(sprintf("Length of repeats: mean = %.1f, sd = %.2f",mean(AllYYResults$Length[IsRepeat]),sd(AllYYResults$Length[IsRepeat])))
print(sprintf("Length of other SNPs: mean = %.1f, sd = %.2f",mean(AllYYResults$Length[-IsRepeat]),sd(AllYYResults$Length[-IsRepeat])))

AllYYResults=AllYYResults[!IsRepeat,]

FewBB=AllYYResults$CountBB<20
print(sprintf("Number with few BB homozygotes %d",sum(FewBB)))
AllYYResults=AllYYResults[!FewBB,]

write.table(AllYYResults,ResultsFile,sep="\t",quote=FALSE,row.names = FALSE)
print("Number of yin yang haplotypes after removing repeats etc.")
nrow(AllYYResults)

print(sprintf("Number of SNPs: mean = %.1f, sd = %.1f",mean(AllYYResults$Number),sd(AllYYResults$Number)))
MaxRow=AllYYResults[AllYYResults$Number==max(AllYYResults$Number),]
print(sprintf("Maximum number of SNPs = %d at %s:%d-%d",MaxRow$Number,MaxRow$Chr,MaxRow$Start,MaxRow$End))
print("Haplotype with most SNPs:")
print(MaxRow)
print(sprintf("Length: mean = %.1f, sd = %.1f",mean(AllYYResults$Length),sd(AllYYResults$Length)))
MaxRow=AllYYResults[AllYYResults$Length==max(AllYYResults$Length),]
print(sprintf("Maximum length = %f at %s:%d-%d (consisting of %d SNPs)",MaxRow$Length,MaxRow$Chr,MaxRow$Start,MaxRow$End,MaxRow$Number))
print(sprintf("Number of haplotypes longer than 100,000 = %d",sum(AllYYResults$Length>100000)))
print("Longest haplotype:")
print(MaxRow)

BackgroundHet=na.omit(select(AllYYResults,"HetMean0","HetMean2"))
print(sprintf("Mean heterozygosity on RR background is %.3f, maximum is %.3f",mean(BackgroundHet$HetMean0),max(BackgroundHet$HetMean0)))
print("Haplotype with maximum heterozygosity on RR background:")
print(AllYYResults[AllYYResults$HetMean0==max(BackgroundHet$HetMean0),])
print(sprintf("Mean heterozygosity on AA background is %.3f, maximum is %.3f",mean(BackgroundHet$HetMean2),max(BackgroundHet$HetMean2)))
print("Haplotype with maximum heterozygosity on AA background:")
print(AllYYResults[AllYYResults$HetMean2==max(BackgroundHet$HetMean2),]) # need to ignore NA
print(sprintf("Number of pairs for which one haplotype has heterozygosity > 0.01 and the other < 0.001 is %d",
	sum((BackgroundHet$HetMean0>0.01&BackgroundHet$HetMean2<0.001) | (BackgroundHet$HetMean2>0.01&BackgroundHet$HetMean0<0.001))))
HetDiff=abs(BackgroundHet$HetMean0-BackgroundHet$HetMean2)
MaxHetDiff=max(HetDiff)
print("Haplotype with maximum difference in background heterozygosity:")
AllYYResults[AllYYResults$HetMean0==BackgroundHet[HetDiff==MaxHetDiff,1] & AllYYResults$HetMean2==BackgroundHet[HetDiff==MaxHetDiff,2],]

print(sprintf("Global Fst has mean %.3f (sd = %.3f, maximum = %.3f)", mean(AllYYResults$Fst), sd(AllYYResults$Fst), max(AllYYResults$Fst)))
print("Haplotype with maximum Fst:")
print(AllYYResults[AllYYResults$Fst==max(AllYYResults$Fst),])
print(sprintf("Global Fst excluding Africa has mean %.3f (sd = %.3f, maximum = %.3f)", mean(AllYYResults$NonAfrFST), sd(AllYYResults$NonAfrFST), max(AllYYResults$NonAfrFST)))
print("Haplotype with maximum Fst excluding Africa:")
print(AllYYResults[AllYYResults$NonAfrFST==max(AllYYResults$NonAfrFST),])
HighFst=select(AllYYResults[AllYYResults$NonAfrFST>0.2,],"Coord","NonAfrFST","AFR","AMR","EAS","EUR","SAS")
print(sprintf("Number of haplotypes with Non-African Fst>0.2 = %d",nrow(HighFst)))
for (c in 3:7) {
	print(sprintf("MAF in %s: mean = %.3f (sd = %3f)",colnames(HighFst)[c],mean(HighFst[,c]),sd(HighFst[,c])))
}
cr=cor.test(AllYYResults$Length,AllYYResults$NonAfrFST)
print("Correlation between non-African Fst and length:")
print(cr)
t=t.test(AllYYResults$Length[AllYYResults$NonAfrFST<=0.2],AllYYResults$Length[AllYYResults$NonAfrFST>0.2])
print("T test of lengths with non-African Fst <=0.2 or >0.2:")
print(t)

print(sprintf("Number (percent) of yin yang haplotypes intercepting at least one gene = %d (%.1f)",
	sum(AllYYResults$GeneNumber>0),sum(AllYYResults$GeneNumber>0)/nrow(AllYYResults)*100))
print("T test of non-African Fst for haplotypes hitting genes or not:")
print (t.test(AllYYResults[AllYYResults$GeneNumber>0,]$NonAfrFST,AllYYResults[AllYYResults$GeneNumber==0,]$NonAfrFST))

# AllYYResults=read.table(ResultsFile,sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="")
YYGenes=unique(unlist(strsplit(AllYYResults$Genes[AllYYResults$GeneNumber>0]," ")))
library(goseq)
AllGenes=read.table(AllGenesFile,header=FALSE,stringsAsFactors=FALSE) 
Lengths=AllGenes[,3]-AllGenes[,2]
genes=as.integer(AllGenes[,4]%in%YYGenes)
names(genes)=AllGenes[,4]
# pwf=nullp(genes,bias.data=Lengths) - needs X and does not work on login node
pwf=nullp(genes,bias.data=Lengths,plot.fit=FALSE)
GO=goseq(pwf,"hg19","geneSymbol")
write.table(GO,GOResultsFile,sep="\t",quote=FALSE,row.names = FALSE)
enriched.GO=GO$category[p.adjust(GO$over_represented_pvalue,method="BH")<.05]
print("Enriched GO terms:")
print(head(enriched.GO))

# Phenos=select(AllYYResults[AllYYResults$Pheno!="",],"Coord","Pheno") does not work - I suspect select function has been replaced
Phenos=AllYYResults[AllYYResults$Pheno!="",][c("Coord","VarName","Pheno")]
colnames(Phenos)=c("Coordinate","SNP","Phenotype")
write.table(Phenos,"YYPhenos.txt",sep="\t",quote=FALSE,row.names = FALSE)


