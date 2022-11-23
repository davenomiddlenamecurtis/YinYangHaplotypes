#!/share/apps/R-4.0.3/bin/Rscript
library(dplyr)
library(tidyr)

# get global Fst
# uses second equation here: https://en.wikipedia.org/wiki/Fixation_index

SummFileTemplate="noSwap.20220727.%s.summ.txt"
MinMAF="0.05"

plink="/share/apps/genomics/plink-1.9/plink"
GenoVcfTemplate="/cluster/project9/bipolargenomes/popGen/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr%s.recalibrated_variants.vcf.gz"

AllChrs=c(c(1:22),"X")
AllSubs="igsr-1000genomes30xOnGrch38.tsv"
				subs=read.table(AllSubs,sep="\t",header=TRUE)
				subs=subs[,c(1,6)]
				colnames(subs)=c("FID","SuperPop")
				AllSuperPops=c("AFR","AMR","EAS","EUR","SAS")
				subs$SuperPop=factor(subs$SuperPop,levels=AllSuperPops)
				subs=subs[order(subs$SuperPop),]
				subs$index=1:nrow(subs) # need this because merge does not keep row order


for (Chr in AllChrs) {
AllResults=data.frame(matrix(nrow=0,ncol=4+length(AllSuperPops)))
colnames(AllResults)[1:4]=c("Chr","Pos","Fst","NonAfrFST")
colnames(AllResults)[5:(4+length(AllSuperPops))]=AllSuperPops
AllResR=1
FstFile=sprintf("chr.%s.Fsts.20220727.txt",Chr)
	if (!file.exists(FstFile)) {
	SummFile=sprintf(SummFileTemplate,Chr)
	Results=read.table(SummFile,header=FALSE)
	colnames(Results)=c("Chr","First","Last","Number","Length")

	for (r in 1: nrow(Results)) {
				StartPos=Results[r,2]
				EndPos=Results[r,2]
				AllResults[AllResR,1]=Chr
				AllResults[AllResR,2]=StartPos
				StartPos=Results[r,2]
				EndPos=Results[r,2]
				VcfName=sprintf(GenoVcfTemplate,Chr)
				RangeVcf=sprintf("TempRange.%s.%d.%d.vcf",Chr,StartPos,EndPos)
				BimName=sprintf("markers.%s.%d.%d.bim",Chr,StartPos,EndPos)
				AwkCmd="awk '{ OFS=\"\t\"; $3=\"Pos\"$2; print $0 }'"
				CommStr=sprintf("if [ ! -e %s -a ! -e %s ] ; then tabix -h %s chr%s:1-1 > %s; tabix %s chr%s:%d-%d | %s >> %s; fi",RangeVcf,BimName,VcfName,Chr,RangeVcf,VcfName,Chr,StartPos,EndPos,AwkCmd,RangeVcf)
				print(CommStr)
				system(CommStr)
				BimName=sprintf("markers.%s.%d.%d.bim",Chr,StartPos,EndPos)
				CommStr=sprintf("if [ ! -e %s ] ; then %s --vcf %s --real-ref-alleles --maf %s --biallelic-only --snps-only  --remove children.txt --make-bed --out markers.%s.%d.%d; fi",BimName,plink,RangeVcf,MinMAF,Chr,StartPos,EndPos)
# means we can be sure we have the right variant, in case more than one at this position
				print(CommStr)
				system(CommStr)
				CommStr=sprintf("rm %s",RangeVcf)
				print(CommStr)
				system(CommStr)
				CommStr=sprintf("cat %s | cut -f 2,5 > refAlleles.%s.%d.%d.txt",BimName,Chr,StartPos,EndPos)
				print(CommStr)
				system(CommStr)
				RawName=sprintf("genotypes.%s.%d.%d.raw",Chr,StartPos,EndPos)
				CommStr=sprintf("if [ ! -e %s ] ; then %s --bfile markers.%s.%d.%d --recodeA --recode-allele refAlleles.%s.%d.%d.txt --out genotypes.%s.%d.%d; fi",RawName,plink,Chr,StartPos,EndPos,Chr,StartPos,EndPos,Chr,StartPos,EndPos)
				print(CommStr)
				system(CommStr)
				genos=read.table(RawName,header=TRUE)
				genos=genos[,-(3:6)]
				colnames(genos)[2]="IDnum"
				genos=merge(subs,genos,by="FID",all=FALSE)
				genos=genos[order(genos$index),]
				genos$IDnum=1:nrow(genos)
# to get overall Fst, use second formula here (based on allele frequencies): https://en.wikipedia.org/wiki/Fixation_index		
				rr=0
				last=""
				SuperPopPos=data.frame(matrix(nrow=length(AllSuperPops),ncol=3))
				colnames(SuperPopPos)=c("SuperPop","N","AC")
				genos$SuperPop=as.character(genos$SuperPop)
				genos=genos[!is.na(genos[,5]),]
				for (rrr in 1:nrow(genos)) {
					if (genos$SuperPop[rrr]!=last) {
						last=genos$SuperPop[rrr]
						rr=rr+1
						SuperPopPos$SuperPop[rr]=genos$SuperPop[rrr]
						SuperPopPos$N[rr]=0
						SuperPopPos$AC[rr]=0
					}
					SuperPopPos$N[rr]=SuperPopPos$N[rr]+2
					SuperPopPos$AC[rr]=SuperPopPos$AC[rr]+genos[rrr,5]
				}
				for (OnlyNonAfr in 0:1)
				{
				TotalN=0
				TotalP=0
				for (rr in (1+OnlyNonAfr):nrow(SuperPopPos)) {
					TotalN=TotalN+SuperPopPos$N[rr]
					TotalP=TotalP+SuperPopPos$AC[rr]
					AllResults[AllResR,4+rr]=SuperPopPos$AC[rr]/SuperPopPos$N[rr]
				}
				SuperPopPos$P=SuperPopPos$AC/SuperPopPos$N
				TotalP=TotalP/TotalN
				S=0
				for (rr in (1+OnlyNonAfr):nrow(SuperPopPos)) {
					S=S+SuperPopPos$N[rr]/TotalN*SuperPopPos$P[rr]*(1-SuperPopPos$P[rr]) # weight by sample size
				}
				Fst=(TotalP*(1-TotalP)-S)/(TotalP*(1-TotalP))
				AllResults[AllResR,3+OnlyNonAfr]=Fst
				}
				AllResR=AllResR+1
				CommStr=sprintf("rm *.%s.%d.%d.*",Chr,StartPos,EndPos)
				system(CommStr)
				}
	write.table(AllResults,FstFile,quote=FALSE,row.names = FALSE)
	}
}

CommStr="head -n1 chr1.Fsts.txt > allFsts.20220705.txt"
system(CommStr)
for (Chr in AllChrs) {
	CommStr=sprintf("tail -n +2 chr%s.Fsts.txt >> allFsts.20220705.txt",Chr)
	system(CommStr)
}

