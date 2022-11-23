# plot genotypes from 1000G
#
# takes options from file provided as first argument as can be seen at the start of the code
# 
# need to have these two files in the local directory:
# igsr-1000genomes30xOnGrch38.tsv - downloaded from 1000g
# children.txt - list of subjects in 1000G who are children, to be excluded, begins:
# HG00405 HG00405
# HG00408 HG00408
# HG00418 HG00418
#
# need to have tabix on the PATH
#
# following two variables must be set to point to the plink executable and the 1000G vcfs
# (note there is a %s in the vcf spec which will be replaced by the chromosome name

plink="/share/apps/genomics/plink-1.9/plink"
GenoVcfTemplate="/cluster/project9/bipolargenomes/popGen/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr%s.recalibrated_variants.vcf.gz"

library(dplyr)
library(tidyr)
library(ggplot2)
wd="C:/Users/dave/OneDrive/sharedseq/popGen"
# setwd(wd)
SelectPos=character()
SelectGeno=character()

SNPsToUse=""

# defaults to be replaced with provided arguments
Chr="3"
HapRange=c(121056668,121325164)
StartPos=HapRange[1]-30000
EndPos=HapRange[2]+30000
MinMAF="0.01"
MinMAF="0.05"
SelectChr="3"
SelectPos[1]="121056668"
SelectGeno[1]="0"
NumToSelect=1

Chr="1"
HapRange=c(113389183,113526463)
StartPos=HapRange[1]-30000
EndPos=HapRange[2]+30000
MinMAF="0.01"
MinMAF="0.05"
SelectChr="1"
SelectPos[1]="113389183"
SelectGeno[1]="0"
NumToSelect=1

# EPAS1
Chr="3"
HapRange=c(45818159,45867532)
StartPos=HapRange[1]
EndPos=HapRange[2]
MinMAF="0.01"
MinMAF="0.05"
SelectChr="3"
SelectPos[1]="45818159"
SelectGeno[1]="0"
NumToSelect=1

ToHighlight=data.frame(matrix(nrow=2,ncol=2))
colnames(ToHighlight)=c("Name","pos")
# ToHighlight$Name[1]="rs10774671"
# ToHighlight$pos[1]=112919388
# ToHighlight$Name[2]="rs15895"
# ToHighlight$pos[2]=113010483

AllSubs="igsr-1000genomes30xOnGrch38.tsv"
# have this local so can run on PC
error=function(ErrStr) {
	print(ErrStr)
	quit("no")
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
	error("Must provide option file name as argument")
}
NumToSelect=0
NumToHighlight=0
TestName=""
KeepVCF=FALSE
KeepFiles=FALSE
UseFixed=TRUE
FixedMAF=0.01
ArgLines=readLines(args[1])
for (l in 1:length(ArgLines)) {
	words=strsplit(ArgLines[l],"\\s+")[[1]]
	if (length(words) >=1) {
	if (words[1]=="--chr") {
		Chr=words[2]
	} else if (words[1]=="--start") {
		StartPos=as.integer(words[2])
	} else if (words[1]=="--end") {
		EndPos=as.integer(words[2])
	} else if (words[1]=="--test-name") {
		TestName=words[2]
	} else if (words[1]=="--keep-vcf") {
		KeepVCF=TRUE
	} else if (words[1]=="--keep-files") {
		KeepFiles=TRUE
	} else if (words[1]=="--not-fixed") { # if selecting by homozygotes, do not include other fixed variants
		KeepFixed=FALSE
	} else if (words[1]=="--maf") {
		MinMAF=words[2]
	} else if (words[1]=="--fixed-maf") {
		FixedMAF=words[2]
	} else if (words[1]=="--snps-to-use") {
		SNPsToUse=words[2]
	} else if (words[1]=="--select") {
		NumToSelect=NumToSelect+1
		SelectPos[NumToSelect]=words[2]
		if (length(words)>2) {
			SelectGeno[NumToSelect]=words[3]
		} else {
			SelectGeno[NumToSelect]="-1"
		}
	} else if (words[1]=="--highlight") {
		if (NumToHighlight>1) {
			error(sprintf("Can highlight maximum of 2 loci: %s",ArgLines[l]))
		}
		NumToHighlight=NumToHighlight+1
		ToHighlight$Name[NumToHighlight]=words[2]
		ToHighlight$pos[NumToHighlight]=as.integer(words[3])
	} else {
		error(sprintf("Unknown option: %s",ArgLines[l]))
	}
	}
}

if (NumToSelect>1) {
	for (s in 1:NumToSelect) {
		if (SelectGeno[s]==-1) {
			error("If selecting on more than one locus must specify selection genotype for each one")
		}
	}
}
SelectChr=Chr

VcfName=sprintf(GenoVcfTemplate,Chr)
if (TestName=="") {
	Root=sprintf("%s.%d.%d",Chr,StartPos,EndPos)
} else {
	Root=sprintf("%s.%s.%d.%d",TestName,Chr,StartPos,EndPos)
}
RangeVcf=sprintf("range.%s.vcf",Root)
BimName=sprintf("markers.%s.bim",Root)
AwkCmd="awk '{ OFS=\"\t\"; $3=\"Pos\"$2; print $0 }'"
if (!file.exists(RangeVcf) & !file.exists(BimName)) {
	CommStr=sprintf("tabix -h %s chr%s:1-1 > %s; tabix %s chr%s:%d-%d | %s >> %s",VcfName,Chr,RangeVcf,VcfName,Chr,StartPos,EndPos,AwkCmd,RangeVcf)
	print(CommStr)
	system(CommStr)
}

if (!file.exists(BimName)) {
	if (SNPsToUse=="") {
		CommStr=sprintf("%s --vcf %s --maf %s --real-ref-alleles --biallelic-only --snps-only  --remove children.txt --make-bed --out markers.%s",plink,RangeVcf,MinMAF,Root)
		print(CommStr)
		system(CommStr)
	} else {
		SNPRange=read.table(SNPsToUse,sep="",header=FALSE,stringsAsFactors=FALSE)
		SNPRange=SNPRange[,1:2]
#		SNPRange[,1]=sprintf("chr%s",SNPRange[,1])
		SNPRange[,3]=SNPRange[,2]
		SNPRange[,4]="ToExtractByPlotGenos"
		write.table(SNPRange,"ToExtractByPlotGenos.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
		CommStr=sprintf("%s --vcf %s --maf %s --real-ref-alleles --biallelic-only --extract range ToExtractByPlotGenos.txt --snps-only  --remove children.txt --make-bed --out markers.%s",plink,RangeVcf,MinMAF,Root)
		print(CommStr)
		system(CommStr)
	}
}
	
	# no need for this to be in a loop, same table applies to all genotypes
	subs=read.table(AllSubs,sep="\t",header=TRUE)
	subs=subs[,c(1,4)]
	colnames(subs)=c("FID","Pop")
	AllPops=c("FIN","CEU","IBS","TSI","GBR","PUR","CLM","PEL","MXL","CHS","CHB","CDX","KHV","JPT","BEB","GIH","PJL","STU","ITU","YRI","LWK","GWD","MSL","ESN","ACB","ASW")
	AfricanPops=c("YRI","LWK","GWD","MSL","ESN","ACB","ASW")
	subs$Pop=factor(subs$Pop,levels=AllPops)
	subs=subs[order(subs$Pop),]
	subs$index=1:nrow(subs) # need this because merge does not keep row order


LastRoot=Root
DefaultGenos=c("0","2")
if (NumToSelect==1 & SelectGeno[1]=="-1" ) {
	NumLoops=2 } else {
	NumLoops=1
}
for (Loop in 1:NumLoops) {
	LastRoot=Root
	if (NumToSelect>0) {
	for (s in 1:NumToSelect) {
		if (SelectGeno[s]==-1) {
			Geno=DefaultGenos[Loop] } else {
			Geno=SelectGeno[s]
		}
		if (TestName=="") {
			SelectRoot=sprintf("%s.%d.%d.%s.%s",Chr,StartPos,EndPos,SelectPos[s],Geno)
			} else {
			SelectRoot=sprintf("%s.%s.%d.%d.%s.%s",TestName,Chr,StartPos,EndPos,SelectPos[s],Geno)
		}
		FilterFile=sprintf("filter.%s.raw",SelectRoot)
		BimName=sprintf("markers.%s.bim",SelectRoot)
		if (!file.exists(FilterFile) & !file.exists(BimName)) {	
			CommStr=sprintf("echo %s %s %s %s > range.txt",Chr,SelectPos,SelectPos,SelectRoot)
			print(CommStr)
			system(CommStr)
			CommStr=sprintf("%s --bfile markers.%s --real-ref-alleles --extract range range.txt --recodeA --out filter.%s",plink,LastRoot,SelectRoot)
			print(CommStr)
			system(CommStr)
		}
		if (!file.exists(BimName)) {
			if (KeepFixed) {
				CommStr=sprintf("%s --bfile markers.%s --real-ref-alleles --make-bed --filter %s %s --mfilter 5 --out markers.%s",plink,LastRoot,FilterFile,Geno,SelectRoot)
			} else {
				CommStr=sprintf("%s --bfile markers.%s --maf %f --real-ref-alleles --make-bed --filter %s %s --mfilter 5 --out markers.%s",plink,LastRoot,FixedMAF,FilterFile,Geno,SelectRoot)
			}
			print(CommStr)
			system(CommStr)
		}
		LastRoot=SelectRoot
	}
	}
	CommStr=sprintf("cat %s | cut -f 2,6 > refAlleles.%s.txt",BimName,LastRoot)
	print(CommStr)
	system(CommStr)
	RawName=sprintf("genotypes.%s.raw",LastRoot)
	if (!file.exists(RawName)) {
		if (KeepFixed) {
			CommStr=sprintf("%s --bfile markers.%s --recodeA --recode-allele refAlleles.%s.txt --out genotypes.%s",plink,LastRoot,LastRoot,LastRoot)
		} else {
			CommStr=sprintf("%s --bfile markers.%s --maf %f --recodeA --recode-allele refAlleles.%s.txt --out genotypes.%s",plink,LastRoot,FixedMAF,LastRoot,LastRoot)
		}
		print(CommStr)
		system(CommStr)
	}
	genos=read.table(RawName,header=TRUE)
	genos=genos[,-(3:6)]
	colnames(genos)[2]="IDnum"
	genos=merge(subs,genos,by="FID",all=FALSE)
	genos=genos[order(genos$index),]
	genos$IDnum=1:nrow(genos)

	rr=1
	last=""
	PopPos=data.frame(matrix(nrow=1,ncol=2))
	colnames(PopPos)=c("Pop","Pos")
	for (r in 1:nrow(genos)) {
		if (genos$Pop[r]!=last) {
			PopPos[rr,1]=toString(genos$Pop[r])
			PopPos[rr,2]=genos$IDnum[r]
			last=genos$Pop[r]
			rr=rr+1
		}
	}

	markers=read.table(BimName,header=FALSE)
	colnames(genos)[5:ncol(genos)]=markers$V4

	toplot=gather(genos,pos,geno,5:ncol(genos))
	toplot$pos=as.numeric(toplot$pos)

	colours=c("blue","plum","red")
	highlightcolours=c("blue4","plum3","red4") # not using these, can change if needed

	toplot$col=colours[toplot$geno+1]
	toplot$col[is.na(toplot$col)]="white"
	for (r in 1:nrow(ToHighlight)) {
		if (ToHighlight$Name[r]!="") {
			toplot$col[toplot$pos==ToHighlight$pos[r]]=highlightcolours[toplot$geno[toplot$pos==ToHighlight$pos[r]]+1]
			highlighted=toplot[toplot$pos==ToHighlight$pos[r],]
			toplot=toplot[toplot$pos!=ToHighlight$pos[r],]
			toplot=rbind(toplot,highlighted) 
		# put these at the bottom so they will be plotted last and not overwritten
		}
	}
	p=ggplot(data=toplot,aes(x=IDnum,y=pos))+geom_point(colour=toplot$col,size=0.05,shape=15) + theme_bw()
	p=p + theme(axis.text.y = element_text(size = 24)) + theme(axis.text.x = element_blank()) + theme(axis.title.x = element_blank())
	p=p+ylab("Position") + theme(axis.title.y = element_text(size = 24))
	lastPos=0
	for (r in 1:nrow(ToHighlight)) {
	# p=p+geom_segment(aes(x = -50, y = ToHighlight$pos[r], xend = 0, yend = ToHighlight$pos[r],colour = "black"))
	# this works differently inside a loop than outside
	# discussion here: https://stackoverflow.com/questions/29037945/using-geom-segment-within-a-for-loop
		p=p+annotate(geom="text",label=ToHighlight$Name[r],x=-130/2000*nrow(genos),y=max(ToHighlight$pos[r],lastPos+2000),size=6)
		lastPos=ToHighlight$pos[r]
	}
	if (NumToHighlight>=1) {
		p=p+geom_segment(aes(x = -50/2000*nrow(genos), y = ToHighlight$pos[1], xend = 0, yend = ToHighlight$pos[1]))
	}
	if (NumToHighlight>1) {
		p=p+geom_segment(aes(x = -50/2000*nrow(genos), y = ToHighlight$pos[2], xend = 0, yend = ToHighlight$pos[2]))
	}
# use vector instead of adding segments in loop
# because loop does not work!!
	for (r in 1:nrow(PopPos)) {
		p=p+annotate(geom="text",label=PopPos$Pop[r],x=PopPos$Pos[r],y=StartPos-5000,size=8)
	# trying to use geom_text() takes forever - I think it is applied to every point
	}
	PngName=sprintf("genos.%s.png",LastRoot)
	options(bitmapType='cairo') # make png without X11
# png(height=12000, width=24000, pointsize=6, file=PngName,res=720)
# 	png(height=9*ceiling((EndPos-StartPos)/100), width=3000+50*nrow(genos), pointsize=6, file=PngName,res=720)
 	png(height=9*ceiling((EndPos-StartPos)/100), width=24000, pointsize=6, file=PngName,res=720)
	print(p)
	dev.off()
}
if (KeepVCF==FALSE) {
	CommStr=sprintf("rm %s",RangeVcf)
	print(CommStr)
	system(CommStr)
}
if (KeepFiles==FALSE) {
	CommStr=sprintf("rm *.%s.*nosex *.%s.*log *.%s.*fam *.%s.*bim *.%s.*bed *.%s.*raw refAlleles.%s.* ToExtractByPlotGenos.txt",Root,Root,Root,Root,Root,Root,Root)
	print(CommStr)
	system(CommStr)
}
