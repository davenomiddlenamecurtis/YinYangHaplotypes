# compare heterozygosity of yin yang haplotype homozygotes in 1000G
#
# takes one arguments: summary file of all yin yang haplotypes with lines:
# chr startPos endPos 
# (if no arguments provided defaults will be used round OAS1)
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

wd="C:/Users/dave/OneDrive/sharedseq/popGen"
# setwd(wd)

SummaryFile="/cluster/project9/bipolargenomes/yinYang/allChrsSummary.20220727.txt"

AllSubs="igsr-1000genomes30xOnGrch38.tsv"
# have this local so can run on PC
error=function(ErrStr) {
	print(Errstr)
	quit("no")
}

TestName=""
MinMAF=0.05
KeepVCF=FALSE
KeepFiles=FALSE
ArgFile="testHWE.arg"

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
	error("Must provide option file name as argument")
}
ArgFile=args[1]

ArgLines=readLines(ArgFile)
for (l in 1:length(ArgLines)) {
	words=strsplit(ArgLines[l],"\\s+")[[1]]
	if (length(words)<1) {
		next
	}
	if (words[1]=="--summary-file") {
		SummaryFile=words[2]
	} else if (words[1]=="--test-name") {
		TestName=words[2]
	} else if (words[1]=="--keep-vcf") {
		KeepVCF=TRUE
	} else if (words[1]=="--keep-files") {
		KeepFiles=TRUE
	} else if (words[1]=="--maf") {
		MinMAF=words[2]
	} else {
		error(sprintf("Unknown option: %s",ArgLines[l]))
	}
}
lines=readLines(SummaryFile)
HetResults=data.frame(matrix(nrow=length(lines),ncol=8))
colnames(HetResults)=c("Coord","Chr","Start","End","HetMean0","HetSD0","HetMean2","HetSD2")
for (ll in 1:length(lines)) {
	words=strsplit(lines[ll],"\\s+")[[1]]
	Chr=words[1]
	StartPos=words[2]
	EndPos=words[3]
	HetResults$Coord[ll]=sprintf("%s:%s",Chr,StartPos)
	HetResults$Chr[ll]=Chr
	HetResults$Start[ll]=StartPos
	HetResults$End[ll]=EndPos
	VcfName=sprintf(GenoVcfTemplate,Chr)
	if (TestName=="") {
		Root=sprintf("HWE.%s.%s.%s",Chr,StartPos,EndPos)
	} else {
		Root=sprintf("%s.%s.%s.%s",TestName,Chr,StartPos,EndPos)
	}
	RangeVcf=sprintf("range.%s.vcf",Root)
	AwkCmd="awk '{ OFS=\"\t\"; $3=\"Pos\"$2; print $0 }'"
	if (!file.exists(RangeVcf)) {
		CommStr=sprintf("tabix -h %s chr%s:1-1 > %s; tabix %s chr%s:%s-%s | %s >> %s",VcfName,Chr,RangeVcf,VcfName,Chr,StartPos,EndPos,AwkCmd,RangeVcf)
		print(CommStr)
		system(CommStr)
	}
	Genos=c("0","2")
	for (g in 1:2) {
		Geno=Genos[g]
		SelectRoot=sprintf("%s.%s",Root,Geno)
		HWEName=sprintf("HWE.%s.txt",SelectRoot)
		if (!file.exists(HWEName)) {
			FilterFile=sprintf("filter.%s.raw",SelectRoot)
			SelectPos=StartPos
			CommStr=sprintf("echo %s %s %s %s > range.txt",Chr,SelectPos,SelectPos,SelectRoot)
			print(CommStr)
			system(CommStr)
			CommStr=sprintf("%s --vcf %s --maf %s --real-ref-alleles --biallelic-only --snps-only  --remove children.txt  --extract range range.txt --recodeA --out filter.%s",plink,RangeVcf,MinMAF,SelectRoot)
			print(CommStr)
			system(CommStr)
			CommStr=sprintf("%s --vcf %s --maf %s --real-ref-alleles --biallelic-only --snps-only  --remove children.txt --filter %s %s --mfilter 5 --hardy",plink,RangeVcf,MinMAF,FilterFile,Geno)
			print(CommStr)
			system(CommStr)
			CommStr=sprintf("tail -n +2 plink.hwe > %s",HWEName)
			print(CommStr)
			system(CommStr)
		}
		if (file.exists("plink.hwe")) {
			HWE=read.table(HWEName,sep="",header=FALSE)
			HetResults[ll,3+g*2]=mean(HWE[,7],na.rm=TRUE)
			HetResults[ll,4+g*2]=sd(HWE[,7],na.rm=TRUE)
		} else { # no valid subjects with this homozygous genotype
			HetResults[ll,3+g*2]=NA
			HetResults[ll,4+g*2]=NA
		}
		if (KeepFiles==FALSE) {
			CommStr=sprintf("rm plink.* *%s*" ,SelectRoot)
			print(CommStr)
			system(CommStr)
		}
	}
	if (KeepVCF==FALSE) {
		CommStr=sprintf("rm %s",RangeVcf)
		print(CommStr)
		system(CommStr)
	}
}

write.table(HetResults,sprintf("HetResults.%s.txt",TestName),quote=FALSE,row.names = FALSE)

 
