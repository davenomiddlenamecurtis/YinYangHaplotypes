#!/share/apps/R-4.0.3/bin/Rscript

AllChr=c(1:22) # we do not have ancestral alleles for X chromosome
# AllChr=c(4)
ResTemplate="/cluster/project9/bipolargenomes/yinYang/chr%s.results.txt"
AllsTemplate="YinYangAlls.chr%s.20221101"
VCFTemplate="/SAN/ugi/Speidel/datasets/hg38/Bergstrom2018HGDP/hgdp.v0.5.archaics.chr%s.vcf.gz"
TsvTemplate="YinYangAlls.chr%s.20221101.txt"
plink="/share/apps/genomics/plink-1.9/plink"
AllResultsFile="/home/rejudcu/yinYang/AllYYResults.20221107.txt"
ValidRegions=read.table(AllResultsFile,header=TRUE,stringsAsFactors=FALSE,sep="\t",quote="")[c("Chr","Start","End")] # removed duplications etc.
SubsToKeep=c("Chimp","ancestor_ensembl_e86","Vindija33.19","AltaiNeandertal","Denisova")
writeLines(SubsToKeep,file("SubsToKeep.txt"))

for (Chr in AllChr) {
	ResFile=sprintf(ResTemplate,Chr)
	VCFFile=sprintf(VCFTemplate,Chr)
	AllsRoot=sprintf(AllsTemplate,Chr)
	Valid=ValidRegions[ValidRegions$Chr==Chr,]
	Results=read.table(ResFile,header=FALSE,stringsAsFactors=FALSE)
	Results[,1]=sprintf("chr%s",Results[,1])
	TabixRegionsFile=sprintf("ForTabix.%s.txt",Chr)
	for (r in 1:nrow(Valid)) {
		ThisHap=Results[Results$V2>=Valid$Start[r]&Results$V2<=Valid$End[r],]
		if (r==1) {
			ForTabix=ThisHap
		} else {
			ForTabix=rbind(ForTabix,ThisHap)
		}
	}
	write.table(ForTabix[,1:2],TabixRegionsFile,quote=FALSE,row.names = FALSE, col.names=FALSE,sep = "\t")
	PlinkFile=sprintf("AllsForPlink.%s.vcf",Chr)
	CommStr=sprintf("tabix %s -h -R %s > %s",VCFFile,TabixRegionsFile,PlinkFile)
	system(CommStr)
	CommStr=sprintf("%s --vcf %s --real-ref-alleles --biallelic-only --snps-only --double-id --keep-fam SubsToKeep.txt --recode transpose --out %s",plink,PlinkFile,AllsRoot)
	system(CommStr)
	JoinedFile=sprintf("AllsJoined.chr%s.txt",Chr)	
	CommStr=sprintf("bash -c 'cut -f 1,2,4,5 %s | sort -k 2b,2 | join - <(cat %s.tped | sort -k 4b,4) -1 2 -2 4 | sort -n > %s'",PlinkFile,AllsRoot,JoinedFile)
	system(CommStr)
	Joined=read.table(JoinedFile,header=FALSE,stringsAsFactors=FALSE)
	Joined[,2]=Joined[,1]
	Joined[,1]=Chr
	Joined=Joined[,-c(5,6,7,10,11,12,13)]
	colnames(Joined)=c("Chr","Pos","REF","ALT","Chimp1","Chimp2","Altai1","Altai2","Den1","Den2")
	TsvFile=sprintf(TsvTemplate,Chr)
	write.table(Joined,TsvFile,quote=FALSE,row.names = FALSE, sep = "\t")
}


