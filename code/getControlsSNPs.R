#!/share/apps/R-4.0.3/bin/Rscript
# select a sample of control SNPs

ResultsFile="/home/rejudcu/yinYang/SupplementaryTable.20221101.txt"
DistTemplate="/cluster/project9/bipolargenomes/yinYang/results.distances.noSwap.20220727.%s.%s.1.250000000.txt"
SampledSNPsTemplate="/cluster/project9/bipolargenomes/yinYang/SampledSNPs.%s.txt"
SampledSNPsVCFTemplate="/cluster/project9/bipolargenomes/yinYang/SampledSNPs.%s.vcf"
SampledSNPsFreqTemplate="/cluster/project9/bipolargenomes/yinYang/SampledSNPs.freq.%s"
ControlSNPsTemplate="Controls.%s.summ.txt"

plink="/share/apps/genomics/plink-1.9/plink"
GenoVcfTemplate="/cluster/project9/bipolargenomes/popGen/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr%s.recalibrated_variants.vcf.gz"

AllChrs=c(c(1:22),"X")
Results=read.table(ResultsFile,sep="\t",header=TRUE,quote="")

error=function(ErrStr) {
	print(Errstr)
	quit("no")
}

for (Chr in AllChrs) {
	SampledSNPsFile=sprintf(SampledSNPsTemplate,Chr)
	if (file.exists(SampledSNPsFile)) {
		next
	}
	DistFile=sprintf(DistTemplate,Chr,Chr)
	Dist=read.table(DistFile,header=FALSE)
	colnames(Dist)[1]="Pos"
	Dist$Chr=sprintf("chr%s",Chr)
	Dist=Dist[c("Chr","Pos")]
	Dist$RowNum=as.numeric(rownames(Dist))
	Dist=Dist[(Dist$RowNum%%100)==0,]
	Dist=Dist[c("Chr","Pos")]
	Dist$Pos2=Dist$Pos
	write.table(Dist,SampledSNPsFile,quote=FALSE,row.names = FALSE, col.names=FALSE, sep = "\t")
}

for (Chr in AllChrs) {
	SampledSNPsVCF=sprintf(SampledSNPsVCFTemplate,Chr)
	if (file.exists(SampledSNPsVCF)) {
		next
	}
	VcfName=sprintf(GenoVcfTemplate,Chr)
	SampledSNPsFile=sprintf(SampledSNPsTemplate,Chr)
	AwkCmd="awk '{ OFS=\"\t\"; $3=\"Pos\"$2; print $0 }'"
	CommStr=sprintf("tabix -h %s chr%s:1-1 > %s; tabix %s -R %s | %s >> %s",VcfName,Chr,SampledSNPsVCF,VcfName,SampledSNPsFile,AwkCmd,SampledSNPsVCF)
	print(CommStr)
	system(CommStr)
}

for (Chr in AllChrs) {
	SampledSNPsFreq=sprintf(SampledSNPsFreqTemplate,Chr)
	if (file.exists(sprintf("%s.frq",SampledSNPsFreq))) {
		next
	}
	CommStr=sprintf("%s --vcf %s --maf 0.1 --real-ref-alleles --biallelic-only --snps-only  --remove children.txt --freq --out %s",plink,SampledSNPsVCF,SampledSNPsFreq)
	print(CommStr)
	system(CommStr)
}

for (Chr in AllChrs) {
	ControlsFile=sprintf(ControlSNPsTemplate,Chr)
	if (file.exists(ControlsFile)) {
		next
	}
	SampledSNPsFreq=sprintf(SampledSNPsFreqTemplate,Chr)
	Freqs=read.table(sprintf("%s.frq",SampledSNPsFreq),header=TRUE,)
	Res=Results[Results$Chr==Chr,]
	Controls=data.frame(matrix(nrow=nrow(Res),ncol=5))
	colnames(Controls)=colnames(Res)[1:5]
	Res$MAF=(Res$CountAA+Res$CountRA/2)/(Res$CountRR+Res$CountRA+Res$CountAA)
	for (r in 1:nrow(Res)) {
		F=Freqs[Freqs$MAF>=Res$MAF[r]-0.02 & Freqs$MAF<=Res$MAF[r]+0.02,]
		if (nrow(F)<1) {
			error(sprintf("No match found for MAF %f",Res$MAF[r]))
		}
		C=F[sample(nrow(F),1),]
		Pos=substring(C$SNP[1],4)
		Controls[r,]=c(Chr,Pos,Pos,0,0)
	}
	write.table(Controls,ControlsFile,quote=FALSE,row.names = FALSE, col.names=FALSE, sep = "\t")
	if (Chr=="1") {
		AllControls=Controls
	} else {
		AllControls=rbind(AllControls,Controls)
	}
}

write.table(AllControls,"allChrsControls.txt",quote=FALSE,row.names = FALSE, col.names=FALSE, sep = "\t")


