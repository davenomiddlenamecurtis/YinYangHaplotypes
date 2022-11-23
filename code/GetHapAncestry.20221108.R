#!/share/apps/R-4.0.3/bin/Rscript

AllChr=c(1:22) # we do not have ancestral alleles for X chromosome
# AllChr=c(4)

TsvTemplate="YinYangAlls.chr%s.20221101.txt"
AncTemplate="Anc.chr%s.20221101.txt"
AncSummaryTemplate="/home/rejudcu/yinYang/Anc.chr%s.summ.20221101.txt"
AllResultsFile="/home/rejudcu/yinYang/AllYYResults.20221107.txt"
ValidRegions=read.table(AllResultsFile,header=TRUE,stringsAsFactors=FALSE,sep="\t",quote="")[c("Chr","Start","End")] # removed duplications etc.
FreqFile="/home/rejudcu/yinYang/allFsts.20220705.txt"
sink("/home/rejudcu/yinYang/GetHapAncestry.20221108.log")
for (Chr in AllChr) {
	TsvFile=sprintf(TsvTemplate,Chr)
	AncFile=sprintf(AncTemplate,Chr)
	Summary=ValidRegions[ValidRegions$Chr==Chr,]
	Alleles=read.table(TsvFile,header=TRUE,stringsAsFactors=FALSE)
	Alleles=Alleles[,-(9:10)] # remove Den
	Alleles[Alleles==0]=NA
	Alleles=na.omit(Alleles)
	Alleles$BothAlt=(Alleles$Chimp1==Alleles$ALT & Alleles$Altai1==Alleles$ALT)*1
	Alleles$BothRef=(Alleles$Chimp1!=Alleles$ALT & Alleles$Altai1!=Alleles$ALT)*1
	Alleles$ChimpAlt=(Alleles$Chimp1==Alleles$ALT & Alleles$Altai1!=Alleles$ALT)*1
	Alleles$NeandAlt=(Alleles$Chimp1!=Alleles$ALT & Alleles$Altai1==Alleles$ALT)*1
	write.table(Alleles,AncFile,quote=FALSE,row.names = FALSE, sep = "\t")
	AncSummary=data.frame(matrix(nrow=nrow(Summary),ncol=12))
	colnames(AncSummary)=c("Chr","Start","End","BothAlt","BothRef","ChimpAlt","NeandAlt","Derived","Ancient","Mixed","Partial","MixedPartial")
	AncSummary[,1:3]=Summary[,1:3]
	for (r in 1:nrow(AncSummary)) {
		ThisHap=Alleles[Alleles$Pos>=AncSummary$Start[r] & Alleles$Pos<=AncSummary$End[r],]
		AncSummary[r,4:7]=colSums(ThisHap[,9:12])
		AncSummary$Derived[r]=(((AncSummary$BothAlt[r]!=0 & AncSummary$BothRef[r]==0) | (AncSummary$BothAlt[r]==0 & AncSummary$BothRef[r]!=0)) & AncSummary$ChimpAlt[r]==0 & AncSummary$NeandAlt[r]==0)*1
		AncSummary$Ancient[r]=((AncSummary$BothAlt[r]==0 & AncSummary$BothRef[r]==0) & ((AncSummary$ChimpAlt[r]!=0 & AncSummary$NeandAlt[r]==0)|(AncSummary$ChimpAlt[r]==0 & AncSummary$NeandAlt[r]!=0)))*1
		AncSummary$Mixed[r]=(AncSummary$BothAlt[r]!=0 & AncSummary$BothRef[r]!=0)*1
		AncSummary$Partial[r]=(AncSummary$Ancient[r]==0 & (AncSummary$ChimpAlt[r]!=0 | AncSummary$NeandAlt[r]!=0))*1
		AncSummary$MixedPartial[r]=(AncSummary$Ancient[r]==0 & (AncSummary$ChimpAlt[r]!=0 & AncSummary$NeandAlt[r]!=0))*1
	}
	write.table(AncSummary,sprintf(AncSummaryTemplate,Chr),quote=FALSE,row.names = FALSE, sep = "\t") # will write combined table later
	if (Chr=="1") {
		AllSummary=AncSummary
	} else {
		AllSummary=rbind(AllSummary,AncSummary)
	}
}
write.table(AllSummary,sprintf(AncSummaryTemplate,"All"),quote=FALSE,row.names = FALSE, sep = "\t")

AllSummary=read.table(sprintf(AncSummaryTemplate,"All"),header=TRUE,stringsAsFactors=FALSE)
AllSummary=AllSummary[rowSums(AllSummary[,4:7])>=5,]
sprintf("Number of haplotypes with at least 5 variants informative for ancestry = %d",nrow(AllSummary))
# write.table(AllSummary,sprintf(AncSummaryTemplate,"All"),quote=FALSE,row.names = FALSE, sep = "\t")
# makes sense to only keep these maybe not
sprintf("Number of haplotypes possibly derived in humans = %d",nrow(AllSummary[AllSummary$Derived==1,]))
sprintf("Number of haplotypes possibly ancestral = %d",nrow(AllSummary[AllSummary$Ancient==1,]))
AllSummary$Length=AllSummary$End+1-AllSummary$Start
sprintf("Compare length of non-derived and derived haplotypes:")
t.test(AllSummary[AllSummary$Derived==0,]$Length,AllSummary[AllSummary$Derived==1,]$Length)
sprintf("Compare length of non-ancient and ancient haplotypes:")
t.test(AllSummary[AllSummary$Ancient==0,]$Length,AllSummary[AllSummary$Ancient==1,]$Length)
Ancient=AllSummary[AllSummary$Ancient==1,]
Freqs=read.table(FreqFile,header=TRUE,stringsAsFactors=FALSE)
Ancient$Coord=sprintf("%s:%s",Ancient$Chr,Ancient$Start)
Freqs$Coord=sprintf("%s:%s",Freqs$Chr,Freqs$Pos)
Ancient=merge(Ancient,Freqs,by="Coord")
Ancient=Ancient[c("Coord","Length","ChimpAlt","NeandAlt","AFR")]
Ancient$AFR=sprintf("%.3f",as.numeric(Ancient$AFR))
Ancient





