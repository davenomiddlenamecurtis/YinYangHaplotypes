#!/bin/bash

# get ref and alt alleles at positions specified as list of chr positions

plink="/share/apps/genomics/plink-1.9/plink"
MAF=0.1

locFile=$1
firstLine=`head -n 1 $locFile`
lastLine=`tail -n 1 $locFile`
first=($firstLine)
last=($lastLine)
tabix -h /cluster/project9/bipolargenomes/popGen/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${first[0]}.recalibrated_variants.vcf.gz chr${first[0]}:${first[1]}-${last[1]} > temp.${first[0]}.${first[1]}.${last[1]}.vcf
rm locList.txt
cat $locFile | while read a b ; do echo $a $b $b extract>>locList.txt; done # use read so no worry about spaces or tabs
$plink --vcf temp.${first[0]}.${first[1]}.${last[1]}.vcf --real-ref-alleles --maf $MAF --biallelic-only --snps-only  --remove children.txt --make-bed --out temp.${first[0]}.${first[1]}.${last[1]} --extract range locList.txt



