#!/bin/bash

# get ref and alt alleles at positions specified as list of chr positions

plink="/share/apps/genomics/plink-1.9/plink"
MAF=0.1

locFile=$1
rm getYYHWE.20220801.txt
cat $locFile | while read a b c
do
  tabix -h /cluster/project9/bipolargenomes/popGen/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr$a.recalibrated_variants.vcf.gz chr$a:$b-$b > getYYHWE.temp.$a.$b.$b.vcf
#  $plink --vcf temp.$a.$b.$b.vcf --real-ref-alleles --maf $MAF --biallelic-only --snps-only  --remove children.txt --make-bed --hardy
  $plink --vcf getYYHWE.temp.$a.$b.$b.vcf --real-ref-alleles --maf $MAF --biallelic-only --snps-only  --remove children.txt --hardy --out getYYHWE
  echo $a $b `tail -n 1 getYYHWE.hwe` >> getYYHWE.20220801.txt
  rm getYYHWE.temp.$a.$b.$b.vcf getYYHWE.hwe
done
