#!/bin/bash

# VEP annotations for lead SNP of each haplotype

plink="/share/apps/genomics/plink-1.9/plink"
MAF=0.05

locFile=allChrsSummary.20220727.txt

if [ ! -e forVEP.vars.20220727.txt ]
then
cat $locFile | while read a b c
do
  tabix -h /cluster/project9/bipolargenomes/popGen/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr$a.recalibrated_variants.vcf.gz chr$a:$b-$b > temp.$a.$b.$b.vcf
  rm plink.bim
  $plink --vcf temp.$a.$b.$b.vcf --real-ref-alleles --maf $MAF --biallelic-only --snps-only  --remove children.txt --make-bed 
  cat plink.bim | awk '{ if ($1!="23") c=$1; else c="X" ; print c,$4,$4,$6"/"$5,"+" }' >> forVEP.vars.20220727.txt
  rm temp.$a.$b.$b.vcf
done
fi

if [ ! -e VEP.annots.20220727.txt ]
then
perl /share/apps/ensembl-vep-97/vep --input_file forVEP.vars.20220727.txt \
  --synonyms ~/vep/chr_synonyms.txt \
  --cache --dir /cluster/project9/bipolargenomes/vepcache --merged --force_overwrite \
  --assembly GRCh38 \
  --fasta /cluster/project9/bipolargenomes/vepcache/homo_sapiens_merged/97_GRCh38 \
  --check_existing --pick \
  --output_file VEP.annots.20220727.txt
fi

tail -n +51 VEP.annots.20220727.txt | cut -f 13 | grep rs > YY.SNPs.20220727.txt
grep -f YY.SNPs.20220727.txt gwas_catalog_v1.0-associations_e106_r2022-07-09.tsv > YY.GWAScat.20220727.tsv


