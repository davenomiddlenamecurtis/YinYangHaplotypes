# YinYangHaplotypes
Code, scripts and results files to investigate yin yang haplotypes in 1000 Genomes data

Coordinates of variants forming yin yang haplotypes:
results/chr*.results.txt
Start and end coordinates of yin yang haplotypes:
results/chr*.summ.txt
Please note that these include the haplotypes which appear to represent repeat sequences. They are excluded by the script CombineYYfiles.20221107.R


Example commands:

plotGenos.20220804.R HighFst.20221115.pga
findYinYang --arg-file allchr12.fyyarg
getDepths --arg-file allDepths.tda
matchGenes  --arg-file matchGenes.20220727.mga

allChrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
for c in $allChrs ;do  ~/yinYang/getDists --arg-file ~/yinYang/anyChr.20220620.args --test-name distances.noSwap.20220727.$c --allow-swap 0 --chr $c; done
for c in $allChrs ; do  ~/yinYang/findYinYang --dist-file-name results.distances.noSwap.20220727.$c.$c.1.250000000.txt --chr $c --test-name noSwap.20220727.$c --min-run 20 --max-dist 0.0015 --separation 10; done
for c in $allChrs; do cat noSwap.20220727.$c.summ.txt >> allChrsSummary.20220727.txt; done
 bash  /home/rejudcu/yinYang/getYYHWE.sh /cluster/project9/bipolargenomes/yinYang/allChrsSummary.20220727.txt
 /share/apps/R-3.6.1/bin/Rscript ~/yinYang/getYinYangFst.20220727.R
 ~/yinYang/matchGenes --arg-file ~/yinYang/matchGenes.20220727.mga
 ~/yinYang/getDepths --arg-file ~/yinYang/allDepths.20220727.tda
 /share/apps/R-3.6.1/bin/Rscript ~/yinYang/CompareHet.20220727.R  ~/yinYang/testHet.20220727.arg
 bash ~/yinYang/annotWithVEP.20220727.sh
 /share/apps/R-3.6.1/bin/Rscript /home/rejudcu/yinYang/getControlsSNPs.R
 /share/apps/R-3.6.1/bin/Rscript ~/yinYang/getYinYangFst.20220727.R Controls

Note that some analyses depend on having first distinguished real YY haplotypes from repeats, which mean reading table or supplementary table:
GetAncAlls.20221101.R:AllResultsFile="/home/rejudcu/yinYang/AllYYResults.20221107.txt"
GetHapAncestry.20221108.R:AllResultsFile="/home/rejudcu/yinYang/AllYYResults.20221107.txt"
getControlsSNPs.R:ResultsFile="/home/rejudcu/yinYang/SupplementaryTable.20221101.txt"

