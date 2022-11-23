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
