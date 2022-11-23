#!/bin/bash

c=$1
cd /scratch0
mkdir $c
cd $c
hostname
pwd
cp /cluster/project9/bipolargenomes/yinYang/children.txt .
~/yinYang/getDists --arg-file ~/yinYang/anyChr.20220620.args --test-name distances.noSwap.20220727.$c --allow-swap 0 --chr $c
ls -lrt
cp res*txt /cluster/project9/bipolargenomes/yinYang
cd ..
rm -r $c
