#!/bin/bash

# Get blacklist peaks and map to hg18 from hg19

module load crossmap

curl "https://docs.google.com/spreadsheets/d/e/2PACX-1vSHJvlVJ_YpcSPE8mAiQNw-EzsO9fMfu_c1IgaUBJ_lrf819nbxr1VgtLidQxi17ODo_fb-Ryjjw3yW/pub?output=tsv" > hg19.blacklist.bed

crossmap bed ~/src/Figure3/hg19ToHg18.over.chain.gz hg19.blacklist.bed hg18.blacklist.bed

echo "$(tail -n +1 hg18.blacklist.bed)" > hg18.blacklist.bed
awk '{print $1"\t"$2"\t"$3"\t"$3-$2"\t42\t+"}' hg18.blacklist.bed > tmp.bed
cat tmp.bed > hg18.blacklist.bed
rm tmp.bed

