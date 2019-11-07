#!/bin/bash

module load bedops
module load crossmap

for i in *.bigWig; do
	~/src/Figure3/bigWigToWig $i $i.wig ; wig2bed --zero-indexed < $i.wig > $i.bed ; crossmap bed ~/src/Figure3/hg19ToHg18.over.chain.gz $i.bed $i.hg18.bed
done
