#!/bin/bash 

# Adjust BED file to conform to BED6 format 

for i in `ls`; do
	awk '{print $1"\t"$2"\t"$3"\t"$3-$2"\t42\t+"}' $i > tmp.bed
	cat tmp.bed > $i
done

rm tmp.bed
