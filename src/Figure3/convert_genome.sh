#!/bin/bash

# Map BED files in current directory to hg18 from hg19

module load crossmap 

for i in `ls *.bed`; do
	suffix_format=(${i//./ })
	new_file_name=$suffix_format"_hg18.bed"
	crossmap bed ~/src/Figure3/hg19ToHg18.over.chain.gz $i $new_file_name
done

