#!/bin/bash

# scATAC-seq data acquisition 

# Download data from figshare 
wget https://ndownloader.figshare.com/files/12853220
mv 12853220 uATAC.B.bed
wget https://ndownloader.figshare.com/files/12853376
mv 12853376 uATAC.T.bed
wget https://ndownloader.figshare.com/files/12853373
mv 12853373 uATAC.CD4.bed
wget https://ndownloader.figshare.com/files/12853442
mv 12853442 uATAC.M.bed

# Split merged files 
awk -F"\t" '{print >"bedfiles/mono_"$8".bed"}' uATAC.M.bed
awk -F"\t" '{print >"bedfiles/cd4_"$8".bed"}' uATAC.CD4.bed
awk -F"\t" '{print >"bedfiles/t_cell_"$8".bed"}' uATAC.T.bed
awk -F"\t" '{print >"bedfiles/b_cell_"$8".bed"}' uATAC.B.bed

# Convert BED files from hg19 to hg18 
cd bedfiles
bash ~/src/Figure3/convert_bed.sh
bash ~/src/Figure3/convert_genome.sh

# Generate pseudo-bulk single cell BED files 
cd ..
cat bedfiles/mono_*_hg18.bed > mono.bed
cat bedfiles/b_cell_*_hg18.bed > b_cell.bed
cat bedfiles/t_cell_*_hg18.bed > t_cell.bed
cat bedfiles/cd4_*_hg18.bed >> t_cell.bed

# Downsample 
shuf -n 12000000 mono.bed > mono.sub.bed
shuf -n 2000000 t_cell.bed > t_cell.sub.bed
shuf -n 7000000 b_cell.bed > b_cell.sub.bed
