#!/bin/bash

# iDNase-seq data acquisition 

cp /data/perriejv/idnase_post_anal/bedfiles/* bedfiles 
cat bedfiles/mono_* > mono.bed
cat bedfiles/B_* > b_cell.bed
cat bedfiles/T_* > t_cell.bed

# Downsample 
shuf -n 12000000 mono.bed > mono.sub.bed
shuf -n 2000000 t_cell.bed > t_cell.sub.bed
shuf -n 7000000 b_cell.bed > b_cell.sub.bed
