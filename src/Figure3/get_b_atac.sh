#!/bin/bash

# bulk ATAC-seq data acquisition

cp /data/perriejv/final_draft_sc_dnase_workflow/b_atac/bedfiles/* bedfiles
cat bedfiles/*_Mono* > mono.bed
cat bedfiles/*_T* > t_cell.bed
cat bedfiles/*_B* > b_cell.bed 

# Downsample 
shuf -n 80000000 mono.bed > mono.sub.bed
shuf -n 550000000 t_cell.bed > t_cell.sub.bed
shuf -n 120000000 b_cell.bed > b_cell.sub.bed
