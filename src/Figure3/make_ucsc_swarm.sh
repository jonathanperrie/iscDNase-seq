#!/bin/bash 
# $1:genome length  $2:file name 


declare -a assay_type=("sc_atac" "i_dnase")
declare -a cell_type=("mono.sub.bed" "t_cell.sub.bed" "b_cell.sub.bed")

echo "#!/bin/bash" > $2

for i in "${assay_type[@]}";
	do
	for j in "${cell_type[@]}";
	do
		name=$(echo "$j" | cut -f 1 -d '.')
		echo 'cd ~/data/Figure3/'$i'; ~/src/Figure3/generateRPBMBasedSummary '$j' '$1' 200 75 1 y '$name'.bedgraph; ~/src/Figure3/setgraph4ucsc '$name'.'$i' 1.1 '$name'.bedgraph '$name'.bedgraph.ucsc'>>$2
	done
done


