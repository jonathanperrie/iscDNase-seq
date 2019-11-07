#!/bin/bash

# Generate swarm file for peak set 
# $1:downsample option $2:file name 

# downsample option 
if [ $1 == 1 ]; then
	declare -a assay_type=("b_atac" "b_dnase" "sc_atac" "i_dnase")
	declare -a cell_type=("mono.sub.bed" "t_cell.sub.bed" "b_cell.sub.bed")
	suffix="sub"
else
	declare -a assay_type=("b_atac" "b_dnase" "sc_atac" "i_dnase")
	declare -a cell_type=("mono.bed" "t_cell.bed" "b_cell.bed")
	suffix="full"
fi

echo "#!/bin/bash" > $2

# iterate through all assays and cell types to generate a specific swarm call
for i in "${assay_type[@]}"
do
	for j in "${cell_type[@]}"
	do
		name=$(echo "$j" | cut -f 1 -d '.')
		echo `bash ~/src/Figure3/macs_format.sh ~/data/Figure3/$i/$j $i"_"$name"_"$suffix` >> $2
	done
done
