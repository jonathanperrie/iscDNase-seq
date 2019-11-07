#!/bin/bash

# $1:path $2:expression path 

echo "#!/bin/bash" > anno.swarm


declare -a rand_type=("" "random_")
declare -a peak_type=("" "inter")
declare -a cell_type=("mono" "t_cell" "b_cell")
declare -a assay_type=("b_atac" "b_dnase")

for i in "${peak_type[@]}"
do
        for j in "${cell_type[@]}"
        do
                for k in "${assay_type[@]}"
                do
                        for l in "${rand_type[@]}"
                        do
                                if [ "$l" == "random_" ]; then 
                                        peak_path="rand_peaks"
                                        if [ "$i" == "inter" ]; then
                                                echo "annotatePeaks.pl "$1"/peaks/"$peak_path"/"$l$i"_"$j"_specific_"$k"_peaks.bed hg18 -gene "$2" > "$l$j"_"$i"_"$k"_anno.txt" >> anno.swarm
                                        else 
                                                echo "annotatePeaks.pl "$1"/peaks/"$peak_path"/"$l$j"_specific_"$k"_peaks.bed hg18 -gene "$2" > "$l$j"_"$i"_"$k"_anno.txt" >> anno.swarm
                                        fi                                        
                                elif [ "$i" == "inter" ]; then
                                        peak_path="inter_peaks"
                                        echo "annotatePeaks.pl "$1"/peaks/"$peak_path"/"$l$j"_specific_"$k"_peaks.bed hg18 -gene "$2" > "$l$j"_"$i"_"$k"_anno.txt" >> anno.swarm
                                else
                                        peak_path=""
                                        echo "annotatePeaks.pl "$1"/peaks/"$peak_path"/"$l$j"_specific_"$k"_peaks.bed hg18 -gene "$2" > "$l$j"_"$i"_"$k"_anno.txt" >> anno.swarm
				fi
                        done 
                done
        done
done
