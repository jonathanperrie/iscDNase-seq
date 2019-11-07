#!/bin/bash

echo "#!/bin/bash" > anno_peaks/anno.swarm

declare -a peak_path=("NA" "rand_peaks" "inter_peaks")

for i in "${peak_path[@]}"
do
	if [ "$i" == "NA" ]; then
		for j in *specific*.bed
		do 
			echo "annotatePeaks.pl ~/data/Figure3/peaks/"$j" hg18 -gene "$1" > "$j".anno" >> anno_peaks/anno.swarm
		done
	else 
        	for j in $i/*specific*.bed
        	do
			name=`echo $j | sed "s/\//_/g"`
        		echo "annotatePeaks.pl ~/data/Figure3/peaks/"$j" hg18 -gene "$1" > "$name".anno" >> anno_peaks/anno.swarm
		done
	fi
done
