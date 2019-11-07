#!/bin/bash

# Get peaks from MACS2 callpeaks 

for i in `ls ../macs/*.xls`
do
	name=$(echo "$i" | cut -f 3 -d '/')
	peak_name=$(echo "$name" | cut -f 1 -d '.')".bed"
	{  awk '{print $1"\t"$2"\t"$3"\t"$4"\t42\t+"}' $i | tail -n +28; } > $peak_name
done
