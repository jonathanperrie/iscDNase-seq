#!/bin/bash

# an individual MACS2 callpeak method call
# $1:file $2:label 

echo "macs2 callpeak -t $1 -f BED -g hs --name $2 --nomodel --call-summits --nolambda"
