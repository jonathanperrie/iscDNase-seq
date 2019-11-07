# iscDNase-seq
Workflow for iscDNase-seq analysis

GEO: GSE136987 (private) 

## Figure 3 plot generation workflow
Assumptions:
* Files are locally accessible 
* Access to Biowolf high performance computing resources
Both of these assumptions will have to be corrected for at some point with analysis conducted starting with GEO over less sophisticated servers. 

### Building tools
``` 
HOME=`pwd`

cd ~/src/Figure3
g++ -o generateRPBMBasedSummary OftenUsedOperatLib.cpp bed.cpp ucsc.cpp TypeDefBase.h Matrix.cpp SequenceTransform.cpp expression.cpp UCSCBEDOperator.cpp generateRPBMBasedSummary_v2.cpp 

if [[ ! -f bigWigToWig ]]; then
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
fi

if [[ ! -f hg19ToHg18.over.chain.gz ]]; then
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz
fi

if [[ ! -f rpkmforgenes.py ]]; then
	wget http://sandberg.cmb.ki.se/media/data/rnaseq/rpkmforgenes.py
fi

chmod 777 generateRPBMBasedSummary bigWigToWig hg19ToHg18.over.chain.gz rpkmforgenes.py
```
### Acquiring data 
```
cd ~/data/Figure3
mkdir -p sc_atac i_dnase b_atac b_dnase cluster
mkdir -p sc_atac/bedfiles i_dnase/bedfiles b_atac/bedfiles b_dnase/bedfiles

cd sc_atac
bash ~/src/Figure3/get_sc_atac.sh 
 
cd ../i_dnase
bash ~/src/Figure3/get_i_dnase.sh 

cd ../b_atac
bash ~/src/Figure3/get_b_atac.sh 

cd ../b_dnase
bash ~/src/Figure3/get_b_dnase.sh
```
### Finding peaks
```
cd ..
mkdir -p macs peaks
cd macs 
bash ~/src/Figure3/make_macs_swarm.sh 0 macs.swarm
swarm -f macs.swarm -g 32 --module macs --gres=lscratch:64

bash ~/src/Figure3/make_macs_swarm.sh 1 macs.sub.swarm
swarm -f macs.sub.swarm -g 32 --module macs --gres=lscratch:64

cd ../peaks 
bash ~/src/Figure3/get_peaks.sh 
bash ~/src/Figure3/get_blacklist.sh 
bash ~/src/Figure3/merge_peaks.sh
```
##### 3d
```
bash ~/src/Figure3/make_ucsc_swarm.sh /data/perriejv/genome_len/hg18_chrlen.txt ucsc.swarm
swarm -f ucsc.swarm -g 32 --gres=lscratch:64

mkdir -p inter_peaks rand_peaks anno_peaks

cd ..
```
### Viz 
```
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate base
conda activate general
```
##### 3a/b
```
Rscript ~/src/Figure3/assay_similarties.R
```
##### 3c
```
Rscript ~/src/Figure3/find_uniq_inter_peaks.R
```
##### 3e
```
# navigate to http://hgdownload.cse.ucsc.edu/goldenpath/hg18/phastCons17way/ and make a links to each chromosome 
mkdir -p phastcons rna_seq 
wget -i filepaths
gunzip *.gz
module load bedtools
python ~/src/Figure3/make_rand_peaks.py
python ~/src/Figure3/find_conserved_profile.py
```
##### 3f
```
# add bigWig signal of all reads for the following cell types: 
# Monocyte: https://www.encodeproject.org/experiments/ENCSR905LVO/, ENCFF706XVQ
# B cell: https://www.encodeproject.org/experiments/ENCSR449GLL/, ENCFF517AYD
# T cell: https://www.encodeproject.org/experiments/ENCSR545MEZ/, ENCFF947QSD
# Mar. 2006 (NCBI36/hg18), Other RefSeq, xenoRefGene
cd ../rna_seq
wget -i filepaths
module load python
python ~/src/Figure3/rpkmforgenes.py -i ENCFF517AYD.bigWig.hg18.bed  ENCFF706XVQ.bigWig.hg18.bed  ENCFF947QSD.bigWig.hg18.bed -readcount -a hgTables.txt -o rna_seq_gene_expression.txt -fulltranscript -rmnameoverlap -p 20

cd ../peaks
bash ~/src/Figure3/make_anno.sh ~/data/Figure3/rna_seq/rna_seq_gene_expression.txt
cd anno_peaks
swarm -f anno.swarm -g 16 -t 2 --module homer
python ~/src/Figure3/gene_expr.py
```
