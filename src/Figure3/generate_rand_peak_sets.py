import numpy as np
import pandas as pd
from pybedtools import BedTool
import random

def calc_gc(peaks,hg18_path):
	"""calculate the GC content in peak set. 
	
	Parameters
	----------
	peaks : peak set
	hg18_path : genome path 

	Returns
	-------
	np.mean(gc) : average GC content over all peaks	
	"""
	peaks=BedTool(peaks).nucleotide_content(hg18_path)
	gc=[]
	for peak in peaks:
		gc+=[float(peak[7])]
	return(np.mean(gc))

def is_valid(cand,bp_df):
	"""check if there is overlap between some candidate peak and the peaks in the blacklist regions.
	
	Parameters
	----------
	cand : candidate peak
	bp_df : blacklist regions as a data frame

	Returns 
	-------
	Bool : Boolean describing if match was valid
	"""
	pos_matches=bp_df[bp_df.chr==cand[0]] 

	if len(pos_matches[(pos_matches.start<cand[1]) & (pos_matches.end>cand[1])]):
		return False
	elif len(pos_matches[(pos_matches.start<cand[2]) & (pos_matches.end>cand[2])]):
		return False
	elif len(pos_matches[(pos_matches.start>cand[1]) & (pos_matches.start<cand[2])]):
		return False
	elif len(pos_matches[(pos_matches.start<cand[1]) & (pos_matches.start>cand[2])]):
		return False
	else:
		return True


def find_new_peaks(peaks,hg18_path,bpath,gl_path,width,threshold):
	"""find a new set of peaks with equal proportions of peaks at each chromosome.
	
	Parameters
	----------
	peaks : reference peak set 
	bpath : blacklist path 
	gl_path : genome length path 

	Returns
	-------
	pseudo_rand_peaks : new set of random peaks 
	"""

	# all proposed peaks should have a start point a standard peak-width away from the end of the chromosome 
	genome_len=read_peaks(gl_path)
	genome_len={x[0]:int(x[1])-width for x in genome_len}

	blacklist_peaks=read_peaks(bpath)
	bp_df=pd.DataFrame(blacklist_peaks)
	bp_df.rename(index=str,columns={0:"chr",1:"start",2:"end"},inplace=True)
	bp_df.start=bp_df.start.astype(int)
	bp_df.end=bp_df.end.astype(int)

	rand_peaks=[]
	
	# generate initial random peak set that has same chromosome proportions 
	for peak in peaks:
		chrom=peak[0]
		start=random.randint(0,genome_len[chrom])
		end=start+width
		while not is_valid((chrom,start,end),bp_df):
			chrom=peak[0]
			start=random.randint(0,genome_len[chrom])
			end=start+width
		rand_peaks+=[[chrom,start,end,'.','.','*']]

	# adjust random peak set so that it's gc-content is close to true peaks threshold 
	true_gc=calc_gc(peaks,hg18_path)
	rand_gc=calc_gc(rand_peaks,hg18_path)

	i=0
	while(np.abs(true_gc-rand_gc)>threshold):
		# reset index when hitting the end
		i=np.random.choice(len(peaks))

		peak=rand_peaks[i]
		rand_gc_peak=float(BedTool([rand_peaks[i]]).nucleotide_content(hg18_path)[0][7])
		chrom=peak[0]

		# evaluate a new candidate peak 
		start=random.randint(0,genome_len[chrom])
		end=start+width
		while not is_valid((chrom,start,end),bp_df):
			chrom=peak[0]
			start=random.randint(0,genome_len[chrom])
			end=start+width
		cand_gc_peak=float(BedTool([[chrom,start,end]]).nucleotide_content(hg18_path)[0][4])

		# if change in GC content will make the random peak set more closely related to the reference peak 
		# set, then go ahead with the update 
		if true_gc>=rand_gc and cand_gc_peak>rand_gc_peak:
			rand_peaks[i]=[chrom,start,end,'.','.','*']
			rand_gc=rand_gc-rand_gc_peak/len(peaks)+cand_gc_peak/len(peaks)
		elif true_gc<rand_gc and cand_gc_peak<rand_gc_peak:
			rand_peaks[i]=[chrom,start,end,'.','.','*']
			rand_gc=rand_gc-rand_gc_peak/len(peaks)+cand_gc_peak/len(peaks)
		else:
			continue
	return rand_peaks

def read_peaks(path):
	"""read in list.
	
	Parameters
	----------
	path : list path
	
	Returns
	-------
	a : list with elements on each line separated
	"""
	with open(path,"rt") as f:
		a=f.readlines()
	a=[x.split() for x in a]
	return a

def find_rand_peaks(peak_path,other_args):
	"""downsample peaks and find a random set of peaks with comparable chromosome proportions and GC content.
	
	Parameters
	---------- 
	peak_path : peak set path
	other_args :  other paths and program specs
		hg18_path : genome path
		bpath : blacklist path
		gl_path : genome length path
		width : 500
		threshold : GC-content proportion 	

	Returns
	-------
	dpeaks : downsampled original peak set 
	rand_peaks : random peak set with similar proportions to downsampled peak set 
	"""
	peaks=read_peaks(peak_path)
	peaks=[x for x in peaks if x[0]!='chrM']

	hg18_path=other_args[0]
	bpath=other_args[1]
	gl_path=other_args[2]
	width=other_args[3]
	threshold=other_args[4]
	 
	rand_peaks=find_new_peaks(peaks,hg18_path,bpath,gl_path,width,threshold)
	print(peak_path)
	return(rand_peaks)
