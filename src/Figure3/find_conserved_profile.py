import glob
from multiprocessing import Pool
import signal
import os
import sys
import fileinput
import numpy as np
from itertools import product
import random
import pandas as pd
from pybedtools import BedTool
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# initial worker
def init_worker():
        signal.signal(signal.SIGINT, signal.SIG_IGN)




def read_chr(fpath):
	"""reading in phastCons conservation genome and binning average conservation scores over at most 100 bp windows.
	
	Parameters
	----------
	fpath : path of chromosome with conservation scores to read in 

	Returns
	-------
	chrom_dicts : average conservation score over at most 100 bp windows in a chromosome 
	"""
	# init dict and indices
	chrom_dicts={}
	start=0
	index=0

	# iterate through chromosome scores 
	for line in fileinput.input(fpath):
		x=line.split()
		
		# if chromosome skips some region, then normalize the previous window (<100 bp) and init new window 	
		if len(x)==4:
			if start in chrom_dicts:
				chrom_dicts[start]/=index
			start=int(x[2].split("=")[1])
			chrom_dicts[start]=0
			index=0

		# if not a black region, then make news windows every 100 locations
		if len(x)==1:
			chrom_dicts[start]+=float(x[0])
			if index==100:
				chrom_dicts[start]/=index
				index=0
				start+=100
				chrom_dicts[start]=0
			index+=1
	
	# track chromosomes that have been binned
	print("%s %d" % (fpath,len(chrom_dicts)))
	return(chrom_dicts)


def find_conserved_profile(peak_key,dicts):
	"""find conservation profile for peak set. 
	
	Parameters
	----------
	peak_key : a specific peak set reference 
	dicts : dict objects for chromosome conservation scores and peak sets
	
	Returns
	-------
	profile : average profile for each peak with some window on both sides
	"""

	# read in peaks
	chrom_dict=dicts[0]
	peak_dict=dicts[1]
	peaks=peak_dict[peak_key]

	# profile around peak +/- 20 windows
	profile=np.zeros(41)

	# iterate through each peak and find region in conservation genome
	for x in peaks:
		# get peak chromosome
		chrom=x[0]
		# get peak mid point
		mid=int((int(x[1])+int(x[2]))/2)
		# find peak chromosome  conservation indices
		keys=np.fromiter(chrom_dict[chrom].keys(),dtype=int)
		# find where peak is located
		matches=np.where(mid>keys)[0]
		if len(matches)==0:
			matches=np.array([0])
		key_idx=max(matches)

		# if the match is near the start or end of the chromosome, adjust the profile index
		index_range=np.arange(max(key_idx-20,0),min(key_idx+21,len(keys)))
		if 0 in index_range:
			adj_idx=20-key_idx
		else:
			adj_idx=0
		
		# build profile
		for idx,i in enumerate(index_range):
			profile[idx+adj_idx]+=chrom_dict[chrom][keys[i]]
	profile/=len(peaks)
	return profile

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

if __name__ == "__main__":
	hg18_path="/fdb/indexes/hg18/hg18.fa"
	bpath="/data/perriejv/v3_idnase/peaks/hg18.blacklist.bed"
	gl_path="/data/perriejv/genome_len/hg18_chrlen.txt"

	# bin chromosome conservation scores into 100 bp windows
	nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))
	p = Pool(nproc, init_worker)
	chrom_paths=glob.glob("phastcons/chr*") 

	try:
		results=p.map(read_chr,chrom_paths)
	except (KeyboardInterrupt,SystemExit):
		p.terminate()
		p.join()
		sys.exit(1)
	else:
		p.close()
		p.join()

	chrom_dict={}
	for i in range(len(results)):
		chrom_dict[chrom_paths[i].split("/")[1]]=results[i]

	peak_paths=glob.glob("peaks/*specific*bed")
	peak_paths+=glob.glob("peaks/inter_peaks/*specific*bed")
	peak_paths+=glob.glob("peaks/rand_peaks/*specific*bed")

	nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))
	p = Pool(nproc, init_worker)

	try:
		peak_results=p.map(read_peaks,peak_paths)
	except (KeyboardInterrupt,SystemExit):
		p.terminate()
		p.join()
		sys.exit(1)
	else:
		p.close()
		p.join()


	peak_dict={}
	for path in peak_paths:
		peak_dict[path.split('/')[-1].split('.bed')[0]+"_"+path.split('/')[-2].split('.bed')[0]]=peak_results[peak_paths.index(path)]

	# match peaks with genome 
	nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))
	p = Pool(nproc, init_worker)

	try:
		profile_results=p.starmap(find_conserved_profile,product(peak_dict.keys(),[(chrom_dict,peak_dict)]))
	except (KeyboardInterrupt,SystemExit):
		p.terminate()
		p.join()
		sys.exit(1)
	else:
		p.close()
		p.join()
	
		# visualize peaks with ATAC and DNase peak sets on same plot for unique and bulk intersecting peaks 
		filenames=list(peak_dict.keys())

		plt.rcParams.update({'font.size': 14})
		for cell_type in ["mono","t_cell","b_cell"]:
			fig = plt.figure(figsize=(4,6))
			ax = fig.add_subplot(111)
			ax.spines['top'].set_color('none')
			ax.spines['bottom'].set_color('none')
			ax.spines['left'].set_color('none')
			ax.spines['right'].set_color('none')
			ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
			count=1
			for assay in ["dnase","atac"]:
#				for portion in [False,True]:
				portion=False
				ax1=fig.add_subplot("21"+str(count))
				for peak in [False,True]:
					index=[x for x in filenames if cell_type in x and assay in x and ("inter" in x)==portion and ("rand" in x)==peak][0]
					y=profile_results[filenames.index(index)]
					if peak==True:
						color="#c39999"
						if assay=="atac" and portion==True:
							ax1.set_xlabel({False:"unique",True:"specific"}[portion])
					elif assay=="atac":
						color="skyblue"
#						if portion==False:
#							ax1.set_ylabel({"dnase":"DNase","atac":"ATAC"}[assay])
#							ax1.set_xlabel({False:"unique",True:"specific"}[portion])
					elif assay=="dnase":
						color="tomato"
						if portion==False:
							pass
#							ax1.set_ylabel({"dnase":"DNase","atac":"ATAC"}[assay])
#						else:
					if assay=="atac":
							dnase_patch = mpatches.Patch(color='tomato', label='iscDNase')
							atac_patch = mpatches.Patch(color='skyblue', label='scATAC')
							intersect_patch = mpatches.Patch(color='#c39999', label='random')
							plt.legend(handles=[dnase_patch,atac_patch,intersect_patch], frameon=False,handlelength=0.7,
								loc='upper center',fontsize=10)
					ax1.plot(y,color=color)
					ax1.set_ylim([0,0.15])
					ax1.set_xticklabels(np.arange(-30,21,10))
					ax1.spines['right'].set_visible(False)
					ax1.spines['top'].set_visible(False)
				print("22"+str(count))
				count+=1
			ax.set_xlabel("Relative peak position (kbp)",labelpad=20)
			ax.set_ylabel("Conservation score (phastCons)",labelpad=25)
			plt.tight_layout()
			plt.savefig("../../Figures/Figure3/"+cell_type+"_cons.png")

