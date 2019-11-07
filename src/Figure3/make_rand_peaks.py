from multiprocessing import Pool
import signal
from itertools import product
import glob
import os
import sys

sys.path.append('../../src/Figure3/')
from generate_rand_peak_sets import * 

peak_paths=glob.glob("peaks/*specific*.bed")+glob.glob("peaks/inter_peaks/*specific*.bed")
bpath="peaks/hg18.blacklist.bed"
hg18_path="/fdb/indexes/hg18/hg18.fa"
# outside scope of program
gl_path="/data/perriejv/genome_len/hg18_chrlen.txt"
width=500
threshold=0.05

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", "2"))
p = Pool(nproc, init_worker)

try:
	rand_peaks=p.starmap(find_rand_peaks,product(peak_paths,[(hg18_path,bpath,gl_path,width,threshold)]))
except (KeyboardInterrupt,SystemExit):
	p.terminate()
	p.join()
	sys.exit(1)
else:
	p.close()
	p.join()

def write_peaks(path,peaks):
	"""read in list.
	
	Parameters
	----------
	path : list path
	
	Returns
	-------
	a : list with elements on each line separated
	"""
	if "inter" in path:
		filename="peaks/rand_peaks/random_inter_"+path.split('/')[-1]
	else:
		filename="peaks/rand_peaks/random_"+path.split('/')[-1]
	with open(filename,"w") as f:
		for p in peaks:
			# add tabs
			p[1]=str(p[1])
			p[2]=str(p[2])
			f.write("\t".join(p)+"\n")


[write_peaks(peak_paths[i],rand_peaks[i]) for i in range(len(rand_peaks))]
