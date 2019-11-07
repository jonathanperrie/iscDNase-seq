import matplotlib
matplotlib.use('agg')

from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import ranksums
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu

# returns list of genes that have expression above 0 for cell type within 10,000 bp from tss
# find gene expression for individual annotated peak set 
# x[-3:-1] refer to gene expression for B cells, Monocytes, and T cells
# x[15] refers to the gene's name
# x[9] refers to TSS distance  
def find_cands(file,offset):
	if type(file) is str:
		with open(file,"rt") as f:
			genes=f.readlines()
		genes=[x.split("\t") for x in genes if len(x.split("\t"))==27]
		genes=[x for x in genes if int(x[-offset])>0]
		return {x[15]:(int(x[-3]),int(x[-2]),int(x[-1])) for x in genes if np.abs(int(x[9].split()[0]))<10000}

	# find genes with expression above 0 in multiple annotated peak sets 
	elif type(file) is list:
		gene_list=[]
		for idx,file in enumerate(file):
			with open(file,"rt") as f:
				genes=f.readlines()
			genes=[x.split("\t") for x in genes if len(x.split("\t"))==27]
			genes=[x for x in genes if int(x[-offset])>0]
			genes=[(x[9],x[15],-1) for x in genes]
			gene_list+=genes
		return {x[1]:int(x[2]) for x in gene_list if np.abs(int(x[0].split()[0]))<10000}

filepaths=glob("*.anno")

rand=False
inter=False
cell_type=""
assay=""

expr_dict={}
offset_dict={"t_cell":1,"b_cell":3,"mono":2}

for file in filepaths:
	if "rand" in file:
		rand=True
	else:
		rand=False
	if "inter" in file:
		inter=True
	else:
		inter=False
	if "t_cell" in file:
		cell_type="t_cell"
	if "b_cell" in file:
		cell_type="b_cell"
	if "mono" in file:
		cell_type="mono"
	if "dnase" in file:
		assay="dnase"
	if "atac" in file:
		assay="atac"
	offset=offset_dict[cell_type]

	if rand:
		cand_name=file
		comp_name=[x for x in filepaths if "rand" not in x and cell_type in x]
		cand=find_cands(cand_name,offset)
		comp=find_cands(comp_name,offset)
		filtered_cand={x:cand[x] for x in cand if x not in comp.keys()}
	else:
		cand_name=file
		comp_name=[x for x in filepaths if "rand" not in x and cell_type in x and assay not in x]
		cand=find_cands(cand_name,offset)
		comp=find_cands(comp_name,offset)
		filtered_cand={x:cand[x] for x in cand if x not in comp.keys()}
	if cell_type not in expr_dict:
		expr_dict[cell_type]={}
	if assay not in expr_dict[cell_type]:
		expr_dict[cell_type][assay]={}
	if (inter,rand) not in expr_dict[cell_type][assay]:
		expr_dict[cell_type][assay][(inter,rand)]=filtered_cand

assay_dict={"dnase":0,"atac":1}
assay_col={"dnase":'tomato',"atac":'skyblue'}
title_dict={"t_cell":"T cell","b_cell":"B cell","mono":"Monocyte"}

plt.rcParams.update({'font.size': 14})

for i in expr_dict:
	fig = plt.figure(figsize=(4, 6))
	ax = fig.add_subplot(111)
	ax.spines['top'].set_color('none')
	ax.spines['bottom'].set_color('none')
	ax.spines['left'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
	count=1
	for j in ["dnase","atac"]:
#		for k in [False,True]:
		ax1=fig.add_subplot(2,1,count)
		assay=[x[3-offset_dict[i]] for x in list(expr_dict[i][j][(False,False)].values())]
		rand_assay=[x[3-offset_dict[i]] for x in list(expr_dict[i][j][(False,True)].values())]

		assay=[x for x in assay if x >100]
		rand_assay=[x for x in rand_assay if x >100]

		bp_dict={1:[x for x in assay],2:[x for x in rand_assay]}


		bp=ax1.boxplot(bp_dict.values(),showfliers=False)
		ax1.plot(np.random.normal(1,0.02,size=len(assay)),[x for x in assay],'.',color=assay_col[j],alpha=0.2)
		ax1.plot(np.random.normal(2,0.02,size=len(rand_assay)),[x for x in rand_assay],'.',color='#c39999',alpha=0.2)

		bp['medians'][0].set_color('navy')
		bp['medians'][1].set_color('navy')

#		ax1.set_xticklabels(["assay","random"])
		ax1.text(1,np.max(ax1.get_ybound())+100,"p-value = "+str(np.around(mannwhitneyu(bp_dict[1],bp_dict[2],alternative="greater").pvalue,5)),fontsize=10)
		count+=1

#		if k==False:
#				ax1.set_ylabel({"dnase":"DNase","atac":"ATAC"}[j])
		if j=="atac":
			ax1.set_xticklabels(["scATAC","random"])
#			ax1.set_xlabel({False:"unique",True:"specific"}[k])
		if j=="dnase":
			ax1.set_xticklabels(["iscDNase","random"])
#		if j=="dnase" and k==True:
#			dnase_patch = mpatches.Patch(color='tomato', label='DNase')
#			atac_patch = mpatches.Patch(color='skyblue', label='ATAC')
#			intersect_patch = mpatches.Patch(color='#c39999', label='random')
#			ax1.legend(handles=[dnase_patch,atac_patch,intersect_patch], frameon=False,handlelength=0.7,loc='upper center')
		ax1.spines['right'].set_visible(False)
		ax1.spines['top'].set_visible(False)
#	ax.set_xlabel("",labelpad=20)
	ax.set_ylabel("Gene expression (RPKM)",labelpad=20)
	plt.tight_layout()
	plt.savefig("../../../../Figures/Figure3/"+i+"_expr.png")





