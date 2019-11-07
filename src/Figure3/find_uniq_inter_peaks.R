# This program finds the intersecting peaks from the single cell and bulk cell assays, unique to bulk 
# and specific to one bulk assay. It can also visualize the counts of each category as a Venn diagram. 

rm(list=ls())

library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg18)
library(BiocParallel)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(pheatmap)
library(plyr)
library(RColorBrewer)
library(colorRamps)
library(ggplot2)
library(plyr)
library(eulerr)
library(stringr)

register(MulticoreParam(workers = future::availableCores()))

# find intersecting peaks unique to single cell and bulk cell and write to file
# @param peaks total set of peaks
# @param cell cell type string for file name 
# @param sc_atac_idx peak indices specific to scATAC-seq
# @param i_dnase_idx peak indices specific to iDNase-seq
# @param b_atac_idx peak indices specific to bulk ATAC-seq
# @param b_dnase_idx peak indices specific to bulk DNase-seq
# @author Jonathan Perrie
find_uniq_peaks<-function(peaks,cell,sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx){
	# unique single cell peaks
	joint_sc_idx<-as.logical(sc_atac_idx+i_dnase_idx)
	i_dnase_uniq_idx<-as.logical(joint_sc_idx-sc_atac_idx)
	sc_atac_uniq_idx<-as.logical(joint_sc_idx-i_dnase_idx)

	# unique bulk peaks
	joint_b_idx<-as.logical(b_atac_idx+b_dnase_idx)
	b_dnase_uniq_idx<-as.logical(joint_b_idx-b_atac_idx)
	b_atac_uniq_idx<-as.logical(joint_b_idx-b_dnase_idx)

	# intersecting single cell and bulk peaks
	joint_dnase_idx<-as.logical(b_dnase_uniq_idx+i_dnase_uniq_idx)
	spc_dnase<-as.logical(i_dnase_uniq_idx-as.logical(joint_dnase_idx-b_dnase_uniq_idx))
	
	joint_atac_idx<-as.logical(b_atac_uniq_idx+sc_atac_uniq_idx)
	spc_atac<-as.logical(sc_atac_uniq_idx-as.logical(joint_atac_idx-b_atac_uniq_idx))

	# write to file
	gr<-peaks[spc_atac]
	gr_df <- data.frame(seqnames=seqnames(gr),
	  starts=start(gr)-1,
	  ends=end(gr),
	  names=c(rep(".", length(gr))),
	  scores=c(rep(".", length(gr))),
	  strands=strand(gr))
	write.table(gr_df, file=paste0(getwd(),"/peaks/",cell,"_specific_b_atac_peaks.bed"), quote=F, sep="\t", row.names=F, col.names=F)

	gr<-peaks[spc_dnase]
	gr_df <- data.frame(seqnames=seqnames(gr),
	  starts=start(gr)-1,
	  ends=end(gr),
	  names=c(rep(".", length(gr))),
	  scores=c(rep(".", length(gr))),
	  strands=strand(gr))
	write.table(gr_df, file=paste0(getwd(),"/peaks/",cell,"_specific_b_dnase_peaks.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}

# plot Venn diagrams of intersecting peak groups
# @param peaks total set of peaks
# @param cell cell type string for file name 
# @param sc_atac_idx peak indices specific to scATAC-seq
# @param i_dnase_idx peak indices specific to iDNase-seq
# @param b_atac_idx peak indices specific to bulk ATAC-seq
# @param b_dnase_idx peak indices specific to bulk DNase-seq
# @author Jonathan Perrie
make_venn<-function(peaks,cell,sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx){
	joint_sc_idx<-as.logical(sc_atac_idx+i_dnase_idx)
	i_dnase_uniq_idx<-as.logical(joint_sc_idx-sc_atac_idx)
	sc_atac_uniq_idx<-as.logical(joint_sc_idx-i_dnase_idx)

	joint_b_idx<-as.logical(b_atac_idx+b_dnase_idx)
	b_dnase_uniq_idx<-as.logical(joint_b_idx-b_atac_idx)
	b_atac_uniq_idx<-as.logical(joint_b_idx-b_dnase_idx)

	joint_dnase_idx<-as.logical(b_dnase_uniq_idx+i_dnase_uniq_idx)
	spc_dnase<-as.logical(i_dnase_uniq_idx-as.logical(joint_dnase_idx-b_dnase_uniq_idx))

	joint_atac_idx<-as.logical(b_atac_uniq_idx+sc_atac_uniq_idx)
	spc_atac<-as.logical(sc_atac_uniq_idx-as.logical(joint_atac_idx-b_atac_uniq_idx))

	# single cell
	v1<-euler(c("iscDNase"=1,"scATAC"=1,"iscDNase&scATAC"=0.5))
	v1$original.values<-c(sum(i_dnase_uniq_idx),sum(sc_atac_uniq_idx),sum(joint_sc_idx-sc_atac_uniq_idx-i_dnase_uniq_idx))
	# bulk cell	
	v2<-euler(c("bulk DNase"=1,"bulk ATAC"=1,"bulk DNase&bulk ATAC"=0.5))
	v2$original.values<-c(sum(b_dnase_uniq_idx),sum(b_atac_uniq_idx),sum(joint_b_idx-b_atac_uniq_idx-b_dnase_uniq_idx))
	# DNase
	v3<-euler(c("iscDNase"=1,"bulk DNase"=1,"iscDNase&bulk DNase"=0.5))
	v3$original.values<-c(sum(i_dnase_uniq_idx-spc_dnase),sum(b_dnase_uniq_idx-spc_dnase),sum(spc_dnase))
	# ATAC
	v4<-euler(c("scATAC"=1,"bulk ATAC"=1,"scATAC&bulk ATAC"=0.5))
	v4$original.values<-c(sum(sc_atac_uniq_idx-spc_atac),sum(b_atac_uniq_idx-spc_atac),sum(spc_atac))

	# plot Venn diagrams
	p<-plot(v1,edges=list(col=c("#ff6347", "#87ceeb", "#c39999"),lwd=5),quantities=list(type=c("counts"),cex=2),labels=NULL,fills= c("white","white","white"),
		main=list(label=paste0(str_pad(names(v1$fitted.values[1]),12,side="right"),"\t\t",str_pad(names(v1$fitted.values[2]),12,side="right"),"\t\n",
		str_pad(paste0(v1$original.values[1]+v1$original.values[3]),12,side="right"),"\t\t",str_pad(paste0(v1$original.values[2]+v1$original.values[3]),12,side="right"),"\t"),cex=2))

	#p<-plot(v1,edges=c("#ff6347", "#87ceeb", "#c39999",lwd=5),quantities=list(cex=2),labels=list(cex=2),fills =c("white","white","white") )
	png(filename=paste0(getwd(),"/../../Figures/Figure3/",cell,"_sc_venn.png"))
	print(p)
	dev.off ();

	p<-plot(v2,edges=list(col=c("#9c1708", "#4c7e92", "#744b4d"),lwd=5),quantities=list(type=c("counts"),cex=2),labels=NULL,fills= c("white","white","white"),
		main=list(label=paste0(str_pad(names(v2$fitted.values[1]),12,side="right"),"\t\t",str_pad(names(v2$fitted.values[2]),12,side="right"),"\t\n",
		str_pad(paste0(v2$original.values[1]+v2$original.values[3]),12,side="right"),"\t\t",str_pad(paste0(v2$original.values[2]+v2$original.values[3]),12,side="right"),"\t"),cex=2))

	#p<-plot(v2,edges=c("#9c1708", "#4c7e92", "#744b4d",lwd=5),quantities=list(cex=2),labels=list(cex=2),fills=c("white","white","white"))
	png(filename=paste0(getwd(),"/../../Figures/Figure3/",cell,"_bulk_venn.png"))
	print(p)
	dev.off ();

	p<-plot(v3,edges=list(col=c("#ff6347", "#9c1708", "#bb4a42"),lwd=5),quantities=list(type=c("counts"),cex=2),labels=NULL,fills= c("white","white","white"),
		main=list(label=paste0(str_pad(names(v3$fitted.values[1]),12,side="right"),"\t\t",str_pad(names(v3$fitted.values[2]),12,side="right"),"\t\n",
		str_pad(paste0(v3$original.values[1]+v3$original.values[3]),12,side="right"),"\t\t",str_pad(paste0(v3$original.values[2]+v3$original.values[3]),12,side="right"),"\t"),cex=2))

	#p<-plot(v3,edges=c("#ff6347", "#9c1708", "#bb4a42",lwd=5),quantities=list(cex=2),labels=list(cex=2),fills =c("white","white","white"))
	png(filename=paste0(getwd(),"/../../Figures/Figure3/",cell,"_dnase_venn.png"))
	print(p)
	dev.off ();


	p<-plot(v4,edges=list(col=c("#87ceeb", "#4c7e92", "#6aa6bf"),lwd=5),quantities=list(type=c("counts"),cex=2),labels=NULL,fills= c("white","white","white"),
		main=list(label=paste0(str_pad(names(v4$fitted.values[1]),12,side="right"),"\t\t",str_pad(names(v4$fitted.values[2]),12,side="right"),"\t\n",
		str_pad(paste0(v4$original.values[1]+v4$original.values[3]),12,side="right"),"\t\t",str_pad(paste0(v4$original.values[2]+v4$original.values[3]),12,side="right"),"\t"),cex=2))

	#p<-plot(v4,edges=c("#87ceeb", "#4c7e92", "#6aa6bf",lwd=5),quantities=list(cex=2),labels=list(cex=2),fills =c("white","white","white"))
	png(filename=paste0(getwd(),"/../../Figures/Figure3/",cell,"_atac_venn.png"))
	print(p)
	dev.off ();
}

# find intersecting peaks unique to single cell and in bulk cell and write to file
# @param peaks total set of peaks
# @param cell cell type string for file name 
# @param sc_atac_idx peak indices specific to scATAC-seq
# @param i_dnase_idx peak indices specific to iDNase-seq
# @param b_atac_idx peak indices specific to bulk ATAC-seq
# @param b_dnase_idx peak indices specific to bulk DNase-seq
# @author Jonathan Perrie
find_uniq_inter_peaks<-function(peaks,cell,sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx){
	joint_sc_idx<-as.logical(sc_atac_idx+i_dnase_idx)
	i_dnase_uniq_idx<-as.logical(joint_sc_idx-sc_atac_idx)
	sc_atac_uniq_idx<-as.logical(joint_sc_idx-i_dnase_idx)

	joint_b_idx<-as.logical(b_atac_idx+b_dnase_idx)
	b_dnase_uniq_idx<-as.logical(b_dnase_idx)
	b_atac_uniq_idx<-as.logical(b_atac_idx)

	joint_dnase_idx<-as.logical(b_dnase_uniq_idx+i_dnase_uniq_idx)
	spc_dnase<-as.logical(i_dnase_uniq_idx-as.logical(joint_dnase_idx-b_dnase_uniq_idx))

	joint_atac_idx<-as.logical(b_atac_uniq_idx+sc_atac_uniq_idx)
	spc_atac<-as.logical(sc_atac_uniq_idx-as.logical(joint_atac_idx-b_atac_uniq_idx))

	gr<-peaks[spc_atac]
	gr_df <- data.frame(seqnames=seqnames(gr),
	  starts=start(gr)-1,
	  ends=end(gr),
	  names=c(rep(".", length(gr))),
	  scores=c(rep(".", length(gr))),
	  strands=strand(gr))
	write.table(gr_df, file=paste0(getwd(),"/peaks/inter_peaks/",cell,"_specific_b_atac_peaks.bed"), quote=F, sep="\t", row.names=F, col.names=F)

	gr<-peaks[spc_dnase]
	gr_df <- data.frame(seqnames=seqnames(gr),
	  starts=start(gr)-1,
	  ends=end(gr),
	  names=c(rep(".", length(gr))),
	  scores=c(rep(".", length(gr))),
	  strands=strand(gr))
	write.table(gr_df, file=paste0(getwd(),"/peaks/inter_peaks/",cell,"_specific_b_dnase_peaks.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}

peak_width=500

# read in blacklist peaks
peakfile<-"peaks/hg18.blacklist.bed"
black_peaks<-getPeaks(peakfile)

# read in scATAC-seq peaks 
sc_atac_peakfile<-"peaks/sc_atac_peaks.bed"
sc_atac_peaks<-sort(resize(setdiff(unique(getPeaks(sc_atac_peakfile)),black_peaks),width=peak_width,fix="center"))

# read in iDNase-seq peaks 
i_dnase_peakfile<-"peaks/i_dnase_peaks.bed"
i_dnase_peaks<-sort(resize(setdiff(unique(getPeaks(i_dnase_peakfile)),black_peaks),width=peak_width,fix="center"))

# read in bulk ATAC-seq peaks
bulk_atac_peakfile<-"peaks/b_atac_peaks.bed"
b_atac_peaks<-sort(resize(setdiff(unique(getPeaks(bulk_atac_peakfile)),black_peaks),width=peak_width,fix="center"))

# read in bulk DNase-seq peaks
bulk_dnase_peakfile<-"peaks/b_dnase_peaks.bed"
b_dnase_peaks<-sort(resize(setdiff(unique(getPeaks(bulk_dnase_peakfile)),black_peaks),width=peak_width,fix="center"))

# read in merged B cell peaks
b_bulk_peakfile<-"peaks/merged_b_cell_peaks.bed"
b_peaks<-sort(resize(setdiff(unique(getPeaks(b_bulk_peakfile)),black_peaks),width=peak_width,fix="center"))

# read in merged T cell peaks
t_bulk_peakfile<-"peaks/merged_t_cell_peaks.bed"
t_peaks<-sort(resize(setdiff(unique(getPeaks(t_bulk_peakfile)),black_peaks),width=peak_width,fix="center"))

# read in merged monocyte peaks 
m_bulk_peakfile<-"peaks/merged_mono_peaks.bed"
m_peaks<-sort(resize(setdiff(unique(getPeaks(m_bulk_peakfile)),black_peaks),width=peak_width,fix="center"))



# find index matches for each cell type and then generate unique peak sets and Venn diagram
# B cells 
sc_atac_idx<-seq(length(b_peaks)) %in% queryHits(findOverlaps(b_peaks,sc_atac_peaks))
i_dnase_idx<-seq(length(b_peaks)) %in% queryHits(findOverlaps(b_peaks,i_dnase_peaks))
b_atac_idx<-seq(length(b_peaks)) %in% queryHits(findOverlaps(b_peaks,b_atac_peaks))
b_dnase_idx<-seq(length(b_peaks)) %in% queryHits(findOverlaps(b_peaks,b_dnase_peaks))
find_uniq_peaks(b_peaks,"b_cell",sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx)
make_venn(b_peaks,"b_cell",sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx)
find_uniq_inter_peaks(b_peaks,"b_cell",sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx)

# T cells
sc_atac_idx<-seq(length(t_peaks)) %in% queryHits(findOverlaps(t_peaks,sc_atac_peaks))
i_dnase_idx<-seq(length(t_peaks)) %in% queryHits(findOverlaps(t_peaks,i_dnase_peaks))
b_atac_idx<-seq(length(t_peaks)) %in% queryHits(findOverlaps(t_peaks,b_atac_peaks))
b_dnase_idx<-seq(length(t_peaks)) %in% queryHits(findOverlaps(t_peaks,b_dnase_peaks))
find_uniq_peaks(t_peaks,"t_cell",sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx)
make_venn(t_peaks,"t_cell",sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx)
find_uniq_inter_peaks(t_peaks,"t_cell",sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx)

# Monocytes 
sc_atac_idx<-seq(length(m_peaks)) %in% queryHits(findOverlaps(m_peaks,sc_atac_peaks))
i_dnase_idx<-seq(length(m_peaks)) %in% queryHits(findOverlaps(m_peaks,i_dnase_peaks))
b_atac_idx<-seq(length(m_peaks)) %in% queryHits(findOverlaps(m_peaks,b_atac_peaks))
b_dnase_idx<-seq(length(m_peaks)) %in% queryHits(findOverlaps(m_peaks,b_dnase_peaks))
find_uniq_peaks(m_peaks,"mono",sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx)
make_venn(m_peaks,"mono",sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx)
find_uniq_inter_peaks(m_peaks,"mono",sc_atac_idx,i_dnase_idx,b_atac_idx,b_dnase_idx)
