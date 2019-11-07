# This program generates plots using the chromVAR package to show how our assay differs from scATAC-seq.

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

set.seed(1)
register(MulticoreParam(workers = future::availableCores()))

# calculate the deviations using chromVAR
# @param bedfiles bedfile file paths 
# @param peaks chromVAR peaks 
# @param motifs reference motif database
# @return chromVARDeviations object 
# @author Jonathan Perrie
calculate_deviations <- function(bedfiles,peaks,motifs){
	# get counts and filter out peaks that are out of range
	frag_counts<-getCounts(bedfiles,peaks,paired=FALSE,format="bed")
	frag_counts<-frag_counts[!(end(peaks) > seqlengths(BSgenome.Hsapiens.UCSC.hg18)[as.character(seqnames(peaks))])&(start(peaks)>0),]
	frag_counts <- addGCBias(frag_counts, genome = BSgenome.Hsapiens.UCSC.hg18)
	counts_filtered <- filterSamples(frag_counts,shiny = FALSE)

	# compute deviations
	counts_filtered <- filterPeaks(counts_filtered, non_overlapping = TRUE,min_fragments_per_peak = 1)
	motif_ix <- matchMotifs(motifs, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg18)
	bg <- getBackgroundPeaks(object = counts_filtered)
	dev <- computeDeviations(object = counts_filtered, annotations = motif_ix,background_peaks = bg)
	return(dev)
}

# merge deviations from two separate chromVAR runs 
# @param dev1 first deviations object 
# @param dev2 second deviations object 
# @return chromVARDeviations object
# @author Jonathan Perrie 
merge_deviations <- function(dev1,dev2) { 
	z<-cbind(assays(dev1)$z,assays(dev2)$z)
	deviations<-cbind(assays(dev1)$deviations,assays(dev2)$deviations)
	depth<-rbind(colData(dev1),colData(dev2))
	tfs<-data.frame(name=rowData(dev1)$name,
		fractionMatches=rowData(dev1)$fractionMatches+rowData(dev2)$fractionMatches,
		fractionBackgroundOverlap=rowData(dev1)$fractionBackgroundOverlap+rowData(dev2)$fractionBackgroundOverlap)
	
	dev<-SummarizedExperiment(assays=list(deviations=deviations,z=z),
		colData=depth,rowData=tfs)
	dev<-as(dev,"chromVARDeviations")
	return(dev)
}
      
# find order of motifs for each assay separately 
# @param dev deviations object
# @param variability variability deviations object
# @return var approximate variability object with no adjusted p-values and a decreasing order
# @author Jonathan Perrie
findMotifOrder <- function(dev,variability){
	zscore<-assays(dev)$z
	zscore[is.na(zscore)]<-0

	var<-data.frame(variability=rowSds(zscore))
	var$name<-variability$name
	rownames(var)<-rownames(variability)
	var$order=-1
	var[order(var$variability,decreasing=TRUE),]$order=seq(dim(var)[1])
	return(var)
}

# find counts for some cluster with a mixture of cell types across assays
# @param cluster_col cluster slice
# @return cell_count counts of each cell type 
# @author Jonathan Perrie
count_cell_freq <- function(cluster_col){
        b_cell=sum(cluster_col[grepl("B cell",names(cluster_col))])
        t_cell=sum(cluster_col[grepl("T cell",names(cluster_col))])
        mono=sum(cluster_col[grepl("Mono",names(cluster_col))])
        cell_count=c(b_cell,t_cell,mono)
        names(cell_count)<-c("B cell cluster","T cell cluster","Monocyte cluster")
        return(cell_count)
}

# adjust ggplot object colours
# @param p ggplot object
# @param palette colours  
# @author jbaums, https://stackoverflow.com/questions/34601194/change-colours-to-defined-palette-for-ggplot-objects/34601726
change_colours <- function(p, palette) {
  n <- nlevels(p$data[[deparse(p$mapping$group)]])
  tryCatch(as.character(palette), 
           error=function(e) stop('palette should be a vector of colours', call.=FALSE))
  if(n > length(palette)) stop('Not enough colours in palette.')
  pal <- function(n) palette[seq_len(n)]
  p + theme(panel.border = element_blank()) + discrete_scale('colour', 'foo', pal) 
}      
 
# cluster cells using top 50 transcription factors
# @param tf_z matrix of accessibility z-scores for cells over motifs 
# @param assay string for writing file and denoting assay type
# @param dev deviations chromVAR data frame 
# @param top_inter_var special case where we have a subset of motifs with variabilities above some threshold     
# @author Jonathan Perrie 
cluster_tfs<-function(tf_z,assay,dev,top_inter_var=NULL,n_clusters){
	var<-findMotifOrder(dev,variability)

	# get top 50 transcription factors
	if (assay=="iDNase"){
		index=grepl("GB",colnames(tf_z))
		top_50_tf<-var[c(order(var$variability,decreasing=TRUE)[1:50]),]
	} else if (assay=="scATAC"){
		index=!grepl("GB",colnames(tf_z))
		top_50_tf<-var[c(order(var$variability,decreasing=TRUE)[1:50]),]
	} else if (!is.null(top_inter_var)){
		top_50_tf<-top_inter_var[c(order(top_inter_var$variability,decreasing=TRUE)[1:50]),]
		assay="shared"
		index=rep(TRUE,dim(tf_z)[2])
	} else {
		index=rep(TRUE,dim(tf_z)[2])
		top_50_tf<-var[c(order(var$variability,decreasing=TRUE)[1:50]),]
	}

	# cluster the transcription factors
	top_50_tf_z<-tf_z[rownames(top_50_tf),index]
	rownames(top_50_tf_z)<-top_50_tf$name
	cluster<-kmeans(t(top_50_tf_z),n_clusters,)$cluster
	cluster<-mapvalues(cluster,from=c(1,2,3),to=c("cluster 1","cluster 2","cluster 3"))

	# rename columns
	cell_names<-toupper(sapply(strsplit(colnames(dev), "_"), "[", 1))
	cell_names[cell_names=="B"]=paste0(assay," B cell")
	cell_names[cell_names=="T"]=paste0(assay," T cell")
	cell_names[cell_names=="MONO"]=paste0(assay," Monocyte")
	cell_names[cell_names=="CD4"]=paste0(assay," T cell")

	# match cluster labels with cell type and write to file
	cell_labels=as.numeric(as.factor(cell_names))
	cluster_row<-data.frame(cell_names)
	cluster_row$cluster<-cluster
	rownames(cluster_row)<-names(cluster)
	colnames(cluster_row)<-c("cell type","cluster")
	cluster_row<-cluster_row[order(cluster_row$cluster),]
	write.table(table(cluster_row),paste0("cluster/",assay,"_cluster_results.txt"))
	return(cluster_row)
}
      
# get the peak and blacklist sets and standardize peak set  
peakfile<-"peaks/merged_sc_peaks.bed"
peaks<-unique(getPeaks(peakfile))
peakfile<-"peaks/hg18.blacklist.bed"
black_peaks<-getPeaks(peakfile)
peaks<-setdiff(peaks,black_peaks)
peaks<-resize(peaks,width=500,fix="center")
peaks<-sort(peaks)

# get list of ATAC  and DNAse file names 
atac_bedfiles<-Sys.glob("sc_atac/bedfiles/*hg18.bed")
dnase_bedfiles<-c(Sys.glob("i_dnase/bedfiles/B_GB*"),
	Sys.glob("i_dnase/bedfiles/T_GB*"),
	Sys.glob("i_dnase/bedfiles/mono_GB*"))

motifs <- getJasparMotifs()     
atac_dev<-calculate_deviations(atac_bedfiles,peaks,motifs)
dnase_dev<-calculate_deviations(dnase_bedfiles,peaks,motifs)

# filter in post 
atac_dev<-atac_dev[,!is.na(colSums(assays(atac_dev)$z))]
dnase_dev<-dnase_dev[,!is.na(colSums(assays(dnase_dev)$z))]

# downsample
n_samples<-210
atac_dev<-atac_dev[,c(sample(colnames(atac_dev)[grepl("t_cell",colnames(atac_dev))],round(n_samples/14*13)),
sample(colnames(atac_dev)[grepl("cd4",colnames(atac_dev))],round(n_samples/14)),
sample(colnames(atac_dev)[grepl("b_cell",colnames(atac_dev))],n_samples),
sample(colnames(atac_dev)[grepl("mono",colnames(atac_dev))],n_samples))]

dnase_dev<-dnase_dev[,c(sample(colnames(dnase_dev)[grepl("mono_GB",colnames(dnase_dev))],n_samples),
			sample(colnames(dnase_dev)[grepl("B_GB",colnames(dnase_dev))],n_samples),
			sample(colnames(dnase_dev)[grepl("T_GB",colnames(dnase_dev))],n_samples))]


dev<-merge_deviations(atac_dev,dnase_dev)

# plot the variability scores for each TF
variability <- computeVariability(dev)
png(filename="../../Figures/Figure3/var.png")
plotVariability(variability, use_plotly = FALSE)
dev.off()

# truncate a copy of the deviation z scores 
tf_z<-assays(dev)$z
tf_z[tf_z>2]=2
tf_z[tf_z<(-2)]=-2

# find top 50 most variable motifs that have a variability of at least 1.05 in both assays
atac_var<-findMotifOrder(atac_dev,variability)
dnase_var<-findMotifOrder(dnase_dev,variability)


threshold=1.05
top_inter_var<-variability[rownames(variability) %in% intersect(rownames(variability[dnase_var$variability>threshold,]),rownames(variability[atac_var$variability>threshold,])),]
top_50_tf<-top_inter_var[c(order(top_inter_var$variability,decreasing=TRUE)[1:50]),]
top_50_tf_z<-tf_z[rownames(top_50_tf),]
rownames(top_50_tf_z)<-top_50_tf$name
# write clusterings to file
n_clusters=3
dnase_clusters <- cluster_tfs(tf_z,"DNase",dnase_dev,top_inter_var=NULL,n_clusters)
atac_clusters <- cluster_tfs(tf_z,"ATAC",atac_dev,top_inter_var=NULL,n_clusters)
shared_clusters <- cluster_tfs(tf_z,"all",dev,top_inter_var=NULL,n_clusters)
top_clusters <- cluster_tfs(tf_z,"all",dev,top_inter_var,n_clusters)

# cluster deviations
cluster<-kmeans(t(top_50_tf_z),n_clusters,)$cluster
cluster<-mapvalues(cluster,from=c(1,2,3),to=c("cluster 1","cluster 2","cluster 3"))

# rename files for presentation
cell_names<-sapply(strsplit(colnames(dev), "_GB"), "[", 1)
cell_names[cell_names=="B"]="iscDNase B cell"
cell_names[cell_names=="T"]="iscDNase T cell"
cell_names[cell_names=="mono"]="iscDNase Monocyte"
cell_names<-sapply(strsplit(cell_names, "_"), "[", 1)
cell_names[cell_names=="b"]="scATAC B cell"
cell_names[cell_names=="t"]="scATAC T cell"
cell_names[cell_names=="mono"]="scATAC Monocyte"
cell_names[cell_names=="cd4"]="scATAC T cell"

# build cluster data structure for pheatmap
cell_labels=as.numeric(as.factor(cell_names))
cluster_row<-data.frame(cell_names)
cluster_row$cluster<-cluster
rownames(cluster_row)<-names(cluster)
colnames(cluster_row)<-c("cell type","cluster")
cluster_row<-cluster_row[order(cluster_row$cluster),]


# match names to colors
cell_cluster_color=c("B cell cluster"="skyblue","T cell cluster"="goldenrod","Monocyte cluster"="tomato")
cluster_color=list("cell type"=c("iscDNase B cell"="skyblue4","scATAC B cell"="skyblue1",
                                                                "iscDNase Monocyte"="tomato4","scATAC Monocyte"="tomato1",
                                                                "iscDNase T cell"="goldenrod4","scATAC T cell"="goldenrod1"),
                                        "cluster"=c("cluster 1"=unname(cell_cluster_color[names(which.max(count_cell_freq(table(cluster_row)[,1])))]),
						"cluster 2"=unname(cell_cluster_color[names(which.max(count_cell_freq(table(cluster_row)[,2])))]),
						"cluster 3"=unname(cell_cluster_color[names(which.max(count_cell_freq(table(cluster_row)[,3])))])),
                                        "atac_cluster"=c("cluster 1"=unname(cell_cluster_color[names(which.max(count_cell_freq(table(atac_clusters)[,1])))]),
                                                "cluster 2"=unname(cell_cluster_color[names(which.max(count_cell_freq(table(atac_clusters)[,2])))]),
                                                "cluster 3"=unname(cell_cluster_color[names(which.max(count_cell_freq(table(atac_clusters)[,3])))])),
                                        "dnase_cluster"=c("cluster 1"=unname(cell_cluster_color[names(which.max(count_cell_freq(table(dnase_clusters)[,1])))]),
                                                "cluster 2"=unname(cell_cluster_color[names(which.max(count_cell_freq(table(dnase_clusters)[,2])))]),
                                                "cluster 3"=unname(cell_cluster_color[names(which.max(count_cell_freq(table(dnase_clusters)[,3])))])),
                                        "assay"=c("iscDNase"="maroon4","scATAC"="darkorange2","both"="forestgreen"))

# only plot specific motif
tf_display_idx<-rownames(top_50_tf_z) %in% c("SPIC","CEBPA","TCF4","TBX2","PAX5","RUNX2") 
tf_rownames<-rownames(top_50_tf_z)[tf_display_idx]
tf_display<-mapvalues(tf_display_idx,from=c(0,1),to=c("","fill"))
tf_display[tf_display=="fill"]=tf_rownames
rownames(top_50_tf_z)=tf_display

# get heatmap from data structure (figure 3A)
png(filename="../../Figures/Figure3/heatmap_select.png",width=797,height=960)
pheatmap(top_50_tf_z[,rownames(cluster_row)],color=colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(8),cluster_cols=FALSE,
	show_colnames=FALSE,annotation_col=cluster_row,clustering_method="ward.D2",clustering_distance_rows = "correlation",annotation_color=cluster_color,
	,fontsize=14)
dev.off();

# set up cluster data structure again, this time using dev directly (requirement for chromvar t-SNE method)
cluster_row<-data.frame(cell_names)
cluster_row$cluster<-cluster
rownames(cluster_row)<-names(cluster)
colnames(cluster_row)<-c("cell_type","cluster")
colData(dev)$cell_type<-as.character(cluster_row$cell_type)
colData(dev)$cluster<-as.factor(cluster_row$cluster)
colData(dev)$single_cluster <- as.factor(rbind(atac_clusters,dnase_clusters)[rownames(colData(dev)),]$cluster)


# truncate deviation z scores 
assays(dev)$z[which(assays(dev)$z<(-2))] = -2
assays(dev)$z[which(assays(dev)$z>(2))] = 2
variability <- computeVariability(dev)
png(filename="../../Figures/Figure3/adjusted_var.png")
plotVariability(variability, use_plotly = FALSE)
dev.off();

# plot t-SNE results for cell type and specific motifs (figure 3B)
tsne_results <- deviationsTsne(dev, threshold = 1.15, perplexity = 30)
tsne_plots <- plotDeviationsTsne(dev, tsne_results, sample_column="cell_type",shiny = FALSE)
p<-change_colours(tsne_plots$cell_type,c("skyblue1","tomato1","goldenrod1","skyblue4","tomato4","goldenrod4"))
png(filename="../../Figures/Figure3/tsne_cell.png")
#p+labs(color="cell type")
p+guides(color=FALSE)
dev.off();

tsne_plots <- plotDeviationsTsne(dev, tsne_results, sample_column="cluster",shiny = FALSE)
png(filename="../../Figures/Figure3/tsne_cluster.png")
p<-change_colours(tsne_plots$cluster,cluster_color$cluster)
#p+labs(color="cluster")
p+guides(color=FALSE)
dev.off();

index=grepl("GB",colnames(tf_z))
tsne_plots <- plotDeviationsTsne(dev[,index], tsne_results[index,], sample_column="single_cluster",shiny = FALSE)
png(filename="../../Figures/Figure3/tsne_dnase_cluster.png")
p<-change_colours(tsne_plots$single_cluster,cluster_color$dnase_cluster)
#p+labs(color="single_cluster")
p+guides(color=FALSE)
dev.off();

index=!grepl("GB",colnames(tf_z))
tsne_plots <- plotDeviationsTsne(dev[,index], tsne_results[index,], sample_column="single_cluster",shiny = FALSE)
png(filename="../../Figures/Figure3/tsne_atac_cluster.png")
p<-change_colours(tsne_plots$single_cluster,cluster_color$atac_cluster)
#p+labs(color="single_cluster")
p+guides(color=FALSE)
dev.off();


tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation_name="PAX5",shiny = FALSE)
png(filename="../../Figures/Figure3/tsne_pax5.png")
tsne_plots[[1]]+guides(color=FALSE)
dev.off();

tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation_name="CEBPA",shiny = FALSE)
png(filename="../../Figures/Figure3/tsne_cebpa.png")
tsne_plots[[1]]+guides(color=FALSE)
dev.off();

tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation_name="SPIC",shiny = FALSE)
png(filename="../../Figures/Figure3/tsne_spic.png")
tsne_plots[[1]] +guides(color=FALSE)
dev.off();

tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation_name="TCF4",shiny = FALSE)
png(filename="../../Figures/Figure3/tsne_tcf4.png")
tsne_plots[[1]]+guides(color=FALSE)
dev.off();

tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation_name="RUNX2",shiny = FALSE)
png(filename="../../Figures/Figure3/tsne_runx2.png")
tsne_plots[[1]]+guides(color=FALSE)
dev.off();

tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation_name="TBX2",shiny = FALSE)
png(filename="../../Figures/Figure3/tsne_tbx2.png")
tsne_plots[[1]]+guides(color=FALSE)
dev.off();

png(filename="../../Figures/Figure3/tsne_scale.png")
tsne_plots[[1]]
dev.off();


