/*
 * UCSCBEDOperator.h
 *
 *  Created on: Sep 28, 2009
 *      Author: hugq
 */

#ifndef UCSCBEDOPERATOR_H_
#define UCSCBEDOPERATOR_H_

#include "bed.h"
#include "ucsc.h"
#include "Matrix.h"
#include "SequenceTransform.h"

struct TMP{
	Str strand;
	Str gid;
	Str chrom;
	set<pair<int,int> > exons;
};

map<Str, set<UCSC, UCSC::sortByTxStart> > UCSCGenesFrmGTF(Str fileName);

vector<pair<double,int> > exp_islandPercent(map<Str, vector<UCSC> > chr_genes_nolap,
	map<Str,double> name_rpkm, GenomePartIdentify regionId,
	Str islandLabel, Str filename, int step = 100);
//assume chr_genes are sorted according to txstat
//FCFS = first come first server. An island only assigned to the region that first call the function.
void assignIsland2UCSC(map<Str, vector<UCSC> >& chr_genes,
	map<Str, set<BED3> >& chr_islands, GenomePartIdentify regionId, Str islandLabel, bool FCFS = false);
//assign source islands to target islands
void assignIsland2Island(map<Str, set<BED3> >& target_islands,
	const map<Str, set<BED3> >& source_islands, Str source_island_label);
//RPBM
int assignReadCount2Island(map<Str, set<BED3> >& target_island,
		Str readf, Str label, int shift =0);

int smoothed_counts_in_islands(map<Str, vector<BED3> >& islands,
		Str readf, Str label, int shift, int upR, int dnR );
//not Dustin's method, the summit should be initialized before use
//upR and DnR are relative to the summit, resolution 1 bp

void assignReadCount5Columns2Island(map<Str, set<BED3> >& target_island,
		Str readf, Str label, int shift =0);
//Filter reads outsides islands
void filterReadsOutOfIsland(const map<Str, set<BED3> >& target_island,
		Str readf, Str outf, int shift =0);

//assume chr_genes are sorted according to txstat
void setReadNumberBoundUCSCOn(map<Str, vector<UCSC> >& chr_genes,
		map<Str, set<BED6, BED6::sortByShifted2Point> > chr_reads,
		GenomePartIdentify regionId, Str islandLabel);

Pa_I_I getIntrestRegion(const UCSC& ucsc, GenomePartIdentify regionId);

map<Str, set<BED4> > getIntrestReions(const map<Str, set<UCSC, UCSC::sortByTxStart> >& chr_genes,
		GenomePartIdentify regionId);

void assignGenomicRegions2Islands(
		map<Str, set<BED4> >& chr_regions,
		const map<Str, set<BED3> >& chr_islands,
		map<Str, map<BED3,set<GenomePartIdentify> > >& island_regions,
		GenomePartIdentify regionId);

void assignUCSC2Islands(  map<Str, set<UCSC, UCSC::sortByTxStart> >& chr_genes,
		map<Str, set<BED3> >& chr_islands );
//use genes from both strands
void generateReadDensityProfile(const map<Str, vector<UCSC> >& chr_genes, Str readf, Str out, int shift =0);
//use genes from both strands
void generateTSSTESReadDensityHeatMap(const map<Str, vector<UCSC> >& chr_genes, Str readf, Str out, int shift =0, int bins = 200 );
//only use genes in positive strands
void generateReadDensityProfile(const map<Str, vector<UCSC> >& chr_genes, Str readf, Str out, bool is5Pr);
vector<double> generateReadDensityProfileOnRegions(map<Str, set<BED3> > regions, Str readf, Str out, int shift,  Str label,int window);
//assume rcds in readf sorted
pair<vector<double>, vector<double> > generateReadDensityProfileOnRegions(map<Str, set<BED4> >& regions, Str readf, Str outf, Str label, int sliceN);
pair<vector<double>, vector<double> > generateReadDensityProfile5ColumnsOnRegions(map<Str, set<BED4> >& regions, Str readf, Str outf, Str label, int sliceN);
void generateRPBMSummaryGraph(Str bed_f, Str chr_len_f, Str out_f, int windowsize, int shift, bool invert = false );
void generateRPBMSummaryGraphStrand(Str bed_f, Str chr_len_f, Str out_f, int windowsize, int shift, Str _strand );
vector<double> readCountRatio(Str read1_f, Str read2_f, map<Str, set<BED3> >& regions, int shift);
std::pair<int,int> readInIsland(map<Str, set<BED3> >& regions, Str readf, int shift );
map<Str, set<BED3> > get_intrest_regions(map<Str, set<UCSC, UCSC::sortByTxStart> > chr_genes,
		Str chr_len_f, Str tag);
map<Str, set<BED3> > notOverlapWithPromoter(map<Str, set<BED3> > pro, map<Str, set<BED3> > beds);
map<Str, set<BED3> > overlapWithPromoter(map<Str, set<BED3> > pro, map<Str, set<BED3> > beds);
void attachNuclsScore2Islands(
		map<Str, set<BED3> >& regions, Str nuc_readf, Str label, int window);
void attachNuclsScore2Islands(
		map<Str, set<BED4> >& regions, Str nuc_readf, Str label, int window);
void attachPairEndTag2Islands( map<Str, set<BED3> >& regions, Str paie_end_f, Str label);

map<Str, SignalPos > scoreWithPWM(M_D pwm, vector<pair<Str, Str> > fast, bool bothStrand = false );
map<Str, set<BED4> > parsSp2Bed4(map<Str, SignalPos > sp, int upL); //sp.first in format chr1:1111-222323

map<Str, set<BED3> > BEDFrmFASTA( Str fileName );

void write4quatiles(vector<pair<double, BED3> > beds, Str tag );

map<Str,set<BED3> > processSICER11( Str f, double fc, double p );
map<Str, set<BED3> > intersect(map<Str, set<BED3> > bed1, map<Str, set<BED3> > bed2);
void generateRPBMSummaryWig(Str bed_f, Str chr_len_f, Str out_f, int windowsize, int shift, Str chr );
void generateTSSTESReadDensity4Profiles(const map<Str, vector<UCSC> >& chr_genes_tmp, Str readf, Str out_f, int shift);
void generalHeatmap(map<Str, vector<BED3> > chr_genes_tmp, Str readf, Str out_f, int upWndN, int dwnWndN, int wndSz, int shift, int bins, Str lg);
void generalHeatmap_PE(map<Str, vector<BED3> > chr_genes_tmp, Str readf, Str out_f, int upWndN, int dwnWndN, int wndSz, int bins);
void generateTSSTESReadDensity4Profiles_Low_preGrp(const map<Str, vector<UCSC> >& chr_genes_tmp, Str readf, Str out_f, int shift);
void generateTSSTESReadDensity4Profiles_Low(const map<Str, vector<UCSC> >& chr_genes_tmp, Str readf, Str out_f, int shift, int bins );

#endif /* UCSCBEDOPERATOR_H_ */
