/*
 * BED.h
 *
 *  Created on: Sep 9, 2009
 *      Author: hugq
 */

#ifndef BED_H_
#define BED_H_

#include "TypeDefBase.h"
#include "OftenUsedOperatLib.h"

struct SignalPos{
	Str seq;
	int relative_start_pos;
	double scr;
	bool inOppositeStrand;
};

class BED3;

class UCSC{
public:
	UCSC(){
	}
	//Read gene record from UCSC gene file
	friend std::istream &operator>>(std::istream &stream, UCSC& ucsc);
	//Output as a format of UCSC gene file
	friend std::ostream &operator<<(std::ostream &stream, const UCSC& ucsc);

	bool operator <( const UCSC & rhs )const{
		return make_pair(txStart, make_pair(txEnd,info)) < make_pair( rhs.txStart,
				make_pair(rhs.txEnd,rhs.info));
	//	return start < rhs.start;
	}//*/
	struct sortByTxStart
	{
	  bool operator()(const UCSC& lhs, const UCSC& rhs) const
	  {
	    return make_pair(lhs.txStart, lhs.name+lhs.strand+lhs.alignID+lhs.proteinID+lhs.info)
	    		< make_pair(rhs.txStart, rhs.name+rhs.strand+rhs.alignID+rhs.proteinID+rhs.info);
	  }
	};
	struct sortByTxStartTxEnd
	{
	  bool operator()(const UCSC& lhs, const UCSC& rhs) const
	  {
	    return //lhs.txEnd < rhs.txEnd;
		make_pair(lhs.txStart,make_pair(lhs.txEnd,lhs.info))
			< make_pair( rhs.txStart, make_pair(rhs.txEnd,rhs.info));;
	  }
	};

	struct sortByTSS
	{
	  bool operator()(const UCSC& lhs, const UCSC& rhs) const
	  {
	    return make_pair(lhs.TSS, lhs.name+lhs.strand) < make_pair(rhs.TSS, rhs.name+rhs.strand);
	  }
	};

	Pa_I_I getFrstIntron() const{
		if(exons.size() > 1){
			if( strand == "+"){
				set<Pa_I_I>::iterator it2 = exons.begin();
				++it2;
				return make_pair(exons.begin()->second, it2->first);
			}
			else{
				set<Pa_I_I>::iterator it2 = exons.end();
				--it2;
				set<Pa_I_I>::iterator it1 = it2;
				--it1;
				return make_pair(it1->second, it2->first);
			}
		}else{
			return make_pair(-1,-1);
		}
	}
	Str name, chrom, strand, proteinID, alignID, alia;
	Str relatives, info;
	map<GenomePartIdentify,set<Str> > genepart_boundby;
	map<GenomePartIdentify, map<Str, int> > genepart_islandLabel_readNumber;
	map<GenomePartIdentify, map<Str, set<BED3> > > islandsOn;//Str label
	int txStart, txEnd, TSS, TSE;
	int cdsStart, cdsEnd;
	set<Pa_I_I> exons;
	struct BOX{
		int start;
		int end;
		Str beg_seq;
		Str end_seq;
		Str strand;
		Str info;
		double rcd;
	};
	typedef BOX EXON;
	typedef BOX INTRON;
	vector<Pa_I_I> exons53, intron53;
	vector<INTRON> INTRON53;
	map<Str, double > id_score;
	double island_score;

	bool ismasked;

};

//Class for bed lines with 3 values: chrom, start and end
class BED3{
public:
	BED3(Str chro, int st, int en): chrom(chro), start(st), end(en){
	}
	BED3(){}
	Pa_I_I getRegion(){
		return make_pair(start,end);
	}
	bool isOverlap(BED3 r){
		return ((end > r.start) && (start <r.end) );
	}
	bool operator <( const BED3 & rhs )const{
		return make_pair(make_pair(chrom,start), make_pair(end, info)) < make_pair(make_pair(rhs.chrom, rhs.start), make_pair(rhs.end,rhs.info));
//		return make_pair(start, make_pair(end, chrom+info)) < make_pair(rhs.start, make_pair(rhs.end,rhs.chrom+rhs.info));
//		return make_pair(start, make_pair(end, chrom)) < make_pair(rhs.start, make_pair(rhs.end,rhs.chrom));
	//	return start < rhs.start;
	}//*/
	bool operator == ( const BED3 & rhs )const{
		return start == rhs.start && chrom == rhs.chrom;
	//	return start < rhs.start;
	}//*/
	//Read gene record from UCSC gene file
	void initFromLineOfFile(std::istream &stream){
		stream>>chrom>>start>>end;
		Str line;
		std::getline( stream, info);
	}
	void shift(int shift){
	}
	Str chrom, info, info1, strand, seq, status;//infor is to store some additional infor
	int start;
	int end, txStart, txEnd;
	int shifted2Point, sumit, sumitTagNum, n, shft, wp, oth, pos_N, neg_N;
	bool used;
	double FDR, FC, pvalue, readInTreat, readInControl, island_score, den, info_p
	, Normalized_Readcount_A, Normalized_Readcount_B;
	map<Str, set<BED3> > islandsOn;
	map<Str, map<BED3,int> > islandsOn2;
	map<GenomePartIdentify, set<UCSC, UCSC::sortByTxStart> > UCSCOn;
	map<Str, double > readCountOn, readDensityOn, id_score;
	map<Str, int > sumitsOn;
	map<Str, vector<double>  > readCountsOn;
	map<Str, pair<int,int>  > sumitPosCountOn;
	map<Str, vector<double> > nucl_scr;
	map<Str, vector<pair<int,int> > > pair_end_mid_len;
	map<Str, vector<double>  > gaussianDensityOn;
	set<Str> infos;
	vector<double> den_array;
	SignalPos sp;
};

class BED4: public BED3{
public:
	BED4(Str chro, int st, int en, Str stran): BED3(chro, st, en), strand(stran){
	}
	BED4(){}
	Pa_I_I getRegion(){
		return make_pair(start,end);
	}
	bool isOverlap(BED3 r){
		return ((end > r.start) && (start <r.end) );
	}
	bool operator <( const BED4 & rhs )const{
		if(start < rhs.start)
			return true;
		else if( start == rhs.start ){
			if( end < rhs.end )
				return true;
		else if(strand < rhs.strand )
			return true;
		}
		return false;
	}//*/
	//Read gene record from UCSC gene file
	void initFromLineOfFile(std::istream &stream){
		stream>>chrom>>start>>end>>strand;
		Str line;
		std::getline( stream, line);
		info = "";
	}
	void shift(int shift){
			shifted2Point = strand == "+" ? start + shift
					: end - shift;
		}
	Str strand;
};
//Class for bed lines with 6 values:  chrom, start, end, name, score, strand
class BED6: public BED3{
public:
	BED6(Str chrom, int start, int end, Str nam, Str sc, Str stran):
		BED3(chrom, start, end), name(nam), score(sc), strand(stran){
	}
	BED6(){}
	struct sortByShifted2Point
	{
		bool operator()(const BED6& lhs, const BED6& rhs) const
		{
			return make_pair(lhs.shifted2Point, lhs.strand) <
				make_pair(rhs.shifted2Point, rhs.strand);
		}
	};
	struct sortByStart
	{
		bool operator()(const BED6& lhs, const BED6& rhs) const
		{
			return make_pair(make_pair(lhs.start, lhs.end), lhs.shifted2Point)
				< make_pair(make_pair(rhs.start, rhs.end), rhs.shifted2Point);
		}
	};
	void shift(int shift){
			shifted2Point = strand == "+" ? start + shift
					: end - shift;
		}
	void initFromLineOfFile(std::istream &stream){
		BED3::initFromLineOfFile(stream);
		stream>>name>>score>>strand;
	}
	Str name, score, strand;
};

template <class T, class S>
void BEDFrmFile( Str fileName, map<Str, vector<T, S> >& chr_bed, int shift = 0){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	int NN = 0;
	while( !in.eof() ){
		T island;
		island.initFromLineOfFile(in);
		island.shift(shift);
		island.island_score = ++NN;
		if( !island.chrom.empty() ){
			chr_bed[island.chrom].push_back(island);
		}
	}
	cout<<"# of records from "<<fileName<<": "<<getSize(chr_bed)<<endl;
	in.close();
}

template <class T, class S>
void BEDFrmFile( Str fileName, map<Str, set<T, S> >& chr_bed, int shift = 0){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	double n = 0;
	//std::map<Str,Str> filter = get_a_b("/home/hugq/workspace/CommonProg/mRNA_rand_list.txt");
	while( !in.eof() ){
		T island;
		island.initFromLineOfFile(in);
		island.shift(shift);
		if( !island.chrom.empty() ){
			++n;
			island.island_score = n;
			//if( filter.find(island.chrom) != filter.end() )
				chr_bed[island.chrom].insert(island);
		}
	}
	cout<<"# of records from \""<<fileName<<"\": "<<n<<endl;
	cout<<"# of records (1 read per position): "<<getSize(chr_bed)<<'('<<100*getSize(chr_bed)/n<<')'<<endl;
	in.close();
}

void BEDFrmSICCERSummary( Str fileName, map<Str, set<BED3> >& chr_bed, int shift = 0);
void BEDFrmPeakSplitter( Str fileName, map<Str, set<BED3> >& chr_bed);
void BEDFrmSICERprobscoreisland( Str fileName, map<Str, set<BED3> >& chr_bed, int shift = 0);
void BEDFrmMACSXLS( Str fileName, map<Str, set<BED3> >& chr_bed, int shift = 0, bool hascontrol = false);
void BEDFrmMACS2XLS( Str fileName, map<Str, set<BED3> >& chr_bed);
void BEDFrmNPS( Str fileName, map<Str, set<BED3> >& chr_bed, double p = 0.01 );
void BEDFrmPIQ( Str fileName, map<Str, set<BED3> >& chr_bed, int len, double purity = 0.7 );
void BEDFrmDanposXlS( Str fileName, map<Str, set<BED3> >& chr_bed );
map<Str, set<BED3> > BEDFrmSICCERDIFF( Str fileName);
map<Str, set<BED3> > BEDFrmSICCERDIFF_2( Str fileName, double fdr, double FC);
map<Str, set<BED3> > BEDFrmFIMO( Str fileName);
map<Str, set<BED3> > BEDFrmtags_around_summit( Str fileName);

pair<map<int,double>,map<int,double> > tagDensityOnRegion
(map<Str, set<BED6> >& chr_tags, map<Str, set<BED4> >& chr_regions, int windowN);
//number of overlap between two sets of islands
std::pair<map<Str, vector<BED3> >,map<Str, vector<BED3> > > getOverlap(map<Str, vector<BED3> >,map<Str, vector<BED3> >);
map<Str, set<BED3> > getOverlap(map<Str, set<BED3> >,map<Str, set<BED3> >);
std::pair<map<Str, vector<BED3> >,map<Str, vector<BED3> > > getOverlap(Str fileName1, Str fileName2);
//islands number is small
bool isOverlap(Pa_I_I region, Str chrom, map<Str, set<BED3> >& islands);

map<Str, set<BED3> > bed42bed3( map<Str, set<BED4> > beds4);
map<Str, set<BED4> > bed32bed4( map<Str, set<BED3> > beds3);

map<Str, set<BED3> > mergeIslandsWithGapLessThan( map<Str, set<BED3> >& beds3, int gap);

map<Str, set<BED3> > getIsolatedIslands( map<Str, set<BED3> >& beds3, int gap);

map<Str, set<BED3> > generateRandomIslands( int island_N, int island_len, map<Str, int> chr_len);

map<Str, set<BED3> > generateRandomIslands( map<Str, set<BED3> > templa , map<Str, int> chr_len);

inline int getlen(map<Str, set<BED3> > beds3){
	int len = 0;
	for( map<Str, set<BED3> >::iterator chrIt = beds3.begin(); chrIt != beds3.end(); ++chrIt ){
		for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt){
			if( bIt->end > bIt->start){
				len += bIt->end - bIt->start;
			}else{
				cout<<"ERR in inline int getlen(map<Str, set<BED3> >& beds3)"<<endl;
				exit(1);
			}
		}
	}
	return len;
}
map<Str, set<BED3> > get_complement_set(map<Str, set<BED3> > pro, Str file );
map<Str, set<BED3> > get_intergenic(map<Str, set<BED3> > pro, map<Str, set<BED3> > genic, Str file );
map<Str, set<BED3> > self_union( map<Str, set<BED3> > beds );
map<Str, set<BED3> > merge(map<Str, set<BED3> > bed1, map<Str, set<BED3> > bed2);
map<Str, set<BED3> > rand_position_in(map<Str, set<BED3> > beds, int island_N, Str chr_len_file, int length );
map<Str, set<BED3> > filter_region_out_of_chr_range(map<Str, set<BED3> > bed, Str file );
BED3 parsStr2Bed3(Str s); //Str like "chr1:22333-44434";
map<Str, set<BED3> > BEDFrmEdgeR(Str fileName, double _FDR, double _log2FC);
map<Str, set<BED3> > PartitionGenomeIntoBins(Str chrlen, int binsz);

#endif /* BED_H_ */
