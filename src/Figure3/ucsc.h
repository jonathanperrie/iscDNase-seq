/*
Authors: Gang-Qing Hu & keji Zhao
Date created: Sep 9 2009
Date modified: -
Contact: Gangqing.Hu@nih.gov
Distributed under GUN GPL license

Description: this class is to deal with UCSC known genes
*/

#ifndef UCSC_H
#define UCSC_H

#include "TypeDefBase.h"
#include "OftenUsedOperatLib.h"
#include "bed.h"

template<class SortCreteria>
void UCSCGenesFrmFile(Str fileName, map<Str, set<UCSC, SortCreteria> >& chr_genes){
	std::ifstream ucscIn(fileName.data());
	if( !ucscIn.good() )
		FileNotFoundERR(fileName);
	Str skip;
	std::getline(ucscIn, skip);//skip the first descriptive line
	int N = 0;
	while( !ucscIn.eof() ){
		UCSC g;
		ucscIn>>g;
		if( !g.name.empty() //not the end of file
				//&& g.chrom.find("random") == std::string::npos// exclude chrom wiht random
				//&& names.find(g.name) == names.end()
				){//reduce redundancy from records with same gene names;
			g.island_score = ++N;
			chr_genes[g.chrom].insert(g);
		}
	}
	ucscIn.close();
	cout<<"# of records from \""<<fileName<<"\": "<<getSize(chr_genes)<<endl;
}

template<class SortCreteria>
void RefSeqGenesFrmFile(Str fileName, map<Str, set<UCSC, SortCreteria> >& chr_genes){
	std::ifstream ucscIn(fileName.data());
	if( !ucscIn.good() )
		FileNotFoundERR(fileName);
	Str skip;
	std::getline(ucscIn, skip);//skip the first descriptive line
	int N = 0;
	while( !ucscIn.eof() ){
		UCSC g;
		ucscIn>>skip>>g;
		if( !g.name.empty() //not the end of file
				//&& g.chrom.find("random") == std::string::npos// exclude chrom wiht random
				//&& names.find(g.name) == names.end()
				){//reduce redundancy from records with same gene names;
			g.island_score = ++N;
			chr_genes[g.chrom].insert(g);
		}
	}
	ucscIn.close();
	cout<<"# of records from \""<<fileName<<"\": "<<getSize(chr_genes)<<endl;
}

map<Str, set<UCSC, UCSC::sortByTxStart> > removeRedundantUCSC(map<Str, set<UCSC, UCSC::sortByTxStart> >& chr_genes);

map<Str, set<UCSC, UCSC::sortByTxStart> > removeRedundantUCSC
(map<Str, set<UCSC, UCSC::sortByTxStart> > chr_genes, Str isoform_ID_f );

struct JUNC{
	Str chrom;
	int left, left_o;
	int right, right_o;
	Str id, a, b, c, d, e, f;
	int left_span;
	int right_span;
	Str strand;
	int count;
};

map<Str, JUNC> get_junction(Str f);
map<Str, JUNC> intersct_junction(map<Str, JUNC> a, map<Str, JUNC> b);
map<Str, JUNC> union_junction(map<Str, JUNC> a, map<Str, JUNC> b);
void write_junction(map<Str, JUNC> a, Str f );

#endif
