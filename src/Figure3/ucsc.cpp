/*
 * UCSC.cpp
 *
 *  Created on: 2009-9-25
 *      Author: appleapple
 */
#include "ucsc.h"

std::istream &operator>>(std::istream &stream, UCSC& ucsc)
{
	Str line;
	std::getline(stream,line);
	std::istringstream buf(line);//ensure one line per gene

	Str exonStarts, exonEnds;
	buf>>ucsc.name>>ucsc.chrom>>ucsc.strand>>ucsc.txStart>>ucsc.txEnd>>ucsc.cdsStart>>ucsc.cdsEnd
	>>exonStarts>>exonStarts>>exonEnds>>ucsc.proteinID>>ucsc.alignID;
//	if( ucsc.alignID.empty())
	//	cout<<ucsc.name<<endl;
	vector<Str> exStarts = split( exonStarts, ',');
	vector<Str> exEnds = split( exonEnds, ',');
	assert( exStarts.size() == exEnds.size() );
	ucsc.TSS = ucsc.strand == "+" ? ucsc.txStart : ucsc.txEnd;
	ucsc.TSE = ucsc.strand == "+" ? ucsc.txEnd : ucsc.txStart;

	for(unsigned int i = 0; i < exStarts.size(); ++i){
		if(exStarts[i].empty() || exEnds[i].empty() )
			continue;
		Pa_I_I x;
		x.first = convertString2<int>(exStarts[i]);
		x.second = convertString2<int>(exEnds[i]);

		ucsc.exons.insert(x);
	}

	ucsc.relatives = ucsc.name;
	ucsc.ismasked = false;
	if( !ucsc.exons.empty() ){
		if( ucsc.strand == "+" ){
			for( set<Pa_I_I>::iterator eIt = ucsc.exons.begin(); eIt != ucsc.exons.end(); ++eIt ){
				ucsc.exons53.push_back(*eIt);
			}
			set<Pa_I_I>::iterator e1It = ucsc.exons.begin();
			set<Pa_I_I>::iterator e2It = e1It;
			++e2It;
			for( ; e2It != ucsc.exons.end(); ++e2It, ++e1It ){
				ucsc.intron53.push_back(make_pair(e1It->second, e2It->first));
			}
		}else{
			if( !ucsc.exons.empty() ){
				set<Pa_I_I>::iterator eIt = ucsc.exons.end();
				do{
					--eIt;
					ucsc.exons53.push_back(*eIt);
				}while( eIt != ucsc.exons.begin() );
			}
			if( ucsc.exons.size() > 1){
				set<Pa_I_I>::iterator e1It = ucsc.exons.end();
				set<Pa_I_I>::iterator e2It = e1It;
				--e2It;
				do{
					--e2It;
					--e1It;
					ucsc.intron53.push_back(make_pair(e2It->second, e1It->first));
				}while( e2It != ucsc.exons.begin() );
			}
		}

	}
	for( unsigned int i = 0; i < ucsc.intron53.size(); ++i ){
		UCSC::INTRON intron;
		intron.start = ucsc.intron53[i].first;
		intron.end = ucsc.intron53[i].second;
		intron.strand = ucsc.strand;
		intron.rcd == 0;
		ucsc.INTRON53.push_back(intron);
	}
	return stream;
}//*/

std::ostream &operator<<(std::ostream &stream, const UCSC& ucsc){
	stream<<ucsc.name<<'\t'<<ucsc.chrom<<'\t'<<ucsc.strand<<'\t'<<ucsc.txStart<<'\t'<<ucsc.txEnd
	<<'\t'<<ucsc.cdsStart<<'\t'<<ucsc.cdsEnd<<'\t'<<ucsc.exons.size()<<'\t';

	for(set<Pa_I_I>::iterator i = ucsc.exons.begin();
			i != ucsc.exons.end(); ++i){
		stream<<i->first<<",";
	}
	stream<<'\t';
	for(set<Pa_I_I>::iterator i = ucsc.exons.begin();
			i != ucsc.exons.end(); ++i){
		stream<<i->second<<",";
	}//
	stream<<'\t'<<ucsc.proteinID<<'\t'<<ucsc.alignID<<"\t"<<ucsc.relatives;
	return stream;
}//*/

map<Str, set<UCSC, UCSC::sortByTxStart> > removeRedundantUCSC(map<Str, set<UCSC, UCSC::sortByTxStart> >& chr_genes){
	int GN = 0;
	map<Str, set<UCSC, UCSC::sortByTxStart> > chr_genes_nolap;
	srand((unsigned)time(0));
	for( map<Str, set<UCSC, UCSC::sortByTxStart> >::iterator it = chr_genes.begin()
			; it != chr_genes.end(); ++it){
		vector<UCSC> genes;
		for( set<UCSC>::iterator i = it->second.begin(); i != it->second.end(); ++i){
			genes.push_back(*i);//sorted
		}
		for( vector<UCSC>::iterator cur = genes.begin(); cur != genes.end(); ++cur){
			if( !cur->ismasked){
				vector<UCSC>::iterator next = cur;
				for( ++next; next != genes.end() && next->txStart < cur->txEnd; ++next ){
					if( (!next->ismasked) /*&& (next->strand == cur->strand)*/ ){
						bool hasExonOverlap = false;
						for( set<Pa_I_I>::iterator it1 = cur->exons.begin ( );
								it1 != cur->exons.end ( ); ++it1 ){
							for( set<Pa_I_I>::iterator it2 = next->exons.begin ( );
								it2 != next->exons.end ( ); ++it2 ){
								if( _min(it1->second,it2->second) > _max(it1->first, it2 ->first) ){
									hasExonOverlap = true;
									break;
								}
							}
						}
						if( hasExonOverlap ){
							bool keepNex = (rand()%2 == 1);
							if( keepNex ){
								cur->ismasked = true;
								next->relatives += "," + cur->relatives;
								break;
							}
							else{
								next->ismasked = true;
								cur->relatives += "," + next->relatives;
							}
						}
					}
				}
				//if not mark as overlap genes after comparing with all following overlapping genes
				if( !cur->ismasked ){
					++GN;
					chr_genes_nolap[cur->chrom].insert(*cur);
				}
			}
		}
	}
//	cout<<"# of UCSC known genes (non-overlap): "<<GN<<endl;
	return chr_genes_nolap;
}


map<Str, set<UCSC, UCSC::sortByTxStart> > removeRedundantUCSC
(map<Str, set<UCSC, UCSC::sortByTxStart> > chr_genes, Str isoform_ID_f ){
	set<Str> ucsc_names;
	for( map<Str, set<UCSC, UCSC::sortByTxStart> >::iterator it = chr_genes.begin(); it != chr_genes.end(); ++it){
		for(set<UCSC, UCSC::sortByTxStart>::iterator gIt = it->second.begin(); gIt != it->second.end(); ++gIt){
			ucsc_names.insert(gIt->name);
		}
	}
	srand((unsigned)time(0));
	map<Str, Str> iso_id = getIsoform_ID(isoform_ID_f);
	map<Str, vector<Str> > id_isos;
	for( map<Str, Str>::iterator it = iso_id.begin(); it != iso_id.end(); ++it ){
		if(ucsc_names.find(it->first) != ucsc_names.end() )
			id_isos[it->second].push_back(it->first);
	}
	set<Str> names;
	for( map<Str, vector<Str> >::iterator it = id_isos.begin(); it != id_isos.end(); ++it ){
		int index = rand()%it->second.size();
		names.insert(it->second[index]);
	}

	map<Str, set<UCSC, UCSC::sortByTxStart> > chr_genes_nolap;
	for( map<Str, set<UCSC, UCSC::sortByTxStart> >::iterator it = chr_genes.begin(); it != chr_genes.end(); ++it){
		for(set<UCSC, UCSC::sortByTxStart>::iterator gIt = it->second.begin(); gIt != it->second.end(); ++gIt){
			if( names.find(gIt->name) != names.end() ){
				chr_genes_nolap[gIt->chrom].insert(*gIt);
			}
		}
	}
	return chr_genes_nolap;
}

map<Str, JUNC> get_junction(Str f){
	Str tmp;
	std::ifstream in(f.data());
	cout<<f<<endl;
	if( !in.good()){
		cout<<f<<"\t not found"<<endl;
		exit(1);
	}
	map<Str, JUNC> id_junc;
	while( !in.eof() ){
		JUNC junc;
		in>>junc.chrom>>junc.left_o>>junc.right_o>>junc.id>>junc.count>>junc.strand>>junc.a>>junc.b>>junc.c>>junc.d>>junc.e;
		getline( in, junc.f );
		if( !junc.chrom.empty() ){
			junc.left_span = atoi(junc.e.substr(0, junc.e.find(",")).data());
			junc.right_span = atoi(junc.e.substr(junc.e.find(",")+1).data());
			if( junc.strand == "+") {
				junc.left = junc.left_o + junc.left_span + 1;
				junc.right = junc.right_o - junc.right_span;
			}else{
				junc.left = junc.left_o + junc.left_span + 1;
				junc.right = junc.right_o - junc.right_span;
			}
			Str id = junc.chrom + "|"+convert2String(junc.left)
			        +"|"+convert2String(junc.right)+"|"+junc.strand
			        +"_"+convert2String(junc.left_span)
			        +"|"+convert2String(junc.right_span);
			if( id_junc.find(id) == id_junc.end() )
				id_junc[id] = junc;
			else{
				cout<<"ERR"<<endl;
			}
		}
	}
	return id_junc;
}

map<Str, JUNC> intersct_junction(map<Str, JUNC> a, map<Str, JUNC> b){
	map<Str, JUNC> inter;
	for( map<Str, JUNC>::iterator i = a.begin(); i != a.end(); ++i ){
		map<Str, JUNC>::iterator j = b.find(i->first);
		if( j != b.end() ){
			JUNC junc = i->second;
			junc.left_span = _max(junc.left_span, j->second.left_span);
			junc.right_span = _max(junc.right_span, j->second.right_span);
			junc.count = junc.count + j->second.count;
			inter[i->first] = junc;
		}
	}
	cout<<"from "<<a.size()<<"  "<<b.size()<<" to "<<inter.size()<<endl;
	return inter;
}

map<Str, JUNC> union_junction(map<Str, JUNC> a, map<Str, JUNC> b){
	map<Str, JUNC> _union;
	for( map<Str, JUNC>::iterator i = a.begin(); i != a.end(); ++i ){
		map<Str, JUNC>::iterator j = b.find(i->first);
		JUNC junc = i->second;
		if( j != b.end() ){
			junc.left_span = _max(junc.left_span, j->second.left_span);
			junc.right_span = _max(junc.right_span, j->second.right_span);
			junc.count = junc.count + j->second.count;
		}
		_union[i->first] = junc;
	}
	for( map<Str, JUNC>::iterator i = b.begin(); i != b.end(); ++i ){
		map<Str, JUNC>::iterator j = a.find(i->first);
		JUNC junc = i->second;
		if( j == a.end() ){
			_union[i->first] = junc;
		}
	}
	cout<<"from "<<a.size()<<"  "<<b.size()<<" to "<<_union.size()<<endl;
	return _union;
}

void write_junction(map<Str, JUNC> a, Str f ){
	std::ofstream out(f.data());
	int n = 0;
	for( map<Str, JUNC>::iterator i = a.begin(); i != a.end(); ++i ){
		out<<i->second.chrom
			<<"\t"<<(i->second.left - i->second.left_span)
			<<"\t"<<(i->second.right + i->second.right_span)
			<<"\tJC_"<<(++n)<<"_"<<i->second.count
			<<"_"<<i->second.left_span<<"_"<<i->second.right_span
			<<"\t"<<i->second.count<<"\t"<<i->second.strand
			<<"\t"<<(i->second.left - i->second.left_span)
			<<"\t"<<(i->second.right + i->second.right_span)
			<<"\t255,0,0 2"
			<<"\t"<<i->second.left_span
			<<","<<i->second.right_span
			<<"\t0,"<<(i->second.right - i->second.left - i->second.right_span)<<endl;//*/
	}
	out.close();
}
