 /*
 * UCSCBEDOperator.cpp
 *
 *  Created on: Sep 28, 2009
 *      Author: hugq
 */

#include "UCSCBEDOperator.h"
#include "OftenUsedOperatLib.h"

void assignIsland2UCSC(map<Str, vector<UCSC> >& chr_genes,
		map<Str, set<BED3> >& chr_islands, GenomePartIdentify regionId, Str islandLabel, bool FCFS){
	int N = 0;
	for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
		for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
			Pa_I_I intrestRegion = getIntrestRegion(*gIt, regionId);

			//see if interested region overlap a island
			set<BED3>& islands = chr_islands[gIt->chrom];
			if( intrestRegion.second - intrestRegion.first > MININTRESTREGIONLEN ){
				if( !islands.empty() ){
					set<BED3>::iterator islandIt = islands.begin();
					for( ; islandIt != islands.end(); ++islandIt){
						if( _min(islandIt->end, intrestRegion.second) > _max(islandIt->start, intrestRegion.first)){
							gIt->genepart_boundby[regionId].insert(islandLabel);
							gIt->islandsOn[regionId][islandLabel].insert(*islandIt);
							++N;
							if(FCFS)
								islands.erase(islandIt);
						}
					}
				}
			}
		}
	}//*/
//	cout<<"# islands assigned: "<<N<<endl;
}

void assignIsland2Island(map<Str, set<BED3> >& target_island,
		const map<Str, set<BED3> >& source_island, Str source_island_label){
	map<Str, vector<BED3> > target_islands = set2vector(target_island);
	map<Str, vector<BED3> > source_islands = set2vector(source_island);
	for( map<Str, vector<BED3> >::iterator chrIt = target_islands.begin(); chrIt != target_islands.end(); ++chrIt ){
		if( !source_islands[chrIt->first].empty() ){
			vector<BED3>::const_iterator sIt = source_islands[chrIt->first].begin();
			for( vector<BED3>::iterator tIt = chrIt->second.begin(); tIt != chrIt->second.end(); ++tIt ){
//				if( tIt->start > 123962000 && tIt->end < 123993820 && tIt->chrom == "chr9" ){
//					cout<<tIt->start<<'\t'<<tIt->end<<endl;
//				}
				while( (sIt->end > tIt->start && sIt != source_islands[chrIt->first].begin())
						|| sIt == source_islands[chrIt->first].end() ){
					--sIt;
				}
				for( ; sIt != source_islands[chrIt->first].end() && (sIt->start < tIt->end); ++sIt){
					if( _min(tIt->end, sIt->end) > _max(tIt->start, sIt->start)){
						tIt->islandsOn[source_island_label].insert(*sIt);
					}
				}
			}
		}
	}
	target_island = vector2set(target_islands);//*/
}

void setReadNumberBoundUCSCOn(map<Str, vector<UCSC> >& chr_genes,
	map<Str, set<BED6, BED6::sortByShifted2Point> > chr_reads,
	GenomePartIdentify regionId, Str islandLabel){
	for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
		set<BED6, BED6::sortByShifted2Point>::const_iterator curRd = chr_reads[chrIt->first].begin();
		for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
			Pa_I_I region = getIntrestRegion(*gIt, regionId);
			if( curRd->shifted2Point >= region.first ){
				while( curRd->shifted2Point >= region.first && curRd != chr_reads[chrIt->first].begin())
					--curRd;
			}
			while( curRd->shifted2Point < region.second && curRd != chr_reads[chrIt->first].end()){
				if( curRd->shifted2Point > region.first ){
					++gIt->genepart_islandLabel_readNumber[regionId][islandLabel];
				}
				++curRd;
			}
		}
	}
}

Pa_I_I getIntrestRegion(const UCSC& ucsc, GenomePartIdentify regionId){
	Pa_I_I intrestRegion;
	if( regionId == PROMOTER )
		return (ucsc.strand == "+")? make_pair(ucsc.TSS - PROMOTERUPLEN, ucsc.TSS + PROMOTERDOWNLEN )
				: make_pair(ucsc.TSS - PROMOTERDOWNLEN, ucsc.TSS + PROMOTERUPLEN ); //for promoter region
	else if( regionId == GENEBODY )
		return ((ucsc.strand == "+")? make_pair(ucsc.txStart + PROMOTERDOWNLEN, ucsc.txEnd )
				: make_pair(ucsc.txStart, ucsc.txEnd - PROMOTERDOWNLEN )); //for gene body
	else if( regionId == TSSTES){
		return make_pair(ucsc.txStart, ucsc.txEnd);
	}else{
		cout<<"regionId should be PROMOTER OR GENEBODY only for this functions\n";
		make_pair(0,0);
	}
	return make_pair(0,0);
}

//note the algorithm will ignore intergenic region at ends of chromosom
map<Str, set<BED4> > getIntrestReions(const map<Str, set<UCSC, UCSC::sortByTxStart> >& chr_genes_tmp, GenomePartIdentify regionId){
	map<Str, vector<UCSC> > chr_genes = set2vector(chr_genes_tmp);
	map<Str, set<BED4> > regions;
	for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin();
			chrIt != chr_genes.end(); ++chrIt){
		vector<UCSC>::iterator gIt = chrIt->second.begin();
		for( ; gIt != chrIt->second.end(); ++gIt  ){
			BED4 region;
			region.chrom = chrIt->first;
			region.info = gIt->name;
			if( regionId == PROMOTER){
				region.strand = gIt->strand;
				if( region.chrom == "+"){
					region.start = gIt->TSS - PROMOTERUPLEN;
					region.end = gIt->TSS + PROMOTERDOWNLEN;
				}else{
					region.start = gIt->TSS - PROMOTERDOWNLEN;
					region.end = gIt->TSS + PROMOTERUPLEN;
				}
				regions[region.chrom].insert(region);
			}
			if( regionId == GENEBODY){
				region.strand = gIt->strand;
				if( gIt->strand == "+"){
					region.start = gIt->TSS+ PROMOTERDOWNLEN+1;
					region.end = gIt->txEnd;
				}else{
					region.start = gIt->txStart;
					region.end = gIt->TSS - PROMOTERDOWNLEN-1;
				}
				if(region.end>region.start){
					regions[region.chrom].insert(region);
				}
			}
		}
	}
	return regions;
}

vector<pair<double,int> > exp_islandPercent(map<Str, vector<UCSC> > chr_genes_nolap,
		map<Str,double> name_rpkm, GenomePartIdentify regionId, Str islandLabel, Str filename, int step){
	vector<pair<double,int> > exp_hasIslandOrNot;
	for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes_nolap.begin(); chrIt != chr_genes_nolap.end(); ++chrIt ){
		for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
			//if( (count( gIt->relatives.begin(), gIt->relatives.end(), ',') == 0) )//if the gene is 1 isform 1 gene
			{
				map<Str,double>::iterator exp = name_rpkm.find(gIt->name);//has a rpkm expression value
				if( exp != name_rpkm.end() ){
					if( gIt->genepart_boundby[regionId].find(islandLabel) != gIt->genepart_boundby[regionId].end())
						exp_hasIslandOrNot.push_back(make_pair(exp->second,1));
					else
						exp_hasIslandOrNot.push_back(make_pair(exp->second,0));
				}//
			}
		}
	}//*/
	std::ofstream out(filename.data());
	out<<"\texp\tp\tcat\n";
	int NN = 0;
	double inactiveGN = 0, inactiveIN = 0;
	sort(exp_hasIslandOrNot.begin(), exp_hasIslandOrNot.end());
	for(unsigned int i = 0; i <  exp_hasIslandOrNot.size();){
		if(  exp_hasIslandOrNot[i].first >0 ){
			double sum = 0;
			double avexp = 0;
			double eN = 0;
			for( int j = 0; j < step; ++j ){
				if( i +j < exp_hasIslandOrNot.size() ){
					avexp +=  exp_hasIslandOrNot[i+j].first;
					sum +=  exp_hasIslandOrNot[i+j].second;
					++eN;
				}
			}
			out<<(++NN)<<'\t'<<(avexp/step)<<'\t'<<100*sum/eN<<'\t'
			<<(islandLabel+genomePartIdentify2Str(regionId))<<endl;
			i += step;
		}else{
			++i;
			++inactiveGN;
			inactiveIN +=  exp_hasIslandOrNot[i].second;
		}
	}//*/
	if( inactiveGN >0 )
		out<<(++NN)<<'\t'<<0<<'\t'<<100*inactiveIN/inactiveGN<<'\t'
		<<(islandLabel+genomePartIdentify2Str(regionId))<<endl;
	out.close();
	return exp_hasIslandOrNot;
}

void assignGenomicRegions2Islands(
		map<Str, set<BED4> >& chr_regions,
		const map<Str, set<BED3> >& chr_islands,
		map<Str, map<BED3,set<GenomePartIdentify> > >& island_regions,
		GenomePartIdentify regionId ){
	for( map<Str, set<BED3> >::const_iterator chrIt = chr_islands.begin();
			chrIt != chr_islands.end(); ++chrIt ){
		set<BED4>& regions = chr_regions[chrIt->first];
		for( set<BED3>::const_iterator irIt = chrIt->second.begin();
				irIt != chrIt->second.end(); ++irIt){
			for(set<BED4>::const_iterator rgIt = regions.begin(); rgIt != regions.end(); ++rgIt ){
				if(_min(rgIt->end, irIt->end) > _max(rgIt->start, irIt->start )){
					island_regions[chrIt->first][*irIt].insert(regionId);
					break;
				}
			}
		}
	}
}

void generateReadDensityProfile(const map<Str, vector<UCSC> >& chr_genes_tmp, Str readf, Str outf, int shift){
	//windows
	int promoterWindow = 500;
	int endWindow = 500;
	//number of windows
	int promoterSliceN = (PROMOTERUPLEN+PROMOTERDOWNLEN)/promoterWindow;
	int endSliceN = PROMOTERUPLEN/endWindow;
	int genebodySliceN = 10;
	//densities
	vector<double> prom_density(promoterSliceN, 0);
	vector<double> body_density(genebodySliceN, 0);
	vector<double> end_density(endSliceN, 0);
	map<Str, vector<UCSC> > chr_genes;
	double GN = 0;
	for( map<Str, vector<UCSC> >::const_iterator chrIt = chr_genes_tmp.begin(); chrIt != chr_genes_tmp.end(); ++chrIt ){
		for( vector<UCSC>::const_iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
			if(gIt->txStart + PROMOTERUPLEN + 100 < gIt->txEnd){
				UCSC ucsc = *gIt;
				ucsc.txStart -= PROMOTERUPLEN;
				ucsc.txEnd += PROMOTERUPLEN;
				chr_genes[chrIt->first].push_back(ucsc);
				++GN;
			}
		}
	}
	cout<<"UCSC gene with gene body > 100 bps:"<<getSize(chr_genes)<<endl;
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end, n= 0, rn = 0;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty()){
			++rn;
			int pos = strand =="+" ? beg + shift : end - shift;
			chr_reads[chrom].push_back(pos);
			if( (++n) % 5000000 == 0 || in.eof()){
				cout<<"# reads in: "<<n<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
						if( *curRd >= gIt->txStart){
							while( *curRd >= gIt->txStart && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						int bodyWindow = (gIt->txEnd - gIt->txStart + 1 - 2*PROMOTERUPLEN - PROMOTERDOWNLEN)/genebodySliceN;
						while( *curRd < gIt->txEnd && curRd != chr_reads[chrIt->first].end()){
							if( *curRd > gIt->txStart ){
								if( gIt->strand == "+"){
									if( *curRd < gIt->txStart + PROMOTERUPLEN + PROMOTERDOWNLEN){
										unsigned int index = (*curRd - gIt->txStart)/promoterWindow;
										if( index < prom_density.size())
											prom_density[index] += 1./promoterWindow;
									}else if(*curRd < gIt->txEnd - PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart-PROMOTERUPLEN-PROMOTERDOWNLEN)/bodyWindow;
										if(index<body_density.size())
											body_density[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (*curRd - gIt->txEnd + PROMOTERUPLEN)/endWindow;
										if( index < end_density.size() )
											end_density[index] += 1./endWindow;
									}
								}else{
									if( *curRd > gIt->txEnd - PROMOTERUPLEN - PROMOTERDOWNLEN){
										unsigned int index = (gIt->txEnd-*curRd)/promoterWindow;
										if( index < prom_density.size())
											prom_density[index] += 1./promoterWindow;
									}else if( *curRd > gIt->txStart + PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd - PROMOTERUPLEN - PROMOTERDOWNLEN-*curRd)/bodyWindow;
										if( index < body_density.size() )
											body_density[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (gIt->txStart+PROMOTERUPLEN-*curRd)/endWindow;
										if( index < end_density.size())
											end_density[index] += 1./endWindow;
									}
								}
							}
							++curRd;
						}
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();

	std::ofstream out(outf.data());
	int NN = 0, i =0;
	out<<"\tn\tden\tgrp"<<endl;
	for( vector<double>::iterator it = prom_density.begin(); it != prom_density.end(); ++it, ++i ){
		out<<(++NN)<<'\t'<<(i-PROMOTERUPLEN/promoterWindow)*promoterWindow<<'\t'<<(*it)*1000000./rn/GN<<"\tpromoter"<<endl;
	}
	int l = (i-PROMOTERUPLEN/promoterWindow + 4)*promoterWindow;
	i = 0;
	for( vector<double>::iterator it = body_density.begin(); it != body_density.end(); ++it, ++i ){
		out<<(++NN)<<'\t'<<l+i*promoterWindow*2<<'\t'<<(*it)*1000000./rn/GN<<"\tgenebody"<<endl;
	}
	l = l+i*promoterWindow*2;
	i = 0;
	for( vector<double>::iterator it = end_density.begin(); it != end_density.end(); ++it, ++i ){
		out<<(++NN)<<'\t'<<(l+i*endWindow)<<'\t'<<(*it)*1000000./rn/GN<<"\tafterTES"<<endl;
	}
	out.close();//*/
}

void generateReadDensityProfile(const map<Str, vector<UCSC> >& chr_genes_tmp, Str readf, Str outf, bool is5Pr){
	//windows
	int promoterWindow = 100;
	int endWindow = 100;
	//number of windows
	int promoterSliceN = PROMOTERUPLEN*2/promoterWindow;
	int endSliceN = PROMOTERUPLEN/endWindow;
	int genebodySliceN = 10;
	//densities
	vector<double> prom_density(promoterSliceN, 0);
	vector<double> body_density(genebodySliceN, 0);
	vector<double> end_density(endSliceN, 0);
	map<Str, vector<UCSC> > chr_genes;
	double GN = 0;
	for( map<Str, vector<UCSC> >::const_iterator chrIt = chr_genes_tmp.begin(); chrIt != chr_genes_tmp.end(); ++chrIt ){
		for( vector<UCSC>::const_iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
			if(gIt->txStart + PROMOTERUPLEN + 1000 < gIt->txEnd)
			{
				UCSC ucsc = *gIt;
				ucsc.txStart -= PROMOTERUPLEN;
				ucsc.txEnd += PROMOTERUPLEN;
				if(ucsc.strand == "+"){
					chr_genes[chrIt->first].push_back(ucsc);
					++GN;
				}
			}
		}
	}
	cout<<"UCSC gene with gene body > 1000 bps:"<<getSize(chr_genes)<<endl;
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end, n= 0, rn = 0;
	int shift = 0;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty() ){
			int pos = strand =="+" ? beg + shift : end - shift;
			if(strand == (is5Pr ? "+":"-")){
				++rn;
				chr_reads[chrom].push_back(pos);
			}
			if( (++n) % 5000000 == 0 || in.eof()){
				cout<<"# reads in: "<<n<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
						if( *curRd >= gIt->txStart){
							while( *curRd >= gIt->txStart && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						int bodyWindow = (gIt->txEnd - gIt->txStart + 1 -3*PROMOTERUPLEN)/genebodySliceN;
						while( *curRd < gIt->txEnd && curRd != chr_reads[chrIt->first].end()){
							if( *curRd > gIt->txStart ){
								if( gIt->strand == "+"){
									if( *curRd < gIt->txStart + 2*PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart)/promoterWindow;
										if( index < prom_density.size())
											prom_density[index] += 1./promoterWindow;
									}else if(*curRd < gIt->txEnd - PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart-2*PROMOTERUPLEN)/bodyWindow;
										if(index<body_density.size())
											body_density[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (*curRd - gIt->txEnd + PROMOTERUPLEN)/endWindow;
										if( index < end_density.size() )
											end_density[index] += 1./endWindow;
									}
								}else{
									if( *curRd > gIt->txEnd - 2*PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd-*curRd)/promoterWindow;
										if( index < prom_density.size())
											prom_density[index] += 1./promoterWindow;
									}else if( *curRd > gIt->txStart + PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd - 2*PROMOTERUPLEN-*curRd)/bodyWindow;
										if( index < body_density.size() )
											body_density[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (gIt->txStart+PROMOTERUPLEN-*curRd)/endWindow;
										if( index < end_density.size())
											end_density[index] += 1./endWindow;
									}
								}
							}
							++curRd;
						}
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();

	std::ofstream out(outf.data());
	int NN = 0, i =0;
	out<<"\tn\tden\tgrp"<<endl;
	for( vector<double>::iterator it = prom_density.begin(); it != prom_density.end(); ++it, ++i ){
		out<<(++NN)<<'\t'<<(i-promoterSliceN/2)*promoterWindow<<'\t'<<(*it)*1000000./rn/GN<<"\tpromoter"<<endl;
	}
	int l = (i-promoterSliceN/2)*promoterWindow;
	i = 0;
	for( vector<double>::iterator it = body_density.begin(); it != body_density.end(); ++it, ++i ){
		out<<(++NN)<<'\t'<<(l+(i+0.5)*PROMOTERUPLEN/2)<<'\t'<<(*it)*1000000./rn/GN<<"\tgenebody"<<endl;
	}
	l = (l+(i)*PROMOTERUPLEN/2);
	i = 0;
	for( vector<double>::iterator it = end_density.begin(); it != end_density.end(); ++it, ++i ){
		out<<(++NN)<<'\t'<<(l+i*endWindow)<<'\t'<<(*it)*1000000./rn/GN<<"\tafterTES"<<endl;
	}
	out.close();//*/
}


void generateTSSTESReadDensityHeatMap(const map<Str, vector<UCSC> >& chr_genes_tmp, Str readf, Str out_f, int shift, int bins){
	//windows
	int promoterWindow = 200;
	int endWindow = 200;
	//number of windows
	int promoterSliceN = PROMOTERUPLEN*2/promoterWindow;
	int endSliceN = PROMOTERUPLEN/endWindow;
	int genebodySliceN = 10;
	//densities
	vector<double> prom_density(promoterSliceN, 0);
	vector<double> body_density(genebodySliceN, 0);
	vector<double> end_density(endSliceN, 0);
	map<Str, vector<UCSC> > chr_genes;
	for( map<Str, vector<UCSC> >::const_iterator chrIt = chr_genes_tmp.begin(); chrIt != chr_genes_tmp.end(); ++chrIt ){
		if( chrIt->first.find("_") == std::string::npos ){
			for( vector<UCSC>::const_iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
				if(gIt->txStart + PROMOTERUPLEN + 1000 < gIt->txEnd)
				{
					UCSC ucsc = *gIt;
					ucsc.txStart -= PROMOTERUPLEN;
					ucsc.txEnd += PROMOTERUPLEN;
					chr_genes[chrIt->first].push_back(ucsc);
				}
			}
		}
	}
	cout<<"# Genes with gene body > 3000 bps:"<<getSize(chr_genes)<<endl;
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end, n= 0, rn = 0;
	vector<pair<double, pair<UCSC, vector<double> > > > indx_bed_den;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty() ){
			int pos = strand =="+" ? beg + shift : end - shift;
			++rn;
			chr_reads[chrom].push_back(pos);
			if( (++n) % 50000000 == 0 || in.eof()){
				//cout<<"# reads in: "<<n<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
						if( *curRd >= gIt->txStart){
							while( *curRd >= gIt->txStart && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						int bodyWindow = (gIt->txEnd - gIt->txStart + 1 -3*PROMOTERUPLEN)/genebodySliceN;
						//cout<<bodyWindow<<endl;
						//exit(1);
						vector<double> pro(promoterSliceN,0.000001);
						vector<double> gb(genebodySliceN,0.000001);
						vector<double> end(endSliceN,0.000001);
						while( *curRd < gIt->txEnd && curRd != chr_reads[chrIt->first].end()){
							if( *curRd > gIt->txStart ){
								if( gIt->strand == "+"){
									if( *curRd < gIt->txStart + 2*PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart)/promoterWindow;
										if( index < prom_density.size())
											pro[index] += 1./promoterWindow;
									}else if(*curRd < gIt->txEnd - PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart-2*PROMOTERUPLEN)/bodyWindow;
									//	cout<<bodyWindow<<endl;exit(1);
										if(index<body_density.size())
											gb[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (*curRd - gIt->txEnd + PROMOTERUPLEN)/endWindow;
										if( index < end_density.size() )
											end[index] += 1./endWindow;
									}//*/
								}else{
									if( *curRd > gIt->txEnd - 2*PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd-*curRd)/promoterWindow;
										if( index < prom_density.size())
											pro[index] += 1./promoterWindow;
									}else if( *curRd > gIt->txStart + PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd - 2*PROMOTERUPLEN-*curRd)/bodyWindow;
										if( index < body_density.size() )
											gb[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (gIt->txStart+PROMOTERUPLEN-*curRd)/endWindow;
										if( index < end_density.size())
											end[index] += 1./endWindow;
									}//*/
								}
							}
							++curRd;
						}
						vector<double> den = pro;
						for( int i = 0; i < gb.size(); ++i)
							den.push_back(gb[i]);
						for( int i = 0; i < end.size(); ++i)
							den.push_back(end[i]);
						indx_bed_den.push_back(make_pair(gIt->island_score, make_pair(*gIt,den)));
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();
	int lib_sz = n;
//	int bins =200;
    cout<<"# reads : "<<lib_sz<<endl;
//    cout<<"gene sz: "<<indx_bed_den.size()<<endl;

    sort( indx_bed_den.begin(), indx_bed_den.end() );

    std::ofstream out(out_f.data());
    for( int i = 0; i < indx_bed_den[0].second.second.size(); ++i ){
    	if( i < promoterSliceN )
    		out<<"\tPro_"<<(-promoterSliceN/2+i)*promoterWindow+promoterWindow/2<<"_(bp)";
    	else if( i < promoterSliceN + genebodySliceN)
    		out<<"\tGenBD_"<<i-promoterSliceN+1;
    	else
		  out<<"\tPosTES_"<<(i-promoterSliceN - genebodySliceN)*endWindow+endWindow/2<<"_(bp)";
    }
    out<<endl;

    int step = indx_bed_den.size()/bins > 1 ?
                    indx_bed_den.size()/bins : 1;
    cout<<"# genes per group"<< step<<endl;
    //      int max = step == 1 ? 1 :  indx_bed_den.size() - step/2;
    for( int i = 0; i < indx_bed_den.size() - step; i += step ){
            vector<double> sum(indx_bed_den[i].second.second.size(), 0);
            for( int j = 0; j < step; ++j ){
                    vector<double> den = indx_bed_den[i+j].second.second;
                    for( int k = 0; k < den.size(); ++k ){
//                            sum[k] += log((den[k])*100000000./lib_sz)/log(2);
                            sum[k] += ((den[k])*1000000./lib_sz);
                    }
            }
            if (step == 1 )
                    out<<indx_bed_den[i].second.first.chrom<<":"
                            <<indx_bed_den[i].second.first.TSS-PROMOTERUPLEN<<"-"
                            <<indx_bed_den[i].second.first.TSS+PROMOTERUPLEN
                            <<"-"<<indx_bed_den[i].second.first.strand;
            for( int k = 0; k < sum.size(); ++k )
                    out<<"\t"<<(sum[k]/step);
            out<<endl;
    }
    out.close();
    cout<<"File saved in "<<out_f<<endl;
}

void generateTSSTESReadDensity4Profiles(const map<Str, vector<UCSC> >& chr_genes_tmp, Str readf, Str out_f, int shift){
	//windows
	int promoterWindow = 100;
	int endWindow = 100;
	//number of windows
	int promoterSliceN = PROMOTERUPLEN*2/promoterWindow;
	int endSliceN = PROMOTERUPLEN*2/endWindow;
	int genebodySliceN = 100;
	//densities
	vector<double> prom_density(promoterSliceN, 0);
	vector<double> body_density(genebodySliceN, 0);
	vector<double> end_density(endSliceN, 0);
	map<Str, vector<UCSC> > chr_genes;
	for( map<Str, vector<UCSC> >::const_iterator chrIt = chr_genes_tmp.begin(); chrIt != chr_genes_tmp.end(); ++chrIt ){
		for( vector<UCSC>::const_iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
			if(gIt->txStart + 2*PROMOTERUPLEN + 5000 < gIt->txEnd)
			{
				UCSC ucsc = *gIt;
				ucsc.txStart -= PROMOTERUPLEN;
				ucsc.txEnd += PROMOTERUPLEN;
				chr_genes[chrIt->first].push_back(ucsc);
			}
		}
	}
	cout<<"UCSC gene with gene body > 5000 bps:"<<getSize(chr_genes)<<endl;
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end, n= 0, rn = 0;
	vector<pair<double, pair<UCSC, vector<double> > > > indx_bed_den;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty() ){
			int pos = strand =="+" ? beg + shift : end - shift;
			++rn;
			chr_reads[chrom].push_back(pos);
			if( (++n) % 50000000 == 0 || in.eof()){
				cout<<"# reads in: "<<n<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
						if( *curRd >= gIt->txStart){
							while( *curRd >= gIt->txStart && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						int bodyWindow = (gIt->txEnd - gIt->txStart + 1 -4*PROMOTERUPLEN)/genebodySliceN;
						//cout<<bodyWindow<<endl;
						//exit(1);
						vector<double> pro(promoterSliceN,0.000001);
						vector<double> gb(genebodySliceN,0.000001);
						vector<double> end(endSliceN,0.000001);
						while( *curRd < gIt->txEnd && curRd != chr_reads[chrIt->first].end()){
							if( *curRd > gIt->txStart ){
								if( gIt->strand == "+"){
									if( *curRd < gIt->txStart + 2*PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart)/promoterWindow;
										if( index < prom_density.size())
											pro[index] += 1./promoterWindow;
									}else if(*curRd < gIt->txEnd - 2*PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart-2*PROMOTERUPLEN)/bodyWindow;
									//	cout<<bodyWindow<<endl;exit(1);
										if(index<body_density.size())
											gb[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (*curRd - gIt->txEnd + 2*PROMOTERUPLEN)/endWindow;
										if( index < end_density.size() )
											end[index] += 1./endWindow;
									}//*/
								}else{
									if( *curRd > gIt->txEnd - 2*PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd-*curRd)/promoterWindow;
										if( index < prom_density.size())
											pro[index] += 1./promoterWindow;
									}else if( *curRd > gIt->txStart + 2*PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd - 2*PROMOTERUPLEN-*curRd)/bodyWindow;
										if( index < body_density.size() )
											gb[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (gIt->txStart+2*PROMOTERUPLEN-*curRd)/endWindow;
										if( index < end_density.size())
											end[index] += 1./endWindow;
									}//*/
								}
							}
							++curRd;
						}
						vector<double> den = pro;
						for( int i = 0; i < gb.size(); ++i)
							den.push_back(gb[i]);
						for( int i = 0; i < end.size(); ++i)
							den.push_back(end[i]);
						indx_bed_den.push_back(make_pair(gIt->island_score, make_pair(*gIt,den)));
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();
	int lib_sz = n;
	int bins =4;
    cout<<"lib sz: "<<lib_sz<<endl;
    cout<<"gene sz: "<<indx_bed_den.size()<<endl;

    sort( indx_bed_den.begin(), indx_bed_den.end() );

    int step = indx_bed_den.size()/bins > 1 ?
                    indx_bed_den.size()/bins : 1;
    cout<<"per bin"<< step<<endl;
    //      int max = step == 1 ? 1 :  indx_bed_den.size() - step/2;
    vector<vector<double> > aves;
    for( int i = 0; i <= indx_bed_den.size() - step; i += step ){
            vector<double> sum(indx_bed_den[i].second.second.size(), 0);
            for( int j = 0; j < step; ++j ){
                    vector<double> den = indx_bed_den[i+j].second.second;
                    for( int k = 0; k < den.size(); ++k ){
//                            sum[k] += log((den[k])*100000000./lib_sz)/log(2);
                            sum[k] += ((den[k])*100000000./lib_sz);
                    }
            }
	    vector<double> tmp;
		for( int k = 0; k < sum.size(); ++k )
                    tmp.push_back((sum[k]/step));
            aves.push_back(tmp);
    }

    std::ofstream out(out_f.data());
  for( int i = 0; i < aves[0].size(); ++i ){
	out<<(i+1);
	for( int j = 0; j < aves.size(); ++j )
		out<<"\t"<<aves[j][i];
	out<<endl;
    }
    out.close();
}


void generateTSSTESReadDensity4Profiles_Low(const map<Str, vector<UCSC> >& chr_genes_tmp, Str readf, Str out_f, int shift, int bins ){
	//windows
	int promoterWindow = 200;
	int endWindow = 200;
	//number of windows
	int promoterSliceN = PROMOTERUPLEN*2/promoterWindow;
	int endSliceN = PROMOTERUPLEN/endWindow;
	int genebodySliceN = 10;
	//densities
	vector<double> prom_density(promoterSliceN, 0);
	vector<double> body_density(genebodySliceN, 0);
	vector<double> end_density(endSliceN, 0);
	map<Str, vector<UCSC> > chr_genes;
	for( map<Str, vector<UCSC> >::const_iterator chrIt = chr_genes_tmp.begin(); chrIt != chr_genes_tmp.end(); ++chrIt ){
		if( chrIt->first.find("_") == std::string::npos ){
			for( vector<UCSC>::const_iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
				if(gIt->txStart + PROMOTERUPLEN + 1000 < gIt->txEnd)
				{
					UCSC ucsc = *gIt;
					ucsc.txStart -= PROMOTERUPLEN;
					ucsc.txEnd += PROMOTERUPLEN;
					chr_genes[chrIt->first].push_back(ucsc);
				}
			}
		}
	}
	cout<<"# genes with gene body > 3K bps:"<<getSize(chr_genes)<<endl;
	cout<<"\tOther genes are excluded from data analysis."<<getSize(chr_genes)<<endl;
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end, n= 0, rn = 0;
	vector<pair<double, pair<UCSC, vector<double> > > > indx_bed_den;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty() ){
			int pos = strand =="+" ? beg + shift : end - shift;
			++rn;
			chr_reads[chrom].push_back(pos);
			if( (++n) % 50000000 == 0 || in.eof()){
				//cout<<"# read in: "<<n<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
					//cout<<chrIt->first<<endl;
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
						if( *curRd >= gIt->txStart){
							while( *curRd >= gIt->txStart && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						int bodyWindow = (gIt->txEnd - gIt->txStart + 1 -3*PROMOTERUPLEN)/genebodySliceN;
						//cout<<bodyWindow<<endl;
						//exit(1);
						vector<double> pro(promoterSliceN,0.000001);
						vector<double> gb(genebodySliceN,0.000001);
						vector<double> end(endSliceN,0.000001);
						while( *curRd < gIt->txEnd && curRd != chr_reads[chrIt->first].end()){
							if( *curRd > gIt->txStart ){
								if( gIt->strand == "+"){
									if( *curRd < gIt->txStart + 2*PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart)/promoterWindow;
										if( index < prom_density.size())
											pro[index] += 1./promoterWindow;
									}else if(*curRd < gIt->txEnd - PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart-2*PROMOTERUPLEN)/bodyWindow;
									//	cout<<bodyWindow<<endl;exit(1);
										if(index<body_density.size())
											gb[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (*curRd - gIt->txEnd + PROMOTERUPLEN)/endWindow;
										if( index < end_density.size() )
											end[index] += 1./endWindow;
									}//*/
								}else{
									if( *curRd > gIt->txEnd - 2*PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd-*curRd)/promoterWindow;
										if( index < prom_density.size())
											pro[index] += 1./promoterWindow;
									}else if( *curRd > gIt->txStart + PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd - 2*PROMOTERUPLEN-*curRd)/bodyWindow;
										if( index < body_density.size() )
											gb[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (gIt->txStart+PROMOTERUPLEN-*curRd)/endWindow;
										if( index < end_density.size())
											end[index] += 1./endWindow;
									}//*/
								}
							}
							++curRd;
						}
						vector<double> den = pro;
						for( int i = 0; i < gb.size(); ++i)
							den.push_back(gb[i]);
						for( int i = 0; i < end.size(); ++i)
							den.push_back(end[i]);
						indx_bed_den.push_back(make_pair(gIt->island_score, make_pair(*gIt,den)));
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();
	int lib_sz = n;
//	int bins =4;
    cout<<"# total read: "<<lib_sz<<endl;

    sort( indx_bed_den.begin(), indx_bed_den.end() );

    int step = indx_bed_den.size()/bins > 1 ?
                    indx_bed_den.size()/bins : 1;
    cout<<"# genes per grp: "<< step<<endl;
    //      int max = step == 1 ? 1 :  indx_bed_den.size() - step/2;
    vector<vector<double> > aves;
    for( int i = 0; i <= indx_bed_den.size() - step/2; i += step ){
            vector<double> sum(indx_bed_den[i].second.second.size(), 0);
            for( int j = 0; j < step; ++j ){
                    vector<double> den = indx_bed_den[i+j].second.second;
                    for( int k = 0; k < den.size(); ++k ){
//                            sum[k] += log((den[k])*100000000./lib_sz)/log(2);
                            sum[k] += ((den[k])*1000000./lib_sz);
                    }
            }
	    vector<double> tmp;
		for( int k = 0; k < sum.size(); ++k )
                    tmp.push_back((sum[k]/step));
            aves.push_back(tmp);
    }

    std::ofstream out(out_f.data());
    out<<"Genomic_region\tLabel";
    for( int j = 0; j < aves.size(); ++j ){
    	out<<"\tgrp"<<j;
	}
	out<<endl;
    for( int i = 0; i < aves[0].size(); ++i ){
      if( i < promoterSliceN )
    	  out<<"Promoter\t"<<(-promoterSliceN/2+i)*promoterWindow+promoterWindow/2<<"_(bp)";
	  else if( i < promoterSliceN + genebodySliceN)
		  out<<"Gene_body\tFraction_"<<i-promoterSliceN+1;
	  else
		  out<<"After_TES\t"<<(i-promoterSliceN - genebodySliceN)*endWindow+endWindow/2<<"_(bp)";
      for( int j = 0; j < aves.size(); ++j )
		  out<<"\t"<<aves[j][i];
	  out<<endl;
    }
    out.close();
    cout<<"File saved in "<<out_f<<endl;
}

vector<double> generateReadDensityProfileOnRegions(map<Str, set<BED3> > regions, Str readf, Str outf, int shift, Str label, int window){
	//windows
	int sliceN = (regions.begin()->second.begin()->end - regions.begin()->second.begin()->start)/window;

	//densities
	vector<double> density(sliceN, 0);
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end, rn = 0;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty()){
			int pos = strand =="+" ? beg + shift : end - shift;
			//if( chrom == "chr1" )
			{
				chr_reads[chrom].push_back(pos);
				++rn;
			}
			if( rn % 5000000 == 0 || in.eof()){
				cout<<"# reads in: "<<rn<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, set<BED3> >::iterator chrIt = regions.begin(); chrIt != regions.end(); ++chrIt ){
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
						if( *curRd >= bIt->start){
							while( *curRd >= bIt->start && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						while( *curRd < bIt->end && curRd != chr_reads[chrIt->first].end()){
							if( *curRd > bIt->start ){
								unsigned int index = (*curRd - bIt->start)/window;
								if( index < density.size() )
									density[index] += 1./window;
								else{
								//	cout<<index<<'\t'<<sliceN<<endl;
								}
							}
							++curRd;
						}
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();

	std::ofstream out(outf.data());
	int NN = 0;
	out<<"\tn\tden\tgrp"<<endl;
	int average_window = 10;
	for( unsigned int j = average_window; j != density.size()-average_window; ++j){
		int pos = (j-sliceN/2)*window;
		double sum = 0;
		for( int k = 0; k < 2*average_window + 1; ++k){
			sum += density[j-average_window+k];
		}
		sum /= (2*average_window + 1);
		out<<(++NN)<<'\t'<<pos<<'\t'<<(sum)*1000000./rn/getSize(regions)<<"\t"<<label<<endl;
//		*it = (*it)*1000000./rn/getSize(regions);
	}
	out.close();//*/
	return density;
}

std::pair<int,int> readInIsland(map<Str, set<BED3> >& regions, Str readf, int shift ){
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	std::pair<int,int> in_all = make_pair(0,0);
	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty()){
			int pos = strand =="+" ? beg + shift : end - shift;
			chr_reads[chrom].push_back(pos);
			if( (++in_all.second) % 5000000 == 0 || in.eof()){
				cout<<"# reads in: "<<in_all.second<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, set<BED3> >::iterator chrIt = regions.begin(); chrIt != regions.end(); ++chrIt ){
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
						if( *curRd >= bIt->start){
							while( *curRd >= bIt->start && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						while( *curRd < bIt->end && curRd != chr_reads[chrIt->first].end()){
							if( *curRd > bIt->start ){
								++in_all.first;
							}
							++curRd;
						}
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();
	return in_all;
}
pair<vector<double>, vector<double> > generateReadDensityProfileOnRegions(
		map<Str, set<BED4> >& regions, Str readf, Str outf, Str label, int sliceN){
	//densities
	vector<double> density5(sliceN, 0);
	vector<double> density3(sliceN, 0);
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	map<Str, vector<BED4> > chr_reads;
	int rn = 0, n =0;
	map<Str, int> ptn_n;
	while( true ){
		BED4 bed;
		in>>bed.chrom>>bed.start>>bed.end>>bed.strand>>bed.strand>>bed.strand;

		chr_reads[bed.chrom].push_back(bed);
		++rn;
		//++ptn_n[bed.chrom];
		if( (++n) % 5000000 == 0 || in.eof()){
			cout<<"# reads in: "<<n<<endl;
			//map reads to region
			for( map<Str, set<BED4> >::iterator chrIt = regions.begin(); chrIt != regions.end(); ++chrIt ){
				if( chr_reads.find(chrIt->first) == chr_reads.end() )
					continue;
				vector<BED4>::const_iterator curRd = chr_reads[chrIt->first].begin();
				for( set<BED4>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
					if( curRd->end >= bIt->start){
						while( curRd->end >= bIt->start && curRd != chr_reads[chrIt->first].begin())
							--curRd;
					}
					int window = (bIt->end - bIt->start)/sliceN;
					while( curRd->start < bIt->end && curRd != chr_reads[chrIt->first].end()){
						if( curRd->end > bIt->start )
						{
							//++ptn_n[bIt->chrom];
							if( bIt->strand == "+" ){
								if( curRd->strand == "+" ){
									if( curRd->start > bIt->start ){
										unsigned int index = (curRd->start - bIt->start)/window;
										if( index < density5.size() ){
											density5[index] += 1./window;
										//	++ptn_n["++"+bIt->chrom];
										}
										else{
											--rn;
										}
									}
								}
								else{
									if( curRd->end < bIt->end ){
										unsigned int index = sliceN - (bIt->end - curRd->end)/window;
										if( index < density3.size() ){
											density3[index] += 1./window;
										//	++ptn_n["+-"+bIt->chrom];
										}
										else
											--rn;
									}
								}//*/
							}
							else{
								if( curRd->strand == "+" ){
									if( curRd->start > bIt->start ){
										unsigned int index = sliceN - (curRd->start - bIt->start)/window;
										if( index < density3.size() ){
											density3[index] += 1./window;
										//	++ptn_n["-+"];
										}
										else
											--rn;
									}
								}else{
									if( curRd->end < bIt->end ){
										unsigned int index = (bIt->end - curRd->end)/window;
										if( index < density5.size() ){
										//	++ptn_n["--"];
											density5[index] += 1./window;
										}
										else
											--rn;
									}
								}
							}//*/
						}
						++curRd;
					}
				}
			}
			chr_reads.clear();
			if( in.eof() && bed.chrom.empty())
				break;
		}
	}
	in.close();

	std::ofstream out((outf+"5").data());
	int NN = 0, i =0;
	out<<"\tn\tden\tgrp"<<endl;
	for( vector<double>::iterator it = density5.begin(); it != density5.end(); ++it, ++i ){
		out<<(++NN)<<'\t'<<(i-sliceN/2)<<'\t'<<(*it)*2000000./rn/getSize(regions)<<"\t"<<label<<'5'<<endl;
		*it = (*it)*2000000./rn/getSize(regions);
	}
	out.close();
	i = 0;
	NN = 0;
	std::ofstream out1((outf+"3").data());
	out1<<"\tn\tden\tgrp"<<endl;
	for( vector<double>::iterator it = density3.begin(); it != density3.end(); ++it, ++i ){
		out1<<(++NN)<<'\t'<<(i-sliceN/2)<<'\t'<<(*it)*2000000./rn/getSize(regions)<<"\t"<<label<<'3'<<endl;
		*it = (*it)*2000000./rn/getSize(regions);
	}
	out1.close();//*/
	for( map<Str, int>::iterator it = ptn_n.begin(); it != ptn_n.end(); ++it ){
		cout<<it->first<<'\t'<<it->second<<endl;
	}
	return make_pair(density5,density3);
}


pair<vector<double>, vector<double> > generateReadDensityProfile5ColumnsOnRegions(map<Str, set<BED4> >& regions, Str readf, Str outf, Str label, int sliceN){
	//densities
	vector<double> density5(sliceN, 0);
	vector<double> density3(sliceN, 0);
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	map<Str, set<BED4> > chr_reads;
	int rn = 0, n =0;
	while( true ){
		BED4 bed;
		in>>bed.chrom>>bed.start>>bed.end>>bed.strand>>bed.strand;

		chr_reads[bed.chrom].insert(bed);
		++rn;
		if( (++n) % 5000000 == 0 || in.eof()){
			cout<<"# reads in: "<<n<<endl;
			//map reads to region
			for( map<Str, set<BED4> >::iterator chrIt = regions.begin(); chrIt != regions.end(); ++chrIt ){
				if( chr_reads.find(chrIt->first) == chr_reads.end() )
					continue;
				set<BED4>::const_iterator curRd = chr_reads[chrIt->first].begin();
				for( set<BED4>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
					if( curRd->end >= bIt->start){
						while( curRd->end >= bIt->start && curRd != chr_reads[chrIt->first].begin())
							--curRd;
					}
					int window = (bIt->end - bIt->start)/sliceN;
					while( curRd->start < bIt->end && curRd != chr_reads[chrIt->first].end()){
						if( curRd->end > bIt->start )
						{
							if( bIt->strand == "+" ){
								if( curRd->strand == "+" ){
									if( curRd->start > bIt->start ){
										unsigned int index = (curRd->start - bIt->start)/window;
										if( index < density5.size() ){
											density5[index] += 1./window;
										}
										else{
											--rn;
										}
									}
								}
								else{
									if( curRd->end < bIt->end ){
										unsigned int index = sliceN - (bIt->end - curRd->end)/window;
										if( index < density3.size() ){
											density3[index] += 1./window;
										}
										else
											--rn;
									}
								}//*/
							}
							else{
								if( curRd->strand == "+" ){
									if( curRd->start > bIt->start ){
										unsigned int index = sliceN - (curRd->start - bIt->start)/window;
										if( index < density3.size() ){
											density3[index] += 1./window;
										}
										else
											--rn;
									}
								}else{
									if( curRd->end < bIt->end ){
										unsigned int index = (bIt->end - curRd->end)/window;
										if( index < density5.size() ){
											density5[index] += 1./window;
										}
										else
											--rn;
									}
								}
							}//*/
						}
						++curRd;
					}
				}
			}
			chr_reads.clear();
			if( (in.eof() && bed.chrom.empty() ))
				break;
		}
	}
	in.close();

	std::ofstream out((outf+"5").data());
	int NN = 0, i =0;
	out<<"\tn\tden\tgrp"<<endl;
	for( vector<double>::iterator it = density5.begin(); it != density5.end(); ++it, ++i ){
		out<<(++NN)<<'\t'<<(i-sliceN/2)*10<<'\t'<<(*it)/rn/getSize(regions)*1000000.<<"\t"<<label<<'5'<<endl;
		*it = (*it)*1000000./rn/getSize(regions);
	}
	out.close();
	i = 0;
	NN = 0;
	std::ofstream out1((outf+"3").data());
	out1<<"\tn\tden\tgrp"<<endl;
	for( vector<double>::iterator it = density3.begin(); it != density3.end(); ++it, ++i ){
		out1<<(++NN)<<'\t'<<(i-sliceN/2)*10<<'\t'<<(*it)/rn/getSize(regions)*1000000.<<"\t"<<label<<'3'<<endl;
		*it = (*it)*1000000./rn/getSize(regions);
	}
	out1.close();//*/
	return make_pair(density5,density3);
}

int assignReadCount2Island(map<Str, set<BED3> >& island,
		Str readf, Str label, int shift){
	cout<<readf<<endl;
	map<Str, vector<BED3> > islands = set2vector(island);
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	for( map<Str, vector<BED3> >::iterator chrIt = islands.begin(); chrIt != islands.end(); ++chrIt ){
		for( vector<BED3>::iterator iIt = chrIt->second.begin(); iIt != chrIt->second.end(); ++iIt ){
			if( iIt->readDensityOn.find(label) == iIt->readDensityOn.end() ){
				iIt->readDensityOn[label] = 0.000001;
				iIt->readCountOn[label] = 0;//0.000001;
			//	iIt->readCountsOn[label] = vector<double>( iIt->end - iIt->start + 1, 0);
			}
		}
	}//*/

	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end, n = 0, rn = 0;
	map<pair<Str,int>, vector<BED3> > rd_bed3;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty()){
			++rn;
			int pos = strand =="+" ? beg + shift : end - shift;
			chr_reads[chrom].push_back(pos);
			if( (++n) % 5000000 == 0 || in.eof()){
				cout<<"# reads in: "<<n<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				int mapN = 0;
				for( map<Str, vector<BED3> >::iterator chrIt = islands.begin(); chrIt != islands.end(); ++chrIt ){
					if( chrIt->first == "chrY" || chr_reads.find(chrIt->first) == chr_reads.end() )
						continue;
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<BED3>::iterator iIt = chrIt->second.begin(); iIt != chrIt->second.end(); ++iIt ){
						//cout<<iIt->chrom<<"\t"<<iIt->start<<"\t"<<iIt->end<<endl;
						if( *curRd >= iIt->start){
							while( *curRd >= iIt->start && curRd != chr_reads[chrIt->first].begin()){
						//		cout<<"!"<<*curRd<<endl;
								--curRd;
							}
						}
						//cout<<iIt->chrom<<"\t"<<iIt->start<<"\t"<<iIt->end<<endl;
						for( ; *curRd < iIt->end && curRd != chr_reads[chrIt->first].end(); ++curRd){
							if( *curRd >= iIt->start){
								iIt->readDensityOn[label] += 1./ (iIt->end-iIt->start);
								iIt->readCountOn[label] += 1;
								/*rd_bed3[make_pair(iIt->chrom, *curRd)].push_back(*iIt);
								if( rd_bed3[make_pair(iIt->chrom, *curRd)].size() > 1){
									for( vector<BED3>::const_iterator iit = rd_bed3[make_pair(iIt->chrom, *curRd)].begin();
											iit != rd_bed3[make_pair(iIt->chrom, *curRd)].end(); ++iit ){
										cout<<iit->chrom<<":"<<iit->start<<'\t'<<iit->end<<'\t'<<*curRd<<endl;
									}
									cout<<"=====\n";
									exit(1);
								}//*/
								int index = (*curRd - iIt->start);
								++mapN;
								//++iIt->readCountsOn[label][index];
							}
						}
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();
	int window = 50;
	//cout<<rn<<endl;
	for( map<Str, vector<BED3> >::iterator chrIt = islands.begin(); chrIt != islands.end(); ++chrIt ){
		for( vector<BED3>::iterator iIt = chrIt->second.begin(); iIt != chrIt->second.end(); ++iIt ){
			if( iIt->readDensityOn.find(label) != iIt->readDensityOn.end() ){
				iIt->readDensityOn[label] = iIt->readDensityOn[label] * 1000000. / rn;//RPBM
			}
			//if( iIt->readCountOn.find(label) != iIt->readCountOn.end() ){
			//	iIt->readCountOn[label] = iIt->readCountOn[label];// * 1000000. / rn;//RPBM
			//}
			/*map<Str, vector<double>  >::iterator it = iIt->readCountsOn.find(label);
			if( it != iIt->readCountsOn.end() ){
				double max = 0;
				int index = -1;
				for( unsigned int i = 0; i < it->second.size() - window; ++i ){
					int sum = 0;
					for( int j = 0; j < window; ++ j)
						sum += it->second[i+j];
					if( sum > max ){
						max = sum;
						index = i;
					}
				}
				iIt->sumitPosCountOn[label].first = iIt->start + index + window/2;
				iIt->sumitPosCountOn[label].second = max;
			}//*/
		}
	}
	island = vector2set(islands);
	//cout<<rn<<endl;
	return rn;
}

int smoothed_counts_in_islands(map<Str, vector<BED3> >& islands,
		Str readf, Str label, int shift, int upR, int dnR ){
	int N = getSize( islands );
	for( map<Str, vector<BED3> >::iterator chrIt = islands.begin(); chrIt != islands.end(); ++chrIt ){
		for( vector<BED3>::iterator iIt = chrIt->second.begin(); iIt != chrIt->second.end(); ++iIt ){
			if( iIt->readCountsOn.find(label+"5") == iIt->readCountsOn.end() ){
				iIt->readCountsOn[label+"5"] = vector<double>(upR + dnR, 0);
			}
			if( iIt->readCountsOn.find(label+"3") == iIt->readCountsOn.end() ){
				iIt->readCountsOn[label+"3"] = vector<double>(upR + dnR, 0);
			}
		}
	}//*/

	std::ifstream in(readf.data());
	if( !in.good() ){
		return -1;
		FileNotFoundERR(readf);
	}

	Str chrom, strand;
	int beg, end, n = 0, rn = 0;
	map<Str, vector<pair<int,Str> > > chr_reads;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty()){
			++rn;
			int pos = strand =="+" ? beg + shift: end - shift;
			chr_reads[chrom].push_back(make_pair(pos,strand));
			if( (++n) % 5000000 == 0 || in.eof() ){
				cout<<"# reads in: "<<n<<endl;
				//sort the reads to speed up
				//cout<<chrom<<":"<<beg<<"-"<<end<<"\t"<<strand<<endl;
				for( map<Str, vector<pair<int,Str> > >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());

				for( map<Str, vector<BED3> >::iterator chrIt = islands.begin(); chrIt != islands.end(); ++chrIt ){
					if( chr_reads.find(chrIt->first) == chr_reads.end() )
						continue;
					vector<pair<int,Str> >::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<BED3>::iterator iIt = chrIt->second.begin(); iIt != chrIt->second.end(); ++iIt ){
						int sumit = iIt->sumit;
						int lft = sumit - upR;
						int rgt = sumit + dnR;
						if(iIt->strand == "+"){
							if( curRd->first >= lft ){
								while( curRd->first >= lft && curRd != chr_reads[chrIt->first].begin())
									--curRd;
							}
							for( ; curRd->first < rgt && curRd != chr_reads[chrIt->first].end(); ++curRd){
								if( curRd->first >= lft ){
									if( curRd->second == "+" )
										++iIt->readCountsOn[label+"5"][curRd->first-lft];
									else
										++iIt->readCountsOn[label+"3"][curRd->first-lft];
								}
							}
						}else{
							if( curRd->first >= lft ){
								while( curRd->first >= lft && curRd != chr_reads[chrIt->first].begin()){
									--curRd;
								}
							}
							for( ; curRd->first <= rgt && curRd != chr_reads[chrIt->first].end(); ++curRd){
								if( curRd->first > lft ){
									if( curRd->second == "+" )
										++iIt->readCountsOn[label+"3"][rgt-curRd->first];
									else
										++iIt->readCountsOn[label+"5"][rgt-curRd->first];
								}
							}
						}
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	/*cout<<"Normlization ...\n";
	//normlize on total tag
	for( map<Str, vector<BED3> >::iterator chrIt = islands.begin(); chrIt != islands.end(); ++chrIt ){
		for( vector<BED3>::iterator iIt = chrIt->second.begin(); iIt != chrIt->second.end(); ++iIt ){
			map<Str, vector<double>  >::iterator l5 = iIt->readCountsOn.find(label+"5");
			if( l5 != iIt->readCountsOn.end() ){
				for( int i = 0; i < l5->second.size(); ++i ){
					l5->second[i] = l5->second[i];
				}
			}
			map<Str, vector<double>  >::iterator l3 = iIt->readCountsOn.find(label+"3");
			if( l3 != iIt->readCountsOn.end() ){
				for( int i = 0; i < l3->second.size(); ++i ){
					l3->second[i] = l3->second[i];
				}
			}
		}
	}
	cout<<"Done.\n";//*/
	return n;
}

void assignReadCount5Columns2Island(map<Str, set<BED3> >& island,
		Str readf, Str label, int shift){
	map<Str, vector<BED3> > islands = set2vector(island);
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	for( map<Str, vector<BED3> >::iterator chrIt = islands.begin(); chrIt != islands.end(); ++chrIt ){
		for( vector<BED3>::iterator iIt = chrIt->second.begin(); iIt != chrIt->second.end(); ++iIt ){
			iIt->readCountsOn[label] = vector<double>( iIt->end- iIt->start + 1,0);
		}
	}//*/
	map<Str, set<int> > chr_reads;
	Str chrom, strand;
	int beg, end, n = 0, rn = 0;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand;
		if( !chrom.empty()){
			++rn;
			int pos = strand =="+" ? beg + shift : end - shift;
			chr_reads[chrom].insert(pos);
			if( (++n) % 5000000 == 0 || in.eof()){
				cout<<"# reads in: "<<n<<endl;
				//sort the reads to speed up
			//	for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
				//	std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, vector<BED3> >::iterator chrIt = islands.begin(); chrIt != islands.end(); ++chrIt ){
					set<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<BED3>::iterator iIt = chrIt->second.begin(); iIt != chrIt->second.end(); ++iIt ){
						if( *curRd >= iIt->start){
							while( *curRd >= iIt->start && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						for( ; *curRd < iIt->end && curRd != chr_reads[chrIt->first].end(); ++curRd){
							if( *curRd >= iIt->start){
								iIt->readDensityOn[label] += 1./ (iIt->end-iIt->start);
								iIt->readCountOn[label] += 1;
								int index = (*curRd - iIt->start);
								++iIt->readCountsOn[label][index];
							}
						}
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();
	int window = 10;
	for( map<Str, vector<BED3> >::iterator chrIt = islands.begin(); chrIt != islands.end(); ++chrIt ){
		for( vector<BED3>::iterator iIt = chrIt->second.begin(); iIt != chrIt->second.end(); ++iIt ){
			if( iIt->readDensityOn.find(label) != iIt->readDensityOn.end() ){
				iIt->readDensityOn[label] = iIt->readDensityOn[label] * 1000000. / rn;//RPBM
			}
			map<Str, vector<double>  >::iterator it = iIt->readCountsOn.find(label);
			if( it != iIt->readCountsOn.end() ){
				int max = 0;
				int index = -1;
				for( unsigned int i = 0; i < it->second.size() - window; ++i ){
					int sum = 0;
					for( int j = 0; j < window; ++ j)
						sum += it->second[i+j];
					if( sum > max ){
						max = sum;
						index = i;
					}
				}
				iIt->sumitPosCountOn[label].first = iIt->start + index + window/2;
				iIt->sumitPosCountOn[label].second = max;
			}
		}
	}//*/
	island = vector2set(islands);
}

//Filter reads outsides islands
void filterReadsOutOfIsland(const map<Str, set<BED3> >& islands,
		Str readf, Str outf, int shift){
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	map<Str, map<int, BED6> > chr_reads;
	Str chrom, a, b, strand;
	int n = 0;
	std::ofstream out(outf.data());
	while( true ){
		BED6 b;
		in>>b.chrom>>b.start>>b.end>>b.name>>b.score>>b.strand;

		int pos = b.strand =="+" ? b.start + shift : b.end - shift;
		chr_reads[b.chrom][pos] = b;
		if( (++n) % 5000000 == 0 || in.eof() ){
			cout<<"# reads in: "<<n<<endl;
			//map reads to region
			for( map<Str, set<BED3> >::const_iterator chrIt = islands.begin(); chrIt != islands.end(); ++chrIt ){
				map<int, BED6>::iterator curRd = chr_reads[chrIt->first].begin();
				for( set<BED3>::const_iterator iIt = chrIt->second.begin(); iIt != chrIt->second.end(); ++iIt ){
					if( curRd->first >= iIt->start){
						while( curRd->first >= iIt->start && curRd != chr_reads[chrIt->first].begin()){
							--curRd;
						}
					}
					for( ; curRd->first < iIt->end && curRd != chr_reads[chrIt->first].end(); ++curRd){
						if( curRd->first >= iIt->start){
							BED6& b = curRd->second;
							out<<b.chrom<<'\t'<<b.start<<'\t'<<b.end<<'\t'<<b.name<<'\t'<<b.score<<'\t'<<b.strand<<endl;
						}
					}
				}
			}
			chr_reads.clear();
			if(b.chrom.empty())
				break;
		}
	}
	in.close();
	out.close();
}

void assignUCSC2Islands( map<Str, set<UCSC, UCSC::sortByTxStart> >& chr_genes,
		map<Str, set<BED3> >& chr_island){
	map<Str, vector<BED3> > chr_islands = set2vector(chr_island);
	for( map<Str, set<UCSC, UCSC::sortByTxStart> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
		vector<BED3>::iterator irIt = chr_islands[chrIt->first].begin();
		for( set<UCSC, UCSC::sortByTxStart>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
			while( irIt->end > gIt->txStart - PROMOTERUPLEN && irIt != chr_islands[chrIt->first].begin() )
				--irIt;
			for( ; irIt->start < gIt->txEnd + PROMOTERUPLEN && irIt != chr_islands[chrIt->first].end(); ++irIt ){
				if( irIt->end > gIt->txStart - PROMOTERUPLEN){
					Pa_I_I region = getIntrestRegion(*gIt, PROMOTER);
					if( _min(region.second, irIt->end) > _max(region.first, irIt->start) ){
						irIt->UCSCOn[PROMOTER].insert(*gIt);
					}else{
						region = getIntrestRegion(*gIt, GENEBODY);
						if( region.second > region.first){
							if( _min(region.second, irIt->end) > _max(region.first, irIt->start) ){
								irIt->UCSCOn[GENEBODY].insert(*gIt);
							}
						}
					}
				}
			}
		}
	}
	chr_island = vector2set(chr_islands);
}

void generateRPBMSummaryGraph(Str bed_f, Str chr_len_f, Str out_f, int windowsize, int shift, bool invert ){
	//Initiate graph file windows;
	cout<<"Window size: "<<windowsize<<endl;
	cout<<"Fragment shift: "<<shift<<endl;
	map<Str, int> chr_len = getChrlen(chr_len_f);
	map<Str, vector<double> > chr_RPBMs;
	for( map<Str, int>::iterator it = chr_len.begin(); it != chr_len.end(); ++it ){
		chr_RPBMs[it->first] = vector<double>( it->second/windowsize, 0);
	}
	//
	std::ifstream in(bed_f.data());
	if( !in.good() )
		FileNotFoundERR(bed_f);
	int rn = 0;
	while( !in.eof() ){
		BED6 b;
		in>>b.chrom>>b.start>>b.end>>b.name>>b.score>>b.strand;
		if( !b.chrom.empty() ){
			int pos = b.strand =="+" ? b.start + shift : b.end - shift;
			map<Str, vector<double> >::iterator chrIt = chr_RPBMs.find(b.chrom);
			if( chrIt != chr_RPBMs.end() ){
				++rn;
				unsigned int index = pos / windowsize;
				if( index < chrIt->second.size()){
					chrIt->second[index] += 1;
				}
			}
		}
	}
	in.close();

	cout<<"# reads: "<<rn<<endl;
	std::ofstream out(out_f.data());
	for( map<Str, vector<double> >::iterator it = chr_RPBMs.begin(); it != chr_RPBMs.end(); ++it ){
		for( unsigned int i = 0; i < it->second.size(); ++i ){
			if( it->second[i] > 0){
				double rpbm = it->second[i]*1000000./windowsize/rn;
				if(invert)
					rpbm *= -1;
				//out<<it->first<<'\t'<<i*windowsize<<'\t'<<((i+1)*windowsize-1)<<'\t'<<rpbm<<endl;
                                out<<it->first<<'\t'<<i*windowsize<<'\t'<<((i+1)*windowsize-1)<<'\t'<<it->second[i]<<endl;
			}
		}
	}
	out.close();
}


void generateRPBMSummaryGraphStrand(Str bed_f, Str chr_len_f, Str out_f, int windowsize, int shift, Str _strand ){
	//Initiate graph file windows;
	cout<<"Window size: "<<windowsize<<endl;
	cout<<"Fragment shift: "<<shift<<endl;
	map<Str, int> chr_len = getChrlen(chr_len_f);
	map<Str, vector<double> > chr_RPBMs;
	for( map<Str, int>::iterator it = chr_len.begin(); it != chr_len.end(); ++it ){
		chr_RPBMs[it->first] = vector<double>( it->second/windowsize, 0);
	}
	//
	std::ifstream in(bed_f.data());
	if( !in.good() )
		FileNotFoundERR(bed_f);
	int rn = 0;
	while( !in.eof() ){
		BED6 b;
		in>>b.chrom>>b.start>>b.end>>b.name>>b.score>>b.strand;
		if( !b.chrom.empty() ){
			int pos = b.strand =="+" ? b.start + shift : b.end - shift;
			map<Str, vector<double> >::iterator chrIt = chr_RPBMs.find(b.chrom);
			if( chrIt != chr_RPBMs.end() && b.strand == _strand ){
				++rn;
				unsigned int index = pos / windowsize;
				if( index < chrIt->second.size()){
					chrIt->second[index] += 1;
				}
			}
		}
	}
	in.close();

	cout<<"# reads: "<<rn<<endl;
	std::ofstream out(out_f.data());
	for( map<Str, vector<double> >::iterator it = chr_RPBMs.begin(); it != chr_RPBMs.end(); ++it ){
		for( unsigned int i = 0; i < it->second.size(); ++i ){
			if( it->second[i] > 0){
				double rpbm = it->second[i]*1000000./windowsize/rn;
				out<<it->first<<'\t'<<i*windowsize<<'\t'<<((i+1)*windowsize-1)<<'\t'<<rpbm<<endl;
			}
		}
	}
	out.close();
}

void generateRPBMSummaryWig(Str bed_f, Str chr_len_f, Str out_f, int windowsize, int shift, Str chr ){
	//Initiate graph file windows;
	cout<<"Window size: "<<windowsize<<endl;
	cout<<"Fragment shift: "<<shift<<endl;
	map<Str, int> chr_len = getChrlen(chr_len_f);
	map<Str, vector<double> > chr_RPBMs;
	for( map<Str, int>::iterator it = chr_len.begin(); it != chr_len.end(); ++it ){
		chr_RPBMs[it->first] = vector<double>( it->second/windowsize, 0);
	}
	//
	std::ifstream in(bed_f.data());
	if( !in.good() )
		FileNotFoundERR(bed_f);
	int rn = 0;
	while( !in.eof() ){
		BED6 b;
		in>>b.chrom>>b.start>>b.end>>b.name>>b.score>>b.strand;
		if( !b.chrom.empty() ){
			int pos = b.strand =="+" ? b.start + shift : b.end - shift;
			map<Str, vector<double> >::iterator chrIt = chr_RPBMs.find(b.chrom);
			if( chrIt != chr_RPBMs.end() ){
				++rn;
				unsigned int index = pos / windowsize;
				if( index < chrIt->second.size()-1){
					chrIt->second[index] += 1;
				}
			}
		}
	}
	in.close();

	cout<<"# reads: "<<rn<<endl;
	std::ofstream out(out_f.data());
	for( map<Str, vector<double> >::iterator it = chr_RPBMs.begin(); it != chr_RPBMs.end(); ++it ){
//		out<<"fixedStep chrom="<<it->first<<" start=1 step="<<windowsize<<" span="<<windowsize<<endl;
		out<<"variableStep chrom="<<it->first<<" span="<<windowsize<<endl;
		if( it->first == chr || chr == "all" ){
			for( unsigned int i = 0; i < it->second.size(); ++i ){
				//if( it->second[i] > 0)
				{
					double rpbm = it->second[i]*1000000./windowsize/rn;
                        	        //out<<int(rpbm)<<endl;
					if( it->second[i] > 0)
						out<<(i*windowsize+windowsize/2)<<"\t"<<rpbm<<endl;
				}
			}
		}
	}
	out.close();
}

map<Str, set<BED3> > get_intrest_regions(map<Str, set<UCSC, UCSC::sortByTxStart> > chr_genes,
		Str chr_len_f, Str tag){
	map<Str, set<BED4> > _promoter = getIntrestReions(chr_genes, PROMOTER);
	map<Str, set<BED3> > promoter = bed42bed3(_promoter);
	map<Str, set<BED4> > _genic = getIntrestReions(chr_genes, GENEBODY);
	map<Str, set<BED3> > genic = bed42bed3(_genic);

	map<Str, set<BED3> > intergenic = get_intergenic(promoter, genic
				, chr_len_f );
	map<Str, set<BED3> > pro_ = self_union(promoter);//fix
	map<Str, set<BED3> > inter_ = self_union(intergenic);//fix

	map<Str, set<BED3> > pro_inter = merge(pro_, inter_);
	map<Str, set<BED3> > genic_ = get_complement_set(pro_inter, chr_len_f);//fix
	if( tag == "pro")
		return pro_;
	else if( tag == "genic")
		return self_union(genic_);
	else if( tag == "inter")
		return inter_;
	else if( tag == "exon"){
		map<Str, set<BED3> > exons;
		for( map<Str, set<UCSC, UCSC::sortByTxStart> >::iterator chrIt = chr_genes.begin();
				chrIt != chr_genes.end(); ++chrIt ){
			for( set<UCSC, UCSC::sortByTxStart>::iterator gIt = chrIt->second.begin();
					gIt != chrIt->second.end(); ++gIt ){
				vector<Pa_I_I>::const_iterator eIt = gIt->exons53.begin();
				for( ; eIt != gIt->exons53.end(); ++eIt ){
					BED3 b;
					b.chrom = chrIt->first;
					b.start = eIt->first;
					b.end = eIt->second;
					exons[b.chrom].insert(b);
				}
			}
		}
		return self_union(exons);
	}else if( tag == "genebody"){
		map<Str, set<BED3> > gb;
		for( map<Str, set<UCSC, UCSC::sortByTxStart> >::iterator chrIt = chr_genes.begin();
				chrIt != chr_genes.end(); ++chrIt ){
			for( set<UCSC, UCSC::sortByTxStart>::iterator gIt = chrIt->second.begin();
					gIt != chrIt->second.end(); ++gIt ){
				BED3 b;
				b.chrom = chrIt->first;
				b.start = gIt->txStart;
				b.end = gIt->txEnd;
				gb[b.chrom].insert(b);
			}
		}
		return self_union(gb);
	}
	else{
		return map<Str, set<BED3> >();
	}
}

map<Str, set<BED3> > notOverlapWithPromoter(map<Str, set<BED3> > pro, map<Str, set<BED3> > beds){
	assignIsland2Island(beds, pro, "pro");

	map<Str, set<BED3> > rst;
	map<Str, vector<BED3> > tmp = set2vector(beds);
	for( map<Str, vector<BED3> >::iterator chrIt = tmp.begin(); chrIt != tmp.end(); ++chrIt ){
		for( vector<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			if( bIt->islandsOn.find("pro") == bIt->islandsOn.end() ){
				rst[bIt->chrom].insert(*bIt);
			}
		}
	}
//	cout<<getSize(rst)<<" out of "<<getSize(beds)<<endl;
	return rst;
}

map<Str, set<BED3> > overlapWithPromoter(map<Str, set<BED3> > pro, map<Str, set<BED3> > beds){
	assignIsland2Island(beds, pro, "pro");

	map<Str, set<BED3> > rst;
	map<Str, vector<BED3> > tmp = set2vector(beds);
	for( map<Str, vector<BED3> >::iterator chrIt = tmp.begin(); chrIt != tmp.end(); ++chrIt ){
		for( vector<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			if( bIt->islandsOn.find("pro") != bIt->islandsOn.end() ){
				bIt->islandsOn.erase(bIt->islandsOn.find("pro"));
				rst[bIt->chrom].insert(*bIt);
			}
		}
	}
//	cout<<getSize(rst)<<" out of "<<getSize(beds)<<endl;
	return rst;
}

void attachNuclsScore2Islands(
		map<Str, set<BED3> >& regions, Str nuc_readf, Str label, int window){
	//initiate nucls scores
	map<Str, vector<BED3> > _regions = set2vector(regions);
	for( map<Str, vector<BED3> >::iterator chrIt = _regions.begin(); chrIt != _regions.end(); ++chrIt ){
		for( vector<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			bIt->nucl_scr[label] = vector<double>((bIt->end - bIt->start)/window, 0);
		}
	}
	//densities
	std::ifstream in(nuc_readf.data());
	if( !in.good() )
		FileNotFoundERR(nuc_readf);
	int n = 0;
	map<Str, map<int,double> > chr_pos_scr;
	while( true ){
		Str chrom;
		int pos1, pos2;
		double scr;
		in>>chrom>>pos1>>pos2>>scr;
		if( !chrom.empty() )
			chr_pos_scr[chrom][(pos1+pos2)/2] = scr;
		if( (++n) % 5000000 == 0 || in.eof()){
			cout<<"# reads in: "<<n<<endl;
			//map reads to region
			for( map<Str, vector<BED3> >::iterator chrIt = _regions.begin(); chrIt != _regions.end(); ++chrIt ){
				if( chr_pos_scr.find(chrIt->first) == chr_pos_scr.end() )
					continue;
				map<int,double>::const_iterator curRd = chr_pos_scr[chrIt->first].begin();
				for( vector<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
					if( curRd->first >= bIt->start){
						while( curRd->first >= bIt->start && curRd != chr_pos_scr[chrIt->first].begin())
							--curRd;
					}
					while( curRd->first < bIt->end && curRd != chr_pos_scr[chrIt->first].end()){
						if( curRd->first > bIt->start )
						{
							unsigned int index = (curRd->first - bIt->start)/window;
							if( index < bIt->nucl_scr[label].size() ){
								bIt->nucl_scr[label][index] += curRd->second;
							}else{
								cout<<"ERR\n";
							}
						}
						++curRd;
					}
				}
			}
			chr_pos_scr.clear();
			if( in.eof() && chrom.empty())
				break;
		}
	}
	in.close();
	regions = vector2set(_regions);
}

void attachNuclsScore2Islands(
		map<Str, set<BED4> >& regions, Str nuc_readf, Str label, int window){
	//initiate nucls scores
	map<Str, vector<BED4> > _regions = set2vector(regions);
	for( map<Str, vector<BED4> >::iterator chrIt = _regions.begin(); chrIt != _regions.end(); ++chrIt ){
		for( vector<BED4>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			bIt->nucl_scr[label] = vector<double>((bIt->end - bIt->start)/window, 0);
		}
	}
	//densities
	std::ifstream in(nuc_readf.data());
	if( !in.good() )
		FileNotFoundERR(nuc_readf);
	int n = 0;
	map<Str, map<int,double> > chr_pos_scr;
	while( true ){
		Str chrom;
		int pos1, pos2;
		double scr;
		in>>chrom>>pos1>>pos2>>scr;
		if( !chrom.empty() )
			chr_pos_scr[chrom][(pos1+pos2)/2] = scr;
		if( (++n) % 5000000 == 0 || in.eof()){
			cout<<"# reads in: "<<n<<endl;
			//map reads to region
			for( map<Str, vector<BED4> >::iterator chrIt = _regions.begin(); chrIt != _regions.end(); ++chrIt ){
				if( chr_pos_scr.find(chrIt->first) == chr_pos_scr.end() )
					continue;
				map<int,double>::const_iterator curRd = chr_pos_scr[chrIt->first].begin();
				for( vector<BED4>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
					if( curRd->first >= bIt->start){
						while( curRd->first >= bIt->start && curRd != chr_pos_scr[chrIt->first].begin())
							--curRd;
					}
					while( curRd->first < bIt->end && curRd != chr_pos_scr[chrIt->first].end()){
						if( curRd->first > bIt->start )
						{
							unsigned int index = bIt->strand == "+" ?
									(curRd->first - bIt->start)/window :
									(bIt->end - curRd->first)/window;
							if( index < bIt->nucl_scr[label].size() ){
								bIt->nucl_scr[label][index] += curRd->second;
							}else{
								cout<<"ERR\n";
							}
						}
						++curRd;
					}
				}
			}
			chr_pos_scr.clear();
			if( in.eof() && chrom.empty())
				break;
		}
	}
	in.close();
	regions = vector2set(_regions);
}


map<Str, SignalPos > scoreWithPWM(M_D pwm, vector<pair<Str, Str> > fast, bool bothStrand ){
	map<Str, SignalPos > id_SignalPos;
	pwm.toBeLoged();
	for( vector<pair<Str, Str> >::iterator it = fast.begin(); it != fast.end(); ++it ){
		Str s = it->second;
		int index = 0;
		double maxS = -1000000;
		for( unsigned int i = 0; i < s.size() - pwm.colum +1; ++i ){
			double scr = 0;
			for( unsigned int k = 0; k < pwm.colum; ++k ){
				scr += pwm(s[i+k]-'0', k);
			}
			if( scr > maxS ){
				maxS = scr;
				index = i;
			}
		}
		Str substr = it->second.substr(index, pwm.colum);
		bool inOpposite = false;
		if( bothStrand ){
			s = getDNAInOppositeStrand(s);
			for( unsigned int i = 0; i < s.size() - pwm.colum + 1; ++i ){
				double scr = 0;
				for( unsigned int k = 0; k < pwm.colum; ++k ){
					scr += pwm(s[i+k]-'0', k);
				}
				if( scr > maxS ){
					maxS = scr;
					index = i;
					inOpposite = true;
				}
			}
			if( inOpposite ){
				substr = s.substr(index, pwm.colum);
				index = it->second.size() - index + 1;
			}
		}
		SignalPos sp;
		sp.relative_start_pos = index;
		sp.seq = SequenceTransform_T::digital2CharSeq(substr);
		sp.scr = maxS;
		sp.inOppositeStrand = inOpposite;
		id_SignalPos[it->first] = sp;
	}
	return id_SignalPos;
}

map<Str, set<BED4> > parsSp2Bed4(map<Str, SignalPos > sp, int upL){
	map<Str, set<BED4> > b4s;
	for( map<Str, SignalPos >::iterator it = sp.begin(); it != sp.end(); ++it ){
		BED3 b3 = parsStr2Bed3(it->first);
		BED4 b4;
		b4.chrom = b3.chrom;
		b4.strand = it->second.inOppositeStrand ? "-" : "+";
		int pos = (b4.strand=="+" ? b3.start + it->second.relative_start_pos + it->second.seq.size()/2
				: b3.start + it->second.relative_start_pos - it->second.seq.size()/2 );
		b4.start = pos - upL;
		b4.end = pos + upL;
		b4.info = it->second.seq;
		b4.sp = it->second;
		b4s[b4.chrom].insert(b4);
//		cout<<it->first<<'\t'<<it->second.relative_start_pos
	//			<<'\t'<<b4.info<<'\t'<<b4.strand<<'\t'<<b4.chrom<<":"<<b4.start<<"-"<<b4.end<<endl;
	}
	return b4s;
}

void attachPairEndTag2Islands( map<Str, set<BED3> >& regions, Str paie_end_f, Str label){
	std::ifstream in(paie_end_f.data());
	if( !in.good() )
		FileNotFoundERR(paie_end_f);
	map<Str, vector<BED3> > _regions = set2vector(regions);
	map<Str, vector<pair<int,int> > > chr_reads;
	Str chrom;
	int beg, end, rn = 0;
	double len = 0;
	while( true ){
		in>>chrom>>beg>>end;
		if( !chrom.empty()){
			chr_reads[chrom].push_back(make_pair((beg+end)/2, end - beg));
			len += end - beg;
			++rn;

			if( rn % 5000000 == 0 || in.eof()){
				cout<<"# reads in: "<<rn<<endl;
				//sort the reads to speed up
				for( map<Str, vector<pair<int,int> > >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, vector<BED3> >::iterator chrIt = _regions.begin(); chrIt != _regions.end(); ++chrIt ){
					vector<pair<int,int> >::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
						if( curRd->first >= bIt->start){
							while( curRd->first >= bIt->start && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						while( curRd->first < bIt->end && curRd != chr_reads[chrIt->first].end()){
							if( curRd->first > bIt->start ){
								bIt->pair_end_mid_len[label].push_back(*curRd);
							}
							++curRd;
						}
					}
				}
				chr_reads.clear();
				if( in.eof() //|| chrom == "chr2"
						)
					break;
			}
		}
	}
	in.close();
	regions = vector2set(_regions);
	cout<<"Average length "<<len/rn<<endl;

}


map<Str, set<BED3> > BEDFrmFASTA( Str fileName ){
	vector<pair<Str, Str> > fasta = readInFast( fileName );

	map<Str, set<BED3> > chr_bed;

	for( vector<pair<Str, Str> >::iterator it = fasta.begin(); it != fasta.end(); ++it ){
		BED3 b = parsStr2Bed3(it->first.substr(1));
		b.seq = it->second;
		chr_bed[b.chrom].insert(b);
	}
	cout<<"# of records from \""<<fileName<<"\": "<<getSize(chr_bed)<<endl;

	return chr_bed;
}

void write4quatiles(vector<pair<double, BED3> > pi_ucsc, Str tag ){
    sort( pi_ucsc.begin(), pi_ucsc.end() );
    unsigned int step = pi_ucsc.size()/4;
    map<int, vector<BED3> > pi_ucscs;
    for( unsigned int i = 0; i < pi_ucsc.size(); i += step ){
            for( unsigned int j = 0; j < step; ++j ){
                    if( i+j < pi_ucsc.size() )
                            pi_ucscs[(i+j)/step].push_back(pi_ucsc[i+j].second);
            }
    }

    for( map<int, vector<BED3> >::iterator it = pi_ucscs.begin(); it != pi_ucscs.end(); ++it ){
            std::ofstream out( Str( tag+convert2String(it->first)).data() );
            for( vector<BED3>::iterator i = it->second.begin(); i != it->second.end(); ++i ){
                     out<<i->chrom<<'\t'<<i->sumit -3000<<'\t'<<i->sumit + 3000<<'\t'<<i->strand<<endl;
             }
             out.close();
     }//*/
}

map<Str,set<BED3> > processSICER11( Str f, double fc, double p ){
	map<Str,set<BED3> > status_beds;
	std::ifstream in(f.data());
	Str line;
	std::getline( in, line);
//	cout<<line<<endl;
//	std::ofstream out(Str(f+".r").data());
//	out<<"\ta\tb\tgrp\n";
	int NN = 0;
	while( !in.eof()){
		BED3 b;
		int Readcount_A, ReadcountB;
		double Normalized_Readcount_A,Normalized_Readcount_B,
		   Fc_A_vs_B, pvalue_A_vs_B, FDR_A_vs_B,
		   Fc_B_vs_A, pvalue_B_vs_A, FDR_B_vs_A;
		in>>b.chrom>>b.start>>b.end>>Readcount_A>>b.Normalized_Readcount_A>>ReadcountB
			>>b.Normalized_Readcount_B>>Fc_A_vs_B>>pvalue_A_vs_B>>FDR_A_vs_B
			>>Fc_B_vs_A>>pvalue_B_vs_A>>FDR_B_vs_A;
		//cout<<b.chrom<<":"<<b.start<<"-"<<b.end<<endl;
		//exit(1);
		if( !b.chrom.empty() ){
			//if( Readcount_A > 5 && ReadcountB > 5)
			{
				if( Fc_A_vs_B > fc && pvalue_A_vs_B < p ){
					b.info = "decr";
					status_beds[b.chrom].insert(b);
//					out<<(++NN)<<'\t'<<Normalized_Readcount_A<<'\t'<<Normalized_Readcount_B<<"\tdecr"<<endl;
				}else if( Fc_B_vs_A > fc && pvalue_B_vs_A < p ){
					b.info = "incr";
					status_beds[b.chrom].insert(b);
//					out<<(++NN)<<'\t'<<Normalized_Readcount_A<<'\t'<<Normalized_Readcount_B<<"\tincr"<<endl;
				}else{
					b.info = "nchg";
					//if( Fc_A_vs_B > 0.8 && Fc_A_vs_B < 1.2 && pvalue_A_vs_B > 0.5 )
						status_beds[b.chrom].insert(b);
					//else
					//	status_beds[b] = "oth";
//					out<<(++NN)<<'\t'<<Normalized_Readcount_A<<'\t'<<Normalized_Readcount_B<<"\tnchg"<<endl;
				}
			}
		}
	}
	return status_beds;
}

map<Str, set<BED3> > intersect(map<Str,set<BED3> > bed1, map<Str,set<BED3> > bed2 ){
        map<Str, set<BED3> > beds = merge( bed1, bed2 ), rst;
        assignIsland2Island(beds, bed1, "bed1");
        assignIsland2Island(beds, bed2, "bed2");

        for( map<Str, set<BED3> >::iterator chrIt = beds.begin(); chrIt != beds.end(); ++chrIt ){
                for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
                        if( bIt->islandsOn.size() == 2 ){
                                rst[bIt->chrom].insert(*bIt);
                        }
                }
        }
        cout<<"\t"<<getSize(bed1)<<"^"<<getSize(bed2)<<"="<<getSize(rst)<<endl;
        return rst;
}

void generalHeatmap(map<Str, vector<BED3> > chr_genes_tmp, Str readf, Str out_f, int upWndN, int dwnWndN, int wndSz, int shift, int bins, Str lg){
	//densities
	vector<double> density(upWndN+dwnWndN, 0);
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);

	map<Str, vector<BED3> > chr_genes;
	for( map<Str, vector<BED3> >::iterator chrIt = chr_genes_tmp.begin(); chrIt != chr_genes_tmp.end(); ++chrIt ){
		for( vector<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			BED3 b = *bIt;
			b.start = bIt->sumit - upWndN*wndSz;
			b.end = bIt->sumit + dwnWndN*wndSz;
			chr_genes[b.chrom].push_back(b);
		}
	}

	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end, n= 0, rn = 0;
	vector<pair<double, pair<BED3, vector<double> > > > indx_bed_den;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty() ){
			int pos = strand =="+" ? beg + shift : end - shift;
			++rn;
			chr_reads[chrom].push_back(pos);
			if( (++n) % 100000000 == 0 || in.eof()){
				cout<<"# reads in: "<<n<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, vector<BED3> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<BED3>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
						if( *curRd >= gIt->start){
							while( *curRd >= gIt->start && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						vector<double> den(upWndN+dwnWndN,0.001);
						while( *curRd < gIt->end && curRd != chr_reads[chrIt->first].end()){
							if( *curRd > gIt->start ){
								if( gIt->strand == "+"){
									unsigned int index = (*curRd - gIt->start)/wndSz;
									if( index < den.size())
										den[index] += 1;
								}else{
									unsigned int index = (gIt->end-*curRd)/wndSz;
									if( index < den.size())
										den[index] += 1;
								}
							}
							++curRd;
						}
						indx_bed_den.push_back(make_pair(gIt->island_score, make_pair(*gIt,den)));
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();
	int lib_sz = n;
    cout<<"lib sz: "<<lib_sz<<endl;
    cout<<"gene sz: "<<indx_bed_den.size()<<endl;

    sort( indx_bed_den.begin(), indx_bed_den.end() );

    std::ofstream out(out_f.data());
    int step = indx_bed_den.size()/bins;
    cout<<"per bin"<< step<<endl;
    //      int max = step == 1 ? 1 :  indx_bed_den.size() - step/2;
    for( int i = 0; i < indx_bed_den.size() - step; i += step ){
            vector<double> sum(indx_bed_den[i].second.second.size(), 0);
            for( int j = 0; j < step; ++j ){
                    vector<double> den = indx_bed_den[i+j].second.second;
                    for( int k = 0; k < den.size(); ++k ){
								sum[k] += den[k];
                    }
            }
            for( int k = 0; k < sum.size(); ++k )
                if( lg == "y" )
                	out<<"\t"<<log(sum[k]*100000000./wndSz/lib_sz/step)/log(2);
                else
                	out<<"\t"<<sum[k]*100000000./wndSz/lib_sz/step;
            out<<endl;
    }
    out.close();
}

void generalHeatmap_PE(map<Str, vector<BED3> > chr_genes_tmp, Str readf, Str out_f, int upWndN, int dwnWndN, int wndSz, int bins){
	//densities
	vector<double> density(upWndN+dwnWndN, 0);
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);

	map<Str, vector<BED3> > chr_genes;
	for( map<Str, vector<BED3> >::iterator chrIt = chr_genes_tmp.begin(); chrIt != chr_genes_tmp.end(); ++chrIt ){
		for( vector<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			BED3 b = *bIt;
			b.start = bIt->sumit - upWndN*wndSz;
			b.end = bIt->sumit + dwnWndN*wndSz;
			chr_genes[b.chrom].push_back(b);
		}
	}

	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end, n= 0, rn = 0;
	vector<pair<double, pair<BED3, vector<double> > > > indx_bed_den;
	while( true ){
		in>>chrom>>beg>>end;
		if( !chrom.empty() ){
			int pos = (beg + end)/2;
			++rn;
			chr_reads[chrom].push_back(pos);
			if( (++n) % 50000000 == 0 || in.eof()){
				cout<<"# reads in: "<<n<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, vector<BED3> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<BED3>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
						if( *curRd >= gIt->start){
							while( *curRd >= gIt->start && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						vector<double> den(upWndN+dwnWndN,0.001);
						while( *curRd < gIt->end && curRd != chr_reads[chrIt->first].end()){
							if( *curRd > gIt->start ){
								if( gIt->strand == "+"){
									unsigned int index = (*curRd - gIt->start)/wndSz;
									if( index < den.size())
										den[index] += 1;
								}else{
									unsigned int index = (gIt->end-*curRd)/wndSz;
									if( index < den.size())
										den[index] += 1;
								}
							}
							++curRd;
						}
						indx_bed_den.push_back(make_pair(gIt->island_score, make_pair(*gIt,den)));
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();
	int lib_sz = n;
    cout<<"lib sz: "<<lib_sz<<endl;
    cout<<"gene sz: "<<indx_bed_den.size()<<endl;

if( bins != 1 ){
    sort( indx_bed_den.begin(), indx_bed_den.end() );

    std::ofstream out(out_f.data());
    int step = indx_bed_den.size()/bins;
    cout<<"per bin"<< step<<endl;
    //      int max = step == 1 ? 1 :  indx_bed_den.size() - step/2;
    for( int i = 0; i < indx_bed_den.size() - step; i += step ){
            vector<double> sum(indx_bed_den[i].second.second.size(), 0);
            for( int j = 0; j < step; ++j ){
                    vector<double> den = indx_bed_den[i+j].second.second;
                    for( int k = 0; k < den.size(); ++k ){
                            sum[k] += den[k];
//                            sum[k] += ((den[k])*100000000./lib_sz);
                    }
            }
            for( int k = 0; k < sum.size(); ++k )
                    out<<"\t"<<(sum[k]*100000000./wndSz/lib_sz/step);
            out<<endl;
    }
    out.close();
	}else{
	    std::ofstream out(out_f.data());
		out<<"\tindex\tden\n";
	    for( int j = 0; j < indx_bed_den[0].second.second.size(); ++j ){
		double sum = 0;
		for( int i = 0; i < indx_bed_den.size(); ++i ){
                    sum += indx_bed_den[i].second.second[j];
                }
		out<<j<<"\t"<<j<<"\t"<<(sum*100000000./wndSz/lib_sz/indx_bed_den.size())<<endl;
            }

	}
}


void generateTSSTESReadDensity4Profiles_Low_preGrp(const map<Str, vector<UCSC> >& chr_genes_tmp, Str readf, Str out_f, int shift){
	//windows
	int promoterWindow = 500;
	int endWindow = 500;
	//number of windows
	int promoterSliceN = PROMOTERUPLEN*2/promoterWindow;
	int endSliceN = PROMOTERUPLEN/endWindow;
	int genebodySliceN = 10;
	//densities
	vector<double> prom_density(promoterSliceN, 0);
	vector<double> body_density(genebodySliceN, 0);
	vector<double> end_density(endSliceN, 0);
	map<Str, vector<UCSC> > chr_genes;
	for( map<Str, vector<UCSC> >::const_iterator chrIt = chr_genes_tmp.begin(); chrIt != chr_genes_tmp.end(); ++chrIt ){
		for( vector<UCSC>::const_iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
			if(gIt->txStart + PROMOTERUPLEN + 1000 < gIt->txEnd)
			{
				UCSC ucsc = *gIt;
				ucsc.txStart -= PROMOTERUPLEN;
				ucsc.txEnd += PROMOTERUPLEN;
				chr_genes[chrIt->first].push_back(ucsc);
			}
		}
	}
	cout<<"UCSC gene with gene body > 6000 bps:"<<getSize(chr_genes)<<endl;
	std::ifstream in(readf.data());
	if( !in.good() )
		FileNotFoundERR(readf);
	map<Str, vector<int> > chr_reads;
	Str chrom, strand;
	int beg, end, n= 0, rn = 0;
	map<double, vector<pair<UCSC, vector<double> > > > indx_bed_den;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty() ){
			int pos = strand =="+" ? beg + shift : end - shift;
			++rn;
			chr_reads[chrom].push_back(pos);
			if( (++n) % 50000000 == 0 || in.eof()){
				cout<<"# reads in: "<<n<<endl;
				//sort the reads to speed up
				for( map<Str, vector<int> >::iterator chrIt = chr_reads.begin(); chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());
				//map reads to region
				for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
					vector<int>::const_iterator curRd = chr_reads[chrIt->first].begin();
					for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
						if( *curRd >= gIt->txStart){
							while( *curRd >= gIt->txStart && curRd != chr_reads[chrIt->first].begin())
								--curRd;
						}
						int bodyWindow = (gIt->txEnd - gIt->txStart + 1 -3*PROMOTERUPLEN)/genebodySliceN;
						//cout<<bodyWindow<<endl;
						//exit(1);
						vector<double> pro(promoterSliceN,0.000001);
						vector<double> gb(genebodySliceN,0.000001);
						vector<double> end(endSliceN,0.000001);
						while( *curRd < gIt->txEnd && curRd != chr_reads[chrIt->first].end()){
							if( *curRd > gIt->txStart ){
								if( gIt->strand == "+"){
									if( *curRd < gIt->txStart + 2*PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart)/promoterWindow;
										if( index < prom_density.size())
											pro[index] += 1./promoterWindow;
									}else if(*curRd < gIt->txEnd - PROMOTERUPLEN){
										unsigned int index = (*curRd - gIt->txStart-2*PROMOTERUPLEN)/bodyWindow;
									//	cout<<bodyWindow<<endl;exit(1);
										if(index<body_density.size())
											gb[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (*curRd - gIt->txEnd + PROMOTERUPLEN)/endWindow;
										if( index < end_density.size() )
											end[index] += 1./endWindow;
									}//*/
								}else{
									if( *curRd > gIt->txEnd - 2*PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd-*curRd)/promoterWindow;
										if( index < prom_density.size())
											pro[index] += 1./promoterWindow;
									}else if( *curRd > gIt->txStart + PROMOTERUPLEN){
										unsigned int index = (gIt->txEnd - 2*PROMOTERUPLEN-*curRd)/bodyWindow;
										if( index < body_density.size() )
											gb[index] += 1./bodyWindow;//gene body
									}else{
										unsigned int index = (gIt->txStart+PROMOTERUPLEN-*curRd)/endWindow;
										if( index < end_density.size())
											end[index] += 1./endWindow;
									}//*/
								}
							}
							++curRd;
						}
						vector<double> den = pro;
						for( int i = 0; i < gb.size(); ++i)
							den.push_back(gb[i]);
						for( int i = 0; i < end.size(); ++i)
							den.push_back(end[i]);
						indx_bed_den[gIt->island_score].push_back(make_pair(*gIt,den));
					}
				}
				chr_reads.clear();
				if( in.eof() )
					break;
			}
		}
	}
	in.close();
	int lib_sz = n;
    vector<pair<double, vector<double> > > aves;
    for( map<double, vector<pair<UCSC, vector<double> > > >::iterator grpIt = indx_bed_den.begin(); grpIt != indx_bed_den.end(); ++grpIt ){
            vector<pair<UCSC, vector<double> > > genes = grpIt->second;
		cout<<"grp "<<grpIt->first<<"\t"<<genes.size()<<endl;
	    vector<double> ave(genes[0].second.size(), 0);
            for( int j = 0; j < genes.size(); ++j ){
                    for( int k = 0; k < genes[j].second.size(); ++k ){
//                            ave[k] += log((genes[j].second[k])*100000000./lib_sz)/log(2);
                            ave[k] += ((genes[j].second[k])*100000000./lib_sz/genes.size());
                    }
            }
		aves.push_back(make_pair(grpIt->first, ave ));
    }

    	std::ofstream out(out_f.data());
  	for( int i = 0; i < aves.size(); ++i )
		out<<"\t"<<aves[i].first;
	out<<endl;

	for( int i = 0; i < aves[0].second.size(); ++i ){
		out<<(i+1);
		for( int j = 0; j < aves.size(); ++j )
			out<<"\t"<<aves[j].second[i];
		out<<endl;
    	}
    	out.close();
}

map<Str, set<UCSC, UCSC::sortByTxStart> > UCSCGenesFrmGTF(Str fileName){
	Str line;
	map<Str, TMP> tid_tscpt;
	std::ifstream gtfIn(fileName.data());
	while( !gtfIn.eof()){
		TMP tmp;
		pair<int,int> exon;
		gtfIn>>tmp.chrom>>line>>line>>exon.first>>exon.second>>line>>tmp.strand>>line>>line;
		if( !tmp.chrom.empty() ){
			if( line != "gene_id" ){
				cout<<"Format ERROR on gene_id"<<endl;
				exit(1);
			}else{
				gtfIn>>tmp.gid;
				gtfIn>>line;
				if( line != "transcript_id" ){
					cout<<"Format ERROR on transcript_id"<<endl;
					exit(1);
				}else{
					Str id;
					gtfIn>>id;
					std::map<Str, TMP>::iterator i = tid_tscpt.find(id);
					if( i == tid_tscpt.end() ){
						tmp.exons.insert(exon);
						tid_tscpt[id] = tmp;
					}else{
						if( i->second.gid != tmp.gid ){
							cout<<"ERR multiple GID"<<endl;
							exit(1);
						}else{
							if( i->second.exons.find(exon) != i->second.exons.end() ){
							//	cout<<"ERR multiple same exon"<<"\t"<<tmp.gid<<"\t"<<exon.first<<"\t"<<exon.second<<endl;
							//	exit(1);
							}else{
								i->second.exons.insert(exon);
							}
						}
					}//*/
				}
			}
		}
		std::getline(gtfIn, line);
	}
	gtfIn.close();

	map<Str, set<UCSC, UCSC::sortByTxStart> > genes;
	for( map<Str, TMP>::iterator tIt = tid_tscpt.begin(); tIt != tid_tscpt.end(); ++tIt ){
		UCSC gene;
		gene.chrom = tIt->second.chrom;
		gene.info = gene.name = tIt->first;
		gene.alia = gene.alignID = tIt->second.gid;
		gene.strand = tIt->second.strand;
		vector<pair<int,int> > exons;
		for( set<pair<int,int> >::iterator eIt = tIt->second.exons.begin(); eIt != tIt->second.exons.end(); ++eIt ){
			exons.push_back(*eIt);
		}
		std::sort(exons.begin(), exons.end());
		gene.txStart = exons[0].first;
		gene.txEnd = exons[exons.size()-1].second;
		gene.exons = tIt->second.exons;
		gene.TSS = gene.strand == "+" ? gene.txStart : gene.txEnd;
		genes[gene.chrom].insert(gene);
//		if( gene.info.find("272329") != std::string::npos )
//							cout<<gene.info<<endl;
	}
	return genes;
}
