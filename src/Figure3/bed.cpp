/*
 * BED.cpp
 *
 *  Created on: 2009-9-26
 *      Author: appleapple
 */
#include "bed.h"

std::pair< map<Str, vector<BED3> >,map<Str, vector<BED3> > >
getOverlap(map<Str, vector<BED3> > chrom_islands1, map<Str, vector<BED3> > chrom_islands2){
	std::pair<map<Str, vector<BED3> >,map<Str, vector<BED3> > > common_unique;
	map<Str, vector<BED3> >::iterator it = chrom_islands1.begin();
	for( ; it != chrom_islands1.end(); ++it ){
		for( vector<BED3>::iterator it1 = it->second.begin();
				it1 != it->second.end(); ++it1 ){
			const vector<BED3>& islands2 = chrom_islands2[it1->chrom];
			vector<BED3>::const_iterator it2 = islands2.begin();
			for( ; it2 != islands2.end(); ++it2 ){
				if( it1->isOverlap(*it2) ){
					common_unique.first[it1->chrom].push_back(*it1);
					break;
				}
			}
			if( it2 == islands2.end() ){
				common_unique.second[it1->chrom].push_back(*it1);
			}
		}
	}
	return common_unique;
}

map<Str, set<BED3> > getOverlap(map<Str, set<BED3> >islands1,map<Str, set<BED3> >islands2){
	map<Str, vector<BED3> > chrom_islands1 =set2vector(islands1);
	map<Str, vector<BED3> > chrom_islands2 =set2vector(islands2);
	map<Str, vector<BED3> > rst1 = getOverlap(chrom_islands1, chrom_islands2).first;
	map<Str, vector<BED3> > rst2 = getOverlap(chrom_islands2, chrom_islands1).first;
	if( getSize(rst1) < getSize(rst2))
		return vector2set(rst1);
	else
		return vector2set(rst2);
}

std::pair<map<Str, vector<BED3> >,map<Str, vector<BED3> > > getOverlap(Str fileName1, Str fileName2){
	map<Str, vector<BED3> > islands1, islands2;

	BEDFrmFile( fileName1, islands1);

//	map<Str, set<BED3> > tmp = vector2set(islands1);
//	tmp = getIsolatedIslands( tmp, 1000);
//	islands1 = set2vector(tmp);

	BEDFrmFile( fileName2, islands2);
	return getOverlap( islands1, islands2 );
}

pair<map<int,double>,map<int,double> > tagDensityOnRegion
(map<Str, set<BED6> >& chr_tags, map<Str, set<BED4> >& chr_regions, int windowN){
	map<int,double> regionTagSameDirections, regionTagDiffDirections;
	int regionN = 0;
	map<Str, set<BED4> >::iterator it = chr_regions.begin();
	for( ; it != chr_regions.end(); ++it ){
		regionN += it->second.size();
		set<BED6>::iterator tagIt = chr_tags[it->first].begin();
		for( set<BED4>::iterator pIt = it->second.begin();
				pIt != it->second.end(); ++pIt ){
			double window = double(pIt->end - pIt->start)/windowN;
			if( window < 1)
				continue;
			//find the first tag that overlaps the BED4 on the left
			while( tagIt->end > pIt->start && tagIt != chr_tags[it->first].begin())
				--tagIt;//if current tag pointer is on the right of BED4, shift it to the right
			while( tagIt->end < pIt->start && tagIt != chr_tags[it->first].end() )
				++tagIt;//move current pointer to overlap tag
			//calculate tag # on the BED4
			Pa_I_I tag = make_pair(tagIt->start, tagIt->end);
			while( (pIt->end > tagIt->start) && (pIt->start < tagIt->end) ){
				if(pIt->strand == "+"){
					if( tagIt->strand == "+"){
						//sense region with sense tag
						if( tagIt->start > pIt->start ){
							regionTagSameDirections[int((tagIt->start - pIt->start) / window)] += 1./window;
						}
					}else{
						//sense region with anti-sence tag
						if( pIt->end > tagIt->end ){
							regionTagDiffDirections[int((tagIt->end - pIt->start) / window)] += 1./window;
						}
					}
				}else{
					if( tagIt->strand == "+"){
						//anti-sense region with sense tag
						if( tagIt->start > pIt->start ){
							regionTagDiffDirections[int((pIt->end - tagIt->start) / window)] += 1./window;
						}
					}else{
						//anti-sense region with anti-sense tag
						if( tagIt->end < pIt->end ){
							regionTagSameDirections[int((pIt->end - tagIt->end) / window)] += 1./window;
						}
					}
				}
				if((++tagIt) == chr_tags[it->first].end())
					break;
				tag = make_pair(tagIt->start, tagIt->end);
			}
		}
	}//*/
	for( map<int,double>::iterator it = regionTagSameDirections.begin(); it != regionTagSameDirections.end(); ++it ){
		it->second = it->second/regionN;
	}

	for( map<int,double>::iterator it = regionTagDiffDirections.begin(); it != regionTagDiffDirections.end(); ++it ){
		it->second = it->second/regionN;
	}
//	cout<<"# region: "<<regionN<<endl;
	return make_pair(regionTagSameDirections,regionTagDiffDirections);
}

bool isOverlap(Pa_I_I region, Str chrom, map<Str, set<BED3> >& islands){
	for( set<BED3>::const_iterator it = islands[chrom].begin(); it != islands[chrom].end(); ++it ){
		if( _max(region.first, it->start) < _min(region.second, it->end) )
			return true;
	}
	return false;
}


map<Str, set<BED3> > bed42bed3( map<Str, set<BED4> > beds4){
	map<Str, set<BED3> > bed3s;
	for( map<Str, set<BED4> >::iterator chrIt = beds4.begin(); chrIt != beds4.end(); ++chrIt){
		for( set<BED4>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			BED4 bed4 = *bIt;
			BED3 bed3;
			bed3.chrom = bed4.chrom;
			bed3.start = bed4.start;
			bed3.end = bed4.end;
			bed3.info = bed4.info;
			bed3.strand = bed4.strand;
			bed3.sp = bed4.sp;
			bed3s[chrIt->first].insert(bed3);
		}
	}
	return bed3s;
}

map<Str, set<BED4> > bed32bed4( map<Str, set<BED3> > beds3){
	map<Str, set<BED4> > bed4s;
	for( map<Str, set<BED3> >::iterator chrIt = beds3.begin(); chrIt != beds3.end(); ++chrIt){
		for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			BED3 bed3 = *bIt;
			BED4 bed4;
			bed4.chrom = bed3.chrom;
			bed4.start = bed3.start;
			bed4.end = bed3.end;
			bed4.info = bed3.info;
			bed4.strand = "+";
			bed4s[chrIt->first].insert(bed4);
		}
	}
	return bed4s;
}

map<Str, set<BED3> > mergeIslandsWithGapLessThan( map<Str, set<BED3> >& beds3, int gap){
	map<Str, set<BED3> > chr_out;
	for( map<Str, set<BED3> >::iterator chrIt = beds3.begin(); chrIt != beds3.end(); ++chrIt){
		for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ){
			set<BED3>::iterator it = bIt;
			BED3 bed = *bIt;
			for( ;it->start - bed.end < gap && it != chrIt->second.end(); ++it ){
				bed.end = it->end;
			}
			chr_out[chrIt->first].insert(bed);
			bIt = it;
		}
	}//*/
	return chr_out;
}

map<Str, set<BED3> > getIsolatedIslands( map<Str, set<BED3> >& beds3, int gap){
	map<Str, set<BED3> > chr_out;
	for( map<Str, set<BED3> >::iterator chrIt = beds3.begin(); chrIt != beds3.end(); ++chrIt){
		set<BED3>::iterator bIt = chrIt->second.begin();
		set<BED3>::iterator posIt = bIt;
		++posIt;
		if( posIt->start - bIt->end > gap )
			chr_out[chrIt->first].insert(*bIt);
		for( ++bIt; bIt != chrIt->second.end(); ++bIt){
			set<BED3>::iterator preIt = bIt;
			set<BED3>::iterator posIt = bIt;
			--preIt;
			++posIt;
			bool isolated = false;
			if( bIt->start - preIt->end > gap ){
				if( posIt != chrIt->second.end()){
					if( posIt->start - bIt->end > gap ){
						isolated = true;
					}
				}else{
					isolated = true;
				}
			}
			if( isolated )
				chr_out[chrIt->first].insert(*bIt);
		}
	}//*/
	return chr_out;
}


map<Str, set<BED3> > generateRandomIslands( map<Str, set<BED3> > templa, map<Str, int> chr_len ){
	srand((unsigned)time(0));
	map<Str, set<BED3> > chr_islands;
	for( map<Str, set<BED3> >::iterator chrIt = templa.begin(); chrIt != templa.end(); ++chrIt){
		int len = chr_len[chrIt->first];
		for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt){
			int islandL = bIt->end - bIt->start;
			int pos = rand() % (len - islandL);
			BED3 bed = *bIt;
			bed.start = pos;
			bed.end = pos + islandL;
			bed.sumit = (bed.start+bed.end)/2;
			chr_islands[bed.chrom].insert(bed);
		}
	}
	return chr_islands;
}

map<Str, set<BED3> > generateRandomIslands( int island_N, int island_len, map<Str, int> chr_len){
	srand((unsigned)time(0));
	map<Str, set<BED3> > chr_islands;

	map<int, pair<Str, int> > index_chr_len;
	int i = 0;
	for(map<Str, int>::iterator it = chr_len.begin(); it != chr_len.end(); ++it, ++i ){
		index_chr_len[i] = make_pair(it->first, it->second);
	}

//	std::ofstream out("len.txt");
//	out<<"\trnd\tn\tgrp\n";
//	int NN = 0;
	int n = 0;
	int chrN = chr_len.size();
	while( n < island_N){
		int chrIndex = rand() % chrN;
		int pos = int(rand() / double(RAND_MAX) * (index_chr_len[chrIndex].second - island_len));
		BED3 b;
		b.start = pos;
		b.end = pos + island_len;
		b.chrom = index_chr_len[chrIndex].first;
		b.sumit = (b.start+b.end)/2;
		/*bool overlap = false;
		map<Str, set<BED3> >::iterator chrIt1 = chr_islands.find(b.chrom);
		for( set<BED3>::iterator bIt1 = chrIt1->second.begin(); bIt1 != chrIt1->second.end(); ++bIt1){
			if(  _min(bIt1->end, b.end) > _max(bIt1->start, b.start) ){
				overlap = true;
				break;
			}
		}
		if( !overlap )//*/
		{
			chr_islands[b.chrom].insert(b);
			++n;
		//	cout<<n<<'\t'<<b.start<<endl;
	//		out<<(++NN)<<'\t'<<rand()<<"\t"<<b.start
		//			<<'\t'<<b.chrom<<endl;
		}
	}
	return chr_islands;
}

void BEDFrmSICCERSummary( Str fileName, map<Str, set<BED3> >& chr_bed, int shift){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	double n = 0;
	while( !in.eof() ){
		BED3 island;
		in>>island.chrom>>island.start>>island.end>>island.readInTreat>>island.readInControl
		>>island.pvalue>>island.FC>>island.FDR;
		island.shift(shift);
		if( !island.chrom.empty()// && island.FDR < 0.001 && island.FC > 1.5
				){
			++n;
			island.island_score = n;
			island.sumit = (island.start+island.end)/2;
			chr_bed[island.chrom].insert(island);
		}
	}
	cout<<"# of records from \""<<fileName<<"\": "<<n<<endl;
//	cout<<"# of records (1 read per position): "<<getSize(chr_bed)<<'('<<100*getSize(chr_bed)/n<<')'<<endl;
	in.close();
}

void BEDFrmPeakSplitter( Str fileName, map<Str, set<BED3> >& chr_bed){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	Str line;
	getline( in, line );
	double n = 0;
	while( !in.eof() ){
		BED3 island;
		in>>island.chrom>>island.start>>island.end>>island.sumitTagNum>>island.sumit;
		if( !island.chrom.empty() ){
			++n;
			chr_bed[island.chrom].insert(island);
		}
	}
	cout<<"# of records from \""<<fileName<<"\": "<<n<<endl;
//	cout<<"# of records (1 read per position): "<<getSize(chr_bed)<<'('<<100*getSize(chr_bed)/n<<')'<<endl;
	in.close();
}

void BEDFrmSICERprobscoreisland( Str fileName, map<Str, set<BED3> >& chr_bed, int shift){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	double n = 0;
	while( !in.eof() ){
		BED3 island;
		in>>island.chrom>>island.start>>island.end>>island.island_score;
		island.sumit = (island.start+island.end)/2;
	//	island.start = island.sumit - 100;
	//	island.end = island.sumit + 100;

		island.shift(shift);
		if( !island.chrom.empty() ){
			++n;
			chr_bed[island.chrom].insert(island);
		}
	}
	cout<<"# of records from \""<<fileName<<"\": "<<n<<endl;
//	cout<<"# of records (1 read per position): "<<getSize(chr_bed)<<'('<<100*getSize(chr_bed)/n<<')'<<endl;
	in.close();
}

map<Str, set<BED3> > BEDFrmtags_around_summit( Str fileName){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	cout<<"read in "<<fileName<<endl;
	map<Str, set<BED3> > chr_bed;
	while( !in.eof() ){
		Str line;
		std::getline(in,line);
		std::istringstream buf(line);//ensure one line per gene
		Str coord;
		buf>>coord;
		if( !coord.empty() ){
			//cout<<coord<<'\t';
			BED3 b = parsStr2Bed3(coord);
			buf>>b.sumit;
			//cout<<b.sumit<<endl;
			while( !buf.eof() ){
				double i;
				buf>>i;
				b.readCountsOn[fileName].push_back(i);
				//cout<<i<<'\t';
			}
			//cout<<endl;
			buf.clear();
			chr_bed[b.chrom].insert(b);
		}
	}
	in.close();
	return chr_bed;
}

void BEDFrmMACSXLS( Str fileName, map<Str, set<BED3> >& chr_bed, int shift, bool hascontrol){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	Str line;
	while( line.find("chr	start	end") == std::string::npos){
		getline(in, line);
	}
	int NN = 0;
	while( !in.eof() ){
		BED3 b;
		in>>b.chrom>>b.start>>b.end>>b.sumit>>b.sumit>>b.sumitTagNum>>b.pvalue>>b.FC;
		b.strand = "+";
//		b.start -= 150;
//		b.end += 150;
		if(hascontrol){
			in>>b.FDR;
		//	cout<<b.start<<'\t'<<b.FDR<<endl;
		}
		if( !b.chrom.empty() ){
			b.sumit += b.start;
			b.shift(shift);
			b.island_score = ++NN;
		//	if( b.pvalue > -10*log(0.0000000001)/log(10))
				chr_bed[b.chrom].insert(b);
		}
	}
	cout<<"# of records from \""<<fileName<<"\": "<<getSize(chr_bed)<<endl;
	in.close();
}

void BEDFrmMACS2XLS( Str fileName, map<Str, set<BED3> >& chr_bed){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	Str line;
	while( line.find("chr") == std::string::npos){
		getline(in, line);
	}
	while( !in.eof() ){
		BED3 b;
		Str tmp;
		in>>b.chrom>>b.start>>b.end>>tmp>>b.sumit>>tmp>>b.pvalue>>b.FC;
		std::getline(in, tmp);
		if( !b.chrom.empty() ){
			chr_bed[b.chrom].insert(b);
		}
	}
	cout<<"# of records from \""<<fileName<<"\": "<<getSize(chr_bed)<<endl;
	in.close();
}

map<Str, set<BED3> > BEDFrmFIMO( Str fileName ){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	Str line;
	getline(in, line);

	map<Str, set<BED3> > chr_bed;
	int nn = 0;
	while( !in.eof() ){
		int a, b;
		double pvalue, qvalue, score;
		Str seq;
		Str coord, strand;
		in>>coord>>coord>>a>>b>>strand>>score>>pvalue>>qvalue>>seq;
		if( !coord.empty() ){
			BED3 b3 = parsStr2Bed3(coord);
			b3.start += a;
			b3.end = b3.start + b - a;
			b3.strand = strand;
			b3.pvalue = pvalue;
			b3.seq = seq;
			b3.island_score = score;
			chr_bed[b3.chrom].insert(b3);
		}
	}
	cout<<"# of records from \""<<fileName<<"\": "<<getSize(chr_bed)<<endl;
	in.close();
	return chr_bed;
}

void BEDFrmNPS( Str fileName, map<Str, set<BED3> >& chr_bed, double p ){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	Str line;
	getline(in, line);
	int nps_l = 0;
	while( !in.eof() ){
		BED3 b;
		in>>b.chrom>>b.start>>b.end
		>>b.sumit>>b.wp
		//>>b.sumitsOn["cd34"]>>b.readDensityOn["cd34"]
		>>b.den>>b.pos_N>>b.neg_N>>line>>b.island_score;
		//b.wp = b.sumit;
		if( !b.chrom.empty() && b.island_score > -10*log(p)/log(10) //&& _abs(b.sumit-b.wp) < 20
		){
			/*if( b.end - b.start < 150 ){
				int mid = (b.start + b.end )/2;
				b.start = mid - 75;
				b.end = mid + 75;
			}//*/
			if( b.chrom != "chrY" )
			{
				chr_bed[b.chrom].insert(b);
				nps_l += b.end - b.start;
			}
		}
		getline( in, line );
	}
	cout<<"# of records from \""<<fileName<<"\": "<<getSize(chr_bed)<<endl;
//	cout<<"NPS length: "<<nps_l<<endl;
	in.close();
}

map<Str, set<BED3> > self_union( map<Str, set<BED3> > beds ){
	map<Str, set<BED3> > unions;
	for(map<Str, set<BED3> >::iterator chrIt = beds.begin(); chrIt != beds.end(); ++chrIt ){
		set<BED3>::iterator bIt = chrIt->second.begin();
		set<BED3>::iterator bIt1 = bIt;
		for( ; bIt != chrIt->second.end(); bIt = bIt1){
			int rightmst = bIt->end;
			for( ++bIt1; bIt1 != chrIt->second.end() && bIt1->start <= rightmst ; ++bIt1 ){
				if( bIt->start == 33605303 ){
					cout<<bIt1->chrom<<":"<<bIt1->start<<"-"<<bIt1->end<<endl;
				}
				if( bIt1->end > rightmst){
					rightmst = bIt1->end;
				}
			}
			BED3 b;
			b.start = bIt->start;
			b.end = rightmst;
			b.chrom = bIt->chrom;
			unions[b.chrom].insert(b);
		}
	}
	return unions;
}

map<Str, set<BED3> > get_intergenic(map<Str, set<BED3> > pro, map<Str, set<BED3> > genic, Str file ){
	pro = merge(pro, genic);
	map<Str, int> chr_len = getChrlen(file);
	map<Str, set<BED3> > inter;
	for( map<Str, set<BED3> >::iterator chrIt = pro.begin(); chrIt != pro.end(); ++chrIt ){
		set<BED3>::iterator bIt1 = chrIt->second.begin();
		set<BED3>::iterator bIt2 = bIt1;
		int lftmst = bIt1->strand == "+" ? bIt1->start - PROMOTERUPLEN: bIt1->start;
		if( lftmst > 0 ){
			BED3 b;
			b.start = 1;
			b.end = lftmst;
			b.chrom = bIt1->chrom;
			inter[b.chrom].insert(b);
		}
		for( ++bIt2; bIt2 != chrIt->second.end(); ++bIt1, ++bIt2 ){
			if( bIt2->start > bIt1->end ){
				BED3 b;
				b.start = bIt1->end;
				b.end = bIt2->start;
				b.chrom = bIt1->chrom;
				inter[b.chrom].insert(b);
			}else{
				cout<<bIt1->chrom<<'\t'<<bIt1->end<<'\t'<<bIt2->start<<endl;
			}
		}
		map<Str, int>::iterator lIt = chr_len.find(chrIt->first);
		if( lIt != chr_len.end() ){
			--bIt2;
			int rgtmst = bIt2->strand == "+" ? bIt2->end: bIt2->end + PROMOTERUPLEN;
			if( rgtmst < lIt->second){
				BED3 b;
				b.start = rgtmst;
				b.end = lIt->second;
				b.chrom = bIt2->chrom;
				inter[b.chrom].insert(b);
			}
		}else{
		//	cout<<chrIt->first<<"ERR 22"<<endl;
		}
	}
	cout<<getSize(inter)<<endl;
	inter = self_union(inter);
	cout<<getSize(inter)<<endl;
	return inter;
}

map<Str, set<BED3> > get_complement_set(map<Str, set<BED3> > bed, Str file ){
	map<Str, set<BED3> > pro = self_union(bed);
	map<Str, int> chr_len = getChrlen(file);
	map<Str, set<BED3> > nonpro;
	for( map<Str, set<BED3> >::iterator chrIt = bed.begin(); chrIt != bed.end(); ++chrIt ){
		set<BED3>::iterator bIt1 = chrIt->second.begin();
		set<BED3>::iterator bIt2 = bIt1;
		int lftmst = bIt1->strand == "+" ? bIt1->start - PROMOTERUPLEN: bIt1->start;
		if( lftmst > 0 ){
			BED3 b;
			b.start = 1;
			b.end = lftmst;
			b.chrom = bIt1->chrom;
			if( b.end - b.start > 3)
				nonpro[b.chrom].insert(b);
		}
		for( ++bIt2; bIt2 != chrIt->second.end(); ++bIt1, ++bIt2 ){
			if( bIt2->start > bIt1->end ){
				BED3 b;
				b.start = bIt1->end;
				b.end = bIt2->start;
				b.chrom = bIt1->chrom;
				if( b.end - b.start > 3)
					nonpro[b.chrom].insert(b);
			}else{
				cout<<bIt1->chrom<<'\t'<<bIt1->end<<'\t'<<bIt2->start<<endl;
			}
		}
		map<Str, int>::iterator lIt = chr_len.find(chrIt->first);
		if( lIt != chr_len.end() ){
			--bIt2;
			int rgtmst = bIt2->strand == "+" ? bIt2->end: bIt2->end + PROMOTERUPLEN;
			if( rgtmst < lIt->second){
				BED3 b;
				b.start = rgtmst;
				b.end = lIt->second;
				b.chrom = bIt2->chrom;
				if( b.end - b.start > 3)
					nonpro[b.chrom].insert(b);
			}
		}else{
		//	cout<<chrIt->first<<"ERR 22"<<endl;
		}
	}
	return nonpro;
}

map<Str, set<BED3> > merge(map<Str, set<BED3> > bed1, map<Str, set<BED3> > bed2){
	for( map<Str, set<BED3> >::iterator chrIt = bed1.begin(); chrIt != bed1.end(); ++chrIt ){
		for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			bed2[bIt->chrom].insert(*bIt);
		}
	}
	return self_union(bed2);
}

//data/phf/hg18/hg18_chrlen.txt
map<Str, set<BED3> > rand_position_in(map<Str, set<BED3> > beds, int island_N, Str chr_len_file, int length ){
	 map<Str, int> chr_len = getChrlen(chr_len_file.data());

	srand((unsigned)time(0));

	vector<pair<Str, int> > index_chr_len;
	int i = 0;
	for(map<Str, int>::iterator it = chr_len.begin(); it != chr_len.end(); ++it, ++i ){
		index_chr_len.push_back(make_pair(it->first, it->second));
	}

	int n = 0;
	int chrN = chr_len.size();
	map<Str, set<BED3> > chr_islands;
	beds = self_union(beds);
	map<Str, vector<BED3> > _beds = set2vector(beds);
	while( n < island_N){
		int chrIndex = rand() % chrN;
		Str chrom = index_chr_len[chrIndex].first;
		if( _beds.find(chrom) != _beds.end() ){
			int rIndex = int(rand()/double(RAND_MAX) * _beds[chrom].size());
			int b_L = _beds[chrom].size();
			if( rIndex < b_L){
				int region_L = _beds[chrom][rIndex].end - _beds[chrom][rIndex].start;
				if( region_L > PROMOTERUPLEN)
				{
					int pos = _beds[chrom][rIndex].start +
						+ int(rand()/double(RAND_MAX) * region_L);
					BED3 b;
					b.start = pos - length/2;
					b.end = pos + length/2;
					b.chrom = chrom;
					chr_islands[b.chrom].insert(b);
					++n;
				}
			}
		}
	}
	return chr_islands;
}

map<Str, set<BED3> > filter_region_out_of_chr_range(map<Str, set<BED3> > beds, Str file ){
	map<Str, set<BED3> > rst;
	map<Str, int> chr_len = getChrlen(file);
	for( map<Str, set<BED3> >::iterator chrIt = beds.begin(); chrIt != beds.end(); ++chrIt ){
		map<Str, int>::iterator maxlenIt = chr_len.find(chrIt->first);
		if( maxlenIt != chr_len.end() ){
			for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
				if( bIt->start > 0 && bIt->end < maxlenIt->second){
					rst[bIt->chrom].insert(*bIt);
				}
			}
		}
	}
	cout<<"From "<<getSize(beds)<<" to "<<getSize(rst)<<endl;
	return rst;
}

BED3 parsStr2Bed3(Str s){
	BED3 b;
	for( int i = 0; i < s.size(); ++i ){
		if( s[i] == ':' || s[i] =='-' || s[i] == '>')
			s[i] = '\t';
	}
	std::istringstream buf(s);//ensure one line per gene
	buf>>b.chrom>>b.start>>b.end;
	std::getline(buf, b.info);
	buf.clear();//*/
	return b;
}

map<Str, set<BED3> > BEDFrmSICCERDIFF( Str fileName ){
	map<Str, set<BED3> > beds;
	std::ifstream in(fileName.data());
	if( !in.good() ){
		cout<<"File not exist: "<<fileName<<endl;
		exit(1);
	}
	Str line;
	getline( in, line );
	while( !in.eof() ){
		BED3 b;
		in>>b.chrom>>b.start>>b.end>>b.info;
		if( !b.info.empty() ){
			int mid = (b.start+b.end)/2;
			//b.start = mid - 5000;
			//b.end = mid + 5000;
			beds[b.chrom].insert(b);
		}
		getline( in, line );
	}
	return beds;
}

map<Str, set<BED3> > BEDFrmSICCERDIFF_2( Str fileName, double fdr, double FC){
	map<Str, set<BED3> > islds;
	double Readcount_A, Normalized_Readcount_A, ReadcountB, Normalized_Readcount_B,
	Fc_A_vs_B, pvalue_A_vs_B, FDR_A_vs_B, Fc_B_vs_A, pvalue_B_vs_A, FDR_B_vs_A;
	Str line;
	std::ifstream in(fileName.data());
	std::getline(in, line);
	while( !in.eof() ){
		BED3 b;
		in>>b.chrom>>b.start>>b.end>>Readcount_A>>Normalized_Readcount_A>>ReadcountB
			>>Normalized_Readcount_B>>Fc_A_vs_B>>pvalue_A_vs_B>>FDR_A_vs_B>>Fc_B_vs_A
			>>pvalue_B_vs_A>>FDR_B_vs_A;
		Str status = "nchg";
		if( Fc_A_vs_B > FC && FDR_A_vs_B < fdr ){
			status = "decr";
		}
		if( Fc_B_vs_A > FC && FDR_B_vs_A < fdr ){
			status = "incr";
		}
		b.FC = Fc_B_vs_A;
		b.info = status;
		if( !b.chrom.empty() )
			islds[b.chrom].insert(b);
	}
	in.close();


	return islds;
}

void BEDFrmPIQ( Str fileName, map<Str, set<BED3> >& chr_bed, int len, double purity ){
	std::ifstream in(fileName.data());
	if( in.good() ){
		Str line;
		while( !in.eof() ){
			std::getline(in, line);
			if( !line.empty() ){
				if( line.find("coord") == std::string::npos ){
					strip(line, '"');
					vector<Str> token = split( line, ',');
					BED3 b;
					b.chrom = token[1];
					b.start = atoi(token[2].data());
					b.end = b.start + len;
					b.sumit = (b.start+b.end)/2;
					b.island_score = atof(token[6].data());
					if( b.island_score > purity )
						chr_bed[b.chrom].insert(b);
				}
			}
		}
	}
	in.close();
}

void BEDFrmDanposXlS( Str fileName, map<Str, set<BED3> >& chr_bed ){
	std::ifstream in(fileName.data());
	if( in.good() ){
		Str line;
		std::getline(in, line);
		while( !in.eof() ){
			BED3 b;
			in>>b.chrom>>b.start>>b.end>>b.sumit>>b.FDR>>b.pvalue;;
			chr_bed[b.chrom].insert(b);
		}
	}
	in.close();
}


map<Str, set<BED3> > BEDFrmEdgeR(Str fileName, double _FDR, double _log2FC){
	map<Str, set<BED3> > isld;
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	Str line;
	getline( in, line );
	int incrN = 0, decrN = 0;
	while( !in.eof() ){
		Str _id;
		double logFC, logConc, pvalue, FDR, LR;
		//in>>_id>>logConc>>logFC>>pvalue>>FDR;
		//in>>_id>>logFC>>logConc>>pvalue>>FDR;//edger3
		in>>_id>>logFC>>logConc>>LR>>pvalue>>FDR;//edger3, GLM
		if( !_id.empty() ){
			Str id, cg("nchg");
			for( int i = 0; i < _id.size(); ++i )
				if(_id[i] !='\"' )
					id += _id[i];
			if( FDR < _FDR ){
				if( logFC > _log2FC ){
					cg = "incr";
					++incrN;
				}
				if( logFC < -_log2FC ){
					cg = "decr";
					++decrN;
				}
			}
			BED3 b = parsStr2Bed3(id);
			b.info = cg;
			b.FDR = FDR;
			b.FC = logFC;
			b.den = logConc;
			isld[b.chrom].insert(b);
		}
	}
	cout<<"\t incr "<<incrN<<"\tdecr "<<decrN<<endl;
	return isld;
}

map<Str, set<BED3> > PartitionGenomeIntoBins(Str chrlen_fn, int binsz){
	map<Str, set<BED3> > bins;
	map<Str,int> chrlen = getChrlen(chrlen_fn);
	for( map<Str,int>::iterator cIt = chrlen.begin(); cIt != chrlen.end(); ++cIt ){
		for( int i = 0; i < cIt->second - binsz; i += binsz ){
			BED3 b;
			b.chrom = cIt->first;
			b.start = i;
			b.end = i + binsz;
			bins[b.chrom].insert(b);
		}
	}
	return bins;
}
