/*
 * HiC.cpp
 *
 *  Created on: Mar 3, 2015
 *      Author: hugq
 */

#include "HiC.h"
#include "OftenUsedOperatLib.h"
#include "UCSCBEDOperator.h"

III::III(Str fn, int bin, int minDis, int maxDis, int minLinkRd){
	III::bin = bin;
	III::minDis = minDis;
	III::maxDis = maxDis;
	III::minLinkRd = minLinkRd;
	totalPETs = 0;
	loadIII(fn);
}

III::III(Str fn, int bin, int minDis, int maxDis, Str stepfnc){
	III::bin = bin;
	III::minDis = minDis;
	III::maxDis = maxDis;
	III::minLinkRd = 1;
	totalPETs = 0;
	loadIII(fn, stepfnc);
}

void III::loadIII(Str fn){
	cout<<"load III  "<<fn<<endl;
	std::ifstream in(fn.data());
	if(  !in.good() )
		cout<<"ERR "<<fn<<endl;

	int totalII = 0, totalCnt = 0;
	while( !in.eof() ){
		Str chr1;
		int pos1, pos2, cnt;
		in>>chr1;
		if( !chr1.empty() && chr1.find("chrUn") == std::string::npos
				&& chr1.find("and") == std::string::npos //&& chr1.find("chrM") == std::string::npos
				&& chr1.find("chrY") == std::string::npos
				//&& chr1 == "chr1"
		){
			in>>pos1>>pos2>>cnt;
			if( cnt >= minLinkRd ){
				if( pos2 - pos1 >= minDis && pos2 - pos1 <= maxDis ){
					++totalII;
					if( totalII%1000000 == 0 )
						cout<<chr1<<"\t"<<pos1<<"\t"<<totalII<<endl;
					totalCnt += cnt;
					pair<int,int> pos = make_pair(pos1/bin,pos2/bin);
					map<Str,map<pair<int,int>, int > >::iterator chrIt = iii.find(chr1);
					if( chrIt != iii.end() ){
						map<pair<int,int>, int >::iterator posIt = chrIt->second.find(pos);
						if( posIt != chrIt->second.end() ){
							posIt->second += cnt;
						}else{
							chrIt->second[pos] = cnt;
						}
					}else{
						iii[chr1][pos] = cnt;
					}
				}
			}
		}else{
			std::getline(in, chr1);
		}
	}
	in.close();
	cout<<"Total II: "<<totalII<<endl;
	cout<<"Total Pet cnt: "<<totalCnt<<endl;
	cout<<"End loading\n";
	totalPETs = totalCnt;
}

void III::loadIII(Str fn, Str stepfnc){
	map<std::pair<int,int>, int> rgn_threshold;
	std::ifstream tIn(stepfnc.data());
	int bsz;
	tIn>>bsz;
	if( bsz != bin){
		std::cout<<"bin sz does not match step function"<<endl;
		exit(1);
	}
	while( !tIn.eof() ){
		std::pair<int,int> rgn;
		int threshold;
		tIn>>rgn.first>>rgn.second>>threshold;
		rgn_threshold[rgn] = threshold;
	}
	tIn.close();

	map<int,int> dis_n;
	effectivePETs = 0;
	int totalII = 0;
	int N = 0;
	cout<<"III "<<fn<<endl;
	cout<<"III stp func "<<stepfnc<<endl;
	cout<<"load III ...\n";
	std::ifstream in(fn.data());
	if(  !in.good() )
		cout<<"ERR "<<fn<<endl;
	while( !in.eof() ){
		Str chr1;
		int pos1, pos2, cnt;
		in>>chr1>>pos1>>pos2>>cnt;
		if( !chr1.empty() //&& chr1 == "chr1"
				){
			int dis = pos2 - pos1;
			if( dis >= minDis && dis <= maxDis ){
				totalPETs += cnt;
				dis /= bin;
				bool hit = false;
				for( map<std::pair<int,int>, int>::iterator rIt = rgn_threshold.begin();
						rIt != rgn_threshold.end(); ++rIt ){
					if( dis > rIt->first.first && dis <= rIt->first.second ){
						if( cnt >= rIt->second ){
							hit = true;
							break;
						}
					}
				}
				if( hit ){
					pair<int,int> pos = make_pair(pos1/bin,pos2/bin);
					++dis_n[pos2/bin-pos1/bin];
					iii[chr1][pos] = cnt;
					++totalII;
					effectivePETs += cnt;
				}
			}
		}
	}
	in.close();
	cout<<"End loading\n";
	cout<<"Total intra-chrom PETs within 2K-2M: "<<totalPETs<<endl;
	cout<<"Total II: "<<totalII<<endl;
	cout<<"Total Pet with interaction: "<<effectivePETs<<endl;
	cout<<"% intra-chrom PETs within interaction %: "<<double(effectivePETs)/totalPETs*100<<endl;
	for( map<int,int>::iterator dIt = dis_n.begin(); dIt != dis_n.end();  ++dIt )
		if( dIt->first <= 20 ){
		//	cout<<"\t"<<dIt->first<<"\t"<<dIt->second<<endl;
		}

}

void III::random_shuffle(std::map<Str, int> chr_len){
//	cout<<"Start shuffle"<<endl;
	srand((unsigned)time(0));
	map<Str,map<pair<int,int>, int > > _iii;
	for( map<Str,map<pair<int,int>, int > >::iterator chrIt = iii.begin(); chrIt != iii.end(); ++chrIt ){
		map<Str, int>::iterator lIt = chr_len.find(chrIt->first);
		if( lIt != chr_len.end() ){
			int chrLen = lIt->second/bin;
			for( map<pair<int,int>, int >::iterator ppIt = chrIt->second.begin(); ppIt != chrIt->second.end(); ++ppIt ){
				for( int t = 0; t < ppIt->second; ++t ){
					int len = (ppIt->first.second - ppIt->first.first);
					if(len < 0 )
						cout<<"ERR !"<<endl;
					bool fit = false;
					while( !fit ){
						int pos = rand()%chrLen;
						pair<int,int> pp = make_pair(pos,pos+len);
						if( pos + len < chrLen){
							map<Str,map<pair<int,int>, int > >::iterator cIt = _iii.find(chrIt->first);
							if( cIt != _iii.end() ){
								map<pair<int,int>, int >::iterator pIt = cIt->second.find(pp);
								if(  pIt == cIt->second.end() ){
									cIt->second[make_pair(pos,pos+len)] += 1;
									fit = true;
								}
							}else{
								_iii[chrIt->first][make_pair(pos,pos+len)] = 1;
								fit = true;
							}
						}
					}
				}
			}
		}
	}
	iii = _iii;
//	cout<<"End shuffle"<<endl;
}

void III::random_shuffle_full(std::map<Str, int> chr_len){
//	cout<<"Start shuffle"<<endl;
//	std::ofstream out("/home/hugq/123.txt");
	srand((unsigned)time(0));
	map<Str,map<pair<int,int>, int > > _iii;
	for( map<Str,map<pair<int,int>, int > >::iterator chrIt = iii.begin(); chrIt != iii.end(); ++chrIt ){
		map<Str, int>::iterator lIt = chr_len.find(chrIt->first);
		if( lIt != chr_len.end() ){
			int chrLen = lIt->second/bin;
			for( map<pair<int,int>, int >::iterator ppIt = chrIt->second.begin(); ppIt != chrIt->second.end(); ++ppIt ){
				for( int t = 0; t < ppIt->second; ++t ){
	                //double pcnt = (rand()%10000/100000.);
					int len = (ppIt->first.second - ppIt->first.first);
	                //len += int((rand()%2*2-1)*pcnt*len+0.5);
					if(len < 0 )
						cout<<"ERR !"<<endl;
					bool fit = false;
					while( !fit ){
						int pos = rand()%(chrLen-len);
						pair<int,int> pp = make_pair(pos,pos+len);
						if( pos + len < chrLen){
							map<Str,map<pair<int,int>, int > >::iterator cIt = _iii.find(chrIt->first);
							if( cIt != _iii.end() ){
								map<pair<int,int>, int >::iterator pIt = cIt->second.find(pp);
								if(  pIt != cIt->second.end() ){
									pIt->second += 1;
								}else{
									cIt->second[make_pair(pos,pos+len)] += 1;
								}
							}else{
								_iii[chrIt->first][make_pair(pos,pos+len)] += 1;
							}
							/*out<<chrIt->first<<":"<<ppIt->first.first
									<<"-"<<ppIt->first.second
									<<"\t"<<chrIt->first<<":"<<pos<<"-"<<(pos+len)
									<<"\t"<<(ppIt->first.second - ppIt->first.first)
									<<"\t"<<len
									<<"\t"<<double(len - (ppIt->first.second - ppIt->first.first))/(ppIt->first.second - ppIt->first.first)<<endl;//*/
							fit = true;
						}
					}
				}
			}
		}
	}
	iii = _iii;
//	cout<<"End shuffle"<<endl;
}

void III::random_shuffle_full_gcmap(std::map<Str, int> chr_len, Str gcfn, Str mpfn ){
	cout<<"Start shuffle"<<endl;
	map<Str,std::vector<pair<double,double> > > chr_gc_mp;
	for( map<Str, int>::iterator cIt = chr_len.begin(); cIt != chr_len.end(); ++cIt ){
		chr_gc_mp[cIt->first] = std::vector<pair<double,double> >(cIt->second/bin,make_pair(-1,-1));
	}
	std::ifstream gcIn(gcfn.data());
	while( !gcIn.eof() ){
		Str chr;
		int start, end;
		double gc;
		gcIn>>chr>>start>>end>>gc;
		if( !chr.empty() ){
			if( chr_gc_mp.find(chr ) != chr_gc_mp.end())
				if( chr_gc_mp[chr][start/bin].first < 0 )
					chr_gc_mp[chr][start/bin].first = gc;
				else
					chr_gc_mp[chr][start/bin].first += gc;
		}
	}
	gcIn.close();

	std::ifstream mpIn(mpfn.data());
	while( !mpIn.eof() ){
		Str chr;
		int start, end;
		double mp;
		mpIn>>chr>>start>>end>>mp;
		if( !chr.empty() ){
			if( chr_gc_mp.find(chr ) != chr_gc_mp.end())
				if( chr_gc_mp[chr][start/bin].second < 0 )
					chr_gc_mp[chr][start/bin].second = mp;
				else
					chr_gc_mp[chr][start/bin].second += mp;
		}
	}
	mpIn.close();

//	std::ofstream out123("/home/hugq/123.txt");
	map<Str,map<pair<int,int>, int > > _iii;
	for( map<Str,map<pair<int,int>, int > >::iterator chrIt = iii.begin(); chrIt != iii.end(); ++chrIt ){
		srand((unsigned)time(0));
//		if( chrIt->first != "chr1")
//			continue;
		int before_N = 0, after_N = 0, sameRegion = 0;
		cout<<chrIt->first<<endl;
		map<Str, int>::iterator lIt = chr_len.find(chrIt->first);
		map<Str,std::vector<pair<double,double> > >::iterator gcmpIt = chr_gc_mp.find(chrIt->first);
		if( lIt != chr_len.end() && gcmpIt != chr_gc_mp.end() ){
			std::vector<pair<double,double> > gcmp = gcmpIt->second;
			int chrLen = lIt->second/bin;
			for( map<pair<int,int>, int >::iterator ppIt = chrIt->second.begin(); ppIt != chrIt->second.end(); ++ppIt ){
				before_N += ppIt->second;
				if( gcmp[ppIt->first.first].first < 0 || gcmp[ppIt->first.first].second < 0
					|| gcmp[ppIt->first.second].first < 0 || gcmp[ppIt->first.second].second < 0
					|| ppIt->first.second >= chrLen )
					continue;
				for( int t = 0; t < ppIt->second; ++t ){
	                int len = (ppIt->first.second - ppIt->first.first);
	                //len += int((rand()%2*2-1)*(rand()%5000/100000.)*len+0.5);
					if(len < 0 )
						cout<<"ERR !"<<endl;
					bool fit = false;
					int tryN = 0;
					while( !fit ){
						int pos = rand()%(chrLen-len);
						pair<int,int> pp = make_pair(pos,pos+len);
						//to ensure sampling from similar gc and mappability
						if( (_abs(gcmp[ppIt->first.first].first-gcmp[pos].first) < 0.02
							&& _abs(gcmp[ppIt->first.first].second-gcmp[pos].second) < 0.05
							&& _abs(gcmp[ppIt->first.second].first-gcmp[pos+len].first) < 0.02
							&& _abs(gcmp[ppIt->first.second].second-gcmp[pos+len].second) < 0.05
							//&& ppIt->first.first != pos
							)
							|| ++tryN>1000000
							)
						{
							if( pos + len < chrLen){
									/*out123<<gcmp[ppIt->first.first].first
									<<"\t"<<gcmp[ppIt->first.first].second
									<<"\t"<<gcmp[ppIt->first.second].first
									<<"\t"<<gcmp[ppIt->first.second].second
									<<"\t"<<gcmp[pos].first
									<<"\t"<<gcmp[pos].second
									<<"\t"<<gcmp[pos+len].first
									<<"\t"<<gcmp[pos+len].second
									<<"\t"<<(ppIt->first.second - ppIt->first.first)
									<<"\t"<<len<<endl;//*/
								map<Str,map<pair<int,int>, int > >::iterator cIt = _iii.find(chrIt->first);
								if( cIt != _iii.end() ){
									map<pair<int,int>, int >::iterator pIt = cIt->second.find(pp);
									if(  pIt != cIt->second.end() ){
										pIt->second += 1;
									}else{
										cIt->second[make_pair(pos,pos+len)] += 1;
									}
								}else{
									_iii[chrIt->first][make_pair(pos,pos+len)] += 1;
								}
								fit = true;
								++after_N;
								if( pos == ppIt->first.first)
									++sameRegion;
							}
						}
					}
				}
			}
		}
		cout<<"\t"<<before_N<<"\t"<<after_N<<"\t"<<sameRegion<<"\t"<<double(sameRegion)/double(after_N)<<endl;
	}
	iii = _iii;
	cout<<"End shuffle \n";
}

/*map<Str, map<int,Str> > loadHiCIsland(Str fn, int bin){
	map<Str, map<int,Str> > chr_pos_target;
	std::cout<<"load "<<fn<<endl;
	int isld_n = 0, bin_n = 0;
	std::ifstream in1(fn.data());
	while( !in1.eof()){
		Str chr;
		int start, end;
		in1>>chr>>start>>end;
		if( !chr.empty() ){
			++isld_n;
			for( int i = start/bin; i <= end/bin; ++i ){
				chr_pos_target[chr][i] = chr+":"+convert2String(start)+"-"+convert2String(end);
				++bin_n;
			}
		}
		std::getline(in1, chr);
	}
	in1.close();
	cout<<"# island "<<isld_n<<endl;
	cout<<"\t# bin "<<bin_n<<endl;
	return chr_pos_target;
}//*/

/*map<Str, map<int,BED3> > BinningIsland(map<Str,set<BED3> > islds, int bin){
	map<Str, map<int,BED3> > chr_pos_target;
	for( map<Str,set<BED3> >::iterator chrIt = islds.begin(); chrIt != islds.end(); ++chrIt ){
		for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			//set this on for analyze bin based interaction
			//this is to ensure that one PETs cant contribute to two adjunt bins
			for( int i = bIt->start/bin; i < bIt->end/bin; ++i )
			//set this on for analysis related to foxp3 KO crisper (Figure 4 of the Nature version)
			//this is to ensure that the me1 region is within assigned bins
			//for( int i = bIt->start/bin; i <= bIt->end/bin; ++i )
			{
				chr_pos_target[bIt->chrom][i] = *bIt;
			}
		}
	}
	return chr_pos_target;
}//*/

/*map<BED3, map<BED3,int> > _assignIsland2IslandThroughIII(map<Str, set<BED3> > target_islands, map<Str, set<BED3> > source_islands, const III& iii){
	map<Str, map<int,BED3> > chr_pos_source = BinningIsland(source_islands, iii.bin);
	map<Str, map<int,BED3> > chr_pos_target = BinningIsland(target_islands, iii.bin);

	map<Str,map<int,map<int, int > > > chr_bin1_bin2_n;
	for( map<Str,map<pair<int,int>, int > >::const_iterator chrIt = iii.iii.begin(); chrIt != iii.iii.end(); ++chrIt ){
		for( map<pair<int,int>, int >::const_iterator pIt = chrIt->second.begin(); pIt != chrIt->second.end(); ++pIt ){
			chr_bin1_bin2_n[chrIt->first][pIt->first.first][pIt->first.second] = pIt->second;
			chr_bin1_bin2_n[chrIt->first][pIt->first.second][pIt->first.first] = pIt->second;
		}
	}

	map<BED3, map<BED3,int> > tgt_src_n;
	for( map<Str, map<int,BED3> >::iterator chrIt1 = chr_pos_source.begin(); chrIt1 != chr_pos_source.end(); ++chrIt1 ){
		map<Str,map<int,map<int, int > > >::iterator cIt = chr_bin1_bin2_n.find(chrIt1->first);
		map<Str, map<int,BED3> >::iterator chrIt2 = chr_pos_target.find(chrIt1->first);
		if(cIt != chr_bin1_bin2_n.end() && chrIt2 != chr_pos_target.end() ){
			for( map<int,BED3>::iterator sIt = chrIt1->second.begin(); sIt != chrIt1->second.end(); ++sIt ){
				map<int,map<int, int > >::iterator lIt1 = cIt->second.find(sIt->first);
				if( lIt1 != cIt->second.end() ){
					for( map<int, int >::iterator lIt2 = lIt1->second.begin(); lIt2 != lIt1->second.end(); ++lIt2 ){
						map<int,BED3>::iterator tIt = chrIt2->second.find(lIt2->first);
						if( tIt != chrIt2->second.end() )
						{
							tgt_src_n[tIt->second][sIt->second] += lIt2->second;
						}
					}
				}
			}
		}
	}
	return tgt_src_n;
}//*/

map<Str, map<int,set<BED3> > > BinningIsland(map<Str,set<BED3> > islds, int bin){
	map<Str, map<int,set<BED3> > > chr_pos_target;
	for( map<Str,set<BED3> >::iterator chrIt = islds.begin(); chrIt != islds.end(); ++chrIt ){
		for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			//set this on for analyze bin based interaction
			//this is to ensure that one PETs cant contribute to two adjunt bins
			int i = bIt->start;
			do
			//set this on for analysis related to foxp3 KO crisper (Figure 4 of the Nature version)
			//this is to ensure that the me1 region is within assigned bins
			//for( int i = bIt->start/bin; i <= bIt->end/bin; ++i )
			{
				chr_pos_target[bIt->chrom][i/bin].insert(*bIt);
				i += bin;
			}while( i < bIt->end );
		}
	}
	return chr_pos_target;
}

map<BED3, map<BED3,int> > _assignIsland2IslandThroughIII(map<Str, set<BED3> > target_islands, map<Str, set<BED3> > source_islands, const III& iii){
	map<Str, map<int,set<BED3> > > chr_pos_source = BinningIsland(source_islands, iii.bin);
	map<Str, map<int,set<BED3> > > chr_pos_target = BinningIsland(target_islands, iii.bin);

	map<Str,map<int,map<int, int > > > chr_bin1_bin2_n;
	for( map<Str,map<pair<int,int>, int > >::const_iterator chrIt = iii.iii.begin(); chrIt != iii.iii.end(); ++chrIt ){
		for( map<pair<int,int>, int >::const_iterator pIt = chrIt->second.begin(); pIt != chrIt->second.end(); ++pIt ){
			chr_bin1_bin2_n[chrIt->first][pIt->first.first][pIt->first.second] = pIt->second;
			chr_bin1_bin2_n[chrIt->first][pIt->first.second][pIt->first.first] = pIt->second;
		}
	}

	map<BED3, map<BED3,int> > tgt_src_n;
	for( map<Str, map<int,set<BED3> > >::iterator chrIt1 = chr_pos_source.begin(); chrIt1 != chr_pos_source.end(); ++chrIt1 ){
		map<Str,map<int,map<int, int > > >::iterator cIt = chr_bin1_bin2_n.find(chrIt1->first);
		map<Str, map<int,set<BED3> > >::iterator chrIt2 = chr_pos_target.find(chrIt1->first);
		if(cIt != chr_bin1_bin2_n.end() && chrIt2 != chr_pos_target.end() ){
			for( map<int,set<BED3> >::iterator sIt = chrIt1->second.begin(); sIt != chrIt1->second.end(); ++sIt ){
				map<int,map<int, int > >::iterator lIt1 = cIt->second.find(sIt->first);
				if( lIt1 != cIt->second.end() ){
					for( map<int, int >::iterator lIt2 = lIt1->second.begin(); lIt2 != lIt1->second.end(); ++lIt2 ){
						map<int,set<BED3> >::iterator tIt = chrIt2->second.find(lIt2->first);
						if( tIt != chrIt2->second.end() )
						{
							for( set<BED3>::const_iterator t1It = tIt->second.begin(); t1It != tIt->second.end(); ++t1It ){
								for( set<BED3>::const_iterator s1It = sIt->second.begin(); s1It != sIt->second.end(); ++s1It )
									tgt_src_n[*t1It][*s1It] += lIt2->second;
							}
						}
					}
				}
			}
		}
	}
	return tgt_src_n;
}

void assignIsland2IslandThroughIII(map<Str, set<BED3> >& target_islands,
	map<Str, set<BED3> > source_islands, const III& iii,
	Str source_island_label){
	map<BED3, map<BED3,int> > tgt_src_n = _assignIsland2IslandThroughIII(target_islands, source_islands, iii);

	map<Str, vector<BED3> > tmp = set2vector(target_islands);
	for( map<Str, vector<BED3> >::iterator chrIt = tmp.begin(); chrIt != tmp.end(); ++chrIt ){
		for( vector<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			map<BED3, map<BED3,int> >::iterator hIt = tgt_src_n.find(*bIt);
			if( hIt != tgt_src_n.end()){
				//if( hIt->second.size() == 1)
				{
					for( map<BED3,int>::iterator it = hIt->second.begin(); it != hIt->second.end(); ++it ){
						bIt->islandsOn2[source_island_label][it->first] = it->second;
					}
				}
			}
		}
	}
	target_islands = vector2set(tmp);//*/
}

void assignIsland2IslandThroughIII2BED12(map<Str, set<BED3> >& target_islands, map<Str, set<BED3> > source_islands, const III& iii, Str fn){
	map<BED3, map<BED3,int> > tgt_src_n = _assignIsland2IslandThroughIII(target_islands, source_islands, iii);

	std::ofstream out(fn.data());
	for( map<BED3, map<BED3,int> >::iterator tIt = tgt_src_n.begin(); tIt != tgt_src_n.end(); ++tIt ){
		for( map<BED3,int>::iterator sIt = tIt->second.begin(); sIt != tIt->second.end(); ++sIt ){
			BED3 s = sIt->first;
			BED3 t = tIt->first;
			if( s.end < t.start ){
				out<<s.chrom<<"\t"<<s.start<<"\t"<<t.end<<"\t"<<s.start<<"-"<<t.end<<";"<<sIt->second
						<<"\t"<<_min(100,sIt->second)/100.*1000<<"\t+\t"<<s.start<<"\t"<<t.end
						<<"\t255,0,0\t2\t"<<(s.end-s.start)<<","<<(t.end-t.start)
						<<"\t"<<s.start-s.start<<","<<t.start-s.start<<endl;
			}else if( t.end < s.start ){
				std::swap(s,t);
				out<<s.chrom<<"\t"<<s.start<<"\t"<<t.end<<"\t"<<s.start<<"-"<<t.end<<";"<<sIt->second
						<<"\t"<<_min(100,sIt->second)/100.*1000<<"\t-\t"<<s.start<<"\t"<<t.end
						<<"\t255,0,0\t2\t"<<(s.end-s.start)<<","<<(t.end-t.start)
						<<"\t"<<s.start-s.start<<","<<t.start-s.start<<endl;
			}
		}
	}
	out.close();
}

TADs::TADs(III iii, Str fn){
	map<Str,set<BED3> > _tads;
	std::ifstream in(fn.data());
	while( !in.eof() ){
		BED3 b;
		in>>b.chrom>>b.start>>b.end;
		if( !b.chrom.empty() ){
			_tads[b.chrom].insert(b);
		}
	}
	in.close();

	map<Str, set<BED3> > iiBins;
    for( map<Str,map<pair<int,int>, int > >::const_iterator chrIt = iii.iii.begin();
                chrIt != iii.iii.end(); ++chrIt ){
        for( map<pair<int,int>, int >::const_iterator pIt = chrIt->second.begin();
        		pIt != chrIt->second.end(); ++pIt ){
        	BED3 b1, b2;
        	b1.chrom = b2.chrom = chrIt->first;
        	b1.start = pIt->first.first*iii.bin;
        	b1.end = b1.start + iii.bin;
        	b2.start = pIt->first.second*iii.bin;
        	b2.end = b2.start + iii.bin;
        	iiBins[b1.chrom].insert(b1);
        	iiBins[b2.chrom].insert(b2);
        }
    }
    int N = 0, NN = 0;
    map<Str, map<Str,Str> > chr_bin_tad;
    assignIsland2Island(iiBins, _tads, "_tads");
	for( map<Str, set<BED3> >::iterator chrIt = iiBins.begin(); chrIt != iiBins.end(); ++chrIt ){
		for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			map<Str,set<BED3> >::const_iterator tIt = bIt->islandsOn.find("_tads");
			if( tIt != bIt->islandsOn.end() ){
				++N;
				if( tIt->second.size() == 1 ){
					chr_bin_tad[chrIt->first][bIt->chrom+":"+convert2String(bIt->start) +"-"+convert2String(bIt->end)]
					    = tIt->second.begin()->chrom+":"+convert2String(tIt->second.begin()->start)+"-"+convert2String(tIt->second.begin()->end);
					//cout<<bIt->chrom+":"+convert2String(bIt->start) +"-"+convert2String(bIt->end)<<endl;
					++NN;
				}
			}
		}
	}
	cout<<NN<<"\t"<<N<<"\t"<<getSize(iiBins)<<endl;

	map<Str,pair<int,int> > tad_ins_ous;
    for( map<Str,map<pair<int,int>, int > >::const_iterator chrIt = iii.iii.begin();
                chrIt != iii.iii.end(); ++chrIt ){
    	map<Str, map<Str,Str> >::const_iterator btIt = chr_bin_tad.find(chrIt->first);
    	if( btIt != chr_bin_tad.end() ){
    		for( map<pair<int,int>, int >::const_iterator pIt = chrIt->second.begin();
					pIt != chrIt->second.end(); ++pIt ){
				BED3 b1, b2;
				b1.chrom = b2.chrom = chrIt->first;
				b1.start = pIt->first.first*iii.bin;
				b1.end = b1.start + iii.bin;
				b2.start = pIt->first.second*iii.bin;
				b2.end = b2.start + iii.bin;
				Str e1 = b1.chrom + ":" + convert2String(b1.start) + "-" + convert2String(b1.end);
				Str e2 = b2.chrom + ":" + convert2String(b2.start) + "-" + convert2String(b2.end);
				map<Str,Str>::const_iterator e1It = btIt->second.find(e1);
				map<Str,Str>::const_iterator e2It = btIt->second.find(e2);
				if( e1It != btIt->second.end() && e2It != btIt->second.end() ){
					if( e1It->second == e2It->second ){
						if( tad_ins_ous.find( e1It->second ) == tad_ins_ous.end() ){
							tad_ins_ous[e1It->second] = make_pair(pIt->second,0);
						}else{
							tad_ins_ous[e1It->second].first += pIt->second;
						}
					}else{
						if( tad_ins_ous.find( e1It->second ) == tad_ins_ous.end() ){
							tad_ins_ous[e1It->second] = make_pair(0,pIt->second);
						}else{
							tad_ins_ous[e1It->second].second += pIt->second;
						}
						if( tad_ins_ous.find( e2It->second ) == tad_ins_ous.end() ){
							tad_ins_ous[e2It->second] = make_pair(0,pIt->second);
						}else{
							tad_ins_ous[e2It->second].second += pIt->second;
						}
					}
				}
			}
    	}
    }

    cout<<tad_ins_ous.size()<<endl;
	for( map<Str, set<BED3> >::iterator chrIt = _tads.begin(); chrIt != _tads.end(); ++chrIt ){
		for( set<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
			TAD tad;
			tad.id = bIt->chrom+":"+convert2String(bIt->start)+"-"+convert2String(bIt->end);
			map<Str,pair<int,int> >::const_iterator it = tad_ins_ous.find(tad.id);
			if( it != tad_ins_ous.end() ){
				tad.insidesPETsN = it->second.first;
				tad.outsidesPETsN = it->second.second;
			}else{
				tad.insidesPETsN = tad.outsidesPETsN = 0;
			}
			tads[chrIt->first].insert(tad);
		}
	}
}

map<Str,set<BED3> > III::getBins(){
	map<Str,set<BED3> > bins;
	for( map<Str,map<pair<int,int>, int > >::const_iterator chrIt = iii.begin(); chrIt != iii.end(); ++chrIt ){
		for( map<pair<int,int>, int>::const_iterator pIt = chrIt->second.begin(); pIt != chrIt->second.end(); ++pIt ){
			BED3 b1, b2;
			b1.chrom = b2.chrom = chrIt->first;
			b1.start = pIt->first.first * bin;
			b1.end = b1.start + bin;
			b2.start = pIt->first.second * bin;
			b2.end = b2.start + bin;
			bins[b1.chrom].insert(b1);
			bins[b2.chrom].insert(b2);
		}
	}
	return bins;
}
