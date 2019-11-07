/*
 * HiC.h
 *
 *  Created on: Mar 3, 2015
 *      Author: hugq
 */

#ifndef HIC_H_
#define HIC_H_

#include "TypeDefBase.h"
#include "OftenUsedOperatLib.h"
#include "bed.h"

//intra-chromosomal interaction index
class III{
public:
	III(Str fn, int bin, int minDis, int maxDis, int minLinkRd);
	III(Str fn, int bin, int minDis, int maxDis, Str stepfnc);
	int bin, minDis, maxDis, minLinkRd, totalPETs, effectivePETs;
	map<Str,map<pair<int,int>, int > > iii;
	void random_shuffle(std::map<Str, int> chr_len);
	void random_shuffle_full(std::map<Str, int> chr_len);
	void random_shuffle_full_gcmap(std::map<Str, int> chr_len, Str gcfn, Str mappbilityfn );
	map<Str,set<BED3> > getBins();
protected:
	void loadIII(Str fn);
	void loadIII(Str fn, Str stepfnc);
};

class TAD{
public:
	Str id;
	int insidesPETsN, outsidesPETsN;
	bool operator <( const TAD & rhs )const{
		return id < rhs.id;
	}//*/
};

class TADs{
public:
	TADs(III iii, Str fn);
	map<Str,set<TAD> > tads;
};

map<Str, map<int,set<BED3> > > BinningIsland(map<Str,set<BED3> > islds, int bin);
void assignIsland2IslandThroughIII(map<Str, set<BED3> >& target_islands, map<Str, set<BED3> > source_islands, const III& iii, Str source_island_label);
void assignIsland2IslandThroughIII2BED12(map<Str, set<BED3> >& target_islands, map<Str, set<BED3> > source_islands, const III& iii, Str fn);
map<BED3, map<BED3,int> > _assignIsland2IslandThroughIII(map<Str, set<BED3> > target_islands, map<Str, set<BED3> > source_islands, const III& iii);
//map<Str,map<pair<int,int>, int > > loadIII(Str fn, int bin, int minDis, int maxDis, int minLinkRd);
//map<Str, map<int,Str> > loadHiCIsland(Str fn, int bin);

#endif /* HIC_H_ */
