/*
 * neuma.h
 *
 *  Created on: Dec 13, 2011
 *      Author: hugq
 */

#ifndef NEUMA_H_
#define NEUMA_H_

#include "TypeDefBase.h"
#include "OftenUsedOperatLib.h"

class gLVKM{
public:
	gLVKM(){}
	bool operator <( const gLVKM & rhs )const{
		return geneid < rhs.geneid;
	}//*/
	bool operator == ( const gLVKM & rhs )const{
		return geneid == rhs.geneid;
	}//*/
	void initFromLineOfFile(std::istream &stream){
		//cout<<geneid<<endl;
		stream>>geneid>>genesymbol>>glvkm;
		stream>>gfvkm>>finalgFVK>>origgFVK>>derivedgFVK
			>>gEUMA>>isoformsN>>measured_isoformsN;
		std::getline( stream, info);
	}
	Str geneid, genesymbol, info;
	double glvkm, gfvkm;
	Str finalgFVK, origgFVK, derivedgFVK, gEUMA, isoformsN, measured_isoformsN;
};

class iLVKM{
public:
	iLVKM(){}
	bool operator <( const iLVKM & rhs )const{
		return isoformid < rhs.isoformid;
	}//*/
	bool operator == ( const iLVKM & rhs )const{
		return isoformid == rhs.isoformid;
	}//*/
	void initFromLineOfFile(std::istream &stream){
		stream>>geneid>>genesymbol>>isoformid>>ilvkm>>ifvkm;
		stream>>finaliFVK>>origiFVK>>iFVKMfactor>>iEUMA>>isoformN>>measured_isoformsN;
		std::getline( stream, info);
	}
	Str geneid, genesymbol, isoformid, info;
	double ilvkm, ifvkm;
	Str finaliFVK, origiFVK, iFVKMfactor, iEUMA, isoformN, measured_isoformsN;
};

class gFVKM{
public:
	gFVKM(){}
	bool operator <( const gFVKM & rhs )const{
		return geneid < rhs.geneid;
	}
	bool operator == ( const gFVKM & rhs )const{
		return geneid == rhs.geneid;
	}
	void initFromLineOfFile(std::istream &stream){
		stream>>geneid>>gfvk>>readcount>>gEUMA;
		std::getline( stream, info);
	}
	Str geneid, info;
	double gfvk, readcount;
	Str gEUMA;
};//*/

class iFVKM{
public:
	iFVKM(){}
	bool operator <( const iFVKM & rhs )const{
		return isoformid < rhs.isoformid;
	}
	bool operator == ( const iFVKM & rhs )const{
		return geneid == rhs.geneid;
	}
	void initFromLineOfFile(std::istream &stream){
		stream>>geneid>>isoformid>>ifvk>>readcount>>iEUMA;
		std::getline( stream, info);
	}
	Str geneid, isoformid;
	double ifvk, readcount;
	Str iEUMA, info;
};//*/

struct ISOFORM{
	Str iEUMA;
	Str ID;
	double readcount;
	double iLVKM;
};

struct GENE{
	Str gEUMA;
	Str ID;
	double readcount;
	double gLVKM;
	vector<ISOFORM> isoforms;
};

class NEUMA{
public:
	NEUMA(){};
	NEUMA(Str fh, Str label);
	map<Str, iLVKM> isoformID_ilvkm;
	map<Str, gLVKM> geneID_glvkm;
	map<Str, iFVKM> isoformID_ifvkm;
	map<Str, gFVKM> geneID_gfvkm;
	Str label;
	map<Str, GENE> genes;
	/*bool operator <( const NEUMA & rhs )const{
		return label < rhs.label;
	}
	bool operator == ( const NEUMA & rhs )const{
		return label == rhs.label;
	}//*/
};

#endif /* NEUMA_H_ */
