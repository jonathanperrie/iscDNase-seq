/*
 * neuma.cpp
 *
 *  Created on: Dec 14, 2011
 *      Author: hugq
 */

#include "neuma.h"

NEUMA::NEUMA(Str fh, Str labl){
	label = labl;
	Str line;
	std::ifstream in((fh+".30.gLVKM").data());
	if( in.good()){
		getline(in, line);//SKIP FIRST LINE
		while( !in.eof( )){
			gLVKM glvkm;
			glvkm.initFromLineOfFile(in);
			if( !glvkm.geneid.empty() && glvkm.gEUMA != "below_cut" ){
				if( geneID_glvkm.find(glvkm.geneid) == geneID_glvkm.end()){
					geneID_glvkm[glvkm.geneid] = glvkm;
				}
				else{
					cout<<"Redundant gene id?\n";
					cout<<glvkm.geneid<<"\t"<<fh<<"30.gLVKM"<<endl;
					exit(1);
				}
			}
		}
	}
	in.close();

	std::ifstream in3((fh+".gFVKM").data());
	if( in3.good()){
		getline(in3, line);//SKIP FIRST LINE
		while( !in3.eof( )){
			gFVKM gfvkm;
			gfvkm.initFromLineOfFile(in3);
			//if( gfvkm.geneid == "99890")
			//	cout<<gfvkm.readcount<<endl;
			if( !gfvkm.geneid.empty() ){
				if( geneID_gfvkm.find(gfvkm.geneid) == geneID_gfvkm.end() )
					geneID_gfvkm[gfvkm.geneid] = gfvkm;
				else{
					cout<<"Redundant gene id?\n";
					cout<<gfvkm.geneid<<"\t"<<fh<<"gFVKM"<<endl;
					exit(1);
				}
			}
		}
	}
	in3.close();

	std::ifstream in1((fh+".30.iLVKM").data());
	if( in1.good() ){
		getline(in1, line);//SKIP FIRST LINE
		while( !in1.eof( )){
			iLVKM ilvkm;
			ilvkm.initFromLineOfFile(in1);
			if( !ilvkm.isoformid.empty() ){
				if( isoformID_ilvkm.find(ilvkm.isoformid) == isoformID_ilvkm.end() )
					isoformID_ilvkm[ilvkm.isoformid] = ilvkm;
				else{
					cout<<"Redundant isoform id?\n";
					cout<<ilvkm.isoformid<<"\t"<<fh<<".30.iLVKM"<<endl;
					exit(1);
				}
			}
		}
	}
	in1.close();

	std::ifstream in2((fh+".iFVKM").data());
	if( in2.good() ){
		getline(in2, line);//SKIP FIRST LINE
		while( !in2.eof( )){
			iFVKM ifvkm;
			ifvkm.initFromLineOfFile(in2);
			if( !ifvkm.isoformid.empty() ){
				if( isoformID_ifvkm.find(ifvkm.isoformid) == isoformID_ifvkm.end() ){
					isoformID_ifvkm[ifvkm.isoformid] = ifvkm;
				}else{
					cout<<"Redundant isoform?\n";
					cout<<ifvkm.isoformid<<"\t"<<fh<<".iFVKM"<<endl;
					exit(1);
				}
			}
		}
	}
	in2.close();

	if( isoformID_ifvkm.empty() ){
		cout<<"isoformID_ifvkm is empty"<<endl;
		cout<<fh<<endl;
		exit(1);
	}else{
		for( map<Str, iFVKM>::iterator it = isoformID_ifvkm.begin(); it != isoformID_ifvkm.end(); ++it ){
			map<Str, iLVKM>::iterator ilIt = isoformID_ilvkm.find(it->second.isoformid);
			if( ilIt != isoformID_ilvkm.end() ){
				ISOFORM isoform;
				isoform.ID = it->second.isoformid;
				isoform.iEUMA = it->second.iEUMA;
				isoform.iLVKM = ilIt->second.ilvkm;
				isoform.readcount = it->second.readcount;

				if( genes.find(it->second.geneid) == genes.end()){
					GENE gene;
					gene.ID = it->second.geneid;
					map<Str, gFVKM>::iterator gfIt = geneID_gfvkm.find(gene.ID);
					map<Str, gLVKM>::iterator glIt = geneID_glvkm.find(gene.ID);
					if( gfIt == geneID_gfvkm.end() ){
					//	cout<<"gene ID missed in geneID_gfvkm "<<gene.ID<<endl;
					}else{
						if( glIt != geneID_glvkm.end() ){
							gene.gEUMA = gfIt->second.gEUMA;
							gene.readcount = gfIt->second.readcount;
							gene.gLVKM = glIt->second.glvkm;
						}
						if( it->second.isoformid.substr(0,2) != "NR")
							gene.isoforms.push_back(isoform);
						genes[gene.ID] = gene;
					}
				}else{
					if( it->second.isoformid.substr(0,2) != "NR")
						genes[it->second.geneid].isoforms.push_back(isoform);
				}
			}else{
				//if( it->second.isoformid.substr(0,2) != "NR")
				//	cout<<"LVKM missing for "<<it->second.isoformid<<endl;
			}
		}
	}

	/*for( map<Str, gLVKM>::iterator it = geneID_glvkm.begin(); it != geneID_glvkm.end(); ++it ){
		if( genes.find(it->first ) == genes.end() ){
			cout<<it->first<<endl;
		}
	}
	exit(1);//*/
}
