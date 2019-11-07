/*
 * Nucs.h
 *
 *  Created on: Aug 23, 2010
 *      Author: hugq
 */

#ifndef NUCS_H_
#define NUCS_H_

#include "TypeDefBase.h"
#include "OftenUsedOperatLib.h"
#include "bed.h"

class Nucs{
public:
	void init_from_single( Str bed_f, int shift, Str strnd = "+-");
	void init_from_double( Str bed_f);
	void gausian_smooth( );
	void gausian_thrshld(int window = 5000);
	Nucs(Str len_f, Str chromsome, double sigma = 15, int window = 200){
		map<Str, int> chr_len = getChrlen(len_f);
		map<Str, int>::iterator it = chr_len.find( chromsome );
		if( it != chr_len.end() ){
			chrlen = it->second;
//			cout<<"chr len: "<<chrlen<<endl;
			tag_raw = vector<double>(chrlen,0);
			//tag_occup = vector<double>(chrlen,0);
			tag_raw_positive = vector<double>(chrlen,0);
			tag_raw_negative = vector<double>(chrlen,0);
//			cout<<"initialization done\n";
			double sum = 0;
			for( map<Str, int>::iterator it = chr_len.begin(); it != chr_len.end(); ++it ){
				sum += it->second;
			}
			factor = chrlen/sum*1000000;
		}else{
			cout<<chromsome <<" has no match in "<<len_f<<endl;
			exit(1);
		}
		dis_p = vector<double>(window/2+1, 0);
		for( int d = 0; d < dis_p.size(); ++d){
			dis_p[d] = exp( -(d*d)/sigma/sigma );
		}
		this->sigma = sigma;
		this->window = window;
		this->chromsome = chromsome;
	}
	void assignDensity2BEDs( map<Str, vector<BED3> >& beds, Str tag);
	void get_graph_file(Str fn, int type = 0);
	void output_nucs_pos(Str fn, int window = 10);
	vector<double> tag_raw, tag_raw_positive, tag_raw_negative, tag_occup;
	vector<double>tag_smooth;
	vector<double> tag_thrshld;
	vector<double> dis_p;
	int window, tag_sum, all_tag_sum, chrlen;
	double sigma, factor;
	Str chromsome;
};

#endif /* NUCS_H_ */
