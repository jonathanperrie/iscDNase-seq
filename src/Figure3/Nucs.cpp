/*
 * Nucs.cpp
 *
 *  Created on: Aug 23, 2010
 *      Author: hugq
 */

#include "Nucs.h"

void Nucs::init_from_single( Str bed_f, int shift, Str strnd){
	cout<<"Read in BED ..."<<endl;
	std::ifstream bIn(bed_f.data());
	if( !bIn.good() ){
		FileNotFoundERR(bed_f);
	}
	cout<<bed_f<<endl;
//	cout<<"shift "<<shift<<endl;
	tag_sum = 0;
	all_tag_sum = 0;
	while( !bIn.eof() ){
		Str chr, s;
		int a, b;
		bIn>>chr>>a>>b>>s>>s>>s;
		++all_tag_sum;
		if( !chr.empty() && chr == chromsome){
			++tag_sum;
			if( strnd == "+" ){
				if( s == "+"){
					if( a+shift < chrlen && a+shift > 0){
						++tag_raw[a+shift];
						++tag_raw_positive[a+shift];
					//	for( int k = 0; k < 2*shift; ++k )
						//	if( a+k < chrlen )
							//	++tag_occup[a+k];
					}
				}
			}else if( strnd == "-" ){
				if( s == "-"){
					if( b-shift > 0 && b-shift <chrlen){
						++tag_raw[b-shift];
						++tag_raw_negative[b-shift];
					//	for( int k = 0; k < 2*shift; ++k )
						//	if( b-k > 0 )
							//	++tag_occup[b-k];
					}
				}
			}else{
				if( s == "+"){
					if( a+shift < chrlen && a+shift > 0){
						++tag_raw[a+shift];
						++tag_raw_positive[a+shift];
					//	for( int k = 0; k < 2*shift; ++k )
						//	if( a+k < chrlen )
							//	++tag_occup[a+k];
					}
				}
				else{
					if( b-shift > 0 && b-shift <chrlen){
						++tag_raw[b-shift];
						++tag_raw_negative[b-shift];
					//	for( int k = 0; k < 2*shift; ++k )
						//	if( b-k > 0 )
							//	++tag_occup[b-k];
					}
				}
			}
		}
	}
	bIn.close();
	factor /= tag_sum;
	cout<<"# tag: "<<tag_sum<<endl;
	cout<<"# tag (all chr): "<<all_tag_sum<<endl;
	cout<<"Begin smooth ...";
	gausian_smooth();
	cout<<" done"<<endl;
}

void Nucs::init_from_double( Str bed_f){
	cout<<"Read in BED ..."<<endl;
	std::ifstream bIn(bed_f.data());
	if( !bIn.good() ){
		FileNotFoundERR(bed_f);
	}
	cout<<bed_f<<endl;
	tag_sum = 0;
	all_tag_sum = 0;
	while( !bIn.eof() ){
		Str chr, s, t;
		int a, b;
		bIn>>chr>>a>>b;
		++all_tag_sum;
		if( !chr.empty() && chr == chromsome){
			++tag_sum;
			int pos = (a+b)/2;
			if( pos < chrlen && pos > 0)
				++tag_raw[pos];
		}
		getline( bIn, t );
	}
	bIn.close();
	factor /= tag_sum;
	cout<<"# tag: "<<tag_sum<<endl;
	cout<<"# tag (all chr): "<<all_tag_sum<<endl;
	cout<<"Begin smooth ...";
	gausian_smooth();
	cout<<" done"<<endl;
}

void Nucs::gausian_smooth(){
	int sz = tag_raw.size();
	int lB = window/2;
	int uB = sz - window/2;
	tag_smooth = vector<double>( sz, 0);
	for( int i = lB; i < uB; ++i){
		if( tag_raw[i] > 0){
			for( int j = -lB; j < lB; ++j ){
				tag_smooth[i+j] += tag_raw[i]*dis_p[_abs(j)];
			}
		}
	}
}


void Nucs::output_nucs_pos(Str fn, int window ){
	int sz = tag_smooth.size();
	double sum = 0;
	for( int i = 0; i < window; ++i ){
		sum += tag_smooth[i] - tag_thrshld[i];
	}
	std::ofstream out(fn.data());
	std::ofstream sub(Str(fn+"_sub").data());
	for( int i = window/2; i < sz - window/2; ){
		int beg = -1, end = -1;
		while( sum <= 0 && i < sz - window/2){
			++i;
			sum += tag_smooth[i+window/2] - tag_thrshld[i+window/2] -
					tag_smooth[i-window/2] + tag_thrshld[i-window/2];
		}
		vector<double> dens;
		if( sum > 0){
			beg = i;
			dens.push_back(sum);
		}

		while( sum > 0 && i < sz - window/2){
			++i;
			sum += tag_smooth[i+window/2] - tag_thrshld[i+window/2] -
					tag_smooth[i-window/2] + tag_thrshld[i-window/2];
			dens.push_back(sum);
		}

		if( sum <= 0 ){
			end = i - 1;
		}
		if( beg > 0 && end > 0 && end > beg){
			out<<chromsome<<'\t'<<beg<<'\t'<<end<<endl;
			for( int i = 0; i < dens.size(); ++ i ){
				bool peak = true;
				for( int j = -20; j < 21; ++j ){
					if( i + j > 0 && i + j < dens.size() ){
						if( dens[i+j] > dens[i] ){
							peak = false;
							break;
						}
					}
				}
				if( peak ){
					sub<<chromsome<<'\t'<<beg+i<<'\t'<<beg+i+1<<endl;
				}
			}
		}
	}
	out.close();
	sub.close();
}
void Nucs::get_graph_file(Str fn, int type ){
	if( type == 1 || type == 2){
		std::ofstream out(fn.data());
		int sz = tag_raw.size();
		int step = 10;
		for( int i = 0; i < sz - step; i += step ){
			double sum1 = 0;
			double sum2 = 0;
			for( int j = 0; j < step; ++j){
				sum1 += tag_smooth[i+j];
				sum2 += tag_raw[i+j];
			}
			if( sum1 > 0 )
				out<<chromsome<<'\t'<<i<<'\t'<<i+step<<'\t'<<(type==2?sum1/step*factor:sum2)<<endl;
		}
		out.close();//*/
	}
}

void Nucs::gausian_thrshld(int window){
	double thi = 1;
	double thr = 0;
	for( vector<double>::iterator i = dis_p.begin(); i != dis_p.end(); ++i ){
		thr += *i;
	}
	thr *= 2*thi;

	double min_thrshld = tag_sum/double(chrlen)*thr;
	int sz = tag_raw.size();
	tag_thrshld = vector<double>( sz, min_thrshld);

	double sum = 0;
	int i = 0;
	for(; i < window; ++i ){
		sum += tag_raw[i];
	}
	i = window/2;
	if( sum/window*thr > min_thrshld )
		tag_thrshld[i] = sum/window*thr;
	for( i = window/2 + 1; i < sz - window/2; ++i){
		sum -= tag_raw[i-window/2];
		sum += tag_raw[i+window/2-1];
		if( sum/window*thr > min_thrshld )
			tag_thrshld[i] = sum/window*thr;
	}
/*	vector<double> tmp = tag_thrshld;
	sort( tmp.begin(), tmp.end ());
	sum = 0;
	for( i = 0; i < tmp.size(); ++i )
		sum += tmp[i];
	int indx = 100;
	cout<<tmp[tmp.size() - indx]<<'\t'<<sum/tmp.size()<<endl;//*/
}

void Nucs::assignDensity2BEDs( map<Str, vector<BED3> >& beds, Str tag){
	for( map<Str, vector<BED3> >::iterator chrIt = beds.begin(); chrIt != beds.end(); ++chrIt ){
		if( chrIt->first == chromsome ){
			for( vector<BED3>::iterator bIt = chrIt->second.begin(); bIt != chrIt->second.end(); ++bIt ){
				bIt->gaussianDensityOn[tag] = vector<double>(bIt->end-bIt->start, 0);
				bIt->gaussianDensityOn[tag+"thrshld"] = vector<double>(bIt->end-bIt->start, 0);
				for( int i = bIt->start; i < bIt->end; ++i){
					bIt->gaussianDensityOn[tag][i-bIt->start] = tag_smooth[i];
					bIt->gaussianDensityOn[tag+"thrshld"][i-bIt->start] = tag_thrshld[i];
				}
				if( bIt->strand == "-" ){
					vector<double> t = bIt->gaussianDensityOn[tag];
					reverse( t.begin(), t.end());
					bIt->gaussianDensityOn[tag] = t;
					t = bIt->gaussianDensityOn[tag+"thrshld"];
					reverse( t.begin(), t.end());
					bIt->gaussianDensityOn[tag+"thrshld"] = t;
				}
			}
		}
	}
}
