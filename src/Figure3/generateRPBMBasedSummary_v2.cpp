//============================================================================
//// Name        : HC2AZ.cpp
//// Author      :
//// Version     :
//// Copyright   : Your copyright notice
//// Description : Hello World in C++, Ansi-style
////============================================================================

#include "UCSCBEDOperator.h"
#include "expression.h"
#include "OftenUsedOperatLib.h"
using namespace std;

int PROMOTERUPLEN = 2500;
int PROMOTERDOWNLEN = 2500;
int MININTRESTREGIONLEN = 0;
int SHIFT = 0;

void generateRPBMSummaryGraph_slide(vector<Str> bed_fs, Str chr_len_f, Str out_f, int windowsize = 5, int shift = 75, int smooth_window_n = 21, bool normalize = true ){
	//Initiate graph file windows;
	cout<<"Window size: "<<windowsize<<endl;
	cout<<"# window for sliding: "<<smooth_window_n<<endl;
	cout<<"Fragment shift: "<<shift<<endl;
	cout<<"Normalize: "<<(normalize?"yes":"no")<<endl;
	map<Str, int> chr_len = getChrlen(chr_len_f);
	map<Str, vector<double> > chr_RPBMs;
	for( map<Str, int>::iterator it = chr_len.begin(); it != chr_len.end(); ++it ){
		chr_RPBMs[it->first] = vector<double>( it->second/windowsize, 0);
	}
	//
	int rn = 0;
	for( int i = 0; i < bed_fs.size(); ++i ){
		std::ifstream in(bed_fs[i].data());
		if( !in.good() )
			FileNotFoundERR(bed_fs[i]);
		cout<<"# reads from "<<bed_fs[i];
		int n =0;
		while( !in.eof() ){
			BED6 b;
			in>>b.chrom>>b.start>>b.end>>b.name>>b.score>>b.strand;
			if( !b.chrom.empty() ){
				int pos = b.strand =="+" ? b.start + shift : b.end - shift;
				map<Str, vector<double> >::iterator chrIt = chr_RPBMs.find(b.chrom);
				if( chrIt != chr_RPBMs.end() ){
					++n;
					unsigned int index = pos / windowsize;
					if( index < chrIt->second.size()){
						chrIt->second[index] += 1;
					}
				}
				if( n % 5000000 == 0){
					cout<<".";
				}
			}
		}
		in.close();
		rn += n;
		cout<<"\t"<<n<<endl;
	}
	cout<<"# total reads: "<<rn<<endl;
	if( normalize ){
		cout<<"normalization ..."<<endl;
		for( map<Str, vector<double> >::iterator it = chr_RPBMs.begin(); it != chr_RPBMs.end(); ++it ){
			for( unsigned int i = 0; i < it->second.size(); ++i ){
				double rpbm = it->second[i]*1000000./windowsize/rn;
				it->second[i] = rpbm;
			}
		}
		cout<<" done!"<<endl;
	}
	//smooth
	map<Str, vector<double> > _chr_RPBMs;
	for( map<Str, int>::iterator it = chr_len.begin(); it != chr_len.end(); ++it ){
		_chr_RPBMs[it->first] = vector<double>( it->second/windowsize, 0);
	}
	cout<<"smooth ...";
	for( map<Str, vector<double> >::iterator it = chr_RPBMs.begin(); it != chr_RPBMs.end(); ++it ){
		cout<<it->first<<" ";
		for( unsigned int i = 0; i < it->second.size() - smooth_window_n; ++i ){
			double sum = 0;
			for( int j = 0; j < smooth_window_n; ++j ){
				sum += it->second[i+j];
			}	
			_chr_RPBMs[it->first][i+smooth_window_n/2] = sum/smooth_window_n;
		}
	}
	cout<<" done!"<<endl;
	std::ofstream out(out_f.data());
	for( map<Str, vector<double> >::iterator it = _chr_RPBMs.begin(); it != _chr_RPBMs.end(); ++it ){
		for( unsigned int i = 0; i < it->second.size(); ++i ){
			if( it->second[i] > 0){
				out<<it->first<<'\t'<<i*windowsize<<'\t'<<((i+1)*windowsize-1)<<'\t'<<it->second[i]<<endl;
			}
		}
	}
	out.close();
}

int main(int argc, char* args[]) {
	if( argc != 8  ){
		cout<<"This program is generate RPBM-based summary-graph file\n"
	    "Usage: generateRPBMBasedSummary_v2 bed6_file chr_length_mapping_file window_size\n"
		" fragment_shift_size slide_window_number to_normalize(y/n) out_put_file\n";
		return 0;
	}

	vector<Str> fs;
	if( Str(args[1]).find(".bed" ) != std::string::npos ){
		fs.push_back(Str(args[1]));
	}else{
		ifstream in(args[1]);
		if( !in.good() )
			FileNotFoundERR(Str(args[1]));
		while( !in.eof() ){
			Str line;
			in>>line;
			if( !line.empty() )
				fs.push_back( line );
		}
		in.close();
	}

	generateRPBMSummaryGraph_slide(fs, args[2], args[7], atoi(args[3]), atoi(args[4]), atoi(args[5]), Str(args[6]) == "y" );

	return 0;
}
