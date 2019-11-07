/*
 * expression.cpp
 *
 *  Created on: Oct 1, 2009
 *      Author: hugq
 */

#include "TypeDefBase.h"
#include "OftenUsedOperatLib.h"
#include "ucsc.h"
#include "expression.h"

map<Str,double> readExpression(Str fileName){
	map<Str,double> name_exp;
	std::ifstream expIn(fileName.data());
	if( !expIn.good() )
		FileNotFoundERR(fileName);
	while( !expIn.eof() ){
		Str name, tmp;
		double exp;
		expIn>>tmp>>exp>>tmp>>tmp>>name;
		if( !name.empty() ){
			name_exp[name] = exp;
		}
		std::getline(expIn,name);
	}
	expIn.close();//*/
	cout<<"# of records from \""<<fileName<<"\": "<<name_exp.size()<<endl;
	return name_exp;
}

map<Str, EDGER> getEdgeRDE(Str fileName, double FDR, double FC){
	map<Str, EDGER> edgers;
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	Str line;
	getline( in, line );
	while( !in.eof() ){
		EDGER e;
		Str id, tmp;
		//in>>id>>e.logConc>>e.logFC>>e.pvalue>>e.FDR;
		in>>id>>e.logFC>>e.logConc>>tmp>>e.pvalue>>e.FDR;
		if( !id.empty() ){
			//e.id = e.id.substr(0, e.id.find("|") );
			for( int i = 0; i < id.size(); ++i )
				if(id[i] !='\"' )
					e.id += id[i];
			//if(e.FDR <=FDR && _abs(e.logFC) >= log(FC)/log(2) )
				edgers[e.id] = e;
		}
	}
	cout<<"# of records from \""<<fileName<<"\": "<<edgers.size()<<endl;
	return edgers;
}

vector<DegSeq> getDegSeq(Str fileName){
	vector<DegSeq> degseqs;
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	Str line;
	getline( in, line );
	while( !in.eof() ){
		DegSeq d;
		in>>d.id>>d.ct1>>d.ct2>>d.logFC>>d.log2FCNorm>>d.zscore>>d.pvalue>>d.qvalue;
		Str line;
		getline( in, line );
		if( !d.id.empty() ){
			d.id = d.id.substr(0, d.id.find("|") );
		//	cout<<d.id<<endl;
			degseqs.push_back(d);
		}
	}
	cout<<"# of records from \""<<fileName<<"\": "<<degseqs.size()<<endl;
	return degseqs;
}

vector<CuffDiff> getCuffDiff(Str fileName){
	vector<CuffDiff> degseqs;
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	Str line;
	getline( in, line );
	while( !in.eof() ){
		CuffDiff d;
		in>>d.id>>d.name>>d.locus>>d.status>>d.FPKM1>>d.FPKM2>>d.logFC
		>>d.test_stat>>d.pvalue>>d.significant;
		getline( in, line );
		if( !d.id.empty() ){
			//d.id = d.id.substr(1, d.id.size() );
//			cout<<d.id<<" "<<d.name<<" "<<d.locus<<" "<<d.status<<" "<<d.FPKM1<<" "<<d.FPKM2<<" "<<d.logFC
	//				<<" "<<d.test_stat<<" "<<d.pvalue<<" "<<d.significant<<endl;
			degseqs.push_back(d);
		}
	}
	cout<<"# of records from \""<<fileName<<"\": "<<degseqs.size()<<endl;
	return degseqs;
}

map<Str,double> readCufflinksGenesExpr(Str fileName){
	map<Str,double> name_exp;
	std::ifstream expIn(fileName.data());
	if( !expIn.good() )
		FileNotFoundERR(fileName);
	Str line;
	getline( expIn, line);
	while( !expIn.eof() ){
		Str id;
		double expr;
		expIn>>id>>line>>line>>line>>line>>expr;
		if( !id.empty() ){
			if( name_exp.find(id) == name_exp.end())
				name_exp[id] = expr;
			else{
//			cout<<"redundant id: "<<id<<endl;
			}
		}
		getline( expIn, line);
	}
	expIn.close();//*/
	cout<<"# of records from \""<<fileName<<"\": "<<name_exp.size()<<endl;
	return name_exp;
}

map<Str,pair<double,Str> > readExpressionFromMicroArr(Str fileName){
	map<Str,pair<double,Str> > name_exp;
	std::ifstream expIn(fileName.data());
	if( !expIn.good() )
		FileNotFoundERR(fileName);
	while( !expIn.eof() ){
		Str name, ptn;
		double exp;
		expIn>>name>>exp>>ptn;
		if( !name.empty() ){
			name_exp[name] = make_pair(exp,ptn);
		}
		std::getline(expIn,name);
	}
	expIn.close();//*/
	cout<<"# of records from \""<<fileName<<"\": "<<name_exp.size()<<endl;
	return name_exp;
}

void RPKMCalculator(map<Str, vector<UCSC> > chr_genes, Str rnaf, Str outf,
		int _5l, int shift  ){
	map<Str, double >name_exon_l;
	map<Str, double> name_N;
	map<Str, int> name_instance;
	//select enxon in the desired region
	for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
		for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
			Str name = gIt->name+"|"+gIt->proteinID+"|"+gIt->alignID+"|"+gIt->strand;
			name_instance[gIt->name] += 1;

			name_N[name] = 0;
			int L = 0;
			for( set<Pa_I_I>::iterator eIt = gIt->exons.begin(); eIt != gIt->exons.end(); ++eIt ){
				L += eIt->second - eIt->first;
			}
			name_exon_l[gIt->name] = L;
			if( L > _5l ){
				int l = 0;
				set<Pa_I_I> exons;
				bool modify = false;
				if( gIt->strand == "+"){
					for( set<Pa_I_I>::iterator eIt = gIt->exons.begin(); eIt != gIt->exons.end(); ++eIt ){
						int preL = l;
						l += eIt->second - eIt->first;
						Pa_I_I exon = *eIt;
						if( !modify ){
							if( l > L - _5l){
								exon.first += L-_5l-preL;
								modify = true;
							}
						}
						if( l > L - _5l){
							exons.insert(exon);
						}
					}
				}else{
					for( set<Pa_I_I>::iterator eIt = gIt->exons.begin(); eIt != gIt->exons.end(); ++eIt ){
						int preL = l;
						l += eIt->second - eIt->first;
						Pa_I_I exon = *eIt;
						if( l <= _5l){
							exons.insert(exon);
						}else{
							exon.second = exon.first+_5l-preL;
							modify = true;
							exons.insert(exon);
							gIt->exons = exons;
							break;
						}
					}
				}
				assert( modify );
				gIt->exons = exons;
			}
			L = 0;
			for( set<Pa_I_I>::iterator eIt = gIt->exons.begin(); eIt != gIt->exons.end(); ++eIt ){
				L += eIt->second - eIt->first;
			}
			name_exon_l[name] = L;//*/
		}
	}

	std::ifstream in(rnaf.data());
	if( !in.good() )
		FileNotFoundERR(rnaf);
	map<Str, vector<long long> > chr_reads;
	Str chrom, strand;
	int beg, end;
	long long n= 0, rn = 0;
	while( true ){
		in>>chrom>>beg>>end>>strand>>strand>>strand;
		if( !chrom.empty()){
			++rn;
			int pos = strand =="+" ? beg + shift : end - shift;
			chr_reads[chrom].push_back(pos);
			if( (++n) % 5000000 == 0 || in.eof()){
				cout<<"# reads in: "<<n<<endl;
				for( map<Str, vector<long long> >::iterator chrIt = chr_reads.begin();
						chrIt != chr_reads.end(); ++chrIt )
					std::sort(chr_reads[chrIt->first].begin(), chr_reads[chrIt->first].end());

				for( map<Str, vector<UCSC> >::iterator chrIt = chr_genes.begin(); chrIt != chr_genes.end(); ++chrIt ){
					if( chr_reads.find(chrIt->first) != chr_reads.end() ){
						vector<long long>::const_iterator curRd = chr_reads[chrIt->first].begin();
						for( vector<UCSC>::iterator gIt = chrIt->second.begin(); gIt != chrIt->second.end(); ++gIt ){
							for( set<Pa_I_I>::iterator eIt = gIt->exons.begin(); eIt != gIt->exons.end(); ++eIt ){
								//cout<<gIt->chrom<<"\t"<<eIt->first<<"\t"<<eIt->second<<endl;
								if( *curRd >= eIt->first ){
									while( *curRd >= eIt->first && curRd != chr_reads[chrIt->first].begin())
										--curRd;
								}
								while( *curRd < eIt->second && curRd != chr_reads[chrIt->first].end()){
									if( *curRd > eIt->first ){
										++name_N[gIt->name+"|"+gIt->proteinID+"|"+gIt->alignID+"|"+gIt->strand];
									}
									++curRd;
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
	in.close();

	std::ofstream out(outf.data());
	std::ofstream out1((outf+".count").data());
	out1<<"ID\tTag_N\tlen\n";
	for( map<Str, double>::iterator it = name_N.begin(); it != name_N.end(); ++it ){
		Str name = it->first.find("|") != std::string::npos ?
				it->first.substr(0, it->first.find("|")): it->first;
		out<<name<<'\t'<<1000000*it->second*1000/name_exon_l[it->first]/rn
			<<'\t'<<name_exon_l[it->first]<<'\t'<<it->second<<'\t'
			<<it->first<<'\t'<<name_instance[name]<<endl;
		out1<<it->first<<'\t'<<it->second<<'\t'<<name_exon_l[it->first]<<endl;
	}
	out.close();
	out1.close();
}


map<Str, set<Str> > getKg_alia(Str f){
	map<Str, set<Str> > kg_alia;
	std::ifstream in(f.data());
	if( !in.good() )
		FileNotFoundERR(f);
	while( !in.eof() ){
		Str kg, alia;
		in>>kg;
		getline(in, alia);
		strip(alia,'\t');
		if( !kg.empty() ){
			kg_alia[kg].insert(alia);
		}
	}
	in.close();
	return kg_alia;
}


map<Str, Str > getKg_entrz(Str f){
	map<Str, Str > kg_entrz;
	std::ifstream in(f.data());
	if( !in.good() )
		FileNotFoundERR(f);
	Str t;
	getline(in, t);
	while( !in.eof() ){
		Str kg, entrz;
		in>>kg;
		getline(in, entrz);
		strip(entrz,'\t');
		if( !kg.empty() ){
			if( kg_entrz.find(kg) == kg_entrz.end() ){
				kg_entrz[kg] = entrz;
			}else{
				cout<<"Multiple entrz id for "<<kg<<endl;
			}
		}
	}
	return kg_entrz;
}

void output_cuflnk(CUFLNK cuflnk, std::ofstream& out){
	out<<cuflnk.chr<<"\t"<<cuflnk.source<<"\t"<<cuflnk.feature
			<<"\t"<<cuflnk.start<<"\t"<<cuflnk.end
			<<"\t"<<cuflnk.score<<"\t"<<cuflnk.strand
			<<"\t"<<cuflnk.frame<<"\t";
	for(map<Str,Str>::iterator it = cuflnk.attr_val.begin();
			it != cuflnk.attr_val.end(); ++it ){
		out<<" "<<it->first<<" \""<<it->second<<"\";";
	}
	out<<endl;
}

map<Str, Str> getPostEdgeR(Str fileName ){
	map<Str, Str > id_status;
	std::ifstream in(fileName.data());
	if( !in.good() ){
		cout<<"File not exist: "<<fileName<<endl;
		exit(1);
	}
	Str line;
	getline( in, line );
	while( !in.eof() ){
		Str id, status;
		in>>id>>status>>status>>status>>status>>status;
		if( !id.empty() ){
			id_status[id] = status;
		}
		getline( in, line );
	}
	return id_status;
}
