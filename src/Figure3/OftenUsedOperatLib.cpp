#include"OftenUsedOperatLib.h"

double GetAver( vector<double>& v )
{
	double aver = 0;
	int i = 0;
	for( ; i < v.size(); ++i )
		aver += v[i];
	return (aver / v.size());
}

double GetDerivate( vector<double>& tmp )
{
	vector<double> v = tmp;
	double aver = GetAver( v );
	double der = 0;
	int i = 0;
	for( ; i < v.size(); ++i )
		der += ( v[i] - aver ) * ( v[i] - aver );
	return sqrt( der / v.size() );
}

vector<Str> split(Str  s, char delim){
	std::vector <Str> ve; //we'll put all of the tokens in here
  	std::string  temp;
	while (s.find(delim, 0) != std::string::npos){ //does the string a comma in it?
		size_t  pos = s.find(delim, 0); //store the position of the delimiter
		temp = s.substr(0, pos);      //get the token
		s.erase(0, pos + 1);          //erase it from the source
		ve.push_back(temp);                //and put it into the array
	}
	ve.push_back(s);           //the last token is all alone
	return ve;
}

vector<Str> split(Str  s, Str delim){
	std::vector <Str> ve; //we'll put all of the tokens in here
  	std::string  temp;
	while (s.find(delim, 0) != std::string::npos){ //does the string a comma in it?
		size_t  pos = s.find(delim, 0); //store the position of the delimiter
		temp = s.substr(0, pos);      //get the token
		s.erase(0, pos + 1 + delim.size());          //erase it from the source
		//if( pos != 0 )
			ve.push_back(temp);                //and put it into the array
	}
	ve.push_back(s);           //the last token is all alone
	return ve;
}

void strip(Str& s, char c){
	Str tmp;
	for( unsigned int i = 0; i < s.size(); ++i){
		if(s[i] != c)
			tmp.push_back(s[i]);
	}
	s = tmp;
}

Str genomePartIdentify2Str(GenomePartIdentify iden){
	switch(iden){
	case PROMOTER:
		return "promoter";
	case GENEBODY:
		return "genebody";
	case INTERGENIC:
		return "intergenic";
	case INTRON:
		return "intron";
	case EXON:
		return "exon";
	default:
		return "ERR";
	}
}


map<Str, int> getChrlen(Str fileName){
	std::ifstream in(fileName.data());
	if( !in.good() )
		FileNotFoundERR(fileName);
	map<Str, int> chr_len;
	while( !in.eof() ){
		Str chr;
		int len;
		in>>chr>>len;
		if( !chr.empty() ){
			chr_len[chr] = len;
		}
	}
	return chr_len;
}


map<Str, Str> getIsoform_ID(Str file){
	map<Str, Str> iso_id;
	std::ifstream in(file.data());
	while( !in.eof() ){
		Str id;
		Str kg;
		in>>id>>kg;
		if( !kg.empty() ){
		//	if( iso_id.find(kg) == iso_id.end() )
				iso_id[kg] = id;
		//	else{
		//		cout<<kg<<endl;
		//	}
		}
	}
	in.close();
	cout<<"# of isoforma that has trascript cluster: "<<iso_id.size()<<endl;
	return iso_id;
}

map<Str, set<Str> > get_kg_alians(Str f ){
	map<Str, set<Str> > kg_alia;
		std::ifstream in4(f.data());
		while( !in4.eof() ){
			Str kg, alia;
			in4>>kg;
			getline(in4, alia);
			strip(alia, '\t');
			strip(alia, '\n');
			if( !kg.empty() ){
				kg_alia[kg].insert(alia);
			}
		}
		cout<<"kg alia mapping: "<<getSize(kg_alia)<<endl;
		return kg_alia;
}
map<Str, Str> get_a_b(Str f){
	map<Str, Str > a_b;
		std::ifstream in4(f.data());
		while( !in4.eof() ){
			Str a, b, tmp;
			in4>>a>>b;
			getline(in4, tmp);
			//strip(b, '\t');
			//strip(b, '\n');
			if( !a.empty() ){
				//if( a_b.find( a ) != a_b.end() )
				//	cout<<"Duplicated key "<<a<<" !"<<endl;
				a_b[a]= b;
			}
		}
		cout<<"a b mapping: "<<a_b.size()<<endl;
		return a_b;
}

map<Str, Str> get_a_after(Str f){
	map<Str, Str > a_b;
		std::ifstream in4(f.data());
		while( !in4.eof() ){
			Str a, b;
			in4>>a;
			getline(in4, b);
			//strip(b, '\t');
			//strip(b, '\n');
			if( !a.empty() ){
				//if( a_b.find( a ) != a_b.end() )
				//	cout<<"Duplicated key "<<a<<" !"<<endl;
				a_b[a]= b;
			}
		}
		cout<<"a b mapping: "<<a_b.size()<<endl;
		return a_b;
}

map<Str, Str> get_b_a(Str f){
	map<Str, Str > b_a;
		std::ifstream in4(f.data());
		while( !in4.eof() ){
			Str a, b;
			in4>>a;
			getline(in4, b);
			strip(b, '\t');
			strip(b, '\n');
			if( !a.empty() ){
				b_a[b]= a;
			}
		}
		cout<<"a b mapping: "<<b_a.size()<<endl;
		return b_a;
}

map<Str, Str> get_post_edgeR(Str f){
	map<Str, Str > a_b;
	std::ifstream in4(f.data());
	Str line;
	getline(in4, line);
	while( !in4.eof() ){
			Str a, b;
			in4>>a>>b>>b>>b>>b>>b;
			if( !a.empty() ){
				a_b[a]= b;
			}
		}
		cout<<"a b mapping: "<<a_b.size()<<endl;
		return a_b;
}

map<Str, Str > get_kg_refseq(Str f ){
	map<Str, Str > kg_ref;
	std::ifstream in4(f.data());
	while( !in4.eof() ){
		Str kg, refseq;
			in4>>kg;
			getline(in4, refseq);
			strip(refseq, '\t');
			strip(refseq, '\n');
			if( !kg.empty() ){
				kg_ref[kg]= refseq;
			}
		}
		cout<<"kg refseq mapping: "<<getSize(kg_ref)<<endl;
		return kg_ref;
}
//contingency table example
/*                         BRG1 biding
 * -----------------------------------------------------
 *						| Presence | Absence |
 * up-regulated gene    |     a    |     b   |  NA =  a+b
 * down-regulated gene  |     c    |     d   |  NB =  c+d
 * -----------------------------------------------------
 * 						|NS = a+c  |NF = b+d | N = a+b+c+d
 * For one from a, b, c, d is less than 7, use Yate's correction
 * kai(1 degree) = N(max(0,|ad-bc|-N/2))^2/NS/NF/NA/NB
 * ortehwise use pearson chi-test
 * From http://en.wikipedia.org/wiki/Yates'_correction_for_continuity
 */
double averageDiffByChiTest(double av_1, int n_1, double av_2, int n_2, bool yatecorection){
	double a = av_1*n_1;
	double b = n_1 - a;
	double c = av_2*n_2;
	double d = n_2 - c;
	double NA = a+b, NB = c+d, NS = a+c, NF = b+d, N = a+b+c+d;
	if( yatecorection ){
		return N*pow((_abs(a*d-b*c)-N/2),2)/NS/NF/NA/NB;
	}else{//pearson
		double a_e = NA*NS/N;
		double b_e = NA*NF/N;
		double c_e = NB*NS/N;
		double d_e = NB*NF/N;
		return pow(a-a_e,2)/a_e + pow(b-b_e,2)/b_e + pow(c-c_e,2)/c_e + pow(d-d_e,2)/d_e;
	}
}

void PAUSE(){
	system("read -p \"Press any key to start backup…\"");
}

double pearson(Ve_D v1, Ve_D v2){
	double v1_a = GetAver(v1);
	double v2_a = GetAver(v2);
	double a = 0, b = 0, c = 0;
	for( int i = 0; i < v1.size(); ++i ){
		a += (v1[i]-v1_a)*(v2[i]-v2_a);
		b += (v1[i]-v1_a)*(v1[i]-v1_a);
		c += (v2[i]-v2_a)*(v2[i]-v2_a);
	}
	if( b == 0 || c == 0)
		return 2;
	return a/sqrt(b)/sqrt(c);
}

map<double, double> get_ecdf(Ve_D v, int bin){
	map<double, double> rst;
	std::sort(v.begin(), v.end());
	int step = v.size()/bin;
	if( step == 0 ){
		cout<<"ecdf errors, samples size < bin size\n";
		exit(1);
	}
	for( int i = 0; i < v.size() - step; i += step ){
		double sum = 0;
		for( int j = 0; j < step; ++j ){
			sum += v[i+j];
		}
		rst[sum/step] = double(i + step)/v.size();
	}
	return rst;
}


map<Str, map<Str,Str> > load_David_GO(Str f){
	std::ifstream in(f.data());
	Str line;
	std::getline(in, line);
	vector<Str> keys = split(line, '\t');
	map<Str, map<Str,Str> > id_key_value;
	while( !in.eof() ){
		getline(in, line);
		if( !line.empty() ){
			vector<Str> values = split(line, '\t');
			if( values.size() != keys.size() ){
				cout<<"DAVID GO format ERR"<<endl;
				exit(1);
			}else{
				map<Str, Str> key_value;
				for( int i = 0; i < keys.size(); ++i ){
					key_value[keys[i]] = values[i];
				}
				if( key_value.find("Category") == key_value.end() ||
					key_value.find("Term") == key_value.end() ){
					cout<<"GO file head line ERR"<<endl;
					exit(1);
				}else{
					id_key_value[key_value["Category"]+"_"+key_value["Term"]] = key_value;
				}
			}
		}
	}
	//cout<<"Load GO #: "<<id_key_value.size()<<endl;
	return id_key_value;
}


int wc(Str f){
	int n = 0;
	Str line;
	std::ifstream in(f.data());
	while( !in.eof()){
		getline(in, line);
		if( !line.empty() )
		{
			++n;
		}
	}
	in.close();
	return n;
}

Str get_gs(Str g){
	Str g1 = g.substr(g.find("|")+1);
	Str g2 = g1.substr(g1.find("|")+1);
	Str g3 = g2.substr(0,g2.find("|"));
	return g3;
}
