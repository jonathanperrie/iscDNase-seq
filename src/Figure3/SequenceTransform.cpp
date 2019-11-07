#include"SequenceTransform.h"
#include"OftenUsedOperatLib.h"

Str SequenceTransform_T::char2DigitalSeq( Str& seq )
{
	int seqLen = seq.size();
	Str tmp;
	tmp.reserve( seqLen );
	int i = 0;
	for( ; i < seqLen; ++i )
		tmp += char2digital( seq[i] );
	return tmp;
}	

char SequenceTransform_T::char2digital( char res )
{
	switch( res)
	{
	case 'A' : return '0';
	case 'C' : return '1';
	case 'G' : return '2';
	case 'T' : return '3';
	case 'N' : return '2';
	case 'X' : return '0';
	case 'H' : return '3';
	case 'M' : return '1';
	case 'K' : return '2';
	case 'D' : return '0';
	case 'R' : return '2';
	case 'Y' : return '3';
	case 'S' : return '1';
	case 'W' : return '0';
	case 'B' : return '1';
	case 'V' : return '2';
	case 'a' : return '0';
	case 'c' : return '1';
	case 'g' : return '2';
	case 't' : return '3';
	case 'n' : return '2';
	case 'x' : return '0';
	case 'h' : return '3';
	case 'm' : return '1';
	case 'k' : return '2';
	case 'd' : return '0';
	case 'r' : return '2';
	case 'y' : return '3';
	case 's' : return '1';
	case 'w' : return '0';
	case 'b' : return '1';
	case 'v' : return '2';
	default  : 
//		assert("Unexpected Character!");
		return '$';
	}
}

void SequenceTransform_T::char2FileDigitalSeq( Str& in, Str& seq )
{
	seq.erase( seq.begin(), seq.end() );
	std::ifstream inFile( in.data() );
	if( !inFile.good() ){
		std::cout<<"file "<<in<<" not found!"<<std::endl;
		exit(1);
	}
	Str firstline;
	std::getline( inFile, firstline );
	if( firstline.find( "gb|" ) != Str::npos )
		firstline= firstline.substr( firstline.find( "gb|" ) + 3);
	if( firstline.find( "bj|" ) != Str::npos )
		firstline= firstline.substr( firstline.find( "bj|" ) + 3);
	if( firstline.find( "emb|" ) != Str::npos )
		firstline= firstline.substr( firstline.find( "emb|" ) + 4);

	seq.reserve( 10000000 );
		
	while( !inFile.eof() )
	{
		char tmpChar;
		inFile>>tmpChar;
		seq += char2digital( tmpChar );
	}
	seq.erase( seq.end() - 1 );
	inFile.close();
}


Str SequenceTransform_T::digital2CharSeq( Str& seq )
{
	Str tmp( seq.size(), '0' );
	int i = 0;
	for( ; i < seq.size(); ++i )
	{
		switch( seq[i] )
		{
		case '0' : tmp[i] = 'A'; break;
		case '1' : tmp[i] = 'C'; break;
		case '2' : tmp[i] = 'G'; break;
		case '3' : tmp[i] = 'T'; break;
		default :
//			assert("Unexpected character!")
			;
		}
	}
	return tmp;
}

Str getDNAInOppositeStrand(Str s){
	std::reverse( s.begin(), s.end() );
	std::for_each( s.begin(), s.end(),
		SequenceTransform_T::ToOppRule_T() );
	return s;
}


vector<pair<Str, Str> > readInFast(Str file){
	std::ifstream in(file.data());
	if( !in.good() ){
		FileNotFoundERR(file);
	}
	vector<pair<Str, Str> > rst;
	Str id, seq;
	Str line;
	while( !in.eof() ){
		std::getline(in, line);
		if( line.find(">") != std::string::npos ){
			//cout<<line<<endl;
			if( !seq.empty() ){
				for( unsigned int j = 0; j < seq.size(); ++j ){
				//	seq[j] = SequenceTransform_T::char2digital( seq[j] );
				}//*/
				rst.push_back(make_pair(id, seq));
			}
			seq.clear();
			id = line;
		}else{
			seq += line;
//			if(seq.size()%1000000==0)
//				cout<<seq.size()<<endl;
		}
	}
	if( !seq.empty() ){
		for( unsigned int j = 0; j < seq.size(); ++j ){
		//	seq[j] = SequenceTransform_T::char2digital( seq[j] );
		}//*/
		rst.push_back(make_pair(id, seq));
	}
	return rst;
}

Str toAminoSeq( const char* seq, int startPos, int endPosition){
	Str tmp;
	for( int i = startPos; i <= endPosition - 3; i+=3 ){
		switch( ((seq[i] - 48)<<4) + ((seq[i + 1] - 48)<<2) + seq[i + 2] - 48 ){
		case 36 : case 37 : case 38 : case 39 :
			tmp += 'A'; break;
		case 57 : case 59 :
			tmp += 'C'; break;
		case 33 : case 35 :
			tmp += 'D'; break;
		case 32 : case 34 :
			tmp += 'E'; break;
		case 61 : case 63 :
			tmp += 'F'; break;
		case 40 : case 41 : case 42 : case 43 :
			tmp += 'G'; break;
		case 17 : case 19 :
			tmp += 'H'; break;
		case 12 : case 13 : case 15 :
			tmp += 'I'; break;
		case 0 : case 2 :
			tmp += 'K'; break;
		case 28 : case 29 : case 30 : case 31 : case 60 : case 62 :
			tmp += 'L'; break;
		case 14 :
			tmp += 'M'; break;
		case 1 : case 3 :
			tmp += 'N'; break;
		case 20 : case 21 : case 22 : case 23 :
			tmp += 'P'; break;
		case 16 : case 18 :
			tmp += 'Q'; break;
		case 8 : case 10 : case 24 : case 25 : case 26 : case 27 :
			tmp += 'R'; break;
		case 9 : case 11 : case 52 : case 53 : case 54 : case 55 :
			tmp += 'S'; break;
		case 4 : case 5 : case 6 : case 7 :
			tmp += 'T'; break;
		case 44 : case 45 : case 46 : case 47 :
			tmp += 'V'; break;
		case 58 :
			tmp += 'W'; break;
		case 49 : case 51 :
			tmp += 'Y'; break;
		case 48 : case 50 : case 56: tmp += 'X'; break;
		default : tmp += '?';
		}
	}
	return tmp;
}

Str get_rev_comp(Str s ){
	Str s_rev;

	for( int j = s.size()-1; j >= 0; --j ){
		if( s[j] == 'A' || s[j] == 'a')
			s_rev += "T";
		else if( s[j] == 'C' || s[j] == 'c' )
			s_rev += "G";
		else if( s[j] == 'G' || s[j] == 'g' )
			s_rev += "C";
		else if( s[j] == 'T' || s[j] == 't' )
			s_rev += "A";
		else
			s_rev += "N";
	}
	return s_rev;

};

Str get_comp(Str s ){
	Str s_rev;

	for( int j = 0; j < s.size(); ++j ){
		if( s[j] == 'A' || s[j] == 'a')
			s_rev += "T";
		else if( s[j] == 'C' || s[j] == 'c' )
			s_rev += "G";
		else if( s[j] == 'G' || s[j] == 'g' )
			s_rev += "C";
		else if( s[j] == 'T' || s[j] == 't' )
			s_rev += "A";
		else
			s_rev += "N";
	}
	return s_rev;
}

Str get_rev_comp_dig(Str s ){
	Str s_rev;

	for( int j = s.size()-1; j >= 0; --j ){
		if( s[j] == '0')
			s_rev += "3";
		else if( s[j] == '1')
			s_rev += "2";
		else if( s[j] == '2')
			s_rev += "1";
		else if( s[j] == '3')
			s_rev += "0";
		else
			s_rev += "-";
	}
	return s_rev;

};
