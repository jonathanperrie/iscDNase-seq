//A set of functions commonly used
#ifndef OFTENUSEDOPERATLIB_H
#define OFTENUSEDOPERATLIB_H

#include"TypeDefBase.h"
#include <stdlib.h>
#include <stdio.h>
#include <cctype>

//split a string by a delimiter. Source code modified from
//http://www.cprogramming.com/faq/cgi-bin/smartfaq.cgi?answer=1057105750&id=1044780608
vector<Str> split(Str  s, char delim);
vector<Str> split(Str  s, Str delim);

void strip(Str& s, char c);

inline Str uppercase(Str s){
	transform(s.begin(), s.end(), s.begin(), toupper);
	return s;
}
//convert a string to a type specified by user. Source code modified from
//http://www.velocityreviews.com/forums/t281186-newbie-string2int.html
template <class T>
T convertString2(const Str& s){
	std::istringstream buf(s);
	T a;
	buf >> a;
	// error checking goes here
	return a;
}
template <class T>
std::string convert2String(const T& i)
{
	std::ostringstream os;
	os << i;
	return(os.str());
}

template <class T>
T _max(T a, T b){
	return a > b ? a:b;
}

template <class T>
T _min(T a, T b){
	return a > b ? b:a;
}

template <class T>
T _abs(T a){
	return a > 0 ? a:0-a;
}

inline void FileNotFoundERR(Str f){
	cout<<"File "<<f<<" not found"<<endl;
	exit(0);
}
//double
double averageDiffByChiTest(double a_av, int a_n, double b_av, int b_n, bool yatecorection = false);

double GetAver( vector<double>& v );
double GetDerivate( vector<double>& tmp );

inline void _pause(){
	std::cout<<"press enter to continue...\n";
	Str line;
	std::getline(std::cin, line);
}

//assume T defines a method named size
template <class T>
int getSize(const map<Str, T >& chr_T ){
	int n = 0;
	for( typename map<Str, T >::const_iterator it = chr_T.begin();
			it != chr_T.end(); ++it)
		n += it->second.size();
	return n;
}

Str genomePartIdentify2Str(GenomePartIdentify iden);

template<typename T>
map<Str, set<T> > vector2set( const map<Str, vector<T> >& a){
	map<Str, set<T> > b;
	for(typename map<Str, vector<T> >::const_iterator chrIt = a.begin(); chrIt != a.end(); ++chrIt ){
		for( typename vector<T>::const_iterator it = chrIt->second.begin(); it != chrIt->second.end(); ++it ){
			b[chrIt->first].insert(*it);
		}
	}
	return b;
}

template<typename T, typename S>
map<Str, set<T, S> > vector2set( const map<Str, vector<T> >& a){
	map<Str, set<T, S> > b;
	for(typename map<Str, vector<T> >::const_iterator chrIt = a.begin(); chrIt != a.end(); ++chrIt ){
		for( typename vector<T>::const_iterator it = chrIt->second.begin(); it != chrIt->second.end(); ++it ){
			b[chrIt->first].insert(*it);
		}
	}
	return b;
}

template<typename T, typename S>
map<Str, vector<T> > set2vector( const map<Str, set<T,S> >& a){
	map<Str, vector<T> > b;
	for(typename map<Str, set<T,S> >::const_iterator chrIt = a.begin(); chrIt != a.end(); ++chrIt ){
		for( typename set<T,S>::const_iterator it = chrIt->second.begin(); it != chrIt->second.end(); ++it ){
			b[chrIt->first].push_back(*it);
		}
	}
	return b;
}

template<typename T>
map<Str, vector<T> > set2vector( const map<Str, set<T> >& a){
	map<Str, vector<T> > b;
	for(typename map<Str, set<T> >::const_iterator chrIt = a.begin(); chrIt != a.end(); ++chrIt ){
		for( typename set<T>::const_iterator it = chrIt->second.begin(); it != chrIt->second.end(); ++it ){
			b[chrIt->first].push_back(*it);
		}
		typedef std::string Str;
	}
	return b;
}

void PAUSE();

map<Str, int> getChrlen(Str file);
map<Str, Str> getIsoform_ID(Str file); //kg 2 isoform_is
map<Str, set<Str> > get_kg_alians(Str f );
map<Str, Str > get_kg_refseq(Str f );
map<Str, Str> get_a_b(Str f);
map<Str, Str> get_a_after(Str f);
map<Str, Str> get_b_a(Str f);
map<Str, Str> get_post_edgeR(Str f);
double pearson(Ve_D v1, Ve_D v2);

map<double, double> get_ecdf(Ve_D v, int bin = 20);
map<Str, map<Str,Str> > load_David_GO(Str f);
int wc(Str f);
Str get_gs(Str g);

#endif
