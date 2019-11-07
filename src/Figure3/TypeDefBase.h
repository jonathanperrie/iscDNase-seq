#ifndef TYPEDEFBASE_H
#define TYPEDEFBASE_H

#include<algorithm>
#include<iomanip>
#include<iostream>
#include<fstream>
#include<map>
#include<set>
#include<string>
#include<sstream>
#include<utility>
#include<vector>

#include"assert.h"
#include"math.h"

using std::vector;
using std::make_pair;
using std::cout;
using std::endl;
using std::setw;
using std::set;
using std::map;
using std::pair;

typedef std::pair<int, int> Pa_I_I;
typedef std::string Str;
typedef std::string string;
typedef std::vector<double> Ve_D;

enum GenomePartIdentify{PROMOTER, GENEBODY, EXON, INTRON, INTERGENIC, TSSTES};
extern int PROMOTERUPLEN, PROMOTERDOWNLEN, MININTRESTREGIONLEN, SHIFT;

#endif //TYPEDEFBASE_H
