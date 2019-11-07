/*
 * Probe.h
 *
 *  Created on: Sep 16, 2009
 *      Author: hugq
 */

#ifndef PROBE_H_
#define PROBE_H_

#include "TypeDefBase.h"

class Probe{
public:
	bool is2UniqueGene;
	Str probeID;
	Str representDBID;
	vector<Str> refIDs;
	friend std::istream &operator>>(std::istream &stream, Probe& p);
	Str apStatus;
	Str rgStatus;
	double foldChange, pvalue;
	vector<double> intensities;
};

std::istream &operator>>(std::istream &stream, Probe& p){
	Str line;
	std::getline(stream, line);
	if( !line.empty()){
		vector<Str> tokens = split(line, '\t');
		p.probeID = tokens[0];
		p.representDBID = tokens[5];
		p.refIDs = split(tokens[6], " /// ");
		p.is2UniqueGene = (p.refIDs.size() == 1);
	}
	return stream;
}

#endif /* PROBE_H_ */
