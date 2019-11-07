/*
 * expression.h
 *
 *  Created on: Oct 1, 2009
 *      Author: hugq
 */

#ifndef EXPRESSION_H_
#define EXPRESSION_H_

struct EDGER{
	Str id;
	double logConc, logFC, pvalue, FDR;
	Str status(double logfc, double adjP){
		if( logFC > logfc && FDR < adjP ){
			return "incr";
		}else if( logFC <-logfc && FDR < adjP ){
			return "decr";
		}else{
			return "nchg";
		}
	}
};

struct DegSeq{
	Str id;
	double ct1, ct2;
	double logFC, log2FCNorm;
	double zscore, pvalue, qvalue;
};

struct CuffDiff{
	Str id, name, locus, status, significant, test_stat;
	double FPKM1, FPKM2, logFC, pvalue;
};

struct CUFLNK{
	Str chr, source, feature, strand, frame;
	int start, end, score;
	map<Str, Str> attr_val;
	vector<CUFLNK> exons;
};

void output_cuflnk(CUFLNK cuflnk, std::ofstream& out);

map<Str, EDGER> getEdgeRDE(Str fileName, double FDR = 100, double FC = 1 );
map<Str, Str> getPostEdgeR(Str fileName );
vector<DegSeq> getDegSeq(Str fileName);
vector<CuffDiff> getCuffDiff(Str fileName);

map<Str,double> readExpression(Str fileName);
map<Str,double> readCufflinksGenesExpr(Str fileName);
map<Str,pair<double,Str> > readExpressionFromMicroArr(Str fileName);
//higher speed if chr_genes are sorted by txStart.
void RPKMCalculator(map<Str, vector<UCSC> > chr_genes, Str rnaf, Str out, int _5l = 5000000, int shift =0);

map<Str, set<Str> > getKg_alia(Str f);
map<Str, Str > getKg_entrz(Str f);
#endif /* EXPRESSION_H_ */
