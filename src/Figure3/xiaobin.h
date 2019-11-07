#ifndef XIAOBIN_H_
#define XIAOBIN_H_

#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>

double gammaln(double xx);
double betaln(double z, double w);
double gammaser(double a, double x);
double gammacf(double a, double x);
double gammap(double a, double x);
double gammaq(double a, double x);
double erff(double x);
double erffc(double x);
double chisquarecdf(double x,double v);
double betacf(double a,double b,double x);
double betai(double a, double b, double x);
double invbetai(double a,double b,double p);
/*
total gene: 9569
set1: 752
set2: 197
overlap: 70

=1-BINOMDIST(70,752,197/9569,TRUE)
xiaobin's reply
The functions are from the numerical recipes.
The pvalue() function is to do the binomial test.

For example, in your case, the function should be called as
pvalue(70,197-70,752.0/(9569-752))
//*/
//n1: overlap
//n2: set2-overlap
//nf: set1/(all-set1)
double pvalue(double n1,double n2,double nf,double fc=1.0);

#endif
