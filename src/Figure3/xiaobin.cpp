/*
 * xiaobin.cpp
 *
 *  Created on: Oct 22, 2012
 *      Author: hugq
 */

#include "xiaobin.h"

double gammaln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
double betaln(double z, double w)
{
	return gammaln(z)+gammaln(w)-gammaln(z+w);
}

double gammaser(double a, double x)
{
	int ITMAX=1000;
	double EPS=1.0e-10;
	if (ITMAX<10*a) ITMAX=floor(10*a);
	int n;
	double sum,del,ap;
	double gln=gammaln(a);
	double gamser;
	if (x <= 0.0)
	{
		gamser=0.0;
		return gamser;
	}
	else
	{
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++)
		{
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del)*ap/(ap-x) < fabs(sum)*EPS)
			{
				gamser=sum*exp(-x+a*log(x)-gln);
				return gamser;
			}
		}
	}
	std::cout<<"a too large, ITMAX too small in routine gser";
	return gamser;
}

double gammacf(double a, double x)
{
	double FPMIN=1.0e-30;
	int ITMAX=10000;
	double EPS=1.0e-10;
	int i;
	double an,b,c,d,del,h;
	double gln=gammaln(a);
	b=x+1.0-a; //Set up for evaluating continued fraction by modified Lentz�s method (�5.2) with b0 = 0.
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) //Iterate to convergence.
	{
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) std::cout<<"a too large, ITMAX too small in gcf"<<std::endl;
	return exp(-x+a*log(x)-gln)*h;
}

double gammap(double a, double x)
{
	if (a <= 0.0)
	{
		std::cout<<("Invalid arguments in routine gammp");
		return -1;
	}
	if (x<=0.0) return 0;
	if (x < (a+1.0)) //Use the series representation.
	{
		return gammaser(a,x);
	}
	else	//Use the continued fraction representation
	{
		return 1.0-gammacf(a,x);
	}
}

double gammaq(double a, double x)
{
	if (a <= 0.0)
	{
		std::cout<<("Invalid arguments in routine gammp");
		return -1;
	}
	if (x<=0.0) return 1.0;
	if (x < (a+1.0)) //Use the series representation.
	{
		return 1.0-gammaser(a,x);
	}
	else	//Use the continued fraction representation
	{
		return gammacf(a,x);
	}
}

double erff(double x)
{
	return x < 0.0 ? -gammap(0.5,x*x) : gammap(0.5,x*x);
}

double erffc(double x)
{
	return x < 0.0 ? 1.0+gammap(0.5,x*x) : gammaq(0.5,x*x);
}

double chisquarecdf(double x,double v)
{
	return gammap(v/2,x/2);
}

double betacf(double a,double b,double x)
{
	int MAXIT=10000;
	double EPS=1.0e-10;
	double FPMIN=1.0e-30;
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;
	qab=a+b; //These q�s will be used in factors that occur
	qap=a+1.0; //in the coefficients (6.4.6).
	qam=a-1.0;
	c=1.0; //First step of Lentz�s method.
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++)
	{
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;// One step (the even one) of the recurrence.
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d; //Next step of the recurrence (the odd one).
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break; //Are we done?
	}
	if (m > MAXIT) std::cout<<("a or b too big, or MAXIT too small in betacf");
	return h;
}

double betai(double a, double b, double x)//Returns the incomplete beta function Ix(a, b).
{
	if (fabs(a)<=1.0e-10) a=1.0e-10;
	if (fabs(b)<=1.0e-10) b=1.0e-10;
	double bt;
	if (x < 0.0 || x > 1.0) std::cout<<("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else //Factors in front of the continued fraction.
	bt=exp(gammaln(a+b)-gammaln(a)-gammaln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0)) //Use continued fraction directly.
	return bt*betacf(a,b,x)/a;
	else //Use continued fraction after making the symreturn
	return 1.0-bt*betacf(b,a,1.0-x)/b;// metry transformation.
}

double invbetai(double a,double b,double p)	//return inverse of the incomplete beta function, using newton method
{
	if (fabs(a)<1.0e-10) return 0;
	if (fabs(b)<1.0e-10) return 1;
	if (fabs(a-1.0)<1.0e-10) return 1.0-pow(1.0-p,1.0/b);
	if (fabs(b-1.0)<1.0e-10) return pow(p,1.0/a);
	double x=(a-1.0)/(a+b-2.0);
	double f=betai(a,b,x);
	while (fabs(f-p)>1.0e-10)
	{
	//	cout<<(f-p)<<endl;
		double q=exp(gammaln(a+b)-gammaln(a)-gammaln(b)+(a-1.0)*log(x)+(b-1.0)*log(1.0-x));
		x=x-(f-p)/q;
		f=betai(a,b,x);
	}
	return x;
}

double pvalue(double n1,double n2,double nf,double fc)
{
	return betai(n1+1,n2,nf*fc/(1.0+nf*fc));
}
