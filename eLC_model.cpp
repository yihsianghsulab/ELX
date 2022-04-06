// xwmvt_model.cpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

//#include "myerror.h"
//#include "rand.h"
//#include "rand-normal-distribution.h"
#include "eLC_model.h"
#include "newmat11/newmatap.h"
#include "newmat11/newmat.h"
#include "newmat11/newmatap.h"
#include "newmat11/newmatio.h"

using namespace std;



int compare_double(const void *n1, const void *n2) {
	
	return *(double *)n1 > *(double *)n2 ? 1 : -1;

}
eLC_model::eLC_model(int ncv, double maxc, double minc, int mC) {
	
	nc = ncv>1 && ncv<100? ncv: 21;
	
	c = new double[nc];
	double d = (maxc - minc)/(nc-1);
	
	int i;
	for (c[0]=minc,i=1; i<nc; i++)
		c[i] = c[i-1] + d;
	
	ncycle = 0;
	minCycle = 1000;
	maxCycle =mC;
    minObs = 50;
	
	s = new double*[nc+2];
	for (i=0; i<nc+2; i++)
		s[i] = new double[maxCycle];
	
}

eLC_model::~eLC_model() {
	
	if (c!=NULL)
		delete[] c;
		
	if (s!=NULL) {
		for (int i=0; i<nc; i++)
			if (s[i] != NULL)
				delete[] s[i];
		delete[] s;
	}
	
}

int eLC_model::qrank(double d, double *t) {
	int pc=0;
	for (int n=0;n<ncycle;n++) {
	    if (d>=t[n]) {
		  	 pc++;
		}
	     //else break;
	     	}
	return pc;
}

int eLC_model::rank(double d, double *t) {
	
	if (d<t[0])
		return 0;
	else if (d>=t[ncycle-1])
		return ncycle;
		
	int k, L1, L2;
    k = 0;
	L1 = 0;
	L2 = ncycle - 1;
	while (L1 < L2 - 1) {
		if (d == t[L2]) {
			k = L2;
			break;
		} 
		else if (d == t[L1]) {
			k = L1;
			break;
		}
		k = (L1 + L2)/2;
		if (d>t[k])
			L1 = k;
		else
			L2 = k;
	}
	while (d>=t[k] && k<ncycle)
		k++;
		
	return k;
}
	
/*
	Z_COR = V_RHO_SQRT*Z
	S = [max(C,Z_COR)]*V_RHO^(-1)*Z_COR	: changed 6/12/2002
*/

// R=V_RHO^(1/2); T=V_RHO^(-1);
static double rand_gaussian(void)
{
    double x = (double)random() / RAND_MAX,
           y = (double)random() / RAND_MAX,
           z = sqrt(-2 * log(x)) * cos(2 * M_PI * y);
    return z;
}
double** eLC_model::null_distribution(int cycles, LowerTriangularMatrix &L, SymmetricMatrix &V) {
	
	    //if (nz<1)
	nz =V.Nrows();
		//cout<<" nz="<<nz<<"\n";
	//if (nz!=L.Ncols() || nz!=R.Nrows())
	//	throw("inversed sqrt of covariance matrix not othorgonal");
	
	ncycle = cycles;
	cout<<" the # permutation is "<<cycles<<"\n";
	ColumnVector Z(nz);
	ColumnVector ZZ(nz);

	ColumnVector D(nz);
	ColumnVector F(nz);
	Matrix S(1,1);
	//double *z = new double[nz];
	//double *d = new double[nz];
	s = new double*[nc];

	
	int i, j, k;
	for (i=0; i<nc; i++)
		s[i] = new double[maxCycle];
	
	for (i=0; i<cycles; i++) {

		for (k=1; k<=nz; k++)
			D(k) = rand_gaussian();
		    Z = L*D;


		for (j=0; j<nc; j++) {
			for (k=1; k<=nz; k++){
				ZZ(k)=fabs(Z(k));
				F(k) = (ZZ(k))>c[j]? (ZZ(k)) : c[j];
			}
			// two side test S=MAX(|Zk|,c)*|Zk|
			S = F.t()*ZZ;
			s[j][i] = S.AsScalar();

			// New test statistics : Product of two weighted Z
			//double SS=F(1)*ZZ(1)*F(2)*ZZ(2);
			//s[j][i] = SS;

     	}

	   double ipg=((i+1.0)*100/(cycles+0.0));
       double prg[5]={20.0,40.0,60.0,80.0,100.0};
       for ( int i=0;i<5;i++){
	   if (fabs(ipg-prg[i])==0.00000)
		  	printf("\r  %2.0f %%Monte Carlo simulation of Null Distribution of eLC test statistics is finished!",ipg);
			//fflush(stdout);
			}
	}
	cout<<"\n";

	return s;

}





// V=covariance matrix of Z
double eLC_model::empirical_p(double z[], SymmetricMatrix &V) {
	
	nz = V.Nrows();
	if (nz!=V.Ncols())
		throw("inversed sqrt of covariance matrix not othorgonal");
	
	double *zs = new double[nc];

	int i=0, j=0, k=0;
	
	// cholesky decomposition for generate multi-normial distribution

	LowerTriangularMatrix L= Cholesky (V);

	//IdentityMatrix T(nz) ;
	
	ColumnVector Z(nz);
	ColumnVector ZZ(nz);
	ColumnVector F(nz);
	Matrix S(1,1);
	
	for (k=1; k<=nz; k++)
		Z(k) = (z[k-1]);
	// two side test S=MAX(|Zk|,c)*|Zk|
	for (j=0; j<nc; j++) {
		for (k=1; k<=nz; k++){
			ZZ(k)=fabs(Z(k));
  			F(k) = (ZZ(k))>c[j]? (ZZ(k)) : c[j];
		}
			S = F.t()*ZZ;
		    zs[j] = S.AsScalar();

		//double SS=F(1)*ZZ(1)*F(2)*ZZ(2);
		//zs[j]=SS;
	}	
	

	double p=0.0,pout=0.0;
	for (ncycle=minCycle; ncycle<=maxCycle; ncycle*=10) {
		cout<< " the cycle running right now is " <<ncycle<<"\n";

		double *sone= new double [ncycle];
		double *stwo= new double [ncycle];
		for (j=0; j<nc; j++) {
				for (i=0; i<ncycle; i++){
					// store the current line
					sone[i] = s[j][i];
				}
				    // sort Snull biggest to smallest on a specific c
				qsort(s[j], ncycle, sizeof(double), compare_double);
				   // rank  every Snull on a specific c
				for (i=0; i<ncycle; i++) {
					k = rank(sone[i], s[j]);
					// record the biggest rank in each replicates
					if (j==0 || k>sone[i])
						stwo[i] = k;
				}
			}
		    // sort all the biggest Snull with respect to each replciates
			qsort(stwo, ncycle, sizeof(double), compare_double);

	for (j=0; j<nc; j++){
			// rank the observed S with respect to each c
			k = rank(zs[j], s[j]);
			if (j==0 || k>i)
				i = k;
		}
	p = qrank(i, stwo);
	pout =1.0 - (p)/(ncycle+1.0);
	delete[] sone;
	delete[] stwo;
	double ptmp=50.0/(ncycle+0.0);
	//cout<<" the p value cut-off of this term is "<<ptmp<<"\n";
	//cout<<" the p value of this term is "<<pout<<"\n";
   if  (pout>=ptmp )
   break;


	}
	delete[] zs;

	return pout;

}


