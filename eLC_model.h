// eLC_model.h
#include "newmat11/newmat.h"
#include "newmat11/newmatap.h"
#include "newmat11/newmatio.h"

class eLC_model {
	int ncycle;
	int minCycle;

	int minObs;
	int nz;
	int nc;
	double *c;
	double **s;
	
	//int compare_double(const void *n1, const void *n2);
	int rank(double d, double *t);
	int qrank(double d, double *t);
	// 6/12/02
	// new definaition of S = Zc*Tao^(-1)*Z, where Zc=max(Z,C), Tao=correlation matrix
	

	public:
	int maxCycle;
	double** null_distribution(int cycles, LowerTriangularMatrix &L, SymmetricMatrix &V);
	eLC_model(int ncv = 51, double maxc = 8, double minc = 1,int mC=10000);
	//void init();
	double empirical_p(double z[], SymmetricMatrix &V);
	~eLC_model();
};

