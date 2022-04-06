// eLCX_main.cpp : Main Programme for eLC package
// update on 6/9/2012:  generating null distribution for all SNPs at the beginning
//  By  Xing Chen S.cD.  email: dr.xingchen@gmail.com


//#include "stdafx.h"
#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
                                     // newmatap.h will get include.h
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
//#include "normbase.h"

//#include "xwxw_mvt.h"
//#include "myutil.h"

#include "eLC_model.h"
//#include "myerror.h"

#include "newmat11/newmat.h"
#include "newmat11/newmatap.h"
#include "newmat11/newmatio.h"
#include "newmat11/include.h"
#include "cdflib.hpp"

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif
using namespace std;
const int kmaxlinesize = 700000;
const int kmaxstrlen = 32;


//float** ReadTable(const char* FileName, int& RowNum, int& ColNum) {
void ReadTable(const char* FileName, int& RowNum, int& ColNum) {
   string line;
   ifstream in(FileName);
    cerr << "read table" << endl;
// Determine number of rows and columns in file.
// Program halts if rows are not of same length.
    while(getline(in,line,'\n')) {
    string segments;
    int ColsPerRow = 0; // Initialize counter.
    stringstream ss;
    ss << line;
    while (getline(ss,segments,'\t')) {
        ColsPerRow++;
    }
    if (RowNum == 0) {
        // Define number of columns in file as the number
        // of columns in the first row.
        ColNum = ColsPerRow;
    } else {
        if (ColsPerRow != ColNum) {
            cerr << "Row " << RowNum << " is not the same length "
            "as row 0." << endl;
            exit(0);
        }
    }
    RowNum++;
 }
    cerr << "row:" << RowNum << " col:" << ColNum << endl;
}

double pearsonCorrelation ( int nsnp, double *x, double *y) {
  int i,n;
  double xbar=0.0, ybar=0.0,sx=0.0,sy=0.0, pcorr=0.0;
  n=nsnp;
  /*compute xbar ybar */
  for(i = 0; i < n; ++i) {
      xbar += x[i];
      ybar += y[i];
   }
   xbar /= n;
   ybar /= n;

   /* compute standard deviation of x*/
   for(i = 0; i < n; ++i) {
      sx += (x[i] - xbar) * (x[i] - xbar);
      sy += (y[i] - ybar) * (y[i] - ybar);
   }
   sx = sqrt(sx /(n - 1));
   sy = sqrt(sy /(n - 1));

   /*compute r, the correlation coefficient between the two arrays */
   for( i = 0; i < n; ++i ) {
      pcorr += (((x[i] - xbar)/sx) * ((y[i] - ybar)/sy));
   }
   pcorr /= (n-1);
   return pcorr;

}


int IntPow(int x, int pow)
{
    int ret = 1;
    while ( pow != 0 )
    {
        if ( (pow & 1) == 1 )
            ret *= x;
        x *= x;
        pow >>= 1;
    }
    return ret;
}


/* usage: eLC -s -e #permutations Input_GWAS_Summary_file Output_eLC_file */
int main(int argc, char *argv[]) {

    int wkST=0,wkEnd=0;
    int nTraits = 0, rTraits;
    //double corr;
    int k, n;
    int nPerm = 0;
    char snpName[kmaxstrlen];
    ifstream InFile;
    ofstream OutFile;
    if (argc < 4) {
            cerr << "usage: " << argv[0]
                    << " -s -e -n #permutations -i Input_GWAS_Summary_file -o Output_eLC_file ; \n"
                    <<" optional : -s starting position ; -e: # of SNPs for analysis \n ";
            return -1;
    }

    int ch=0;
    char *infn,*otfn;
    //string *infn, *otfn;
    int RowNum = 0;
    int ColNum = 0;
    while((ch=getopt(argc,argv,"i:o:s:e:n:"))!=-1){
          switch(ch)
          {

            case 's':  wkST=atoi(optarg);break;
            case 'e':  wkEnd=atoi(optarg);break;
            case 'n':  nPerm=atoi(optarg);break;

            case 'i':  infn=optarg;
                        break;

            case 'o': otfn=optarg;
                      break;
            }
    }
    
    cerr << "input:" << infn << endl;
    cerr << "ooutput:" << otfn << endl;
    cerr << "nPerm:" << nPerm << endl;
      InFile.open(infn);
      //fprintf(stderr,"in: %s\n",infn);

      OutFile.open(otfn);
    /*
    if (InFile == NULL) {
            printf("Error: can't access Univariate_GWAS_Summary_InputFile\n");
              return 1;
              }
    if (OutFile == NULL) {
              printf("Error: can't create output file for writing\n");
              return 1;
              }
    */
    ReadTable(infn,RowNum,ColNum);


    //float beta[nTraits], se[nTraits], z[nTraits], pval[nTraits],a,w[nTraits];


    /*} else {
        //InFile = fopen(argv[2], "r");
        //OutFile = fopen(argv[3], "w+");

      //  InFile.open(argv[2]);
       // OutFile.open(argv[3]);


          //sscanf(argv[1],"%d",&nPerm);
    }*/



    int nperm = IntPow(10, nPerm);
    cerr << "nperm:" << nperm << endl;
    //skip the header
    //char mystring [10000];
    //fgets(mystring , 10000 , InFile);



    nTraits=ColNum-1;
    int nSNPs=RowNum-1;
    
    fprintf(stderr," Read in from %s : %d SNPS, %d traits;\n Now Caculating Correlation Matrix\n",infn,nSNPs,nTraits);
    //cout <<" There are totally "<< nSNPs<<" SNPs and "<< nTraits <<" Traits for analysis\n";


    // read all the z values into vector for calculate Pearson correlations


    rTraits = nTraits;
    //char *trait = new char[nTraits];
    //int *nperson = new int[nTraits];
    //double *beta = new double[nTraits];
    //double *se = new double[nTraits];

    //double *pval = new double[nTraits];

    SymmetricMatrix U(nTraits);

    eLC_model *xwmvtc;




    // Reposition to start of stream buffer.

    // InFile.clear();
    //InFile.seekg (0, ios::beg);


    //skip the header
    string dummyLine;
    getline(InFile, dummyLine);

    string Input;
    vector< vector< double > > iData;
    while ( getline( InFile, Input, '\n' ) )
       {
           // Input contains entire line
        stringstream Parse(Input);

        string Value;
        vector<double> Line;
        while ( getline(Parse,Value,'\t') ) {
            stringstream fs(Value);
            double f=0.0;
            fs >> f;
            Line.push_back( f );
        }
           if ( Line.size() != ColNum )
           {
             std::cout << "Line " << iData.size() + 1 << " expected " <<
             ColNum << " values, received " << Line.size() << " values, aborting." << endl;
             return 1;
           }

           iData.push_back(Line);
       }
       /* // Check whether Data read in >>>>> Successfully!
       cout << "Data:\n";

        for ( size_t i = 0; i < iData.size(); ++i )
        {
            for ( size_t j = 0; j < iData[i].size(); ++j )
                std::cout << iData[i][j] << " ";
            std::cout << "\n";
        }
        +/*/
      // At this point, all read values are in our vector of vector of floats named Data.
       // calculate Pearson Correlation for each pairs ,skip the header and first column


    for (int j = 1; j <= nTraits; j++) {
        for (k = 1; k <= nTraits; k++)
            U(j, k) = 0.0;
    }


    for ( int j = 1; j <= nTraits; j++ ) {
        for ( int k = 1; k <= nTraits; k++ ) {

              double *tmpx= new double [nSNPs];
              double *tmpy= new double [nSNPs];
              for (int i=0; i< nSNPs; i++) {
                   cout<<iData[i][j]<<"\t"<<iData[i][k]<<"\n";
                    tmpx[i]=(iData[i][j]);
                    tmpy[i]=(iData[i][k]);
                  //cout<<tmpx[i]<<"\t"<<tmpy[i]<<"\n";
              }

               double pcr = pearsonCorrelation(nSNPs,tmpx,tmpy);
               U(j,k)=pcr;
               if (j<k){
                   fprintf(stderr," the pair-wise correlation between %d th trait and %d th trait is %.5lf\n",j,k,U(j,k));
               }
               delete[] tmpx;
               delete[] tmpy;
               //cout << " The correlation between the pairwise Univariate GWAS results is "<< U(j,k) << "\n";
                  //         std::cout << "\n";
           }

       }

    InFile.close();
    OutFile.close();
    vector< vector< double > > anull;
    iData.swap(anull);


    /// initiate the null distribution for all
    xwmvtc = new eLC_model ( 26,6,1,nperm);

    SymmetricMatrix T(rTraits);

    for (int j = 1; j <= rTraits; j++) {
        for (k = 1; k <= rTraits; k++) {
                    T(j, k) = 0.0;
                }
            }
    for (int j = 1; j <= rTraits; j++) {
        for (k = 1; k <= rTraits; k++) {
                    T(j, k) = U(j, k);
                }
            }


    LowerTriangularMatrix L= Cholesky (T);


    xwmvtc->null_distribution(nperm ,L,U);



    //  Read in File Line by Line for OB,dLC,eLC

    // Reposition to start of stream buffer.

        // InFile.clear();
        //InFile.seekg (0, ios::beg);

        FILE *Infile;
        FILE *Outfile;

        int stRow=1,endRow=RowNum;

        if (wkST!=0||wkEnd!=0){
            stRow=wkST;
            endRow=wkEnd;
        }

        Infile = fopen(infn, "r");
        Outfile = fopen(otfn, "w+");


        // output header

        fprintf(Outfile, "SNP\t OB\t OB.Pval\t dLC\t dLC.Pval\t eLC.Pval\n");

        char mystring[10000];
        //skip the header and rows
        for ( int nr=0;nr<stRow;nr++)
        fgets (mystring , 10000 , Infile);

    //double *w= new double[nTraits];

    //double a,az,p;

    //calculate the correlation matrix U

    fprintf(stderr," Starting eLC analyzing from %dth SNP......\n",stRow);

    for( int nr=1;nr<endRow;nr++) {
        double *z = new double[nTraits];
        for (int j = 0; j < rTraits; j++) {
                    z[j] = 0.0;
                }

           fscanf(Infile, "%s", snpName);
        for (int i = 0; i < nTraits; i++) {
            fscanf(Infile,"%lf",&z[i] );
            }
            fscanf(Infile,"\n");
        double *zr = new double[rTraits];
        for (int j = 0; j < rTraits; j++) {
            zr[j] = 0.0;
        }

        for (int j = 0; j < rTraits; j++) {
            zr[j] = z[j];
        }

        //Using summary data,

        Matrix V(1, 1), V2(1, 1);
        V = 0.0;
        V2 = 0.0;
        ColumnVector Wd(rTraits) ,W(rTraits),E(rTraits), Wz(rTraits);


         // OB
         double OBz=0.0,OBp=0.0,OBtmp=0.0;
         for (n=0; n<rTraits; n++)
             E(n+1) = 1.0;

         Wd=(E.t()*T.i()*E);
         W = T.i()*E*Wd.i();
         V = W.t()*T*W;

         //printf("V=%.10lf\n",V.AsScalar());

         for (OBz=0.0,n=0; n<rTraits; n++) {
             OBz += W(n+1)*zr[n];
             //cout << W(n+1) << "\n";
         }

         OBz /= sqrt(V.AsScalar());
         OBtmp=abs(OBz);

         //cout<< OBz<<"\n";

         double mean=0;
         double pz=-1.0;
         double qz=-1.0;
         double sd=1.0;
         int status;
         int which=1;
         double bound;


         cdfnor(&which,&pz, &qz, &OBtmp, &mean, &sd,&status, &bound);
         OBp=2*qz;
         //cout << OBz << "\t" << OBp << "\n";



        // dLC

        for (n = 0; n < rTraits; n++)
            Wz(n + 1) = z[n];

        V2 = Wz.t() * T.i() * Wz;
        double Sdc=0.0;
        double Pdc=0.0;
        Sdc=V2.AsScalar();
        // calculate the p-value

          double df=nTraits;
          double pdc1=-1.0;
          double qdc=-1.0;

          cdfchi( &which, &pdc1, &qdc, &Sdc, &df, &status, &bound );
          Pdc=qdc;
          // cout << " The bivariate result for SNP\t"<<snpName<<"\tOB value is "<< OBz
          //    << "\tequal to p-value "<< OBp<< " and Chi_Beta is\t" << Sdc << "and Chisq_Z is " << Pdc<< "\n";




          /// eLC

        double eCHI_p;
         eCHI_p = xwmvtc->empirical_p(z, T);

        //cout << " The eLC for SNP\t " << snpName<< "\t equal to \t" << eCHI_p<< "\n";

        //fprintf(OutFile,"%s\t %lf\t %F\n",snpName,chi1,chi2);
        /*    printf("V=%.10lf\n",V.AsScalar());

         for (a=0.0,n=0; n<rTraits; n++) {
         a += W(n+1)*beta[n];
         //cout << W(n+1) << "\n";
         }

         //a /= sqrt(V.AsScalar());
         //cout<< a<<"\n";
         //fprintf(OutFile,"%s\t %lf\n",snpName,a);
         // printf("%s\t %lf\t \n",snpName,a);



         for (a=0.0,n=0; n<nTraits; n++) {
         for (w[n]=0.0,k=0; k<nTraits; k++)
         w[n] += U(n+1,k+1);
         printf("w[n]=%lf\n",w[n]);
         a += w[n];
         printf("a=%lf",a);
         }
         // now a = 1/var(comb)


         double s;
         for (s=0.0,n=0; n<nTraits; n++) {
         w[n] /= a;
         printf("w[n]=%lf\n",w[n]);
         s += w[n]*beta[n];
         printf("s=%lf\n",s);

         // cout << w[n] << "\t";
         }
         s *= sqrt(a);
         cout << snpName<<"\t"<<s << "\n";

         fprintf(OutFile,"%s %lf\n",snpName,s);
         */

        //           empirical method

        /* beta covariance matrix

       int ii, jj;
        SymmetricMatrix C(nTraits);
            for (ii=1; ii<=nTraits; ii++)
                for (jj=1; jj<=nTraits; jj++)
                    C(ii,jj) =U(ii,jj)/sqrt(U(ii,ii)*U(jj,jj));

        double emp_p;
        xwmvt = new xwmvt_model;
        emp_p = xwmvt->empirical_p(zr, C);

        // z covariance matrix ; 3 traits


        double emp_pz;
        SymmetricMatrix Uz(nTraits);
        Uz(1, 1) = 1.22;
        Uz(2, 2) = 1.17;
        Uz(3, 3) = 1.19;
        Uz(2, 1) = Uz(1, 2) = -0.446;
        Uz(3, 1) = Uz(1, 3) = -0.315;
        Uz(2, 3) = Uz(3, 2) = 0.332;
        // HDL LDL
        Uz(1,1)=1.19;
        Uz(2,2)=1.22;
        Uz(2,1)=Uz(1,2)=-0.315;

        //HDL TG
        //Uz(1,1)=1.22;
         Uz(2,2)=1.17;
         Uz(2,1)=Uz(1,2)=-0.447;
        // LDL TG
        // Uz(1,1)=1.19;
         Uz(2,2)=1.17;
         Uz(2,1)=Uz(1,2)=0.3314;
        xwmvtz = new xwmvt_model;
        emp_pz = xwmvtz->empirical_p(zr, Uz);
        */

        fprintf(Outfile, "%s\t %lf\t  %-G\t  %lf\t %-G\t %-G\t\n", snpName,OBz, OBp,Sdc,Pdc,eCHI_p);

        //free pointer;
        delete[] z;
    }
    delete xwmvtc;
    int endpos= stRow+endRow-1;
    fprintf(stderr," eLC analysis finished at the %dth SNP!\n",endpos);
    fclose(Infile);
    fclose(Outfile);
    return 0;
}


