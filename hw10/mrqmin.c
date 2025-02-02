
#define NRANSI
#include "nr.h"
#include "nrutil.h"

void mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **covar, float **alpha, float *chisq,
	void (*funcs)(float, float [], float *, float [], int, struct Point []), float *alamda, struct Point p[])
{
	//printf("lamda0: %12.6f\n", *alamda);
	void covsrt(float **covar, int ma, int ia[], int mfit);
	void gaussj(float **a, int n, float **b, int m);
	void mrqcof(float x[], float y[], float sig[], int ndata, float a[],
		int ia[], int ma, float **alpha, float beta[], float *chisq,
		void (*funcs)(float, float [], float *, float [], int, struct Point[]), struct Point p[]);
	int j,k,l;
	static int mfit;
	static float ochisq,*atry,*beta,*da,**oneda;
	if (*alamda < 0.0) {
		atry=vector(1,ma);
		beta=vector(1,ma);
		da=vector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]) mfit++;
		oneda=matrix(1,mfit,1,1);
		*alamda=0.001;
		mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs, p);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(*alamda));
		oneda[j][1]=beta[j];
		//printf("%12.9f check \n", covar[j][j]);
	}
	gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		covsrt(covar,ma,ia,mfit);
		covsrt(alpha,ma,ia,mfit);
		free_matrix(oneda,1,mfit,1,1);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		free_vector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
	mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs, p);
	printf("ochi : %12.9f chisq : %12.9f", ochisq, *chisq);
		
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		// printf("jere");
		// printf("ochi : %12.9f", ochisq);
		 for (j=1;j<=mfit;j++) {
			for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
			beta[j]=da[j];
			printf("\n%12.9f check \n", covar[j][j]);
		}
		for (l=1;l<=ma;l++) a[l]=atry[l];
	} else {
		//printf("lamda1: %12.6f\n", *alamda);
		*alamda *= 10.0;
		*chisq=ochisq;
		//printf("lamda2: %12.6f\n", *alamda);
	}
}
#undef NRANSI
