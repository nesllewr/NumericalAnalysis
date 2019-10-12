
#define NRANSI
#include "nrutil.h"


void lubksb(float **a, int n, int *indx, float b[])
{
	int i,ii=0,ip,j;
	float sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


void mprove(float **a, float **alud, int n, int indx[], float b[], float x[])
{
	void lubksb(float **a, int n, int *indx, float b[]);
	int j,i;
	double sdp;
	float *r;

	r=vector(1,n);
	for (i=1;i<=n;i++) {
		sdp = -b[i];
		for (j=1;j<=n;j++) sdp += a[i][j]*x[j];
		r[i]=sdp;
	}
	lubksb(alud,n,indx,r);
	for (i=1;i<=n;i++) x[i] -= r[i];
	free_vector(r,1,n);
}
#undef NRANSI
