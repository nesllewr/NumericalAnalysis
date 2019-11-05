/* Driver for routine gasdev */

#include <stdio.h>
#include <time.h>
// #include <stdlib.h>

#define NRANSI
#include "nr.h"

#define N 100
#define NOVER2 (N/2)
#define NPTS 100
#define ISCAL 1000
#define LLEN 60

int main(void)
{
	char words[LLEN+1];
	int i,j,k,klim,dist[N+1];
	float dd,x;

	float cdf =0.0;
	// srand(time(NULL));
	// long idum = rand();
	long idum=(int)time(NULL);

	for (j=0;j<=N;j++) dist[j]=0;
	for (i=1;i<=NPTS;i++) {
		x = 5*ran1(&idum)-3;
		j=(int)(20*x+60);
		if ((j >= 0) && (j < N)) ++dist[j];
	}

	printf("Uniform distribution of %6d points\n",NPTS);
	printf ("%5s %10s %10s %9s\n","x","p(x)","cdf(x)","graph:");
	for (j=0;j<=N;j++) {

		dd=(float)dist[j]/NPTS;
		for (k=1;k<=LLEN;k++) words[k]=' ';
		klim=(int) (ISCAL*dd);
		if (klim > LLEN)  klim=LLEN;
		for (k=1;k<=klim;k++) words[k]='*';
		cdf = cdf+dd;
		printf("%8.4f %8.4f %8.4f ",-3.0+j*0.05, dd,cdf);
		for (k=1;k<=LLEN;k++) printf("%c",words[k]);
		printf("\n");
	}


	for (j=0;j<=N;j++) dist[j]=0;
	for (i=1;i<=NPTS;i++) {
		x=0.25*N*(1.5*(gasdev(&idum)));
		j=(int)(x > 0 ? x+0.5 : x -0.5	);
		if ((j >= -NOVER2) && (j <= NOVER2)) ++dist[j+NOVER2];
	}

	printf("Normally distributed deviate of %6d points\n",NPTS);
	printf ("%5s %10s %10s %9s\n","x","p(x)","cdf(x)","graph:");
	
	cdf=0.0;
	for (j=0;j<=N;j++) {
		dd=(float) dist[j]/NPTS;
		for (k=1;k<=LLEN;k++) words[k]=' ';
		klim=(int) (ISCAL*dd);
		if (klim > LLEN)  klim=LLEN;
		for (k=1;k<=klim;k++) words[k]='*';
		cdf = cdf + dd;
		printf("%8.4f %8.4f %8.4f ",-3.0+j*0.05, dd, cdf);
		for (k=1;k<=LLEN;k++) printf("%c",words[k]);
		printf("\n");
	}
	return 0;
}
#undef NRANSI