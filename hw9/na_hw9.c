/* Driver for routine gaussj */

#include <stdio.h>
#include <stdlib.h>
#define NRANSI

#include "nr.h"
#include "nrutil.h"

#define NP 20
#define MP 20
#define MAXSTR 80
	


int main(void)
{
	float x, y, xp, yp;
	float xsum, ysum, x2sum, y2sum, xysum, nsum;
	float xpsum, xpxsum, xpysum;
	float ypsum, ypxsum, ypysum;
	float **a, **ai, **b, **r, **sol;
	FILE *fp;
	char* filename = "fitdata3.dat";
	a = matrix(1, NP, 1, NP);
	ai = matrix(1, NP, 1, NP);
	b = matrix(1, NP, 1, MP);
	r = matrix(1, NP, 1, MP);
	sol = matrix(1, NP, 1, MP);
	//initialize
	nsum = xsum = ysum = x2sum = y2sum = xysum = 0.0;
	xpsum = xpxsum = xpysum = ypsum = ypxsum = ypysum = 0.0;

	fp = fopen(filename, "r");
	
	if (fp == NULL)
	{
		printf("Data '%s' file is not found\n", filename);
		return 0;
	}
	printf("file name : %s\n", filename);
	while (!feof(fp)) {
		fscanf(fp, "%f %f %f %f", &x, &y, &xp, &yp);
		nsum += 1.0;
		xsum += x;
		ysum += y;
		x2sum += x*x;
		y2sum += y*y;
		xysum += x*y;
		xpsum += xp;
		xpxsum += xp*x;
		xpysum += xp*y;
		ypsum += yp;
		ypxsum += yp*x;
		ypysum += yp*y;
	}

	a[1][1] = x2sum;
	a[1][2] = xysum;
	a[1][3] = xsum;
	a[2][1] = xysum;
	a[2][2] = y2sum;
	a[2][3] = ysum;
	a[3][1] = xsum;
	a[3][2] = ysum;
	a[3][3] = nsum;
	b[1][1] = xpxsum;
	b[2][1] = xpysum;
	b[3][1] = xpsum;

	for (int l = 1; l <= 3; l++) {
		for (int k = 1; k <= 3; k++) ai[k][l] = a[k][l];
		for (int k = 1; k <= 1; k++) r[l][k] = b[l][k];
	}
	printf("\nOriginal matrix a : \n");
	for (int k = 1; k <= 3; k++) {
		for (int l = 1; l <= 3; l++) printf("%12.6f\t", a[k][l]);
		printf("\n");
	}

	gaussj(ai, 3, r, 3);
	printf("\nInverse of matrix a : \n");
	for (int k = 1; k <= 3; k++) {
		for (int l = 1; l <= 3; l++) printf("%12.6f", ai[k][l]);
		printf("\n");
	}

	printf("\nsolving the equation\n");
	for (int l = 1; l <= 3; l++){
		sol[l][1] = 0.0;
		for (int j = 1; j <= 3; j++){
			sol[l][1] += ai[l][j] * b[j][1];
		}
		printf("a%d is %12.6f\n", l,sol[l][1]);
	}

	b[1][1] = ypxsum;
	b[2][1] = ypysum;
	b[3][1] = ypsum;
	

	for (int l = 1; l <= 3; l++){
		sol[l][1] = 0.0;
		for (int j = 1; j <= 3; j++){
			sol[l][1] += ai[l][j] * b[j][1];
		}
		printf("a%d is %12.6f\n", l+3, sol[l][1]);
	}
	
	fclose(fp);
	
	free_matrix(r,1,NP,1,MP);
	free_matrix(b,1,NP,1,MP);
	
	free_matrix(ai,1,NP,1,NP);
	free_matrix(a,1,NP,1,NP);
	
	return 0;
}
#undef NRANSI