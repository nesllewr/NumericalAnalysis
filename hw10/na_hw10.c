/*
Levenberg-Marquardt Method - Algorithm
1. 	Guess a_cur
2. 	compute x^2(a_cur)
	x^2(a) = Sigma((y - f(x, a) / sigma)^2)
3. 	pick a modest lambda, say lambda = 0.001
4. 	solve Hessian' * delta_a = -gradient_x^2(a_cur) where Hessian' = Hessian + lambda * I
5. 	if x^2(a_cur + delta_a) >= x^2(a_cur)
		increase lambda by a factor of 10 and go to 4
	else
		decrease lambda by a factor of 10 and update a_cur = a_cur + delta_a and go to 4
*/

/* Driver for routine mrqmin */
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#define NRANSI
#include <math.h>
#include "nr.h"
#include "nrutil.h"

#define NPT 180
#define MA 8
#define SPREAD 0.001
#define MAX_SIZE 10000

static int counter = 0;

void func(float x, float a[], float *y, float dyda[], int na, struct Point p[]);
int main(void)
{
	struct Point p[MAX_SIZE];
	long idum=(-911);
	int i,*ia,iter,itst,j,k,mfit=MA;
	float alamda,chisq,ochisq,*x,*y,*sig,**covar,**alpha;
	float *x2;

	static float a[MA+1] = {1.0,5.0,2.0,3.0,2.0,5.0,3.0,4.0};
	static float gues[MA+1] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

	char* filename = "good1.dat";
	FILE *fp;

	if ((fp = fopen(filename, "r")) == NULL) {
		nrerror("Data file not found\n");
	}


	ia=ivector(1,MA);
	x=vector(1,NPT);
	x2 = vector(1,NPT);
	y=vector(1,NPT);
	sig=vector(1,NPT);
	covar=matrix(1,MA,1,MA);
	alpha=matrix(1,MA,1,MA);

	while (!feof(fp)) {
		fscanf(fp, "%f %f %f %f", &p[counter].x, &p[counter].y, &p[counter].x_p, &p[counter].y_p);
		counter++;
	}



	/* In x,y */
	for (i=1;i<=NPT;i++) {
		float divide = a[7] * p[i].x + a[8] * p[i].y + 1;
		float up1 = a[1] * p[i].x + a[2] * p[i].y + a[3];
		float up2 = a[4] * p[i].x + a[5] * p[i].y + a[6];
		x[i]=p[i].x + p[i].y;
		//x2[i] = p[i].y;
		y[i] += ((up1/divide)+ (up2/divide));
		y[i] *= (1.0+SPREAD*gasdev(&idum));
		sig[i] = SPREAD*y[i];
	}
	//printf("where");
	

	for (i=1;i<=MA;i++) ia[i]=1;
	for (i=1;i<=MA;i++) a[i]=gues[i];
	for (iter=1;iter<=2;iter++) {
		alamda = -1.0;
		//printf("alamda before mrqmin: %12.6f\n", alamda);
		mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,func,&alamda, p);
		k=1;
		itst=0;
		for (;;) {
			printf("\n%s %2d %17s %10.4f %10s %12.6f\n","Iteration #",k,
				"chi-squared:",chisq,"alamda:", alamda);
			printf("%8s %8s %8s %8s %8s %8s %8s %8s\n",
				"a[1]","a[2]","a[3]","a[4]","a[5]","a[6]","a[7]","a[8]");
			for (i=1;i<=MA;i++) printf("%9.4f",a[i]);
			printf("\n");
			ochisq=chisq;
			mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,func,&alamda, p);
			k++;
			
			if (chisq > ochisq)
				itst=0;
			else if (fabs(ochisq-chisq) < 0.1)
				itst++;
			if (itst < 4) continue;
			alamda=0.0;
			mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,func,&alamda, p);
			printf("%12.9f", covar[1][1]);
			printf("\nUncertainties:\n");
			printf("\n a :\n");
			for (i=1;i<=8;i++) printf("%f\t", sqrt(covar[i][i]));
			printf("\n\n");
			// float divY = x[1] / 10000;
			// float divX = x[1] - ((int)divY) * 10000;
			// float tmpY = (a[1] * divX + a[2] * divY + a[3]) / (a[4] * divX + a[5] * divY + 1);
			// printf("x : %f , y: %f \n", x[1], y[1]);
			// printf("function(x): %f\n", tmpY);
			break;
		}
		if (iter == 1) {
			printf("press return to continue with constraint\n");
			(void) getchar();
			printf("holding a[4] and a[5] constant\n");
			for (j=1;j<=MA;j++) a[j] += 0.1;
			a[4]=0.1;
			ia[4]=0;
			a[5]=0.1;
			ia[5]=0;
		}
	}
	free_matrix(alpha,1,MA,1,MA);
	free_matrix(covar,1,MA,1,MA);
	free_vector(sig,1,NPT);
	free_vector(y,1,NPT);
	free_vector(x,1,NPT);
	free_ivector(ia,1,MA);
	getchar();
	return 0;
}
#undef NRANSI


void func(float x, float a[], float *y, float dyda[], int na,struct Point p[])
{
	int i;
	float fac, ex, arg;

	*y=0.0;
	for ( i = 1; i <= na; i++ ) {
		dyda[i] = 0.0;
	}
	for ( i = 0; i < counter; i++) {
		float divide = a[7] * p[i].x + a[8] * p[i].y + 1;
		float up1 = a[1] * p[i].x + a[2] * p[i].y + a[3];
		float up2 = a[4] * p[i].x + a[5] * p[i].y + a[6];
		// partial deriviation
		//printf("%d %12.9f %12.9f %12.9f %12.9f %12.9f\n", i, a[1],a[2],up1,up2,divide);
		*y += (up1 / divide )+(up2 / divide);

		//printf("%12.9f\n", *y);
		dyda[1] += -2 * p[i].x * ( -up1 + p[i].x_p * divide ) / (sqrt(divide));
		dyda[2] += -2 * p[i].y * ( -up1 + p[i].x_p * divide ) / ( sqrt(divide));
		dyda[3] += 2 * ( up1 - p[i].x_p * divide ) / ( sqrt(divide));
		dyda[4] += -2 * p[i].x * ( -up2 + p[i].y_p * divide ) / ( sqrt(divide) );
		dyda[5] += -2 * p[i].y * ( -up2 + p[i].y_p * divide ) / ( sqrt(divide) );
		dyda[6] += 2 * ( up2 - p[i].y_p * divide ) / ( sqrt(divide));
		dyda[7] += ( 2 * p[i].x * up1 * ( p[i].x_p - up1 / divide ) ) / ( sqrt(divide) ) + ( 2 * p[i].x * up2 * ( p[i].y_p - up2 / divide ) ) / ( sqrt(divide) );
		dyda[8] += ( 2 * p[i].y * up1 * ( p[i].x_p - up1 / divide ) ) / ( sqrt(divide) ) + ( 2 * p[i].y * up2 * ( p[i].y_p - up2 / divide ) ) / ( sqrt(divide));
	}
}