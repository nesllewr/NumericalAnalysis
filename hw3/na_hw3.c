#include <stdio.h>
#include <string.h>
// #include <stdlib.h>
#include <math.h>

#include "nr.h"
#include "nrutil.h"

#define  E 2.71828182845904523536
#define PI 3.14159265358979323846
#define NBMAX 20

float fx(float x) {
	return bessj0(x);
}

float fx1(float x) {
	return 10 * pow(E, -x) * sin(2 * PI * x) - 2;
}
float fx2(float x) {
	return pow(x, 2) - 2 * x * pow(E, -x) + pow(E, -2 * x);
}
float fx3(float x) {
	return cos(x + sqrt(2)) + x * (x / 2 + sqrt(2));
}
float fx4(float x) {
	return 5 * pow(x, 5) - 5 * x +2;
}
void func(float x, float *y, float *dy) {
	*y = bessj0(x);
	*dy = -bessj1(x);
}

void func1(float x, float *y, float *dy) {
	*y = fx1(x);
	*dy = 10 * pow(E, -x) * (2 * PI * cos(2 * PI * x) - sin(2 * PI * x));
}
void func2(float x, float *y, float *dy) {
	*y = fx2(x);
	*dy = 2 * pow(E, (-2 * x))*(pow(E, x) + 1) * (pow(E, x) * x - 1);
}
void func3(float x, float *y, float *dy) {
	*y = fx3(x);
	*dy = sqrt(2) + x - sin(sqrt(2) + x);
}
void func4(float x, float *y, float *dy) {
	*y = fx4(x);
	*dy = 25 * pow(x, 4) - 5;
}

int main() {

	const float X1 = 1.0;
	const float X2 = 10.0;
	const float XACC = 1e-6;

	int i;

	printf("Bracketing method\n");

	int num_roots = 100;
	float *xb1 = vector(1, NBMAX);
	float *xb2 = vector(1, NBMAX);

	zbrak(fx, X1, X2, 1e4, xb1, xb2, &num_roots);
	
	printf("Root%s found in %d interval%s: \n",
		num_roots > 1 ? "s" : "",
		num_roots,
		num_roots > 1 ? "s" : "");
	for (i = 1; i < num_roots+1; i++)
		printf("%d\t: [%.10f, %.10f]\n", i, xb1[i], xb2[i]);
	

	printf("Bisection method\n");
	for (i = 1; i < num_roots+1; i++)
		printf("%d\t: %.10f\n", i, rtbis(fx, xb1[i], xb2[i], XACC));
	

	printf("Linear interpolation\n");
	for (i = 1; i < num_roots+1; i++)
		printf("%d\t: %.10f\n", i, rtflsp(fx, xb1[i], xb2[i], XACC));
	

	printf("Secant method\n");
	for (i = 1; i < num_roots+1; i++)
		printf("%d\t: %.10f\n", i, rtsec(fx, xb1[i], xb2[i], XACC));
	

	printf("Newton-Raphson\n");
	for (i = 1; i < num_roots+1; i++)
		printf("%d\t: %.10f\n", i, rtnewt(func, xb1[i], xb2[i], XACC));
	

	printf("Newton with bracketing\n");
	for (i = 1; i < num_roots+1; i++)
		printf("%d\t: %.10f\n", i, rtsafe(func, xb1[i], xb2[i], XACC));
	

	
	printf("\n\n\nusing the routine of rtsafe.c\n");
	
	printf("problem 1\n");
	zbrak(fx1, 0.1, 1.0, 1e4, xb1, xb2, &num_roots);
	for ( i = 1; i < num_roots+1 ; i++) {
		float xacc = (1.0e-6)*(xb1[i] + xb2[i]) / 2.0;
		printf("root %d : %f\n", i, rtsafe(func1, xb1[i], xb2[i], xacc));

	}

	printf("\nproblem 2\n");
	zbrak(fx2, 0.0, 1.0, 1e4, xb1, xb2, &num_roots);
	for (i = 1; i < num_roots + 1; i++) {
		float xacc = (1.0e-6)*(xb1[i] + xb2[i]) / 2.0;
		printf("root %d : %f\n", i, rtsafe(func2, xb1[i], xb2[i], xacc));
	}

	printf("\nproblem 3\n");
	zbrak(fx3, -2.0, -1.0, 1e4, xb1, xb2, &num_roots);
	for (i = 1; i < num_roots + 1; i++) {
		float xacc = (1.0e-6)*(xb1[i] + xb2[i]) / 2.0;
		printf("root %d : %f\n", i, rtsafe(func3, xb1[i], xb2[i], xacc));
	}

	printf("\nproblem 4\n");
	zbrak(fx4, -10.0, 2.0, 1e4, xb1, xb2, &num_roots);
	for (i = 1; i < num_roots + 1; i++) {
		float xacc = (1.0e-6)*(xb1[i] + xb2[i]) / 2.0;
		printf("root %d : %f\n", i, rtsafe(func4, xb1[i], xb2[i], xacc));

	}


	printf("\n");
	free_vector(xb1, 1, NBMAX);
	free_vector(xb2, 1, NBMAX);

	return 0;
}