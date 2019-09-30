#include <stdio.h>
#include <string.h>
#include <math.h>

#include "nr.h"
#include "nrutil.h"


#define PI 3.14159265358979323846

#define NBMAX 20
#define JMAX 40

float fx(float R) {
    return exp(-0.005 * R) * cos(sqrt(2000 - R * R * 0.01) * 0.05) - 0.01;
}

float fx1(float x){
	return 0.0885 * PI * pow(x*x + 0.81,1.5) -x;
}
float fx2(float T){
	return -1.2 + 0.99403 + 1.671e-4 * T + 9.7215e-8 * pow(T, 2) - 9.5838e-11 * pow(T, 3) + 1.9520e-14 * pow(T, 4);
}

void func(float R, float *y, float *dy) {
    *y = fx(R);
    *dy = exp(-0.005 * R) * (sin(sqrt(2000 - R * R * 0.01) * 0.05) / sqrt(2000 - R * R * 0.01) * 0.05 - cos(sqrt(2000 - R * R * 0.01) * 0.05)) * 0.005;
}

void func1(float x, float *y, float *dy){
	*y = fx1(x);
	*dy = 0.0885*PI*pow(x*x+0.81,0.5)*3*x-1;
}
void func2(float T, float *y, float *dy){
	*y = fx2(T);
	*dy = 1.671e-4 + 9.7215e-8 * 2 * T - 9.5838e-11 * 3 * pow(T, 2) + 1.9520e-14 * 4 * pow(T, 3);
}

float relative_error(float x, float y){
	return (x - y) / y;
}

int main(int argc, char *argv[]){

	const float x1 = 0.0;
	const float x2 = 400.0;
	float root1, root2;

	printf("R.E. of 1e-4\n");

	root1 = rtbis(fx, x1, x2, 1e-4);
	printf("Bisection :  %f\n", root1);
	printf("relative_error : %f\n",fabs(relative_error(root1,328)));


	root1 = rtflsp(fx, x1, x2, 1e-4);
	printf("Linear interpolation :  %f\n", root1);
	printf("relative_error : %f\n",fabs(relative_error(root1,328)));

	root1 = rtsec(fx, x1, x2, 1e-4);	
	printf("Secant method :  %f\n", root1);
	printf("relative_error : %f\n",fabs(relative_error(root1,328)));

	//root1 = rtnewt(func, x1, x2, 1e-4);		
	//printf("Newton-Raphson :  %f\n", rtnewt(func, x1,x2, 1e-4));
	//printf("relative_error : %f\n",fabs(relative_error(root1,328)));

	root1 = rtsafe(func, x1, x2, 1e-4);	
	printf("Newton with bracketing :  %f\n", root1);
	printf("relative_error : %f\n",fabs(relative_error(root1,328)));

	

	printf("\nR.E. of 1e-6\n");

	root2 = rtbis(fx, x1, x2, 1e-6);
	printf("Bisection :  %f\n", root2);
	printf("relative_error : %f\n",fabs(relative_error(root2,328)));


	root2 = rtflsp(fx, x1, x2, 1e-6);
	printf("Linear interpolation :  %f\n", root2);
	printf("relative_error : %f\n",fabs(relative_error(root2,328)));

	// root2 = rtsec(fx, x1, x2, 1e-6);	
	// printf("Secant method :  %f\n", root2);
	// printf("relative_error : %f\n",fabs(relative_error(root2,328)));

	//root2 = rtnewt(func, x1, x2, 1e-6);		
	//printf("Newton-Raphson :  %f\n", root2);
	//printf("relative_error : %f\n",fabs(relative_error(root2,328)));

	root2 = rtsafe(func, x1, x2, 1e-6);	
	printf("Newton with bracketing :  %f\n", root2);
	printf("relative_error : %f\n",fabs(relative_error(root1,328)));


	printf("\nProblem 8.32\n");

	float root3;
	const float x3 = 0.2;
	const float x4 = 0.25;

	root3 = rtbis(fx1, x3, x4, 1e-4);
	printf("Bisection :  %f\n", root3);
	//printf("relative_error : %f\n",fabs(relative_error(root3,0.221)));


	root3 = rtflsp(fx1, x3, x4, 1e-4);
	printf("Linear interpolation :  %f\n", root3);
	//printf("relative_error : %f\n",fabs(relative_error(root3,0.221)));

	root3 = rtsec(fx1, x3, x4, 1e-4);	
	printf("Secant method :  %f\n", root3);
	//printf("relative_error : %f\n",fabs(relative_error(root3,0.221)));

	// root3 = rtnewt(func1, x3, x4, 1e-4);		
	// printf("Newton-Raphson :  %f\n", root3);
	// printf("relative_error : %f\n",fabs(relative_error(root3,0.221)));

	root3 = rtsafe(func1, x3, x4, 1e-4);	
	printf("Newton with bracketing :  %f\n", root3);
	//printf("relative_error : %f\n",fabs(relative_error(root3,0.221)));

	printf("\nProblem 8.32\n");

	float root4;
	const float x5 = 1000;
	const float x6 = 1500;
	root4 = rtbis(fx2, x5, x6, 1e-4);
	printf("Bisection :  %f\n", root4);
	//printf("relative_error : %f\n",fabs(relative_error(root4,1126)));


	root4 = rtflsp(fx2, x5, x6, 1e-4);
	printf("Linear interpolation :  %f\n", root4);
	//printf("relative_error : %f\n",fabs(relative_error(root4,1126)));

	root4 = rtsec(fx2, x5, x6, 1e-4);	
	printf("Secant method :  %f\n", root4);
	//printf("relative_error : %f\n",fabs(relative_error(root4,1126)));

	root4 = rtnewt(func2, x5, x6, 1e-4);		
	printf("Newton-Raphson :  %f\n", root4);
	//printf("relative_error : %f\n",fabs(relative_error(root4,1126)));

	root4 = rtsafe(func2, x5, x6, 1e-4);	
	printf("Newton with bracketing :  %f\n", root4);
	//printf("relative_error : %f\n",fabs(relative_error(root4,1126)));

	


	return 0;
}