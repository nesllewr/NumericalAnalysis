
#include <math.h>

void fgauss(float x, float a[], float *y, float dyda[], int na)
{
	int i;
	float fac,ex,arg;

	*y=0.0;
	for (i=1;i<=na-1;i+=3) {
		arg=(x-a[i+1])/a[i+2];
		ex=exp(-arg*arg);
		fac=a[i]*ex*2.0*arg;
		*y += a[i]*ex;
		dyda[i]=ex;
		dyda[i+1]=fac/a[i+2];
		dyda[i+2]=fac*arg/a[i+2];
	}

	



	// float divisor = (a[7] * d.x + a[8] * d.y + 1);

 //    *y1 = (a[1] * d.x + a[2] * d.y + a[3]) / divisor;
 //    dyda[1] = d.x / divisor;
 //    dyda[2] = d.y / divisor;
 //    dyda[3] = 1 / divisor;
 //    dyda[7] = -d.x * *y1 / divisor;
 //    dyda[8] = -d.y * *y1 / divisor;


 //    *y2 = (a[4] * d.x + a[5] * d.y + a[6]) / divisor;
 //    dyda[4] = d.x / divisor;
 //    dyda[5] = d.y / divisor;
 //    dyda[6] = 1 / divisor;
 //    dyda[7] = -d.x * *y2 / divisor;
 //    dyda[8] = -d.y * *y2 / divisor;
}
