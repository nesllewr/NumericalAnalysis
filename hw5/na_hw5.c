#include <stdio.h>
//#include <stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define NP 10
#define MAXSTR 80

int main(int argc, char *argv[]){

    int i, j, n, m;

    float **a, *b,*sol;
    a = matrix(1,NP,1,NP);
    b = vector(1, NP); 
    sol = vector(1, NP); // A * x = b

    
    float **ainverse, **x;
    ainverse = matrix(1,NP,1,NP);
    x =  matrix(1,NP,1,NP);

    float **origin;
    origin = matrix(1,NP,1,NP);

    float *w, **u, **v;
    w = vector(1,NP);
    u = matrix(1,NP,1,NP);
    v = matrix(1,NP,1,NP);
    
    FILE *fp;
  
    if ((fp = fopen(argv[1],"r")) == NULL)
        nrerror("There is no such a file\n");
   
    fscanf(fp,"%d %d",&n,&m);

    for (i=1; i<n+1; i++) {
        for (j=1; j<n+1 ; j++) {
            fscanf(fp,"%f",&a[i][j]);
            ainverse[i][j] = a[i][j];
            origin[i][j] = a[i][j];
            u[i][j] = a[i][j];
        }
    }

    for(i=1; i<n+1; i++) {
        fscanf(fp,"%f", &b[i]);        
    }

////////////////////    
    printf("Solution using Gauss-Jordan elimination\n");

    gaussj(ainverse,n,x,m);
   
    float sum;
    for(i=1; i < n+1 ;i++){
        sum=0;
        for(j=1;j<n+1; j++){
            sum += ainverse[i][j] * b[j];
        }
        sol[i]=sum; 
    }

    for(i=1; i<n+1; i++) {
        printf("%12.6f ", sol[i]);
    }
   /////////////////// 
    printf("\nSolution using LU Decomposition\n");

    int *indx;
    indx=ivector(1,NP);
    float d;


    ludcmp(a,n,indx,&d);

    for(i =1; i < n+1 ;i++)
        d *=a[i][i];

    for(i=1; i<n+1; i++) {
        sol[i] = b[i];
    }

    lubksb(a,n,indx,sol);

    for(i=1; i<n+1; i++) {
        printf("%12.6f ", sol[i]);
    }

    //////////////////
     printf("\nSolution using Singular Value Decomposition\n");

    for(i=1; i<n+1; i++) {
        sol[i] = b[i];
    }

    float wmax = 0.0, wmin;
    for(j =1; j< n+1;j++){
        if(w[j] > wmax) wmax = w[j];
    }

    svdcmp(u,m,n,w,v);
    wmin = wmax*1.0e-6;
    for(i=1; i < n+1 ;i++){
        if(w[j] < wmin) w[j]=0.0;
    }
    svbksb(u,w,v,n,n,b,sol);
    
    for(i=1; i<n+1; i++) {
        printf("%12.6f ", sol[i]);
    }

    printf("\n");

/////////////////
    printf("\nIterative impovement of mprove\n");

    for (i=1;i< n+1 ;i++) {
        sol[i] *= (1.0 + 0.01);
        printf("%12.6f ", sol[i]);
    }
    printf("\nAfter mprove\n");
    mprove(origin, a,n,indx,b,sol);

    for(i=1; i<n+1; i++) {
        printf("%12.6f ", sol[i]);
    }
// //////////////
    printf("\nDeterminant of A : %12.6f\n",d);


    
    fclose(fp);

    free_vector(b,1,NP);
    free_vector(sol,1,NP);
    free_vector(w,1,NP);

    free_ivector(indx,1,NP);

    free_matrix(x,1,NP,1,NP);
    free_matrix(ainverse,1,NP,1,NP);
    free_matrix(a,1,NP,1,NP);
    free_matrix(v,1,NP,1,NP);
    free_matrix(u,1,NP,1,NP);
    free_matrix(origin,1,NP,1,NP);

   

    return 0;
}
#undef NRANSI
