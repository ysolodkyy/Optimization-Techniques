/*
 ============================================================================
 Name        : Final_Newtons.c
 Author      : Yevgen Solodkyy
 	 	 	 TUESDAY 6:20p
 Description : Final_SteepestDescent NO Line Search
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//double Lambda(double *p,double x,double y);

void gaussElim(double *A,double *b,int n,double *x);

void NewtonMult( double tolerance, int max_iters, double x0, double y0, double *x_min, double *y_min, double *f_min);

double f(double x, double y);

double f_x(double x, double y);

double f_xx(double x, double y);

double f_xy(double x, double y);

double f_y(double x, double y);

double f_yx(double x, double y);

double f_yy(double x, double y);



#define maxn 3
#define TOLERANCE 0.001
#define NUMiterS 3000

int main(void)
{
	double X_min, Y_min, F_min,X0, Y0;

	printf("enter X0\n");
	scanf("%lf", &X0);
	printf("enter Y0\n");
	scanf("%lf", &Y0);

	printf("\nX0= %.2f, \nY0= %.2f, \n\n", X0, Y0);


	NewtonMult( TOLERANCE, NUMiterS, X0, Y0, &X_min, &Y_min, &F_min);

	printf("\nX_min= %.3f, Y_min= %.3f, F_min= %.3f;",X_min, Y_min, F_min);

	return EXIT_SUCCESS;
}



void NewtonMult( double tolerance, int max_iters, double x0, double y0, double *x_min, double *y_min, double *f_min)

{
	int iter =0; // =1

	double x = x0, y = y0;

	double p[2]={0,0}, g[2]={0,0},g_neg[2]={0,0};


	double LAMBDA;

	double H[4]; // [row][column]

	// variables to calculate the norm :
	double NORM_g, g_SQRD;

	//  g = del_f:

	g[0]= f_x(x,y);
	g[1]= f_y(x,y);

	printf("\ng[0]=%f",g[0]);
	printf("\ng[1]=%f",g[1]);


	while (iter<=max_iters)
	{

		// 1st row
		H[0]=f_xx(x,y);
		H[1]=f_yx(x,y);

		// 2nd row
		H[2]=f_xy(x,y);
		H[3]=f_yy(x,y);

		printf("\nsecond derivatives\n");

		for (int i =0; i<5; i++)
		{
			printf("\n\nH[%d]=%f",i, H[i]);
		}

		g_neg[0] = -g[0];
		g_neg[1] = -g[1];

		// solve Hp=-g; // mistake somewhere here?

		printf("\n solve values of p:");

		gaussElim(H,g_neg,2,p); // changed 3 to 2

		// check the values of p:
		printf("\n\np[0]=%f",p[0]);
		printf("\np[1]=%f",p[1]);

		LAMBDA =1;

		// update variables:
		x = x+LAMBDA*p[0];
		y = y+LAMBDA*p[1];

		// g = del_f:
		g[0]= f_x(x,y);
		g[1]= f_y(x,y);

		printf("\nx= %.2f, \ny= %.2f\n\n", x, y);

		// calculate the NORM_g:

		g_SQRD=g[0]*g[0]+g[1]*g[1];//0;

		NORM_g = sqrt(g_SQRD);

		printf("NORM_g = %f\n",NORM_g);

		printf("\nEND LOOP\n");

		if(NORM_g<tolerance)
		{
			break;
		}
		else
		{
			iter++;
		}

	}

	*x_min = x;
	*y_min = y;

	*f_min = f(x,y);

	printf("\nTOT_ITER=%d", iter);
}

// Gauss Elimination

void gaussElim(double *A,double *b,int n,double *x)
{
	double dtemp, mult, pivot;
	int i, j, k, npivots, pivotrow;

	npivots = 0;
	for(i=0; i<(n-1); i++) {
		pivot = fabs(A[n*i+i]);
		pivotrow = i;
		for(j=(i+1); j<n; j++) {
			dtemp = fabs(A[n*j+i]);
			if(dtemp>pivot) {
				pivot = dtemp;
				pivotrow = j;
			}
		}
		if(fabs(pivot)<1.0e-10) break;
		if(pivotrow!=i) {
			for(j=i; j<n; j++) {
				dtemp = A[n*i+j];
				A[n*i+j] = A[n*pivotrow+j];
				A[n*pivotrow+j] = dtemp;
			}
			dtemp = b[i];
			b[i] = b[pivotrow];
			b[pivotrow] = dtemp;
			npivots++;
		}
		for(j=(i+1); j<n; j++) {
			mult = -A[n*j+i]/A[n*i+i];
			for(k=i; k<n; k++) A[n*j+k] += mult*A[n*i+k];
			b[j] += mult*b[i];
		}
	}

	for(i=(n-1); i>=0; i--) {
		x[i] = b[i];
		for(j=(i+1); j<n; j++) x[i] = x[i] - A[n*i+j]*x[j];
		x[i] /= A[n*i+i];
	}
}

/////




////


////

double f(double x, double y)
{
	return (1.0-x)*(1.0-x)+105.0*(y-x*x)*(y-x*x);
}

//FIRST DERIVATIVES

double f_x(double x, double y)
{
//	return 424.0*x*x*x-4*x*(105*y+1);
	return 420.0*x*x*x - 420.0*x*y + 2*x - 2.0;
}


double f_y(double x, double y)
{
	return 210.0*y - 210.0*x*x;
}


// SECOND DERIVATIVES -- all are constant

//X-

double f_xx(double x, double y)
{
//	return 1272*x*x-4*(105*y+1);//0; //can I do this? yes, I can.
	return 1260*x*x - 420.0*y + 2.0;
}

double f_xy(double x, double y)
{
	return -420.0*x;
}

//Y-


double f_yx(double x, double y)
{
	return -420.0*x;
}

double f_yy(double x, double y)
{
	return 210;
}
