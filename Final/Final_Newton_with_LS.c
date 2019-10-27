/*
 ============================================================================
 Name        : Final_Newtons.c
 Author      : Yevgen Solodkyy
				TUESDAY 6:20p
 Description : Final_SteepestDescent With Line Search
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void gaussElim(double *A,double *b,int n,double *x);

double quadCubicLS (double *p, double x0, double y0, double f_0, double alpha, int maxncuts);

void NewtonMult( double tolerance, int max_iters, double x0, double y0, double *x_min, double *y_min, double *f_min);

double f(double x, double y);

double f_x(double x, double y);

double f_xx(double x, double y);

double f_xy(double x, double y);

double f_y(double x, double y);

double f_yx(double x, double y);

double f_yy(double x, double y);



#define maxn 2//3
#define TOLERANCE 0.001
#define NUMiterS 3000
#define ALPHA 0.0001
#define MAXNCUTS 10


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
	int iter =0;

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

		LAMBDA = quadCubicLS (p,x, y, f(x,y), ALPHA, MAXNCUTS);

		printf("\nlambda = %lf", LAMBDA);

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

	int num_iter = iter;
	*x_min = x;
	*y_min = y;

	*f_min = f(x,y);

	printf("\nNUM_ITER=%d", num_iter);
}
/***********************************/
double quadCubicLS (double *p, double x0, double y0, double f_0, double alpha, int maxncuts)
{

	double del_f[2] = {0,0};

	double x = x0, y= y0;

	double ll =0.1, ul = 0.5;

	double lambda, lambda_prev,lambda2_prev, f_lambda, f_lambda_prev,f_lambda2_prev, lambda_hat, a, b;

	int cut = 1, success = 0;


	//ok we have p[3]

	del_f[0]= f_x(x,y);
	del_f[1]= f_y(x,y);

	//first directional derivative

	// consider taking this f_p or use *p as an input

	double f_p = (del_f[0]*p[0]+del_f[1]*p[1]); // according to the notes no need to divide by (sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]));

	// somehow we need to get delta_x, which I think is actually delta_p[3] = -f_p/f_pp.

	// delta_p is really what's called delta_x in the notes. what's the definition of delta_x ???

	double delta_x = p[0];
	double delta_y = p[1];

	/**************************************/

	// x0 is an input, so change every x0 to x.

	// for calculating a & b
	double A[2][2], B[2]; // [row][column]


	lambda = 1;

	lambda_prev = lambda;

	f_lambda = f(x+lambda*delta_x, y+lambda*delta_y); // I am keeping the function as a function of f(x,y,z)

	f_lambda_prev = f_lambda;

	while(cut <= maxncuts)
	{
		if(f_lambda <= f_0 + alpha*lambda*f_p)
		{
			success = 1;
			break;
		}
		if(cut ==1)
		{
			lambda = - f_p/(2*(f_lambda - f_p-f_0));
		}
		else
		{
			// calculate a & b:

			B[0] = f_lambda_prev-f_p*lambda_prev-f_0;
			B[1] = f_lambda2_prev-f_p*lambda2_prev-f_0;


			// A[row][column]
			A[0][0] = 1/pow(lambda_prev,2);
			A[0][1] = -1/pow(lambda2_prev,2);

			A[1][0] = -lambda2_prev/pow(lambda_prev,2);
			A[1][1] = lambda_prev/pow(lambda2_prev,2);


			a = (A[0][0]*B[0]+A[0][1]*B[1])/(lambda_prev - lambda2_prev);
			b = (A[1][0]*B[0]+A[1][1]*B[1])/(lambda_prev - lambda2_prev);


			// finish calculating a & b;
			lambda  = (-b + sqrt(b*b -3*a*f_p))/(3*a);
		}

		if(lambda > ul*lambda_prev)
		{
			lambda = ul*lambda_prev;
		}
		if(lambda < ll*lambda_prev)
		{
			lambda = ll*lambda_prev;
		}

		f_lambda = f(x+lambda*delta_x, y+lambda*delta_y);
		lambda2_prev = lambda_prev;
		lambda_prev = lambda;

		f_lambda2_prev = f_lambda_prev;
		f_lambda_prev = f_lambda;

		cut++;
	}
	if (success == 1)
	{
		lambda_hat = lambda;
	}
	else
	{
		lambda_hat = 0;
	}

	printf("\nncuts=%d", cut);

	return lambda_hat;
}

/***********************************/


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
