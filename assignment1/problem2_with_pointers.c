/*
 ============================================================================
 Name        : problem2.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double f(double x);

void goldenSearch(double a, double b, double c, float tol, int maxniter, double *x_min, double *f_min);


#define NUM_ITERATIONS 21
#define TOLERANCE 0.0001


void main(void)
{

	double x_a=0, x_c=10;

	double x_b=x_a+0.38197*(x_c-x_a);

	double X_minimum, f_minimum;

	goldenSearch(x_a, x_b, x_c, TOLERANCE, NUM_ITERATIONS, &X_minimum, &f_minimum);

	printf("\nx_min= %f, f_min= %f", X_minimum, f_minimum);
}

void goldenSearch(double a, double b, double c, float tol, int maxniter, double *x_min, double *f_min)
{
	int niter=1;

	while(niter<=maxniter)
	{
		if(abs(c-a)<tol)
		{
			break;
		}
		double d = b+0.38197*(c-b);
		if(f(d)<f(b))
		{
			a=b;
			b=d;
		}
		else
		{
			c=a;
			a=d;
		}
		niter++;
	}
	//printf("\nx_min= %f, f_min= %f", b, f(b));
	*x_min=b;
	*f_min=f(b);
}

double f(double x)
{
	return 2*exp(-0.5*x) + 0.09375*x*x*x - 1.125*x*x + 3.375*x + 2;
}
