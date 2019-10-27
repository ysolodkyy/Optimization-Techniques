/*
 ============================================================================
 Name        : problem1.c
 Author      :
 Version     :
 Copyright   : Your copyright notice
 Description : finite difference method for finding delta_x
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>

void NewtonMin(double x0, double tol, int maxniter, double *x_min, double *f_min);

double f(double x);
double diff1_f (double x);
double diff2_f (double x);

#define NUM_ITERATIONS 7
#define TOLERANCE 0.0001
#define E 0.01 // epsilon

void main(void)
{

	double x0=0.0;
	double X_minimum, f_minimum;

	NewtonMin(x0,TOLERANCE,NUM_ITERATIONS, &X_minimum, &f_minimum);

	printf("\nx_min= %f, f_min= %f", X_minimum, f_minimum);

}

// define the NewtonMin function

void NewtonMin(double x0, double tol, int maxniter, double *x_min, double *f_min)
{
	int iter =1;
	double x=x0;
	double x_prev, delta_x;

	while(iter <= maxniter)
	{
		x_prev = x;
		// finite difference delta_x obviates need for calculating derivatives

		delta_x = -diff1_f(x)/diff2_f(x);

		x=x+delta_x;

		if(0) // not yet clear what the unacceptable value of x is.
		{
			// something here
		}

		if(fabs(x-x_prev)<tol)
		{
			break;
		}
		else
		{
			iter++;
		}
	}
	printf("\niterations=%d",iter);

	*x_min=x;
	*f_min=f(x);
}


// define the input function

double f(double x)
{
	return 20*exp(-0.5*x) + 0.09375*x*x*x - 1.125*x*x + 3.375*x + 2;
}

double diff1_f (double x)
{
	return (f(x+E)-f(x-E))/(2*E);
}


double diff2_f (double x)
{
	return (f(x+E)-2*f(x)+f(x-E))/(E*E);
}
