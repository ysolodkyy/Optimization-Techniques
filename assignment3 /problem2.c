/*
 ============================================================================
 Name        : problem2.c
 Author      :
 Version     :
 Copyright   : Your copyright notice
 Description : Newton Secant method
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "math.h"

double f(double x);
double diff1_f (double x);
double diff2_f (double x);


void NewtonSec(double x0, double tol, int maxniter, double *x_min, double *f_min);

#define NUM_ITERATIONS 21
#define TOLERANCE 0.0001
#define E 0.01 // epsilon

int main(void)
{

	double x0=0.0;
	double X_minimum, f_minimum;

	NewtonSec(x0,TOLERANCE,NUM_ITERATIONS, &X_minimum, &f_minimum);

	printf("\nx_min= %f, f_min= %f", X_minimum, f_minimum);
return 0;
}




void NewtonSec(double x0, double tol, int maxniter, double *x_min, double *f_min)
{
	int iter =1;
	double x=x0;
	double x_prev, f_curr, f_prev, f2_prev;
	double delta_x, delta_x_prev, delta_x2_prev, numer, denom;

	x_prev=x;

	f_curr = f(x);
	f_prev=f_curr;

	delta_x=0.0;
	delta_x_prev=delta_x;

	while(iter <= maxniter)
	{
		delta_x2_prev=delta_x_prev;
		delta_x_prev=delta_x;

		if(iter<3) // should be 2, he said
		{
			delta_x = -diff1_f(x)/diff2_f(x);

			x_prev =x;
			x=x+delta_x;

			//printf("delta_x = %f", delta_x);
		}
		else
		{
				numer = (delta_x*delta_x)*f2_prev-(delta_x_prev+delta_x2_prev)*(delta_x_prev+delta_x2_prev)*f_prev+delta_x2_prev*(2*delta_x_prev+delta_x2_prev)*f_curr;

				denom = 2*delta_x_prev*f2_prev-2*(delta_x_prev+delta_x2_prev)*f_prev+2*delta_x2_prev*f_curr;

				delta_x=-numer/denom;

				x_prev =x;
				x=x+delta_x;
		}

		if(0) // not yet clear what the unacceptable value of x is.
		{
			// something here
		}

		if(fabs(x-x_prev)<tol)
		{
			//printf("\nBreak condition: x,  x-x_prev = %f", x-x_prev);
			break;
		}
		else
		{
			f2_prev = f_prev;
			f_prev = f_curr;
			f_curr = f(x);
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
	return -10*exp(-0.5*x) + 3*0.09375*x*x - 2*1.125*x + 3.375;
}


double diff2_f (double x)
{
	return 5*exp(-0.5*x) + 6*0.09375*x - 2*1.125;
}
