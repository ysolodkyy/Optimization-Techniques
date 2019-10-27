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

void goldenSearch(double a, double b, double c, float tol, int maxniter);


#define NUM_ITERATIONS 21
#define TOLERANCE 0.0001

//static double out[2]; // consider having this guy inside
// create a struct type to output values of the function x_min & f(x_min).


void main(void)
{

	double x_a=0, x_c=10; // declare the boundary values

	double x_b=x_a+0.38197*(x_c-x_a);// declare b;


	goldenSearch(x_a, x_b, x_c, TOLERANCE, NUM_ITERATIONS);

	//printf("\nx_min= %f, f_min= %f", out[0], out[1]);
}

void goldenSearch(double a, double b, double c, float tol, int maxniter)
{
	int niter=1;

	//static double out[2]; // consider having this guy inside

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
	printf("\nx_min= %f, f_min= %f", b, f(b));
/*
	out[0]=b;
	out[1]=f(b);
*/
}

double f(double x)
{
	return 2*exp(-0.5*x) + 0.09375*x*x*x - 1.125*x*x + 3.375*x + 2;
}
