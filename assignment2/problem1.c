/*
 ============================================================================
 Name        : problem1.c
 Author      :
 Version     :
 Copyright   : Your copyright notice
 Description : Newton Root
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double f(double x);
double f_prime(double x);
double NewtonRoot (double x0, double tol, int maxniter);


#define NUM_ITERATIONS 100
#define TOLERANCE 0.0001

void main(void)
{

double x0 = 6.1;

printf("\nx_zero=%f",NewtonRoot (x0, TOLERANCE, NUM_ITERATIONS));

}

// declare NewtonRoot function

double NewtonRoot (double x0, double tol, int maxniter)
{
	int iter =1;
	double x = x0; // I think I should be able to remove this step
	double x_prev, delta_x;

	while(iter <= maxniter)
	{
		x_prev = x;

		// defensive coding block...
		if(f_prime(x)==0)
		{
			printf("f_prime(x)==0: undefined output");
			break; // probably should break to some useful spot... in real life abort would kill the patient anyway!!!
		} // end defensive coding block;

		else // what would happen if I didn't have the else here?
		{
		delta_x = -f(x)/f_prime(x); // might want to create a variable double f_prime so I dont call the function out twice?
		x = x+delta_x;
		}

		if(abs(x-x_prev)<tol)
		break; // otherwise it will keep iterating the loop
		else iter++;
	}
printf("\nnumber of iterations= %d",iter);
return x;

}

// input declare functions;


double f(double x)
{
	return -0.09375*x*x*x + 1.125*x*x - 3.375*x + 6;
}

double f_prime (double x)
{
	return  -3*0.09375*x*x + 2.25*x - 3.375;
}
