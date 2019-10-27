/*
 ============================================================================
 Name        : problem1.c
 Author      : Yevgen Solodkyy
 Version     :
 Copyright   : Your copyright notice
 Description : Bysector Method
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#define ITERATIONS 21

double f(double x);
double bisection(double a, double b, double f_a, double f_b, double tol, int maxniter);

void main(void)
{

	double b=10, a=0; // the input variables
	double tol = 0.0001;

printf("x_zero=%f", bisection(a,b,f(a),f(b),tol,ITERATIONS)); // using 7 as the number of iterations for now.

}


// function declarations:

double f(double x)
{
	double y = -0.09375*x*x*x + 1.125*x*x - 3.375*x + 6;
	return y;
}


// declare the bisector function

double bisection(double a, double b, double f_a, double f_b, double tol, int maxniter)

{
int niter = 1;
double c;

while(niter<=maxniter)
{
 c = (a+b)/2;

double f_c=f(c); // will need to call out the other function from here.

if((b-a) < tol)
break;

if((f_c >0 && f_a >0) || (f_c <0 && f_a <0)) // check if the sign of the functions is the same
{
a = c;
f_a = f_c;
}
else
{
b = c;
f_b = f_c;
}

niter++;

}
double x_zero;

x_zero=c;

return x_zero;
}
