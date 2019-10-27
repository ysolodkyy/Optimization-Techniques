/*
 ============================================================================
 Name        : Final_SteepestDescent_with_LS.c
 Author      : Yevgen Solodkyy
 Version     : TUESDAY 6:20p.
 Copyright   : Final_Newton_with_LS
 Description : Final_SteepestDescent_with_LS
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double Lambda(double *p,double x,double y);
double GradNorm(double x, double y);

double quadCubicLS (double *p, double x0, double y0, double f_0, double alpha, int maxncuts);

void SteepestDescent( double tolerance, int max_iters, double x0, double y0, double *x_min, double *y_min, double *f_min);

double f(double x, double y);

double f_x(double x, double y);

double f_xx(double x, double y);

double f_xy(double x, double y);

double f_y(double x, double y);

double f_yx(double x, double y);

double f_yy(double x, double y);



#define maxn 2//3
#define TOLERANCE 0.001
#define MAX_ITERS 3000
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


	SteepestDescent( TOLERANCE, MAX_ITERS, X0, Y0, &X_min, &Y_min, &F_min);

	printf("\nX_min= %.3f, Y_min= %.3f, F_min= %.3f;",X_min, Y_min, F_min);

	return 0;
}

void SteepestDescent( double tolerance, int max_steps, double x0, double y0, double *x_min, double *y_min, double *f_min)

{
	int dimensions = 2; // specific to this problem
	int step =0;

	double x = x0, y = y0;

	double p[dimensions];

	double LAMBDA;

	//double g_SQRD, NORM_g;
	//double g[2];

	while (step<=max_steps)
	{

		p[0]=-f_x(x,y);
		p[1]=-f_y(x,y);


		printf("\nnum_step=%d", step);
		for(int index =0; index <2; index ++)
		{

			printf("\np[%d]=%.2f", index, p[index]);
		}

		LAMBDA = quadCubicLS (p,x, y, f(x,y), ALPHA, MAXNCUTS); //Lambda(p,x,y);


		printf("\np[0] = %lf",p[0]);
		printf("\np[1] = %lf",p[1]);

		printf("\nLAMBDA = %lf",LAMBDA);



		x = x+LAMBDA*p[0];
		y = y+LAMBDA*p[1];





		if(GradNorm(x,y)< tolerance)
		{
			break;
		}
		else
		{
			step++;
		}

	}

	*x_min = x;
	*y_min = y;

	*f_min = f(x,y);

	printf("\nnum_step=%d", step);
}


/////////BEGIN BICUBIC LINE SEARCH //////////



double quadCubicLS (double *p, double x0, double y0, double f_0, double alpha, int maxncuts)
{

	double del_f[2] = {0,0};

	double x = x0, y= y0;

	double ll =0.1, ul = 0.5;

	double lambda, lambda_prev,lambda2_prev, f_lambda, f_lambda_prev,f_lambda2_prev, lambda_hat, a, b;

	int cut = 1, success = 0;

	del_f[0]= f_x(x,y);
	del_f[1]= f_y(x,y);

	//first directional derivative

	double f_p = (del_f[0]*p[0]+del_f[1]*p[1]); // according to the notes no need to divide by (sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]));

	// somehow we need to get delta_x, which I think is actually delta_p[3] = -f_p/f_pp.

	// delta_p is really what's called delta_x in the notes. what's the definition of delta_x ???

	double delta_x = p[0];
	double delta_y = p[1];

	// ---------------------------------------//

	// x0 is an input, so change every x0 to x.

	// for calculating a & b
	double A[2][2], B[2]; // [row][column]


	lambda = 1;

	lambda_prev = lambda;

	f_lambda = f(x+lambda*delta_x, y+lambda*delta_y);

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

		f_lambda = f(x+lambda*delta_x, y+lambda*delta_y); // incorrectly had alpha instead of lambda
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


//// end  BICUBIC LINE SEARCH  ////////
/*
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

		f_lambda = f(x+alpha*delta_x, y+alpha*delta_y);
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



*/


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


/////////


/*************************************/
double GradNorm(double x, double y)
{
	double del_f[2];
	double del_f_SQRD, NORM_del_f;

	printf("\n calculate Grad_f with x=%lf, y=%lf",x,y);

	del_f[0]=f_x(x,y);
	del_f[1]=f_y(x,y);

	del_f_SQRD=del_f[0]*del_f[0]+del_f[1]*del_f[1];//0;

	NORM_del_f = sqrt(del_f_SQRD);

	printf("\nNORM_del_f = %f\n",NORM_del_f);
	return NORM_del_f;
}
/*************************************/



double Lambda(double *p,double x,double y)
{

	double LAMBDA, lambda_denominator, lambda_numerator;

	double del_f[2];
	double del2_f[2][2]; // [row][column]

	// this is an intermediate step vector for del2(f(x_k))Pk
	double del2_Pk[2];

	// begin calculating LAMBDA
	del_f[0]=f_x(x,y);
	del_f[1]=f_y(x,y);

	lambda_numerator = -1*(p[0]*del_f[0]+p[1]*del_f[1]);

	//printf("\nlambda_numerator=%f",lambda_numerator);

	// (2) calculate the denominator for lambda:

	// del2_f[row][column]:

	del2_f[0][0]=f_xx(x,y);
	del2_f[1][0]=f_xy(x,y);

	del2_f[0][1]=f_yx(x,y);
	del2_f[1][1]=f_yy(x,y);

	// intermediate matrix multiplication step: del2_Pk[row]:=del2_f[row][column]*p[row]:

	del2_Pk[0]	=	p[0]*del2_f[0][0] + p[1]*del2_f[0][1];
	del2_Pk[1]	=	p[0]*del2_f[1][0] + p[1]*del2_f[1][1];

	// so then p[column]*del2_Pk[row]:

	lambda_denominator =p[0]*del2_Pk[0]+p[1]*del2_Pk[1];


	if(lambda_denominator==0)
	{
		printf("\n !! WARNING !! \n\n\tlambda_denominator==0; set LAMBDA=1");
		LAMBDA=1;//break;
	}
	else
	{
		LAMBDA = lambda_numerator/lambda_denominator;
	}

	printf("\nlambda=%.3f\n",LAMBDA);

	// END calculating LAMBDA
	return LAMBDA;
}


