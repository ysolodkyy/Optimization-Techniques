/*
 ============================================================================
 Name        : Final_SteepestDescent.c
 Author      : Yevgen Solodkyy
 Version     :	TUESDAY 6:20p
 Copyright   : STEEPESET DESCENT
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//double Lambda(double *p,double x,double y);
double GradNorm(double x, double y);
//void gaussElim(double *A,double *b,int n,double *x);

void SteepestDescent( double tolerance, int max_iters, double x0, double y0, double *x_min, double *y_min, double *f_min);

double f(double x, double y);

double f_x(double x, double y);

double f_xx(double x, double y);

double f_xy(double x, double y);

double f_y(double x, double y);

double f_yx(double x, double y);

double f_yy(double x, double y);


#define maxn 3
#define TOLERANCE 0.001
#define MAX_ITERS 3000

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


	while (step<=max_steps)
	{

		p[0]=-f_x(x,y);
		p[1]=-f_y(x,y);


		//printf("\nnum_step=%d", step);
		for(int index =0; index <2; index ++)
		{

			printf("\np[%d]=%.2f", index, p[index]);
		}


		LAMBDA =1;

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

	printf("\nTOT_STEP=%d", step);
}

/*************************************/



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

	del_f[0]=f_x(x,y);
	del_f[1]=f_y(x,y);

	del_f_SQRD=del_f[0]*del_f[0]+del_f[1]*del_f[1];//0;

	NORM_del_f = sqrt(del_f_SQRD);

	printf("\nNORM_del_f = %f\n",NORM_del_f);
	return NORM_del_f;
}
/*************************************/
