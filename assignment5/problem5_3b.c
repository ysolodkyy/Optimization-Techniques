/*
 ============================================================================
 Name        : problem5_3b.c
 Author      : YEVGEN SOLODKYY HW_5
 Description : ConjGrad
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void ConjGrad( double tolerance, int max_steps, double x0, double y0, double z0, double *x_min, double *y_min, double *z_min, double *f_min);

double f(double x, double y, double z);

double f_x(double x, double y, double z);

double f_xx(double x, double y, double z);

double f_xy(double x, double y, double z);

double f_xz(double x, double y, double z);

double f_y(double x, double y, double z);

double f_yx(double x, double y, double z);

double f_yy(double x, double y, double z);

double f_yz(double x, double y, double z);

double f_z(double x, double y, double z);

double f_zx(double x, double y, double z);

double f_zy(double x, double y, double z);

double f_zz(double x, double y, double z);



#define TOLERANCE 0.00001
#define NUMSTEPS 15

int main(void)
{
	double X_min, Y_min, Z_min, F_min,X0, Y0, Z0;

	printf("enter X0\n");
	scanf("%f", &X0);
	printf("enter Y0\n");
    scanf("%f", &Y0);
	printf("enter Z0\n");
    scanf("%f", &Z0);


	printf("\nX0= %.2f, \nY0= %.2f, \nZ0= %.2f\n\n", X0, Y0, Z0);


	ConjGrad( TOLERANCE, NUMSTEPS, X0, Y0, Z0, &X_min, &Y_min, &Z_min,&F_min);

	printf("\nX_min= %.3f, Y_min= %.3f, Z_min= %.3f, F_min= %.3f;",X_min, Y_min, Z_min, F_min);

	return EXIT_SUCCESS;
}



void ConjGrad( double tolerance, int max_steps, double x0, double y0, double z0, double *x_min, double *y_min, double *z_min, double *f_min)

{

	//int dimensions = 3; // specific to this problem
	int step =1;

	double x = x0, y = y0, z=z0;

	// try this and then see if it can be reset to all zeros
	double p[3]={0,0,0}, g[3]={0,0,0};

	// for calculating LAMBDA

	double LAMBDA, lambda_denominator, lambda_numerator, Beta, Beta_numerator, Beta_denominator;
	double del_f[3];
	double del2_f[3][3], H[3][3]; // [row][column]


	// this is an intermediate step vector for del2(f(x_k))Pk
	double del2_Pk[3], Hp[3];

	// variables to calculate the norm :
	double NORM_g, g_SQRD;


	g[0]= f_x(x,y,z);
	g[1]= f_y(x,y,z);
	g[2]= f_z(x,y,z);


	printf("\ng[0]=%f",g[0]);
	printf("\ng[1]=%f",g[1]);
	printf("\ng[2]=%f",g[2]);

	while (step<=max_steps)
	{
		printf("\nnum_step=%d", step);

		del_f[0]=f_x(x,y,z);
		del_f[1]=f_y(x,y,z);
		del_f[2]=f_z(x,y,z);


		// H[row][column]: Because all second derivatives are constant, H[][] will be constant throughout

		H[0][0]=f_xx(x,y,z);
		H[1][0]=f_xy(x,y,z);
		H[2][0]=f_xz(x,y,z);

		H[0][1]=f_yx(x,y,z);
		H[1][1]=f_yy(x,y,z);
		H[2][1]=f_yz(x,y,z);

		H[0][2]=f_zx(x,y,z);
		H[1][2]=f_zy(x,y,z);
		H[2][2]=f_zz(x,y,z);

		for(int i=0;i<=2;i++)
		{
			for(int j=0;j<=2;j++)
			{
				printf("\nH[%d][%d]=%f",i,j,H[i][j]);
			}
			printf("\n");
		}

		if(step ==1)
		{
			p[0]=-g[0];
			p[1]=-g[1];
			p[2]=-g[2];
		}
		else
		{
			// calculate Beta
			printf("\nCALCULATE Beta:");

			printf("\n\np[0]=%f",p[0]);
			printf("\np[1]=%f",p[1]);
			printf("\np[2]=%f",p[2]);
			//intermediate matrix multiplication step: Hp[row]:=H[row][column]*p[row]:

			Hp[0]	=	p[0]*H[0][0] + p[1]*H[0][1] + p[2]*H[0][2];
			Hp[1]	=	p[0]*H[1][0] + p[1]*H[1][1] + p[2]*H[1][2];
			Hp[2]	=	p[2]*H[2][0] + p[2]*H[2][1] + p[2]*H[2][2];

			printf("\n\nHp[0]=%f",Hp[0]);
			printf("\nHp[1]=%f",Hp[1]);
			printf("\nHp[2]=%f",Hp[2]);

			// so then p[column]*Hp[row]:

			printf("\n\ng[0]=%f",g[0]);
			printf("\ng[1]=%f",g[1]);
			printf("\ng[2]=%f",g[2]);

			Beta_numerator =g[0]*Hp[0]+g[1]*Hp[1]+g[2]*Hp[2];
			printf("\nBeta_numerator=%f",Beta_numerator);



			printf("\n\np[0]=%f",p[0]);
			printf("\np[1]=%f",p[1]);
			printf("\np[2]=%f",p[2]);

			Beta_denominator =p[0]*Hp[0]+p[1]*Hp[1]+p[2]*Hp[2];
			printf("\nBeta_denominator=%f",Beta_denominator);

			Beta = Beta_numerator/Beta_denominator;

			printf("\nBeta=%f",Beta);
			// update p[]


			p[0]=-g[0]+ Beta*p[0];
			p[1]=-g[1]+ Beta*p[1];
			p[2]=-g[2]+ Beta*p[2];
		}

		printf("\nCALCULATE LAMBDA:");

		lambda_numerator = -1*(p[0]*del_f[0]+p[1]*del_f[1]+p[2]*del_f[2]);

		// (2) calculate the denominator for lambda:

		// del2_f[row][column]:

		del2_f[0][0]=f_xx(x,y,z);
		del2_f[1][0]=f_xy(x,y,z);
		del2_f[2][0]=f_xz(x,y,z);

		del2_f[0][1]=f_yx(x,y,z);
		del2_f[1][1]=f_yy(x,y,z);
		del2_f[2][1]=f_yz(x,y,z);

		del2_f[0][2]=f_zx(x,y,z);
		del2_f[1][2]=f_zy(x,y,z);
		del2_f[2][2]=f_zz(x,y,z);

		// intermediate matrix multiplication step: del2_Pk[row]:=del2_f[row][column]*p[row]:

		del2_Pk[0]	=	p[0]*del2_f[0][0] + p[1]*del2_f[0][1] + p[2]*del2_f[0][2];
		del2_Pk[1]	=	p[0]*del2_f[1][0] + p[1]*del2_f[1][1] + p[2]*del2_f[1][2];
		del2_Pk[2]	=	p[2]*del2_f[2][0] + p[2]*del2_f[2][1] + p[2]*del2_f[2][2];

		// so then p[column]*del2_Pk[row]:

		lambda_denominator =p[0]*del2_Pk[0]+p[1]*del2_Pk[1]+p[2]*del2_Pk[2];

		// OPTIONAL:
		if(lambda_denominator==0)
		{
			printf("\n !! ERROR !! \n\n\tlambda_denominator==0");
			LAMBDA =1; // break;
		}

		else
		{
			LAMBDA = lambda_numerator/lambda_denominator;
		}

		printf("\nlambda=%.3f\n",LAMBDA);

		// END calculating LAMBDA

		// update variables:
		x = x+LAMBDA*p[0];
		y = y+LAMBDA*p[1];
		z = z+LAMBDA*p[2];

		//
		g[0]= f_x(x,y,z);
		g[1]= f_y(x,y,z);
		g[2]= f_z(x,y,z);

		//
		printf("\nx= %.2f, \ny= %.2f, \nz= %.2f\n\n", x, y, z);

		// calculate the NORM_g

		g_SQRD=g[0]*g[0]+g[1]*g[1]+g[2]*g[2];//0;

		NORM_g = sqrt(g_SQRD);

		printf("NORM_g = %f\n",NORM_g);

		printf("\nEND LOOP\n");

		if(NORM_g<tolerance)
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
	*z_min = z;
	*f_min = f(x,y,z);

	printf("\nnum_step=%d", step);
}


double f(double x, double y, double z)
{
	return 1.5*x*x + 0.5*y*y + 1.5*z*z + 1.0*x*y - 11.0*x-y-9.0*z+26.0;
}

//First Derivatives

double f_x(double x, double y, double z)
{
	return 3.0*x+z-11.0;
}


double f_y(double x, double y, double z)
{
	return y - 1.0;
}


double f_z(double x, double y, double z)
{
	return 3.0*z + x - 9.0;
}


// Second Derivatives


double f_xx(double x, double y, double z)
{
	return 3.0; //can I do this? yes, I can.
}

double f_xy(double x, double y, double z)
{
	return 0;
}

double f_xz(double x, double y, double z)
{
	return 1.0;
}


double f_yx(double x, double y, double z)
{
	return 0;
}

double f_yy(double x, double y, double z)
{
	return 1.0;
}

double f_yz(double x, double y, double z)
{
	return 0;
}


double f_zx(double x, double y, double z)
{
	return 1.0;
}

double f_zy(double x, double y, double z)
{
	return 0;
}

double f_zz(double x, double y, double z)
{
	return 3.0;
}
