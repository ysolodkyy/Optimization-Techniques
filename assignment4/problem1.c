#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define TOLERANCE 0.00001

void sort(double x1, double y1, double x2, double y2, double x3, double y3, double *Wx, double *Wy, double *f_W, double *Gx, double *Gy, double *f_G, double *Bx, double *By, double *f_B );

double f(double x, double y);

void nelderMead(double x1, double y1, double x2, double y2, double x3, double y3, double tol, double *x_min, double *y_min, double *f_min);


int main(void)
{

	double x1=0, y1=0, x2=1.2, y2=0, x3=0, y3=0.8;



	double x_min, y_min, f_min;

	nelderMead(x1, y1, x2, y2, x3, y3, TOLERANCE, &x_min, &y_min, &f_min);


	printf("\nx_min = %f, y_min = %f, f_min = %f.",x_min, y_min, f_min);

	return EXIT_SUCCESS;
}


// define the sort function

void sort(double x1, double y1, double x2, double y2, double x3, double y3, double *Wx, double *Wy, double *f_W, double *Gx, double *Gy, double *f_G, double *Bx, double *By, double *f_B )
{

	// this function is super ugly, but I  want to see if the structure works

	double f_array1[3], f_array2[3], f_array3[3];

	double B[2], G[2], W[2]; // function coordinate arrays/vectors.

	//double f_B, f_G, f_W; //f_B == lowest value of f, f_G middle value of f, f_W highest value of f;

	f_array1[0]=x1;
	f_array1[1]=y1;

	f_array1[2]=f(x1,y1);

	f_array2[0]=x2;
	f_array2[1]=y2;
	f_array2[2]=f(x2,y2);

	f_array3[0]=x3;
	f_array3[1]=y3;
	f_array3[2]=f(x3,y3);

	// BEGIN SORTING :

	if(f_array1[2]>f_array2[2] && f_array2[2]>f_array3[2])
	{
		*f_W=f_array1[2];
		W[0]=f_array1[0];
		W[1]=f_array1[1];

		*f_G=f_array2[2];
		G[0]=f_array2[0];
		G[1]=f_array2[1];

		*f_B=f_array3[2];
		B[0]=f_array3[0];
		B[1]=f_array3[1];
	}
	else if(f_array1[2]>f_array3[2] && f_array3[2]>f_array2[2])
	{
		*f_W=f_array1[2];
		W[0]=f_array1[0];
		W[1]=f_array1[1];

		*f_G=f_array3[2];
		G[0]=f_array3[0];
		G[1]=f_array3[1];

		*f_B=f_array2[2];
		B[0]=f_array2[0];
		B[1]=f_array2[1];
	}

	else if(f_array2[2]>f_array1[2] && f_array1[2]>f_array3[2])
	{
		*f_W=f_array2[2];
		W[0]=f_array2[0];
		W[1]=f_array2[1];

		*f_G=f_array1[2];
		G[0]=f_array1[0];
		G[1]=f_array1[1];

		*f_B=f_array3[2];
		B[0]=f_array3[0];
		B[1]=f_array3[1];
	}
	else if(f_array2[2]>f_array3[2] && f_array3[2]>f_array1[2])
	{
		*f_W=f_array2[2];
		W[0]=f_array2[0];
		W[1]=f_array2[1];

		*f_G=f_array3[2];
		G[0]=f_array3[0];
		G[1]=f_array3[1];

		*f_B=f_array1[2];
		B[0]=f_array1[0];
		B[1]=f_array1[1];
	}
	else if(f_array3[2]>f_array1[2] && f_array1[2]>f_array2[2])
	{
		*f_W=f_array3[2];
		W[0]=f_array3[0];
		W[1]=f_array3[1];


		*f_G=f_array1[2];
		G[0]=f_array1[0];
		G[1]=f_array1[1];

		*f_B=f_array2[2];
		B[0]=f_array2[0];
		B[1]=f_array2[1];
	}
	else if(f_array3[2]>f_array2[2] && f_array2[2]>f_array1[2])
	{
		*f_W=f_array3[2];
		W[0]=f_array3[0];
		W[1]=f_array3[1];

		*f_G=f_array2[2];
		G[0]=f_array2[0];
		G[1]=f_array2[1];

		*f_B=f_array1[2];
		B[0]=f_array1[0];
		B[1]=f_array1[1];
	}

	// END SORTING

	*Wx=W[0];
	*Wy=W[1];

	*Gx=G[0];
	*Gy=G[1];

	*Bx=B[0];
	*By=B[1];

}


// define f(x):


double f(double x, double y)
{
	return x*x - 4.0*x + y*y - y - x*y + 8.0;
}


// the beast:


void nelderMead(double x1, double y1, double x2, double y2, double x3, double y3, double tol, double *x_min, double *y_min, double *f_min)
{
	int count=0;

	double Wx, Wy, Gx, Gy, Bx, By, f_W, f_G, f_B;//f_B == lowest value of f, f_G middle value of f, f_W highest value of f;
	sort( x1, y1,  x2, y2, x3, y3, &Wx, &Wy, &f_W, &Gx, &Gy, &f_G, &Bx, &By, &f_B );


	printf("\nf_W=%f,f_G=%f,f_B=%f\n", f_W,f_G,f_B); // this checks out


	double f_bar;


	f_bar = (f_W+f_G+f_B)/3.0; // there is some kind of issue here.

	printf("\nf_bar=%f \n",f_bar);


	double ERR;

	ERR = sqrt(((f_W-f_bar)*(f_W-f_bar)+(f_G-f_bar)*(f_G-f_bar)+(f_B-f_bar)*(f_B-f_bar))/3.0);


	double Mx, My, Rx, Ry, Ex, Ey, Cix, Ciy, Cox, Coy, f_Ci, f_Co, f_C, Cx, Cy, Sx, Sy;
	double f_R,f_E;// R[] & E[] are arrays to represent points
	// start using Mx, My, Rx, Ry, Ex, Ey, Cix, Ciy, Cox, Coy;


	printf("\nerror=%f, before looping \n", ERR);


	while(ERR>tol)
	{
		printf("\nI"); // I
		Mx=(Bx+Gx)/2;
		My=(By+Gy)/2;

		Rx = 2*Mx-Wx;
		Ry = 2*My-Wy;

		f_R = f(Rx,Ry);

		if(f_B<f_R && f_R<f_G) // II
		{
			printf("\nII");
			Wx=Gx;
			Wy=Gy;
			f_W=f_G;

			Gx=Rx;
			Gy=Ry;
			f_G=f_R;
		}
		else		// III
		{
			printf("\nIII");

			if(f_R<f_B) // IV.
			{
				printf("\nIV");
				Ex=2*Rx-Mx;
				Ey=2*Ry-My;

				f_E = f(Ex,Ey);

				if(f_E<f_R) // V
				{
					printf("\nV");
					Wx=Gx;
					Wy=Gy;
					f_W=f_G;

					Gx=Bx;
					Gy=By;
					f_G=f_B;

					Bx=Ex;
					By=Ey;
					f_B=f_E;
				}
				else // VI
				{
					printf("\nVI");
					Wx=Gx;
					Wy=Gy;
					f_W=f_G;

					Gx=Bx;
					Gy=By;
					f_G=f_B;

					Bx=Rx;
					By=Ry;
					f_B=f_R;
				}
			}

			else // VII
			{
				printf("\nVII");
				Cix = (Mx+Wx) /2;
				Ciy = (My+Wy) /2;
				f_Ci = f(Cix,Ciy);

				Cox = (Mx+Rx) /2;
				Coy = (My+Ry) /2;
				f_Co = f(Cox,Coy);

				if(f_Ci<f_Co) // IIX
				{
					printf("\nIIX");
					Cx = Cix;
					Cy = Ciy;
					f_C = f_Ci;
				}
				else // IX
				{
					printf("\nIX");
					Cx = Cox;
					Cy = Coy;
					f_C = f_Co;
				}
				if(f_C<f_G) // X
				{
						// something is wrong in this part of the code
					printf("\nX");
					sort(Bx, By, Gx, Gy, Cx, Cy, &Wx, &Wy, &f_W, &Gx, &Gy, &f_G, &Bx, &By, &f_B);
				}
				else // XI
				{
					printf("\nXI");
					Sx = (Wx+Bx)/2;
					Sy = (Wy+By)/2;

					sort(Bx, By, Sx, Sy, Mx, My, &Wx, &Wy, &f_W, &Gx, &Gy, &f_G, &Bx, &By, &f_B );

				}
			}

		}

		// troubleshooting -- making it go through 19 loops


		printf("\nBx=%f,By=%f,f_B=%f\n", Bx,By,f_B);
		printf("\nGx=%f,Gy=%f,f_G=%f\n", Gx,Gy,f_G);
		printf("\nWx=%f,Wy=%f,f_W=%f\n", Wx,Wy,f_W);

		f_bar = (f_W+f_G+f_B)/3.0; // re-evaluate f_bar

		ERR = sqrt(((f_W-f_bar)*(f_W-f_bar)+(f_G-f_bar)*(f_G-f_bar)+(f_B-f_bar)*(f_B-f_bar))/3.0);

		printf("\nerror=%f\n", ERR);
		printf("\n");




		++count;

		// end troubleshooting


	}
		printf("\ncount=%d", count);
		// the end:
		*y_min=By;
		*x_min=Bx;
		*f_min=f_B;

}
