/*
 ============================================================================
 Name        : test_sort.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Test the sort function 
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void sort(double x1, double y1, double x2, double y2, double x3, double y3, double *Wx, double *Wy, double *f_W, double *Gx, double *Gy, double *f_G, double *Bx, double *By, double *f_B );

double f(double x, double y);




int main(void)
{

	double x1=0, y1=0, x2=1.2, y2=0, x3=0, y3=0.8;

	// inside nedlerMead{

	double Wx, Wy, Gx, Gy, Bx, By, f_W, f_G, f_B;

	sort( x1, y1,  x2, y2, x3, y3, &Wx, &Wy, &f_W, &Gx, &Gy, &f_G, &Bx, &By, &f_B );

	//printf("\nmax f_W= %f,\nmid f_G=%f,\nmin f_B=%f", f_W,f_G,f_B);
	printf("\nWx=%.2f,\nWy=%.2f,\nf_W=%.2f,\nGx=%.2f,\nGy=%.2f,\nf_G=%.2f,\nBx=%.2f,\nBy=%.2f,\nf_B=%.2f,\n", Wx, Wy, f_W, Gx, Gy, f_G, Bx, By, f_B);

	// TEST going in a circle
	//sort(Bx, By, Wx, Wy, Gx, Gy, &Wx, &Wy, &f_W, &Gx, &Gy, &f_G, &Bx, &By, &f_B );
	sort( x2, y2, x1, y1,  x3, y3, &Wx, &Wy, &f_W, &Gx, &Gy, &f_G, &Bx, &By, &f_B );

	printf("\nSECOND RUN:");


	printf("\nWx=%.2f,\nWy=%.2f,\nf_W=%.2f,\nGx=%.2f,\nGy=%.2f,\nf_G=%.2f,\nBx=%.2f,\nBy=%.2f,\nf_B=%.2f,\n", Wx, Wy, f_W, Gx, Gy, f_G, Bx, By, f_B);
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



/////////////


// define f(x):


double f(double x, double y)
{
	return x*x - 4*x + y*y - y - x*y + 8;
}


//////////
