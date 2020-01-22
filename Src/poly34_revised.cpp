// poly34.cpp : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
// Thanks to Alexandr Rakhmanin <rakhmanin (at) gmail.com>
// public domain

// The following two functions were added by Sawada & Minkov
// SolveP3_revised: after SolveP3 of Khasin
// SolveP4_revised: after SolveP3 of Khasin
//
// (c) Group of "Mathematical and Computational Psychology"
// https://social.hse.ru/en/psy/mcp/
//
// Please read the LICENSE and NO WARRANTY statement in :
// MinkovSawada2020_License.txt
//
// Jun/16/2018: Retrieved from http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// Jun/17/2018: Revised by T. Sawada and V. Minkov
// Jan/22/2020: Uploaded to GitHub by T. Sawada


#include "pch.h" // MS Visual Studio C++



#include "pch.h"
#include <math.h>

#include "poly34.h"		   // solution of cubic and quartic equation
#define	TwoPi  6.28318530717958648
const double eps = 1e-14;


// x - array of size n
// For sorting roots x: x[0] <= x[1] <= ...
// Also used for a substitute of dblSort3
void poly34_BubbleSort_revised(double *x, int n)
{
	for (int i = 0; i<n - 1; i++)
		for (int j = n - 1; j>i; j--)
			if (x[j]<x[j - 1])
			{
				double t = x[j];
				x[j] = x[j - 1];
				x[j - 1] = t;
			}
}

//---------------------------------------------------------------------------
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ｱ i*x[2], return 1
int SolveP3_revised(double *x, double a, double b, double c) {	// solve cubic equation x^3 + a*x^2 + b*x + c = 0
	double a2 = a*a;
	double q = (a2 - 3 * b) / 9;
	double r = (a*(2 * a2 - 9 * b) + 27 * c) / 54;
	// equation x^3 + q*x + r = 0
	double r2 = r*r;
	double q3 = q*q*q;
	double A, B;
	//if (0<=q3 && r2 <= q3+eps) {//<<-- FIXED! (Jul/2018)
	if (r2 <= q3)
	{
		double t = r / sqrt(q3);
		if (t<-1) t = -1;
		if (t> 1) t = 1;
		t = acos(t);
		a /= 3; q = -2 * sqrt(q);
		x[0] = q*cos(t / 3) - a;
		x[1] = q*cos((t + TwoPi) / 3) - a;
		x[2] = q*cos((t - TwoPi) / 3) - a;
		return(3);
	}
	else {
		//A =-pow(fabs(r)+sqrt(r2-q3),1./3); 
		A = -root3(fabs(r) + sqrt(r2 - q3));
		if (r<0) A = -A;
		B = A == 0 ? 0 : B = q / A;

		a /= 3;
		x[0] = (A + B) - a;
		x[1] = -0.5*(A + B) - a;
		x[2] = 0.5*sqrt(3.)*(A - B);
		if (fabs(x[2])<eps) { x[2] = x[1]; return(2); }
		return(1);
	}
}// SolveP3(double *x,double a,double b,double c) {	
 //---------------------------------------------------------------------------

int   SolveP4De_revised(double *x, double b, double c, double d)	// solve equation x^4 + b*x^2 + c*x + d
{
	//if( c==0 ) return SolveP4Bi(x,b,d); // After that, c!=0
	if (fabs(c)<1e-14*(fabs(b) + fabs(d))) return SolveP4Bi(x, b, d); // After that, c!=0

	double xt3[3] = {};
	int res3 = SolveP3_revised(xt3, 2 * b, b*b - 4 * d, -c*c); // solve resolvent //(Wikipedia 1a)
													   // by Viet theorem:  x1*x2*x3=-c*c not equals to 0, so x1!=0, x2!=0, x3!=0
	if (res3>1)	// 2 or 3 real roots // 3 real roots
	{
		//dblSort3(x[0], x[1], x[2]);	// sort roots to x[0] <= x[1] <= x[2]
		poly34_BubbleSort_revised(xt3, 3); // xt3[0] <= xt3[1] <= xt3[2]

		// Note: xt3[0]*xt3[1]*xt3[2]= c*c > 0
		//(3a) All xt3[0-2] are >0 (standard)
		//(3b) Only xt3[2] is >0

		if (xt3[0] > 0) // (3a) All xt3[0-2] are >0 (standard)
		{
			double sz1 = sqrt(xt3[0]);
			double sz2 = sqrt(xt3[1]);
			double sz3 = sqrt(xt3[2]);
			if (c < 0) sz3 = -sz3;
			x[0] = (-sz1 - sz2 - sz3) / 2;
			x[1] = (-sz1 + sz2 + sz3) / 2;
			x[2] = (+sz1 - sz2 + sz3) / 2;
			x[3] = (+sz1 + sz2 - sz3) / 2;
			return 4;
		} // if( xt3[0] > 0) // (3a) All xt3[0-2] are >0 (standard) 
		else // if (xt3[2] > 0) // (3b) Only xt3[2] is >0 // Ferrari's solution
		{
			double a2 = -sqrt(xt3[2]);
			double b2_1 = b / 2 + xt3[2] / 2;
			double b2_2 = c / 2 / sqrt(xt3[2]);
			double xt2_1[2] = {};
			double xt2_2[2] = {};
			int res2_1 = SolveP2(xt2_1, +a2, b2_1 + b2_2); //
			int res2_2 = SolveP2(xt2_2, -a2, b2_1 - b2_2); //

			if (res2_1>0)
			{
				x[0] = xt2_1[0]; //1st real root
				x[1] = xt2_1[1]; //2nd real root
				x[2] = xt2_2[0]; //3rd real root or real part of pair of complex roots
				x[3] = xt2_2[1]; //4th real root or imaginary part of pair of complex roots
			}
			else if (res2_2>0)
			{
				x[0] = xt2_2[0]; //1st real root
				x[1] = xt2_2[1]; //2nd real root
				x[2] = xt2_1[0]; //3rd real root or real part of pair of complex roots
				x[3] = xt2_1[1]; //4th real root or imaginary part of pair of complex roots
			}
			else
			{
				x[0] = xt2_1[0]; //Real part of pair of complex roots
				x[1] = xt2_1[1]; //Imaginary part of pair of complex roots
				x[2] = xt2_2[0]; //Real part of pair of complex roots
				x[3] = xt2_2[1]; //Imaginary part of pair of complex roots
			}
			return res2_1 + res2_2; //0, 2 (0+2 or 2+0), or 4 (2+2)
		} // else // if (xt3[2] > 0) // (3b) Only xt3[2] is >0

	} // if( res3>1 ) // 2 or 3 real roots
	else // 1 real root
	{
		// now resoventa have 1 real and pair of compex roots
		// xt3[0] - real root, and xt3[0]>0, 
		// xt3[1] +- i*xt3[2] - complex roots, 
		// x[0] must be >=0. But one times x[0]=~ 1e-17, so:
		if (xt3[0] < 0)
			xt3[0] = 0;
		double sz1 = sqrt(xt3[0]);
		double szr, szi;
		CSqrt(xt3[1], xt3[2], szr, szi);  // (szr+i*szi)^2 = xt3[1]+i*xt3[2]
		if (c < 0) sz1 = -sz1;
		x[0] = -sz1 / 2 - szr;			// 1st real root
		x[1] = -sz1 / 2 + szr;			// 2nd real root
		x[2] = sz1 / 2;
		x[3] = szi;
		return 2;
	} // else // 1 real root
} // SolveP4De_revised(double *x, double b, double c, double d)	// solve equation x^4 + b*x^2 + c*x + d
  //-----------------------------------------------------------------------------
double N4Step_revised(double x, double a, double b, double c, double d)	// one Newton step for x^4 + a*x^3 + b*x^2 + c*x + d
{
	double fx = (((x + a)*x + b)*x + c)*x + d;	// f(x)
	double fxs = ((4 * x + 3 * a)*x + 2 * b)*x + c;	// f'(x)
	
	double nx = x - fx / fxs;
	double fnx = (((nx + a)*nx + b)*nx + c)*nx + d;	// f(nx)
	if (fabs(fx) < fabs(fnx))
		return x;
	else
		return nx;
}
//-----------------------------------------------------------------------------
// x - array of size 4
// return 4: 4 real roots x[0], x[1], x[2], x[3]
// return 2: 2 real roots x[0], x[1] and complex x[2] +- i*x[3], 
// return 0: two pair of complex roots: x[0] +- i*x[1],  x[2] +- i*x[3], 
int   SolveP4_revised(double *x, double a, double b, double c, double d) {	// solve equation x^4 + a*x^3 + b*x^2 + c*x + d by Dekart-Euler method
																	// move to a=0:
	double d1 = d + 0.25*a*(0.25*b*a - 3. / 64 * a*a*a - c);
	double c1 = c + 0.5*a*(0.25*a*a - b);
	double b1 = b - 0.375*a*a;
	int res = SolveP4De_revised(x, b1, c1, d1);

	if (res == 4)
	{
		x[0] -= a / 4;
		x[1] -= a / 4;
		x[2] -= a / 4;
		x[3] -= a / 4;
	}
	else if (res == 2)
	{
		x[0] -= a / 4;
		x[1] -= a / 4;
		x[2] -= a / 4;
	}
	else
	{
		x[0] -= a / 4;
		x[2] -= a / 4;
	}

	// one Newton step for each real root:
	if (res>0)
	{
		x[0] = N4Step_revised(x[0], a, b, c, d);
		x[1] = N4Step_revised(x[1], a, b, c, d);
	}
	if (res>2)
	{
		x[2] = N4Step_revised(x[2], a, b, c, d);
		x[3] = N4Step_revised(x[3], a, b, c, d);
	}

	poly34_BubbleSort_revised(x, res);

	return res;
}

