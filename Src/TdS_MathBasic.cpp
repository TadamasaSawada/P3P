// (c) Group of "Mathematical and Computational Psychology"
// https://social.hse.ru/en/psy/mcp/
//
// Please read the LICENSE and NO WARRANTY statement in :
// MinkovSawada2020_License.txt
//
// Jan/22/2020: Uploaded to GitHub by T. Sawada


#include "pch.h" // MS Visual Studio C++
#include <math.h>
#include <cmath>
#include "TdS_Math.h"
#include "poly34.h"     // solution of cubic and quartic equation

const double eps = 1e-14;


double TdS_RandomNumber(double start, double end)
{
	return (double)rand()*(end - start) / RAND_MAX + start;
}

void TdS_RandomNumbers(double start, double end, int n_v, double *v)
{
	for (int i = 0; i < n_v; i++)
		v[i] = TdS_RandomNumber(start, end);
}

void TdS_BubbleSort(double *x, int n) // Sort: x[0] <= x[1] <= ...
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

void TdS_BubbleSort(double *x, int *index, int n) // Sort: x[0] <= x[1] <= ...
{
	for (int i = 0; i < n; i++)
		index[i] = i;

	for (int i = 0; i < n - 1; i++)
		for (int j = n - 1; j > i; j--)
			if (x[index[j]] < x[index[j - 1]])
			{
				int t = index[j];
				index[j] = index[j - 1];
				index[j - 1] = t;
			}
}

int	TdS_Min(int num_v, int *v, int *i_min)
{
	int temp_min = v[0];
	int temp_i_min = 0;

	for (int i = 1; i < num_v; i++)
		if (v[i] < temp_min)
		{
			temp_min = v[i];
			temp_i_min = i;
		}
	if (i_min != NULL) *i_min = temp_i_min;
	return temp_min;
}

double	TdS_Min(int num_v, double *v, int *i_min)
{
	double temp_min = v[0];
	int temp_i_min = 0;

	for (int i = 1; i < num_v; i++)
		if (v[i] < temp_min)
		{
			temp_min = v[i];
			temp_i_min = i;
		}
	if (i_min != NULL) *i_min = temp_i_min;
	return temp_min;
}

int	TdS_Max(int num_v, int *v, int *i_max)
{
	int temp_max = v[0];
	int temp_i_max = 0;

	for (int i = 1; i < num_v; i++)
		if (v[i] > temp_max)
		{
			temp_max = v[i];
			temp_i_max = i;
		}
	if (i_max != NULL) *i_max = temp_i_max;
	return temp_max;
}

double	TdS_Max(int num_v, double *v, int *i_max)
{
	double temp_max = v[0];
	int temp_i_max = 0;

	for (int i = 1; i < num_v; i++)
		if (v[i] > temp_max)
		{
			temp_max = v[i];
			temp_i_max = i;
		}
	if (i_max != NULL) *i_max = temp_i_max;
	return temp_max;
}

int TdS_QuatricEquation(double *x, double coef4, double coef3, double coef2, double coef1, double coef0)
{
	int res = -1;
	if (fabs(coef4) < 1e-10) // 3rd-polynomial or less
	{
		if (fabs(coef3) < 1e-10) // 2nd-polynomial or less
		{
			if (fabs(coef2) < 1e-10) // 1st-polynomial or less
			{
				if (fabs(coef1) < 1e-10)
					return -1; // error

							   // 1st-polynomial
				x[0] = -coef0 / coef1;
				return 1;
			}

			// 2nd-polynomial
			double b2 = coef1 / coef2;
			double c2 = coef0 / coef2;
			double determinant = b2*b2 - 4 * c2;
			if (fabs(determinant) < 1e-10)
			{
				x[0] = -b2 / 2;
				return 1;
			}
			else if (determinant > 0)
			{
				x[0] = (-b2 - sqrt(determinant)) / 2;
				x[1] = (-b2 + sqrt(determinant)) / 2;
				return 2;
			}
			else
				return 0;
		}

		// 3rd-polynomial
		res = SolveP3_revised(x, coef2 / coef3, coef1 / coef3, coef0 / coef3);
		TdS_BubbleSort(x, res);
	}
	else
		// 4th-polynomial
		res = SolveP4_revised(x, coef3 / coef4, coef2 / coef4, coef1 / coef4, coef0 / coef4);

	// Remove identical roots of 3rd- and 4th-polynomials
	int i, j;
	if (res > 1)
	{
		for (i = 1, j = 1; i < res; i++)
			if (fabs(x[i - 1] - x[i]) > eps) x[j++] = x[i];
		res = j;
	}

	return res;
}

int TdS_SimultaneousEquations3(double *x, double a1, double a2, double a3, double a0, double b1, double b2, double b3, double b0, double c1, double c2, double c3, double c0)
{
	//After http://cplusplus.happycodings.com/mathematics/code5.html

	//a1*x + a2*y + a3*z + a0 = 0
	//b1*x + b2*y + b3*z + b0 = 0
	//c1*x + c2*y + c3*z + c0 = 0

	double denominator = (a1*b2*c3 + a2*c1*b3 + a3*b1*c2) - (a1*b3*c2 + a2*b1*c3 + a3*b2*c1);
	if (fabs(denominator) < eps)
		return -1; // error

	x[0] = ((a2*c3*b0 + a3*b2*c0 + a0*b3*c2) - (a2*b3*c0 + a3*c2*b0 + a0*b2*c3)) / denominator;
	x[1] = ((a1*b3*c0 + a3*c1*b0 + a0*b1*c3) - (a1*c3*b0 + a3*b1*c0 + a0*b3*c1)) / denominator;
	x[2] = ((a1*c2*b0 + a2*b1*c0 + a0*b2*c1) - (a1*b2*c0 + a2*c1*b0 + a0*b1*c2)) / denominator;

	return 0; // fine
}