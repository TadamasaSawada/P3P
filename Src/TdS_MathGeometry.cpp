// (c) Group of "Mathematical and Computational Psychology"
// https://social.hse.ru/en/psy/mcp/
//
// Please read the LICENSE and NO WARRANTY statement in :
// MinkovSawada2020_License.txt
//
// Jan/22/2020: Uploaded to GitHub by T. Sawada


#include "pch.h" // MS Visual Studio C++
#include <math.h>
//#include <cmath>
#include "TdS_Math.h"

const double eps = 1e-14;

/*2D*/
void TdS_Vector2D_Normalize(double *x, double *y)
{
	double l = sqrt( (*x)*(*x) + (*y)*(*y) );
	*x /= l;
	*y /= l;
}

void TdS_Vector2D_Normalize(double *v)
{
	double l = sqrt((v[0])*(v[0]) + (v[1])*(v[1]) );
	v[0] /= l;
	v[1] /= l;
}

double TdS_ACosSin(double rCos, double rSin, bool deg1_rad0)
{
	double radius, radAsin, radAcos, radResult;
	radius = sqrt(rCos*rCos + rSin*rSin);
	if (radius < eps) return 0;
	if (rSin>radius) rSin = radius; if (rSin<-radius) rSin = -radius;
	if (rCos>radius) rCos = radius; if (rCos<-radius) rCos = -radius;
	radAsin = asin(rSin / radius);
	radAcos = acos(rCos / radius);
	if (radAsin>0)	radResult = (radAcos>Pi / 2) ? (Pi - radAsin) : (radAsin);
	else			radResult = (radAcos>Pi / 2) ? (Pi - radAsin) : (radAsin + PiTwo);

	if(deg1_rad0)
		return radResult/Pi*180.0;
	else
		return radResult;
}

/*3D*/
void TdS_Vector3D_Plus(double *vx, double *va, double *vb)
{
	vx[0] = vb[0] + va[0];
	vx[1] = vb[1] + va[1];
	vx[2] = vb[2] + va[2];
}

void TdS_Vector3D_Minus(double *vx, double *va, double *vb, bool bma1_amb0)
{
	if (bma1_amb0) // b-a
	{
		vx[0] = vb[0] - va[0];
		vx[1] = vb[1] - va[1];
		vx[2] = vb[2] - va[2];
	}
	else // a-b
	{
		vx[0] = va[0] - vb[0];
		vx[1] = va[1] - vb[1];
		vx[2] = va[2] - vb[2];
	}
}

void	TdS_Vector3D_Scale(double *vx, double scale)
{
	vx[0] *= scale;
	vx[1] *= scale;
	vx[2] *= scale;
}

double TdS_Vector3D_Length(double *v)
{
	return sqrt((v[0])*(v[0]) + (v[1])*(v[1]) + (v[2])*(v[2]));
}

void TdS_Vector3D_Normalize(double *vx)
{
	double l = TdS_Vector3D_Length(vx);
	if(l>eps)
	{
		vx[0] /= l;
		vx[1] /= l;
		vx[2] /= l;
	}
}

double TdS_Vector3D_DotProduct(double *va, double *vb)
{
	return va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
}

void TdS_Vector3D_CrossProduct(double *vx, double *va, double *vb)
{
	vx[0] = va[1]*vb[2]-va[2]*vb[1];
	vx[1] = va[2]*vb[0]-va[0]*vb[2];
	vx[2] = va[0]*vb[1]-va[1]*vb[0];
}

double TdS_Vector3D_Angle(double *va, double *vb, double *vn, bool deg1_rad0)
{	
	double dotProduct = TdS_Vector3D_DotProduct(va, vb);
	double la = TdS_Vector3D_Length(va);
	double lb = TdS_Vector3D_Length(vb);
	if (la < eps || lb < eps)
		return -1;

	if (vn != NULL)
	{
		TdS_Vector3D_CrossProduct(vn, va, vb);
		TdS_Vector3D_Normalize(vn);
	}

	double cosTheta = dotProduct / (la*lb);
	if (cosTheta> 1) cosTheta =  1;
	if (cosTheta<-1) cosTheta = -1;

	if (deg1_rad0)
		return acos(cosTheta)/Pi*180.0;
	else
		return acos(cosTheta);
}

double TdS_Point3D_Distance(double *pa, double *pb)
{
	double ppa[3] = { pa[0],pa[1],pa[2] };
	double ppb[3] = { pb[0],pb[1],pb[2] };

	return sqrt( pow(pa[0]-pb[0],2) + pow(pa[1]-pb[1],2) + pow(pa[2]-pb[2],2) );
}

double	TdS_Point3D_DistanceManhattan(double *pa, double *pb)
{
	return fabs(pa[0]-pb[0]) + fabs(pa[1]-pb[1]) + fabs(pa[2]-pb[2]);
}

double TdS_Point3D_Angle(double *po, double *pa, double *pb, double *vn, bool deg1_rad0)
{
	double va[3] = {};
	double vb[3] = {};
	TdS_Vector3D_Minus(va, po, pa, 1);
	TdS_Vector3D_Minus(vb, po, pb, 1);
	return TdS_Vector3D_Angle(va, vb, vn, deg1_rad0);
}

double TdS_Line3D_Intersection(double *pxa, double *pxb, double *va, double *pa, double *vb, double *pb)
{
	//After http://www.sousakuba.com/Programming/gs_two_lines_intersect.html
	
	double van[3] = { va[0],va[1],va[2] };
	double vbn[3] = { vb[0],vb[1],vb[2] };
	TdS_Vector3D_Normalize(van);
	TdS_Vector3D_Normalize(vbn);
	if (TdS_Point3D_Distance(van, vbn) < eps)
		return -1; // error: parallel

	double work1 = TdS_Vector3D_DotProduct(van, vbn);
	double work2 = 1 - work1*work1;

	if (work2 < eps)
		return -1; // error: parallel

	double v_papb[3] = {};
	TdS_Vector3D_Minus(v_papb, pa, pb, 1); // v_papb = pb - pa

	double d1 = (TdS_Vector3D_DotProduct(v_papb, van) - work1 * TdS_Vector3D_DotProduct(v_papb, vbn)) / work2;
	double d2 = (work1 * TdS_Vector3D_DotProduct(v_papb, van) - TdS_Vector3D_DotProduct(v_papb, vbn)) / work2;

	pxa[0] = pa[0] + d1 * van[0];
	pxa[1] = pa[1] + d1 * van[1];
	pxa[2] = pa[2] + d1 * van[2];

	pxb[0] = pb[0] + d2 * vbn[0];
	pxb[1] = pb[1] + d2 * vbn[1];
	pxb[2] = pb[2] + d2 * vbn[2];

	return TdS_Point3D_Distance(pxa, pxb);
}

double TdS_Line3D_Intersection(double *pxa, double *pxb, double *lva, double *lvb, double *va, double *pa, double *vb, double *pb)
{
	//After http://www.sousakuba.com/Programming/gs_two_lines_intersect.html

	double pan[3] = { pa[0],pa[1],pa[2] };
	double pbn[3] = { pb[0],pb[1],pb[2] };



	double van[3] = { va[0],va[1],va[2] };
	double vbn[3] = { vb[0],vb[1],vb[2] };
	TdS_Vector3D_Normalize(van);
	TdS_Vector3D_Normalize(vbn);
	if (TdS_Point3D_Distance(van, vbn) < eps)
		return -1; // error: parallel

	double work1 = TdS_Vector3D_DotProduct(van, vbn);
	double work2 = 1 - work1 * work1;

	if (work2 < eps)
		return -1; // error: parallel

	double v_papb[3] = {};
	TdS_Vector3D_Minus(v_papb, pa, pb, 1); // v_papb = pb - pa

	*lva = (TdS_Vector3D_DotProduct(v_papb, van) - work1 * TdS_Vector3D_DotProduct(v_papb, vbn)) / work2;
	*lvb = (work1 * TdS_Vector3D_DotProduct(v_papb, van) - TdS_Vector3D_DotProduct(v_papb, vbn)) / work2;

	pxa[0] = pa[0] + (*lva) * van[0];
	pxa[1] = pa[1] + (*lva) * van[1];
	pxa[2] = pa[2] + (*lva) * van[2];

	pxb[0] = pb[0] + (*lvb) * vbn[0];
	pxb[1] = pb[1] + (*lvb) * vbn[1];
	pxb[2] = pb[2] + (*lvb) * vbn[2];

	return TdS_Point3D_Distance(pxa, pxb);
}

int TdS_Vector3D_From3Vectors3Angles(double *vx, double *va, double *vb, double *vc, double angA, double angB, double angC)
{
	double a0 = cos(angA*Pi / 180.0);
	double b0 = cos(angB*Pi / 180.0);
	double c0 = cos(angC*Pi / 180.0);
	int err = TdS_SimultaneousEquations3(vx, va[0],va[1],va[2],a0, vb[0],vb[1],vb[2],b0, vc[0],vc[1],vc[2],c0);
	return err;
}

int TdS_Tetrahedron3D_CartesianCoordinates(double *eye, double *pA, double *pB, double *pC, double angA, double angB, double angC, double angAB, double angBC, double angCA, double lAB, double lEA, double lEB, double lEC)
{
	char x = 0, y = 1, z = 2;

	if ( fabs(angA+angB+angC - 180) > eps*1000 )
		return -1;

	angA  *= Pi / 180.0;
	angB  *= Pi / 180.0;

	//Set bottom triangle on the XY-plane
	pA[x] = 0;   pA[y] = 0; pA[z] = 0;
	pB[x] = lAB; pB[y] = 0; pB[z] = 0;

	pC[x] = lAB * tan(angB) / (tan(angA)+tan(angB));
	pC[y] = fabs(pC[x] * tan(angA));
	pC[z] = 0;

	//Compute eye position
	double lEA2 = lEA*lEA;
	double lEB2 = lEB*lEB;
	double lEC2 = lEC*lEC;
	eye[x] = (lEA2 - lEB2 + pB[x]*pB[x]) / (2.0*pB[x]);
	eye[y] = (lEA2 - lEC2 + pC[x]*pC[x] + pC[y]*pC[y] -2*eye[x]*pC[x]) / (2.0*pC[y]);
	eye[z] = sqrt(lEA2 - eye[x]*eye[x] - eye[y]*eye[y]);
	
	if (lEA2 < eye[x] * eye[x] + eye[y] * eye[y])
		return -2;

	//Check
	double test_angAB = TdS_Point3D_Angle(eye, pA, pB, NULL, 1);
	double test_angBC = TdS_Point3D_Angle(eye, pB, pC, NULL, 1);
	double test_angCA = TdS_Point3D_Angle(eye, pC, pA, NULL, 1);
	double test_angDiff = fabs(test_angAB-angAB) + fabs(test_angBC-angBC) + fabs(test_angCA-angCA);

	if ( test_angDiff > 1)
		return -3;
	else
		return 1;

}