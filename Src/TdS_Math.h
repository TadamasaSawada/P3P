// (c) Group of "Mathematical and Computational Psychology"
// https://social.hse.ru/en/psy/mcp/
//
// Please read the LICENSE and NO WARRANTY statement in :
// MinkovSawada2020_License.txt
//
// Jan/22/2020: Uploaded to GitHub by T. Sawada


#pragma once

#define	TwoPi   6.28318530717958648
#define	Pi		3.14159265358979323
#define PiTwo   6.28318530717958648
#define PiHalf	1.57079632679489662


/*Basic*/
double	TdS_RandomNumber(double start, double end);
void	TdS_RandomNumbers(double start, double end, int n_v, double *v);
int		TdS_Min(int num_v, int *v, int *i_min);
double	TdS_Min(int num_v, double *v, int *i_min);
int		TdS_Max(int num_v, int *v, int *i_max);
double	TdS_Max(int num_v, double *v, int *i_max);
void	TdS_BubbleSort(double *x, int n);
void	TdS_BubbleSort(double *x, int *index, int n);
int		TdS_QuatricEquation(double *x, double coef4, double coef3, double coef2, double coef1, double coef0);
int		TdS_SimultaneousEquations3(double *x, double a1, double a2, double a3, double a0, double b1, double b2, double b3, double b0, double c1, double c2, double c3, double c0);

/*Geometry 2D*/
void	TdS_Vector2D_Normalize(double *x, double *y);
void	TdS_Vector2D_Normalize(double *v);
double	TdS_ACosSin(double rCos, double rSin, bool deg1_rad0);


/*Geometry 3D*/
void	TdS_Vector3D_Plus(double *vx, double *va, double *vb);
void	TdS_Vector3D_Minus(double *vx, double *va, double *vb, bool bma1_amb0);
void	TdS_Vector3D_Scale(double *vx, double scale);
double	TdS_Vector3D_Length(double *vx);
void	TdS_Vector3D_Normalize(double *vx);
double	TdS_Vector3D_DotProduct(double *va, double *vb);
void	TdS_Vector3D_CrossProduct(double *vx, double *va, double *vb);
double	TdS_Vector3D_Angle(double *va, double *vb, double *vn, bool deg1_rad0);

double	TdS_Point3D_Distance(double *pa, double *pb);
double	TdS_Point3D_DistanceManhattan(double *pa, double *pb);
double	TdS_Point3D_Angle(double *po, double *pa, double *pb, double *vn, bool deg1_rad0);

double	TdS_Line3D_Intersection(double *pxa, double *pxb, double *va, double *pa, double *vb, double *pb);
double TdS_Line3D_Intersection(double *pxa, double *pxb, double *lva, double *lvb, double *va, double *pa, double *vb, double *pb);

  //3P3 problem
int		TdS_Vector3D_From3Vectors3Angles(double *vx, double *va, double *vb, double *vc, double angA, double angB, double angC);
int		TdS_Tetrahedron3D_CartesianCoordinates(double *eye, double *pA, double *pB, double *pC, double angA, double angB, double angC, double angAB, double angBC, double angCA, double lAB, double lEA, double lEB, double lEC);
