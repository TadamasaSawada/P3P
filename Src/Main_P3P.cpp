// (c) Group of "Mathematical and Computational Psychology"
// https://social.hse.ru/en/psy/mcp/
//
// Please read the LICENSE and NO WARRANTY statement in :
// MinkovSawada2020_License.txt
//
// Jan/22/2020: Uploaded to GitHub by T. Sawada


#include "pch.h" // MS Visual Studio C++
#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <ctime>
#include <random>

#include <windows.h> // time
#pragma comment(lib, "Kernel32.lib") // time

#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/calib3d.hpp>


//#if NDEBUG 
//#pragma comment(lib, "opencv_world411.lib")
//#else
//#pragma comment(lib, "opencv_world411d.lib")
//#endif


#include "TdS_Math.h"
#include "P3P_Banno.h"
#include "P3P_MinkovSawada.h"


using namespace std;
using namespace Eigen;

#define pi 3.14159265358979323846264338327950288
#define EPS (1e-14)
#define TriangleNp 3
#define SAVE_SIMULATION true
#define NUM_TRIALS  1000000
#define NUM_SESSIONS 1

#define ERR_THRESHOLD0 0.01 // Compared with the original scene
#define ERR_THRESHOLD1 0.01 // Error of visual angles from the original retinal image
#define ERR_THRESHOLD2 0.01 // Duplicate

#define MAX_ECCENTRICITY 45 // 45 or 85 deg

int ix = 0, iy = 1, iz = 2;

double subfunction_length(double x, double y, double z);
void subfunction_norm_interpretations(int num_interpretations, double(&interpretations)[10][3]);
void subfunction_sort_interpretations(int num_interpretations, double(&interpretations)[10][3]);
double subfunction_compare_angles(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double ang12);
double subfunction_compare_interpretations(int ia, int ib, double(&interpretationsA)[10][3], double(&interpretationsB)[10][3]);

int main()
{
	srand((unsigned)time(NULL));

	unsigned int seed = (unsigned)time(NULL);
	srand(seed);

	cout << "seed = " << seed << endl << endl;


	int numPoints = TriangleNp;
	int minZ = 1;

	double points_x[TriangleNp], points_y[TriangleNp], points_z[TriangleNp], points_d[TriangleNp];
	double point0[3], point1[3], point2[3];
	double eye[3] = { 0,0,0 };
	double triangle_rad[3], triangle_deg[3];
	double triangle_base, triangle_apex[2];
	double results_Original[10][3]; // [10] just for consistency
	double retina_angles_deg[3];
	double retina_xy[TriangleNp][3];//z = 1
	ofstream fout;
	char fout_name[100];
	
	Matrix3d matrix_retina; // columns = points, rows = xyz
	Matrix3d matrix_shape; //  columns = points, rows = xyz
	Matrix3d matrix_results;
	std::vector<Matrix3d> matrix_Rot;
	std::vector<Vector3d> vector_Tran;
	Matrix3d matrix_Tran;
	CGroebner P3P_Banno2018;

	int ocv_flags;
	vector<cv::Mat> ocv_rvecs, ocv_tvecs;
	cv::Mat cameraMatrix(3, 3, cv::DataType<double>::type);
	cv::setIdentity(cameraMatrix);
	cv::Mat distCoeffs(4, 1, cv::DataType<double>::type);
	distCoeffs.at<double>(0) = 0;
	distCoeffs.at<double>(1) = 0;
	distCoeffs.at<double>(2) = 0;
	distCoeffs.at<double>(3) = 0;

	int num_MinkovSawada;
	int num_Banno2018;
	int num_OCV_GaoEtal2003;
	int num_OCV_KeRoumeliotis2017;

	double results_MinkovSawada[10][3];
	double results_Banno2018[10][3];
	double results_CV_GaoEtal2003[10][3];
	double results_CV_KeRoumeliotis2017[10][3];

	double err1_Banno2018[10];
	double err1_CV_GaoEtal2003[10];
	double err1_CV_KeRoumeliotis2017[10];
	double err2_Banno2018[10];
	double err2_CV_GaoEtal2003[10];
	double err2_CV_KeRoumeliotis2017[10];
	bool err3_Banno2018[10];
	bool err3_CV_GaoEtal2003[10];
	bool err3_CV_KeRoumeliotis2017[10];

	double time_MinkovSawada = 0; //msec
	double time_Banno2018 = 0; //msec
	double time_GaoEtal2003 = 0; //msec
	double time_KeRoumeliotis2017 = 0; //msec

	int summary_MinkovSawada[8] = {};		//0:miss, 1-4:false-alarm, 5-7:less/equal/more than MinkovSawada
	int summary_Banno2018[8] = {};			//0:miss, 1-4:false-alarm, 5-7:less/equal/more than MinkovSawada
	int summary_GaoEtal2003[8] = {};		//0:miss, 1-4:false-alarm, 5-7:less/equal/more than MinkovSawada
	int summary_KeRoumeliotis2017[8] = {};	//0:miss, 1-4:false-alarm, 5-7:less/equal/more than MinkovSawada

	LARGE_INTEGER start_time, end_time, ticks_time;
	QueryPerformanceFrequency(&ticks_time);

	if (SAVE_SIMULATION)
	{
		sprintf_s(fout_name, 100, "Output_P3P.txt");
		fout.open(fout_name, ios::out | ios::app);
		fout << endl;
		fout << "Seed = " << seed << endl;
		fout << "#Repeat = " << NUM_TRIALS << endl;
		fout << "Max eccentricity = " << MAX_ECCENTRICITY << endl;
		fout << "Visual angle threshold = " << ERR_THRESHOLD1 << endl;
		fout << "Threshold for comparison with the original = " << ERR_THRESHOLD0 << endl;
		fout << "Threshold for detecting duplicate = " << ERR_THRESHOLD2 << endl;
		fout.close();
	}

	int session = 0;

start_session:

	time_MinkovSawada = 0;
	time_Banno2018 = 0;
	time_GaoEtal2003 = 0;
	time_KeRoumeliotis2017 = 0;

	for (int i = 0; i < 8; i++) //0:miss, 1-4:false-alarm, 5-7:less/equal/more than MinkovSawada
	{
		summary_MinkovSawada[i] = 0;
		summary_Banno2018[i] = 0;
		summary_GaoEtal2003[i] = 0;
		summary_KeRoumeliotis2017[i] = 0;
	}


	for (int trial = 0; trial < NUM_TRIALS; trial++)
	{

		/*Generate a random triangle*/
		double imgvar = 0;
		for (int p = 0; p < numPoints; p++)
		{
			points_z[p] = TdS_RandomNumber(minZ, minZ * 100);

			double eccentricity = TdS_RandomNumber(0, MAX_ECCENTRICITY) * pi / 180.0;
			double ori_xy = TdS_RandomNumber(0, 360) * pi / 180.0;
			points_x[p] = points_z[p] * tan(eccentricity) * cos(ori_xy);
			points_y[p] = points_z[p] * tan(eccentricity) * sin(ori_xy);

			points_d[p] = subfunction_length(points_x[p], points_y[p], points_z[p]);

			retina_xy[p][ix] = points_x[p] / points_z[p];
			retina_xy[p][iy] = points_y[p] / points_z[p];
			retina_xy[p][iz] = 1;

			results_Original[0][p] = points_d[p];
		}

		subfunction_norm_interpretations(1, results_Original);


	/*Minkov & Sawada*/
		point0[ix] = points_x[0]; point0[iy] = points_y[0]; point0[iz] = points_z[0];
		point1[ix] = points_x[1]; point1[iy] = points_y[1]; point1[iz] = points_z[1];
		point2[ix] = points_x[2]; point2[iy] = points_y[2]; point2[iz] = points_z[2];

		triangle_rad[0] = TdS_Point3D_Angle(point0, point2, point1, NULL, 0);//201
		triangle_rad[1] = TdS_Point3D_Angle(point1, point0, point2, NULL, 0);//012
		triangle_rad[2] = TdS_Point3D_Angle(point2, point1, point0, NULL, 0);//120

		triangle_deg[0] = triangle_rad[0] * 180 / pi;
		triangle_deg[1] = triangle_rad[1] * 180 / pi;
		triangle_deg[2] = triangle_rad[2] * 180 / pi;

		retina_angles_deg[0] = TdS_Point3D_Angle(eye, point0, point1, NULL, 1);//01
		retina_angles_deg[1] = TdS_Point3D_Angle(eye, point1, point2, NULL, 1);//12
		retina_angles_deg[2] = TdS_Point3D_Angle(eye, point2, point0, NULL, 1);//20

		num_MinkovSawada = -9999;
		QueryPerformanceCounter(&start_time);
		num_MinkovSawada = P3P_MinkovSawadaB(triangle_deg[0], triangle_deg[1], retina_angles_deg[0], retina_angles_deg[1], retina_angles_deg[2], results_MinkovSawada);
		QueryPerformanceCounter(&end_time);

		subfunction_norm_interpretations(num_MinkovSawada, results_MinkovSawada);
		subfunction_sort_interpretations(num_MinkovSawada, results_MinkovSawada);
		
		int count_miss = 1;
		int count_error = 0; // not used for MinkovSawada
		int count_duplicate = 0;
		int count_sign = 0;
		for (int i = 0; i < num_MinkovSawada; i++)
		{
			if (subfunction_compare_interpretations(0, i, results_Original, results_MinkovSawada) < ERR_THRESHOLD0)
				count_miss = 0;
		}
		//cout << endl;
		time_MinkovSawada += (end_time.QuadPart - start_time.QuadPart) * 1000.0 / ticks_time.QuadPart;
		summary_MinkovSawada[0] += count_miss;


		/*Banno (2018)*/
		triangle_base = TdS_Point3D_Distance(point0, point1); //== Rab in Minkov & Sawada
		triangle_apex[ix] = triangle_base * tan(triangle_rad[1]) / (tan(triangle_rad[0]) + tan(triangle_rad[1]));
		triangle_apex[iy] = triangle_apex[ix] * tan(triangle_rad[0]);

		matrix_retina << retina_xy[0][ix], retina_xy[1][ix], retina_xy[2][ix], // columns = points, rows = xyz
			retina_xy[0][iy], retina_xy[1][iy], retina_xy[2][iy],
			retina_xy[0][iz], retina_xy[1][iz], retina_xy[2][iz];

		matrix_shape << 0, triangle_base, triangle_apex[ix], // columns = points, rows = xyz
			0, 0, triangle_apex[iy],
			0, 0, 0;

		num_Banno2018 = -9999;
		QueryPerformanceCounter(&start_time);
		num_Banno2018 = P3P_Banno2018.CalculatePosition(matrix_shape, matrix_retina, matrix_Rot, vector_Tran);
		QueryPerformanceCounter(&end_time);

		for (int i = 0; i < num_Banno2018; i++)
		{
			matrix_Tran << vector_Tran[i], vector_Tran[i], vector_Tran[i];
			matrix_results = matrix_Rot[i].transpose() * (matrix_shape - matrix_Tran); // columns = points, rows = xyz

			results_Banno2018[i][0] = subfunction_length(matrix_results(ix, 0), matrix_results(iy, 0), matrix_results(iz, 0));
			results_Banno2018[i][1] = subfunction_length(matrix_results(ix, 1), matrix_results(iy, 1), matrix_results(iz, 1));
			results_Banno2018[i][2] = subfunction_length(matrix_results(ix, 2), matrix_results(iy, 2), matrix_results(iz, 2));

			double tempx0 = retina_xy[0][ix] - matrix_results(ix, 0) / matrix_results(iz, 0);
			double tempy0 = retina_xy[0][iy] - matrix_results(iy, 0) / matrix_results(iz, 0);
			double tempx1 = retina_xy[1][ix] - matrix_results(ix, 1) / matrix_results(iz, 1);
			double tempy1 = retina_xy[1][iy] - matrix_results(iy, 1) / matrix_results(iz, 1);
			double tempx2 = retina_xy[2][ix] - matrix_results(ix, 2) / matrix_results(iz, 2);
			double tempy2 = retina_xy[2][iy] - matrix_results(iy, 2) / matrix_results(iz, 2);
			err1_Banno2018[i] = tempx0 * tempx0 + tempy0 * tempy0 + tempx1 * tempx1 + tempy1 * tempy1 + tempx2 * tempx2 + tempy2 * tempy2;

			double errAB = subfunction_compare_angles(	matrix_results(ix, 0), matrix_results(iy, 0), matrix_results(iz, 0),
														matrix_results(ix, 1), matrix_results(iy, 1), matrix_results(iz, 1),
														retina_angles_deg[0]);
			double errBC = subfunction_compare_angles(	matrix_results(ix, 1), matrix_results(iy, 1), matrix_results(iz, 1),
														matrix_results(ix, 2), matrix_results(iy, 2), matrix_results(iz, 2),
														retina_angles_deg[1]);
			double errCA = subfunction_compare_angles(	matrix_results(ix, 2), matrix_results(iy, 2), matrix_results(iz, 2),
														matrix_results(ix, 0), matrix_results(iy, 0), matrix_results(iz, 0),
														retina_angles_deg[2]);
			err2_Banno2018[i] = errAB + errBC + errCA;

			if (matrix_results(iz, 0) < 0 || matrix_results(iz, 1) < 0 || matrix_results(iz, 2) < 0)
				err3_Banno2018[i] = true;
			else
				err3_Banno2018[i] = false;
		}

		subfunction_norm_interpretations(num_Banno2018, results_Banno2018);
		subfunction_sort_interpretations(num_Banno2018, results_Banno2018);
		
		count_miss = 1;
		count_error = 0;
		count_duplicate = 0;
		count_sign = 0;
		for (int i = 0; i < num_Banno2018; i++)
		{
			if (err2_Banno2018[i] > ERR_THRESHOLD1)
			{
				count_error++;
				results_Banno2018[i][0] = -TdS_RandomNumber(1000, 2000);
				results_Banno2018[i][1] = -TdS_RandomNumber(1000, 2000);
				results_Banno2018[i][2] = -TdS_RandomNumber(1000, 2000);
			}
			else if (err3_Banno2018[i])
			{
				count_sign++;
				results_Banno2018[i][0] = -TdS_RandomNumber(4000, 5000);
				results_Banno2018[i][1] = -TdS_RandomNumber(4000, 5000);
				results_Banno2018[i][2] = -TdS_RandomNumber(4000, 5000);
			}
			else if (subfunction_compare_interpretations(0, i, results_Original, results_Banno2018) < ERR_THRESHOLD0)
				count_miss = 0;
			if (i != 0) if (subfunction_compare_interpretations(i, i - 1, results_Banno2018, results_Banno2018) < ERR_THRESHOLD2)
				count_duplicate++;
		}
		
		time_Banno2018 += (end_time.QuadPart - start_time.QuadPart) * 1000.0 / ticks_time.QuadPart;
		summary_Banno2018[0] += count_miss;
		if (count_error > 0)
			summary_Banno2018[count_error]++;
		if (num_MinkovSawada > (num_Banno2018 - count_error - count_duplicate - count_sign)) //5: less than MinkovSawada
			summary_Banno2018[5]++;
		else if (num_MinkovSawada == (num_Banno2018 - count_error - count_duplicate - count_sign)) //6: equal to MinkovSawada
			summary_Banno2018[6]++;
		else //7: more than MinkovSawada
			summary_Banno2018[7]++;


		/*OpenCV*/
		std::vector<cv::Point2f> imagePoints;
		std::vector<cv::Point3f> objectPoints;

		objectPoints.push_back(cv::Point3f(0, 0, 0)); //point0
		objectPoints.push_back(cv::Point3f(triangle_base, 0, 0)); //point1
		objectPoints.push_back(cv::Point3f(triangle_apex[ix], triangle_apex[iy], 0)); //point2

		imagePoints.push_back(cv::Point2f(retina_xy[0][ix], retina_xy[0][iy])); //point0
		imagePoints.push_back(cv::Point2f(retina_xy[1][ix], retina_xy[1][iy])); //point1
		imagePoints.push_back(cv::Point2f(retina_xy[2][ix], retina_xy[2][iy])); //point2
		//n.b. retina_xy[0~2][iz] = 1;


	/*OpenCV: Gao, Hou, Tang, & Chang (2003)*/
		ocv_flags = cv::SolvePnPMethod::SOLVEPNP_P3P; // Gao, Hou, Tang, & Chang (2003) 
		num_OCV_GaoEtal2003 = -9999;
		QueryPerformanceCounter(&start_time);
		num_OCV_GaoEtal2003 = cv::solveP3P(objectPoints, imagePoints, cameraMatrix, distCoeffs, ocv_rvecs, ocv_tvecs, ocv_flags);
		QueryPerformanceCounter(&end_time);

		for (int i = 0; i < num_OCV_GaoEtal2003; i++)
		{
			err1_CV_GaoEtal2003[i] = 0;

			cv::Mat rvec;
			cv::Rodrigues(ocv_rvecs[i], rvec);
			cv::Mat tvec = ocv_tvecs[i];
			double tempPoints[3][3];
			for (int p = 0; p < numPoints; p++)
			{
				cv::Matx31d tempPt0(objectPoints[p].x, objectPoints[p].y, objectPoints[p].z);
				cv::Mat tempPt1 = rvec * tempPt0 + tvec;
				results_CV_GaoEtal2003[i][p] = subfunction_length(tempPt1.at<double>(ix, 0), tempPt1.at<double>(iy, 0), tempPt1.at<double>(iz, 0));

				double tempxp = retina_xy[p][ix] - tempPt1.at<double>(ix, 0) / tempPt1.at<double>(iz, 0);
				double tempyp = retina_xy[p][iy] - tempPt1.at<double>(iy, 0) / tempPt1.at<double>(iz, 0);
				err1_CV_GaoEtal2003[i] += tempxp * tempxp + tempyp * tempyp;

				tempPoints[p][ix] = tempPt1.at<double>(ix, 0);
				tempPoints[p][iy] = tempPt1.at<double>(iy, 0);
				tempPoints[p][iz] = tempPt1.at<double>(iz, 0);
			}

			double errAB = subfunction_compare_angles(	tempPoints[0][ix], tempPoints[0][iy], tempPoints[0][iz],
														tempPoints[1][ix], tempPoints[1][iy], tempPoints[1][iz],
														retina_angles_deg[0]);
			double errBC = subfunction_compare_angles(	tempPoints[1][ix], tempPoints[1][iy], tempPoints[1][iz],
														tempPoints[2][ix], tempPoints[2][iy], tempPoints[2][iz],
														retina_angles_deg[1]);
			double errCA = subfunction_compare_angles(	tempPoints[2][ix], tempPoints[2][iy], tempPoints[2][iz],
														tempPoints[0][ix], tempPoints[0][iy], tempPoints[0][iz],
														retina_angles_deg[2]);
			err2_CV_GaoEtal2003[i] = errAB + errBC + errCA;

			if (tempPoints[0][iz] < 0 || tempPoints[1][iz] < 0 || tempPoints[2][iz] < 0)
				err3_CV_GaoEtal2003[i] = true;
			else
				err3_CV_GaoEtal2003[i] = false;
		}

		subfunction_norm_interpretations(num_OCV_GaoEtal2003, results_CV_GaoEtal2003);
		subfunction_sort_interpretations(num_OCV_GaoEtal2003, results_CV_GaoEtal2003);
		
		count_miss = 1;
		count_error = 0;
		count_duplicate = 0;
		count_sign = 0;
		for (int i = 0; i < num_OCV_GaoEtal2003; i++)
		{
			if (err2_CV_GaoEtal2003[i] > ERR_THRESHOLD1)
			{
				count_error++;
				results_CV_GaoEtal2003[i][0] = -TdS_RandomNumber(1000, 2000);
				results_CV_GaoEtal2003[i][1] = -TdS_RandomNumber(1000, 2000);
				results_CV_GaoEtal2003[i][2] = -TdS_RandomNumber(1000, 2000);
			}
			else if (err3_CV_GaoEtal2003[i])
			{
				count_sign++;
				results_CV_GaoEtal2003[i][0] = -TdS_RandomNumber(4000, 5000);
				results_CV_GaoEtal2003[i][1] = -TdS_RandomNumber(4000, 5000);
				results_CV_GaoEtal2003[i][2] = -TdS_RandomNumber(4000, 5000);
			}
			else if (subfunction_compare_interpretations(0, i, results_Original, results_CV_GaoEtal2003) < ERR_THRESHOLD0)
				count_miss = 0;
			if (i != 0) if (subfunction_compare_interpretations(i, i - 1, results_CV_GaoEtal2003, results_CV_GaoEtal2003) < ERR_THRESHOLD2)
				count_duplicate++;
		}
		
		time_GaoEtal2003 += (end_time.QuadPart - start_time.QuadPart) * 1000.0 / ticks_time.QuadPart;
		summary_GaoEtal2003[0] += count_miss;
		if (count_error > 0)
			summary_GaoEtal2003[count_error]++;
		if (num_MinkovSawada > (num_OCV_GaoEtal2003 - count_error - count_duplicate - count_sign)) //5: less than MinkovSawada
			summary_GaoEtal2003[5]++;
		else if (num_MinkovSawada == (num_OCV_GaoEtal2003 - count_error - count_duplicate - count_sign)) //6: equal to MinkovSawada
			summary_GaoEtal2003[6]++;
		else //7: more than MinkovSawada
			summary_GaoEtal2003[7]++;


		/*OpenCV: Ke & Roumeliotis (2017)*/
		ocv_flags = cv::SolvePnPMethod::SOLVEPNP_AP3P; // Ke & Roumeliotis (2017)
		num_OCV_KeRoumeliotis2017 = -9999;
		QueryPerformanceCounter(&start_time);
		num_OCV_KeRoumeliotis2017 = cv::solveP3P(objectPoints, imagePoints, cameraMatrix, distCoeffs, ocv_rvecs, ocv_tvecs, ocv_flags);
		QueryPerformanceCounter(&end_time);

		for (int i = 0; i < num_OCV_KeRoumeliotis2017; i++)
		{
			err1_CV_KeRoumeliotis2017[i] = 0;

			cv::Mat rvec;
			cv::Rodrigues(ocv_rvecs[i], rvec);
			cv::Mat tvec = ocv_tvecs[i];
			double tempPoints[3][3];
			for (int p = 0; p < numPoints; p++)
			{
				cv::Matx31d tempPt0(objectPoints[p].x, objectPoints[p].y, objectPoints[p].z);
				cv::Mat tempPt1 = rvec * tempPt0 + tvec;
				results_CV_KeRoumeliotis2017[i][p] = subfunction_length(tempPt1.at<double>(ix, 0), tempPt1.at<double>(iy, 0), tempPt1.at<double>(iz, 0));

				double tempxp = retina_xy[p][ix] - tempPt1.at<double>(ix, 0) / tempPt1.at<double>(iz, 0);
				double tempyp = retina_xy[p][iy] - tempPt1.at<double>(iy, 0) / tempPt1.at<double>(iz, 0);
				err1_CV_KeRoumeliotis2017[i] += tempxp * tempxp + tempyp * tempyp;

				tempPoints[p][ix] = tempPt1.at<double>(ix, 0);
				tempPoints[p][iy] = tempPt1.at<double>(iy, 0);
				tempPoints[p][iz] = tempPt1.at<double>(iz, 0);
			}

			double errAB = subfunction_compare_angles(tempPoints[0][ix], tempPoints[0][iy], tempPoints[0][iz],
				tempPoints[1][ix], tempPoints[1][iy], tempPoints[1][iz], retina_angles_deg[0]);
			double errBC = subfunction_compare_angles(tempPoints[1][ix], tempPoints[1][iy], tempPoints[1][iz],
				tempPoints[2][ix], tempPoints[2][iy], tempPoints[2][iz], retina_angles_deg[1]);
			double errCA = subfunction_compare_angles(tempPoints[2][ix], tempPoints[2][iy], tempPoints[2][iz],
				tempPoints[0][ix], tempPoints[0][iy], tempPoints[0][iz], retina_angles_deg[2]);
			err2_CV_KeRoumeliotis2017[i] = errAB + errBC + errCA;

			if (tempPoints[0][iz] < 0 || tempPoints[1][iz] < 0 || tempPoints[2][iz] < 0)
				err3_CV_KeRoumeliotis2017[i] = true;
			else
				err3_CV_KeRoumeliotis2017[i] = false;
		}

		subfunction_norm_interpretations(num_OCV_KeRoumeliotis2017, results_CV_KeRoumeliotis2017);
		subfunction_sort_interpretations(num_OCV_KeRoumeliotis2017, results_CV_KeRoumeliotis2017);
		
		count_miss = 1;
		count_error = 0;
		count_duplicate = 0;
		count_sign = 0;
		for (int i = 0; i < num_OCV_KeRoumeliotis2017; i++)
		{
			if (err2_CV_KeRoumeliotis2017[i] > ERR_THRESHOLD1)
			{
				count_error++;
				results_CV_KeRoumeliotis2017[i][0] = -TdS_RandomNumber(1000, 2000);
				results_CV_KeRoumeliotis2017[i][1] = -TdS_RandomNumber(1000, 2000);
				results_CV_KeRoumeliotis2017[i][2] = -TdS_RandomNumber(1000, 2000);
			}
			else if (err3_CV_KeRoumeliotis2017[i])
			{
				count_sign++;
				results_CV_KeRoumeliotis2017[i][0] = -TdS_RandomNumber(4000, 5000);
				results_CV_KeRoumeliotis2017[i][1] = -TdS_RandomNumber(4000, 5000);
				results_CV_KeRoumeliotis2017[i][2] = -TdS_RandomNumber(4000, 5000);
			}
			else if (subfunction_compare_interpretations(0, i, results_Original, results_CV_KeRoumeliotis2017) < ERR_THRESHOLD0)
				count_miss = 0;
			if (i != 0) if (subfunction_compare_interpretations(i, i - 1, results_CV_KeRoumeliotis2017, results_CV_KeRoumeliotis2017) < ERR_THRESHOLD2)
				count_duplicate++;
		}
		
		time_KeRoumeliotis2017 += (end_time.QuadPart - start_time.QuadPart) * 1000.0 / ticks_time.QuadPart;
		summary_KeRoumeliotis2017[0] += count_miss;
		if (count_error > 0)
			summary_KeRoumeliotis2017[count_error]++;
		if (num_MinkovSawada > (num_OCV_KeRoumeliotis2017 - count_error - count_duplicate - count_sign)) //5: less than MinkovSawada
			summary_KeRoumeliotis2017[5]++;
		else if (num_MinkovSawada == (num_OCV_KeRoumeliotis2017 - count_error - count_duplicate - count_sign)) //6: equal to MinkovSawada
			summary_KeRoumeliotis2017[6]++;
		else //7: more than MinkovSawada
			summary_KeRoumeliotis2017[7]++;
	} // trial


/*Output*/
	fout.open(fout_name, ios::out | ios::app);

	fout << endl;
	fout << "Name" << "\t" << "msec" << "\t" << "miss" << "\t";
	fout << "false-alarm-1" << "\t" << "false-alarm-2" << "\t" << "false-alarm-3" << "\t" << "false-alarm-4" << "\t";
	fout << "less than MinkovSawada" << "\t" << "equal to MinkovSawada" << "\t" << "more than MinkovSawada" << endl;

	fout << "MinkovSawada:" << "\t" << time_MinkovSawada << "\t" << summary_MinkovSawada[0] << endl;

	fout << "Banno (2018):" << "\t" << time_Banno2018 << "\t" << summary_Banno2018[0] << "\t";;
	fout << summary_Banno2018[1] << "\t" << summary_Banno2018[2] << "\t" << summary_Banno2018[3] << "\t" << summary_Banno2018[4] << "\t";
	fout << summary_Banno2018[5] << "\t" << summary_Banno2018[6] << "\t" << summary_Banno2018[7] << endl;

	fout << "Gao et al. (2003):" << "\t" << time_GaoEtal2003 << "\t" << summary_GaoEtal2003[0] << "\t";
	fout << summary_GaoEtal2003[1] << "\t" << summary_GaoEtal2003[2] << "\t" << summary_GaoEtal2003[3] << "\t" << summary_GaoEtal2003[4] << "\t";
	fout << summary_GaoEtal2003[5] << "\t" << summary_GaoEtal2003[6] << "\t" << summary_GaoEtal2003[7] << endl;

	fout << "Ke & Roumeliotis (2017):" << "\t" << time_KeRoumeliotis2017 << "\t" << summary_KeRoumeliotis2017[0] << "\t";
	fout << summary_KeRoumeliotis2017[1] << "\t" << summary_KeRoumeliotis2017[2] << "\t" << summary_KeRoumeliotis2017[3] << "\t" << summary_KeRoumeliotis2017[4] << "\t";
	fout << summary_KeRoumeliotis2017[5] << "\t" << summary_KeRoumeliotis2017[6] << "\t" << summary_KeRoumeliotis2017[7] << endl;

	fout.close();

	session++;
	if (session < NUM_SESSIONS)
		goto start_session;

	int dummy;
	cout << "Type any number to finish" << endl;
	cin >> dummy;


	return 0;
}

double subfunction_length(double x, double y, double z)
{
	return sqrt(x*x + y * y + z * z);
}

void subfunction_norm_interpretations(int num_interpretations, double(&interpretations)[10][3])
{
	double length;
	int point0 = 0, point1 = 1, point2 = 2;
	for (int i = 0; i < num_interpretations; i++)
	{
		length = subfunction_length(interpretations[i][point0], interpretations[i][point1], interpretations[i][point2]);
		interpretations[i][point0] /= length;
		interpretations[i][point1] /= length;
		interpretations[i][point2] /= length;
	}
}

void subfunction_sort_interpretations(int num_interpretations, double(&interpretations)[10][3])
{
	double sortedarray[10];
	int index[10];


	for (int i = 0; i < num_interpretations; i++)
		sortedarray[i] = interpretations[i][0];
	TdS_BubbleSort(sortedarray, index, num_interpretations);

	double copy_interpretation[10][3];
	for (int i = 0; i < num_interpretations; i++)
	{
		copy_interpretation[i][0] = interpretations[i][0];
		copy_interpretation[i][1] = interpretations[i][1];
		copy_interpretation[i][2] = interpretations[i][2];
	}

	for (int i = 0; i < num_interpretations; i++)
	{
		interpretations[i][0] = copy_interpretation[index[i]][0];
		interpretations[i][1] = copy_interpretation[index[i]][1];
		interpretations[i][2] = copy_interpretation[index[i]][2];
	}
}

double subfunction_compare_angles(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double ang12)
{
	double v1[3] = { p1x, p1y, p1z };
	double v2[3] = { p2x, p2y, p2z };

	if (p1z < 0) { v1[ix] *= -1; v1[iy] *= -1; v1[iz] *= -1; }
	if (p2z < 0) { v2[ix] *= -1; v2[iy] *= -1; v2[iz] *= -1; }


	double ang12_recovered = TdS_Vector3D_Angle(v1, v2, NULL, 1);

	return pow(ang12 - ang12_recovered, 2);

}

double subfunction_compare_interpretations(int ia, int ib, double(&interpretationsA)[10][3], double(&interpretationsB)[10][3])
{
	double diff = 0;
	for (int ixyz = 0; ixyz < 3; ixyz++)
	{
		diff += pow(interpretationsA[ia][ixyz] - interpretationsB[ib][ixyz], 2);
	}

	return diff;
}