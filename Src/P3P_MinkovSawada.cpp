// (c) Group of "Mathematical and Computational Psychology"
// https://social.hse.ru/en/psy/mcp/
//
// Please read the LICENSE and NO WARRANTY statement in :
// MinkovSawada2020_License.txt
//
// Jan/22/2020: Upload to GitHub by T. Sawada


#include "pch.h" // MS Visual Studio C++
#include "P3P_MinkovSawada.h"

using namespace std;

//const double EPS = 1e-14;
#define pi 3.14159265358979323846264338327950288
#define EPS (1e-14)
#define MaxNp 100


int P3P_MinkovSawadaA(double angleA, double angleB, double angleAB, double angleBC, double angleCA, double(&interpretations)[10][3])
{
	int numResults = 0;
	//double results[4][3];
	//vector<vector<double>> result;

	// Valid tetrahedron
	if ((angleAB <= 0) || (180 <= angleAB) || (angleBC <= 0) || (180 <= angleBC) || (angleCA <= 0) || (180 <= angleCA))
		return -1;
	if ((angleAB + angleBC < angleCA) || (angleBC + angleCA < angleAB) || (angleCA + angleAB < angleBC))
		return -1;
	if ((angleA <= 0) || (180 <= angleA) || (angleB <= 0) || (180 <= angleB))
		return -1;

	// Valid triangle
	if (180 <= angleA + angleB)
		return -1;
	if (360 <= angleAB + angleBC + angleCA)
		return -1;

	double angleC = 180 - (angleA + angleB);

	// Valid triangle and valid image but invalid tetrahedron
	if ((180 - angleAB) + (180 - angleBC) < angleB)
		return 0;
	if ((180 - angleBC) + (180 - angleCA) < angleC)
		return 0;
	if ((180 - angleCA) + (180 - angleAB) < angleA)
		return 0;

	double cosAB = cos(pi*angleAB / 180.0);
	double cosBC = cos(pi*angleBC / 180.0);
	double cosCA = cos(pi*angleCA / 180.0);
	double cosA = cos(pi*angleA / 180.0);
	double sinA = sin(pi*angleA / 180.0);
	double cosB = cos(pi*angleB / 180.0);
	double sinB = sin(pi*angleB / 180.0);

	double Rab = 1;
	double Rbc = Rab / (cosB + sinB * cosA / sinA);
	double Rca = Rab / (cosA + sinA * cosB / sinB);

	double K1 = ((pow(Rbc, 2)) / (pow(Rca, 2)));
	double K2 = ((pow(Rbc, 2) / pow(Rab, 2)));
	double K1K2 = K1 * K2;
	double cos2CA = cosCA * cosCA;
	double cos2BC = cosBC * cosBC;


	double G4 = pow(K1K2 - K1 - K2, 2) - 4 * K1K2*cos2BC;
	double G3 = 4 * (K1K2 - K1 - K2)*(K2 - K1K2)*cosAB
		+ 4 * K1*cosBC*((K1K2 + K2 - K1)*cosCA + 2 * K2*cosAB*cosBC);
	double G2 = pow(2 * (K2 - K1K2)*cosAB, 2)
		+ 2 * (K1K2 + K1 - K2)*(K1K2 - K1 - K2)
		+ 4 * K1*((K1 - K2)*cos2BC + (K1 - K1K2)*cos2CA - 2 * (K2 + K1K2)*cosAB*cosCA*cosBC);
	double G1 = 4 * (K1K2 + K1 - K2)*(K2 - K1K2)*cosAB
		+ 4 * K1*((K1K2 - K1 + K2)*cosCA*cosBC + 2 * K1K2*cosAB*cos2CA);
	double G0 = pow(K1K2 + K1 - K2, 2) - 4 * K1*K1K2*cos2CA;

	double x[4] = {};
	int numPolySolutions;
	numPolySolutions = TdS_QuatricEquation(x, G4, G3, G2, G1, G0);
	if (numPolySolutions <= 0)
		return 0;


	for (int i = 0; i < numPolySolutions; i++)
	{
		if (x[i] <= 0) continue;

		if (pow(x[i], 2) - 2 * x[i] * cosAB + 1 <= 0) continue;

		double a = Rab / sqrt(pow(x[i], 2) - 2 * x[i] * cosAB + 1);  // A24
		double b = x[i] * a;
		// a > 0 and b > 0 because Rab > 0

		if (angleA == angleBC && a < 1e-10) continue;
		if (angleB == angleCA && b < 1e-10) continue;

		double m1 = 1 - K1;
		double p1 = 2 * (K1 * cosCA - x[i] * cosBC);
		double q1 = (pow(x[i], 2) - K1);

		double m2 = 1;
		double p2 = 2 * (-x[i] * cosBC);
		double q2 = pow(x[i], 2) * (1 - K2) + 2 * x[i] * K2 * cosAB - K2;

		double mqmq = (m1 * q2 - m2 * q1);

		if (fabs(mqmq) > 1e-10)
		{
			double c0 = a * (p2*q1 - p1 * q2) / mqmq;

			if (c0 > 0)
			{
				if (angleC == angleAB && c0 < 1e-10) continue;
				interpretations[numResults][0] = a;
				interpretations[numResults][1] = b;
				interpretations[numResults][2] = c0;
				numResults++;
			}
		}
		else
		{
			double y0 = cosCA * cosCA + (Rca*Rca - a * a) / (a*a);
			if (y0 < 0) continue;

			double c1 = a * (cosCA + sqrt(y0)); // A27
			double c2 = a * (cosCA - sqrt(y0)); // A27

			if (c1 > 0)
			{
				if (angleC == angleAB && c1 < 1e-10) continue;
				if (fabs(b*b + c1 * c1 - 2 * b*c1*cosBC - Rbc * Rbc) > 1e-10) continue; // A3
				interpretations[numResults][0] = a;
				interpretations[numResults][1] = b;
				interpretations[numResults][2] = c1;
				numResults++;
			}

			if (c2 > 0 && sqrt(y0) != 0)
			{
				if (angleC == angleAB && c2 < 1e-10) continue;
				if (fabs(b*b + c2 * c2 - 2 * b*c2*cosBC - Rbc * Rbc) > 1e-10) continue; // A3
				interpretations[numResults][0] = a;
				interpretations[numResults][1] = b;
				interpretations[numResults][2] = c2;
				numResults++;
			}
		}
	}

	return numResults;
}


int P3P_MinkovSawadaB(double angleA, double angleB, double angleAB, double angleBC, double angleCA, double(&interpretations)[10][3])
{
	interpretations[0][0] = 0; interpretations[0][1] = 0; interpretations[0][2] = 0;
	interpretations[1][0] = 0; interpretations[1][1] = 0; interpretations[1][2] = 0;
	interpretations[2][0] = 0; interpretations[2][1] = 0; interpretations[2][2] = 0;
	interpretations[3][0] = 0; interpretations[3][1] = 0; interpretations[3][2] = 0;

	// Valid tetrahedron
	if ((angleAB <= 0) || (180 <= angleAB) || (angleBC <= 0) || (180 <= angleBC) || (angleCA <= 0) || (180 <= angleCA))
		return -1;
	if ((angleAB + angleBC < angleCA) || (angleBC + angleCA < angleAB) || (angleCA + angleAB < angleBC))
		return -2;
	if ((angleA <= 0) || (180 <= angleA) || (angleB <= 0) || (180 <= angleB))
		return -3;

	// Valid triangle
	if (180 <= angleA + angleB)
		return -4;
	if (360 <= angleAB + angleBC + angleCA)
		return -5;

	double angleC = 180 - (angleA + angleB);

	// Valid triangle and valid image but invalid tetrahedron
	//  #interpretations = 0
	if ((180 - angleAB) + (180 - angleBC) < angleB)
		return 0;
	if ((180 - angleBC) + (180 - angleCA) < angleC)
		return 0;
	if ((180 - angleCA) + (180 - angleAB) < angleA)
		return 0;





	// New Jun/28/2018
	double angAB, angBC, angCA, angA, angB, angC;
	double tempA[3] = { angleA, angleB, angleC };
	double tempB[3] = { angleB, angleC, angleA };
	double tempC[3] = { angleC, angleA, angleB };
	double tempAB[3] = { angleAB, angleBC, angleCA };
	double tempBC[3] = { angleBC, angleCA, angleAB };
	double tempCA[3] = { angleCA, angleAB, angleBC };
	double interpretationsX3[3][8][3] = {};
	int numInterpretationsX3[3] = {};

	for (int ang = 0; ang < 3; ang++)
	{
		angA = tempA[ang];  angB = tempB[ang];  angC = tempC[ang];
		angAB = tempAB[ang]; angBC = tempBC[ang]; angCA = tempCA[ang];

		double cosAB = cos(pi*angAB / 180.0);
		double cosBC = cos(pi*angBC / 180.0);
		double cosCA = cos(pi*angCA / 180.0);
		double cosA = cos(pi*angA / 180.0);
		double sinA = sin(pi*angA / 180.0);
		double cosB = cos(pi*angB / 180.0);
		double sinB = sin(pi*angB / 180.0);

		double Rab = 100;
		double Rbc = Rab / (cosB + sinB * cosA / sinA);
		double Rca = Rab / (cosA + sinA * cosB / sinB);

		//double tanABx, tanBCx, tanCAx;
		//double ratio_ang_R[3], ratio_max;
		//int i_max;
		//tanABx = (angAB<=84 ? tan(pi*angAB / 180.0) : 10);
		//tanBCx = (angBC<=84 ? tan(pi*angBC / 180.0) : 10);
		//tanCAx = (angCA<=84 ? tan(pi*angCA / 180.0) : 10);
		//ratio_ang_R[0] = Rab / tanABx;
		//ratio_ang_R[1] = Rbc / tanBCx;
		//ratio_ang_R[2] = Rca / tanCAx;

		//ratio_max = TdS_Max(3, ratio_ang_R, &i_max);

		//Rab /= ratio_max/100.0;
		//Rbc /= ratio_max/100.0;
		//Rca /= ratio_max/100.0;

		double K1 = ((pow(Rbc, 2)) / (pow(Rca, 2)));
		double K2 = ((pow(Rbc, 2) / pow(Rab, 2)));
		double K1K2 = K1 * K2;
		double cos2CA = cosCA * cosCA;
		double cos2BC = cosBC * cosBC;

		double G4 = pow(K1K2 - K1 - K2, 2) - 4 * K1K2*cos2BC;
		double G3 = 4 * (K1K2 - K1 - K2)*(K2 - K1K2)*cosAB
			+ 4 * K1*cosBC*((K1K2 + K2 - K1)*cosCA + 2 * K2*cosAB*cosBC);
		double G2 = pow(2 * (K2 - K1K2)*cosAB, 2)
			+ 2 * (K1K2 + K1 - K2)*(K1K2 - K1 - K2)
			+ 4 * K1*((K1 - K2)*cos2BC + (K1 - K1K2)*cos2CA - 2 * (K2 + K1K2)*cosAB*cosCA*cosBC);
		double G1 = 4 * (K1K2 + K1 - K2)*(K2 - K1K2)*cosAB
			+ 4 * K1*((K1K2 - K1 + K2)*cosCA*cosBC + 2 * K1K2*cosAB*cos2CA);
		double G0 = pow(K1K2 + K1 - K2, 2) - 4 * K1*K1K2*cos2CA;

		double x[4] = {};
		int numPolySolutions;
		numPolySolutions = TdS_QuatricEquation(x, G4, G3, G2, G1, G0);

		if (numPolySolutions <= 0)
			continue;

		for (int i = 0; i < numPolySolutions; i++)
		{
			if (x[i] <= 0) continue;

			if (pow(x[i], 2) - 2 * x[i] * cosAB + 1 <= 0) continue;

			double a = Rab / sqrt(pow(x[i], 2) - 2 * x[i] * cosAB + 1);  // A24
			double b = x[i] * a;
			// a > 0 and b > 0 because Rab > 0

			if (angA == angBC && a < 1e-10) continue;
			if (angB == angCA && b < 1e-10) continue;

			double y0 = cosCA * cosCA + (Rca*Rca - a * a) / (a*a); // Equation A?
			if (y0 < 0) continue;


			double m1 = 1 - K1;
			double p1 = 2 * (K1 * cosCA - x[i] * cosBC);
			double q1 = (pow(x[i], 2) - K1);
			double m2 = 1;
			double p2 = 2 * (-x[i] * cosBC);
			double q2 = pow(x[i], 2) * (1 - K2) + 2 * x[i] * K2 * cosAB - K2;
			double mqmq = (m1 * q2 - m2 * q1);

			y0 += 0;
			if (false) // (fabs(mqmq) > 1e-10) // August/2018
			{
				double c0 = a * (p2*q1 - p1 * q2) / mqmq;

				if (c0 > 0)
				{
					if (angC == angAB && c0 < 1e-10) continue;

					interpretationsX3[ang][numInterpretationsX3[ang]][0] = a;
					interpretationsX3[ang][numInterpretationsX3[ang]][1] = b;
					interpretationsX3[ang][numInterpretationsX3[ang]][2] = c0;
					numInterpretationsX3[ang]++;
				}
			}
			else
			{
				double c1 = a * (cosCA + sqrt(y0)); // A27
				double c2 = a * (cosCA - sqrt(y0)); // A27
				double criterion1 = fabs(b*b + c1 * c1 - 2 * b*c1*cosBC - Rbc * Rbc); // A3
				double criterion2 = fabs(b*b + c2 * c2 - 2 * b*c2*cosBC - Rbc * Rbc); // A3

				if (c1 > 0 && (angC == angAB && c1 < 1e-10) == false && criterion1 <= 1e-4*Rab)
				{
					//if (angC == angAB && c1 < 1e-10) continue;
					//if (fabs(b*b + c1*c1 - 2 * b*c1*cosBC - Rbc*Rbc) > 1e-10) continue; // A3

					interpretationsX3[ang][numInterpretationsX3[ang]][0] = a;
					interpretationsX3[ang][numInterpretationsX3[ang]][1] = b;
					interpretationsX3[ang][numInterpretationsX3[ang]][2] = c1;
					numInterpretationsX3[ang]++;
				}

				if (c2 > 0 && (angC == angAB && c2 < 1e-10) == false && criterion2 <= 1e-4*Rab)
				{
					//if (angC == angAB && c2 < 1e-10) continue;
					//if (fabs(b*b + c2*c2 - 2 * b*c2*cosBC - Rbc*Rbc) > 1e-10) continue; // A3

					interpretationsX3[ang][numInterpretationsX3[ang]][0] = a;
					interpretationsX3[ang][numInterpretationsX3[ang]][1] = b;
					interpretationsX3[ang][numInterpretationsX3[ang]][2] = c2;
					numInterpretationsX3[ang]++;
				}
			}
		} // for (int i = 0; i < numPolySolutions; i++)

		// Check redundancy
		int check_redundancy[8] = {};
		double diff;
		for (int i = 0; i < numInterpretationsX3[ang] - 1; i++)
		{
			if (check_redundancy[i] == 1) continue;
			for (int j = i + 1; j < numInterpretationsX3[ang]; j++)
			{
				diff = pow(interpretationsX3[ang][i][0] - interpretationsX3[ang][j][0], 2.0);
				diff += pow(interpretationsX3[ang][i][1] - interpretationsX3[ang][j][1], 2.0);
				diff += pow(interpretationsX3[ang][i][2] - interpretationsX3[ang][j][2], 2.0);

				if (diff < 1e-4*Rab)
					check_redundancy[j] = 1;
			}
		}
		int k = 0;
		for (int i = 0; i < numInterpretationsX3[ang]; i++)
		{
			if (check_redundancy[i] == 1)
				continue;
			else
			{
				interpretationsX3[ang][k][0] = interpretationsX3[ang][i][0];
				interpretationsX3[ang][k][1] = interpretationsX3[ang][i][1];
				interpretationsX3[ang][k][2] = interpretationsX3[ang][i][2];
				k++;
			}
		}
		numInterpretationsX3[ang] = k;

		// Examine results
		double la, lb, lc;
		double reangV, reangT1, reangT2;
		int validity[8] = {};
		double check3D1a, check3D1b, check3D2;
		for (int i = 0; i < numInterpretationsX3[ang]; i++)
		{
			la = interpretationsX3[ang][i][0];
			lb = interpretationsX3[ang][i][1];
			lc = interpretationsX3[ang][i][2];

			// (a,b)
			reangV = acos((la*la + lb * lb - Rab * Rab) / (2 * la*lb))  * 180.0 / pi;
			reangT1 = acos((Rab*Rab + la * la - lb * lb) / (2 * Rab*la)) * 180.0 / pi;
			reangT2 = acos((Rab*Rab + lb * lb - la * la) / (2 * Rab*lb)) * 180.0 / pi;
			check3D1a = fabs(reangV - angAB);
			check3D1b = angAB;
			check3D2 = 180 - reangV - reangT1 - reangT2;
			if (fabs(reangV - angAB) > 0.01 && fabs(reangV - angAB) > angAB*0.01)
			{
				//											printf("Err1: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
				validity[i] = 1;
			}
			//if (fabs(180 - reangV-reangT1-reangT2) > 0.01)
			//{
			//	printf("Err2: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
			//	validity[i] = 1;
			//}

			// (b,c)
			reangV = acos((lb*lb + lc * lc - Rbc * Rbc) / (2 * lb*lc))  * 180.0 / pi;
			reangT1 = acos((Rbc*Rbc + lb * lb - lc * lc) / (2 * Rbc*lb)) * 180.0 / pi;
			reangT2 = acos((Rbc*Rbc + lc * lc - lb * lb) / (2 * Rbc*lc)) * 180.0 / pi;
			check3D1a = fabs(reangV - angBC);
			check3D1b = angBC;
			check3D2 = 180 - reangV - reangT1 - reangT2;
			if (fabs(reangV - angBC) > 0.01 && fabs(reangV - angBC) > angBC*0.01)
			{
				//											printf("Err1: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
				validity[i] = 1;
			}
			//if (fabs(180 - reangV - reangT1 - reangT2) > 0.01)
			//{
			//	printf("Err2: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
			//	validity[i] = 1;
			//}

			// (c,a)
			reangV = acos((lc*lc + la * la - Rca * Rca) / (2 * lc*la))  * 180.0 / pi;
			reangT1 = acos((Rca*Rca + lc * lc - la * la) / (2 * Rca*lc)) * 180.0 / pi;
			reangT2 = acos((Rca*Rca + la * la - lc * lc) / (2 * Rca*la)) * 180.0 / pi;
			check3D1a = fabs(reangV - angCA);
			check3D1b = angCA;
			check3D2 = 180 - reangV - reangT1 - reangT2;
			if (fabs(reangV - angCA) > 0.01 && fabs(reangV - angCA) > angCA*0.01)
			{
				//											printf("Err1: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
				validity[i] = 1;
			}
			//if (fabs(180 - reangV-reangT1-reangT2) > 0.01)
			//{
			//	printf("Err2: %f, %f, %f\n", check3D1a, check3D1b, check3D2);
			//	validity[i] = 1;
			//}
		}
		k = 0;
		for (int i = 0; i < numInterpretationsX3[ang]; i++)
		{
			if (validity[i] == 1)
				continue;
			interpretationsX3[ang][k][0] = interpretationsX3[ang][i][0];
			interpretationsX3[ang][k][1] = interpretationsX3[ang][i][1];
			interpretationsX3[ang][k][2] = interpretationsX3[ang][i][2];
			k++;
		}
		numInterpretationsX3[ang] = k;

		if (4 < numInterpretationsX3[ang])
		{
			//printf("Err3 (4<num): %f, %f, %f\n", check3D1a, check3D1b, check3D2);
			//numInterpretationsX3[ang] = 4;

			double diff2[10][10] = {}; // 10 is large enough
			for (int i = 0; i < numInterpretationsX3[ang] - 1; i++)
				for (int j = i + 1; j < numInterpretationsX3[ang]; j++)
				{
					diff2[i][j]  = pow(interpretationsX3[ang][i][0] - interpretationsX3[ang][j][0], 2.0);
					diff2[i][j] += pow(interpretationsX3[ang][i][1] - interpretationsX3[ang][j][1], 2.0);
					diff2[i][j] += pow(interpretationsX3[ang][i][2] - interpretationsX3[ang][j][2], 2.0);
				}

			for (int n = 4; n < numInterpretationsX3[ang]; n++)
			{
				double diff2min = 999998;
				int i_diff2min;
				int j_diff2min;
				for (int i = 0; i < numInterpretationsX3[ang] - 1; i++)
					for (int j = i + 1; j < numInterpretationsX3[ang]; j++)
						if(diff2[i][j]<diff2min)
						{
							i_diff2min = i;
							j_diff2min = j;
						}
				diff2[i_diff2min][j_diff2min] = 999999;
				interpretationsX3[ang][j_diff2min][0] = 999999;
				interpretationsX3[ang][j_diff2min][1] = 999999;
				interpretationsX3[ang][j_diff2min][2] = 999999;
			}

			int temp_i = 0;
			for (int i = 0; i < numInterpretationsX3[ang]; i++)
			{
				if (interpretationsX3[ang][i][0] > 999998 && interpretationsX3[ang][i][1] > 999998 && interpretationsX3[ang][i][2] > 999998)
					continue;
				interpretationsX3[ang][temp_i][0] = interpretationsX3[ang][i][0];
				interpretationsX3[ang][temp_i][1] = interpretationsX3[ang][i][1];
				interpretationsX3[ang][temp_i][2] = interpretationsX3[ang][i][2];
				temp_i++;
			}
			numInterpretationsX3[ang] = 4;
		}


		// Sort results
		double vj0, vj1, v_swap[3];
		for (int i = 0; i < numInterpretationsX3[ang] - 1; i++)
			for (int j = numInterpretationsX3[ang] - 1; j > i; j--)
			{
				vj0 = interpretationsX3[ang][j][0] * 10000 + interpretationsX3[ang][j][1] * 100 + interpretationsX3[ang][j][2];
				vj1 = interpretationsX3[ang][j - 1][0] * 10000 + interpretationsX3[ang][j - 1][1] * 100 + interpretationsX3[ang][j - 1][2];

				if (vj0 < vj1)
				{
					for (int k = 0; k < 3; k++) v_swap[k] = interpretationsX3[ang][j][k];
					for (int k = 0; k < 3; k++) interpretationsX3[ang][j][k] = interpretationsX3[ang][j - 1][k];
					for (int k = 0; k < 3; k++) interpretationsX3[ang][j - 1][k] = v_swap[k];
				}
			}

		double standardization = 0;
		if (ang == 0) standardization = Rab;
		if (ang == 1) standardization = Rca;
		if (ang == 2) standardization = Rbc;
		//August 11
		//for (int i = 0; i < numInterpretationsX3[ang]; i++) for (int k = 0; k < 3; k++)
		//	interpretationsX3[ang][i][k] = interpretationsX3[ang][i][0] / standardization;

		for (int i = 0; i < numInterpretationsX3[ang]; i++)
		{
			double temp_interpretations[3];
			if (ang == 0)
			{
				temp_interpretations[0] = interpretationsX3[ang][i][0] / Rab;
				temp_interpretations[1] = interpretationsX3[ang][i][1] / Rab;
				temp_interpretations[2] = interpretationsX3[ang][i][2] / Rab;
			}
			else if (ang == 1)
			{
				temp_interpretations[0] = interpretationsX3[ang][i][2] / Rca;
				temp_interpretations[1] = interpretationsX3[ang][i][0] / Rca;
				temp_interpretations[2] = interpretationsX3[ang][i][1] / Rca;
			}
			else //ang == 2
			{
				temp_interpretations[0] = interpretationsX3[ang][i][1] / Rbc;
				temp_interpretations[1] = interpretationsX3[ang][i][2] / Rbc;
				temp_interpretations[2] = interpretationsX3[ang][i][0] / Rbc;
			}
			for (int k = 0; k < 3; k++)
				interpretationsX3[ang][i][k] = temp_interpretations[k];
		}
	} // for (int ang = 0; ang < 3; ang++)


	// Copy data
	int numInterpretations = 0;
	int check_consistency[3][4] = {};
	for (int ang = 0; ang < 3; ang++) for (int i = 0; i < numInterpretationsX3[ang]; i++) check_consistency[ang][i] = 1;
	//ang = 0
	for (int i = 0; i < numInterpretationsX3[0]; i++)
	{
		for (int abc = 0; abc < 3; abc++) interpretations[i][abc] = interpretationsX3[0][i][abc];
		check_consistency[0][i] = 0;
		numInterpretations++;

		//ang = 0 vs ang = 1+2
		for (int ang = 1; ang < 3; ang++) for (int j = 0; j < numInterpretationsX3[ang]; j++)
		{
			double diff;
			diff = pow(interpretationsX3[0][i][0] - interpretationsX3[ang][j][0], 2.0);
			diff += pow(interpretationsX3[0][i][1] - interpretationsX3[ang][j][1], 2.0);
			diff += pow(interpretationsX3[0][i][2] - interpretationsX3[ang][j][2], 2.0);
			if (diff < 1e-2)
				check_consistency[ang][j] = 0;
		}
	}
	//ang = 1
	for (int j = 0; j < numInterpretationsX3[1]; j++)
	{
		if (check_consistency[1][j] == 0)
			continue;
		for (int abc = 0; abc < 3; abc++) interpretations[numInterpretations][abc] = interpretationsX3[1][j][abc];
		check_consistency[1][j] = 0;
		numInterpretations++;

		//ang = 1 vs ang = 2
		for (int k = 0; k < numInterpretationsX3[2]; k++)
		{
			double diff;
			diff = pow(interpretationsX3[1][j][0] - interpretationsX3[2][k][0], 2.0);
			diff += pow(interpretationsX3[1][j][1] - interpretationsX3[2][k][1], 2.0);
			diff += pow(interpretationsX3[1][j][2] - interpretationsX3[2][k][2], 2.0);
			if (diff < 1e-2)
				check_consistency[2][k] = 0;
		}
	}
	//ang = 2
	for (int k = 0; k < numInterpretationsX3[2]; k++)
	{
		if (check_consistency[2][k] == 0)
			continue;
		for (int abc = 0; abc < 3; abc++) interpretations[numInterpretations][abc] = interpretationsX3[2][k][abc];
		numInterpretations++;
	}

	//int total_consistency = 0;
	//for (int ang = 0; ang < 3; ang++) for (int i = 0; i < 4; i++) total_consistency += check_consistency[ang][i];
	//if (total_consistency > 0)
	//	total_consistency = 0;

	return numInterpretations;
}