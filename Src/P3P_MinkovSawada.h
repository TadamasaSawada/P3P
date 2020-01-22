// (c) Group of "Mathematical and Computational Psychology"
// https://social.hse.ru/en/psy/mcp/
//
// Please read the LICENSE and NO WARRANTY statement in :
// MinkovSawada2020_License.txt
//
// Jan/22/2020: Uploaded to GitHub by T. Sawada

#pragma once

#include <cmath>
#include "TdS_Math.h"

#include <iostream>
//#include <vector>
//#include <iterator>
//#include <algorithm>
//#include <fstream>
//#include <sstream>
//#include <ctime>
//#include <random>


int P3P_MinkovSawadaA(double angleA, double angleB, double angleAB, double angleBC, double angleCA, double(&interpretations)[10][3]);//After Fischler & Bolles (1981)
int P3P_MinkovSawadaB(double angleA, double angleB, double angleAB, double angleBC, double angleCA, double(&interpretations)[10][3]);//Revised for precision+accuracy