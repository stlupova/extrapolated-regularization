#include <iostream>
#include <vector>
#include <cmath>

#include "surface.h"

using namespace std;

void eval_Stokes_SL_offSurf_3del(double h, double DEL1, double DEL2, double DEL3,
				 int N_quad, const vector<Surf_point>& Surfc,
				 int N_target, const vector<Target_point>& Target,
				 vector<double>& SL1,
				 vector<double>& SL2,
				 vector<double>& SL3, 
				 int use_subtraction);

void eval_Stokes_DL_offSurf_3del(double h, double DEL1, double DEL2, double DEL3,
				 int N_quad, const vector<Surf_point>& Surfc,
				 int N_target, const vector<Target_point>& Target,
				 vector<double>& DL1,
				 vector<double>& DL2,
				 vector<double>& DL3);
