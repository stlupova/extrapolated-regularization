#include <iostream>
#include <vector>
#include <iomanip>

#include "stokesIntegrals.h" 

using namespace std;

//*********************************************************************
// evaluates stokes single layer potential (Surfc.f is density)

void eval_Stokes_SL_offSurf_3del(double h, double DEL1, double DEL2, double DEL3,
				 int N_quad, const vector<Surf_point>& Surfc,
				 int N_target, const vector<Target_point>& Target,
				 vector<double>& SL1,
				 vector<double>& SL2,
				 vector<double>& SL3, 
				 int use_subtraction) {

  vector<double> pt(3,0), dx(3,0);
  double f0_dot_n0;
  double EightPI = 1.0/(8.0*PI);
  
  initialize_vector(3*N_target, SL1);
  initialize_vector(3*N_target, SL2);
  initialize_vector(3*N_target, SL3);

  for (int i=0; i<N_target; i++) {    

    pt[0] = Target[i].x;
    pt[1] = Target[i].y;
    pt[2] = Target[i].z;

    if (use_subtraction == 1) {
      f0_dot_n0 = dot_product(Target[i].nrst_f, Target[i].nrst_Nrml);
    }
    
    for (int j=0; j<N_quad; j++) {


      for (int k=0; k<3; k++)  dx[k] = pt[k] - Surfc[j].x[k];      

      double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      
      vector<double> f = Surfc[j].f;

      if (use_subtraction == 1) {
	for (int k=0; k<3; k++)  f[k] -= f0_dot_n0 * Surfc[j].Nrml[k];
      }
      
      double f_dot_dx = dot_product(f, dx);

      double H1_d1, H2_d1, H1_d2, H2_d2, H1_d3, H2_d3;
      H_gauss_opt2(r, DEL1, H1_d1, H2_d1); 
      H_gauss_opt2(r, DEL2, H1_d2, H2_d2); 
      H_gauss_opt2(r, DEL3, H1_d3, H2_d3); 

      for (int k=0; k<3; k++) {
	double m1 = f[k] * Surfc[j].Area;
	double m2 = dx[k] * f_dot_dx * Surfc[j].Area;
	
	SL1[3*i+k] += m1 * H1_d1 + m2 * H2_d1;
	SL2[3*i+k] += m1 * H1_d2 + m2 * H2_d2;
	SL3[3*i+k] += m1 * H1_d3 + m2 * H2_d3;
      }

    }
    for (int k=0; k<3; k++) {
      SL1[3*i+k] *= EightPI;
      SL2[3*i+k] *= EightPI;
      SL3[3*i+k] *= EightPI;
    }
  }
}

//*********************************************************************
// evaluates stokes double layer potential (Surfc.g is density)

void eval_Stokes_DL_offSurf_3del(double h, double DEL1, double DEL2, double DEL3,
				 int N_quad, const vector<Surf_point>& Surfc,
				 int N_target, const vector<Target_point>& Target,
				 vector<double>& DL1,
				 vector<double>& DL2,
				 vector<double>& DL3) {
  
  vector<double> pt(3,0), dx(3,0), x_hat(3,0);
  double EightPI6 = -6.0/(8.0*PI);

  initialize_vector(3*N_target, DL1);
  initialize_vector(3*N_target, DL2);
  initialize_vector(3*N_target, DL3);
  
  for (int i=0; i<N_target; i++) {    

    double u = 0.0, v = 0.0, w = 0.0;
    pt[0] = Target[i].x;
    pt[1] = Target[i].y;
    pt[2] = Target[i].z;
        
    vector<double> x0(3,0), n0(3,0), g0(3,0);

    x0[0] = Target[i].nrst_x;
    x0[1] = Target[i].nrst_y;
    x0[2] = Target[i].nrst_z;
    
    n0[0] = Target[i].nrst_Nrml[0];
    n0[1] = Target[i].nrst_Nrml[1];
    n0[2] = Target[i].nrst_Nrml[2];

    g0[0] = Target[i].nrst_g[0];
    g0[1] = Target[i].nrst_g[1];
    g0[2] = Target[i].nrst_g[2];

    double b = Target[i].b;
    
    for (int j=0; j<N_quad; j++) {
      
      for (int k=0; k<3; k++) {
	dx[k] = pt[k] - Surfc[j].x[k];
	x_hat[k] = Surfc[j].x[k] - x0[k];
      }
      
      double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      double sig_sq = dot_product(x_hat, x_hat) - 2.0*b* dot_product(x_hat, n0);
      
      vector<double> g = Surfc[j].g;
      for (int k=0; k<3; k++)  g[k] -= g0[k];

      double g_dot_n0 = dot_product(g, n0);
      double g_dot_x_hat = dot_product(g, x_hat);
      double n0_dot_n = dot_product(n0, Surfc[j].Nrml);
      double x_hat_dot_n = dot_product(x_hat, Surfc[j].Nrml);

      double c_m1a = g_dot_n0 * n0_dot_n;
      double c_m1b = g_dot_x_hat * n0_dot_n  +  g_dot_n0 * x_hat_dot_n;

      double c_m2a = g_dot_x_hat * x_hat_dot_n;
      double c_m2b = b* (g_dot_x_hat * n0_dot_n  +  g_dot_n0 * x_hat_dot_n);

      double s2_d1, s3_d1, s2_d2, s3_d2, s2_d3, s3_d3;
      Stresslet_orig(r, DEL1, s2_d1, s3_d1);
      Stresslet_orig(r, DEL2, s2_d2, s3_d2);
      Stresslet_orig(r, DEL3, s2_d3, s3_d3);
      
      vector<double> m1(3,0), m2(3,0);

      for (int k=0; k<3; k++) {

	m1[k] = dx[k] * c_m1a - n0[k] * c_m1b;

	m2[k] = ( -sig_sq * m1[k] + dx[k] * c_m2a + x_hat[k] * c_m2b) * Surfc[j].Area;

	m1[k] *= Surfc[j].Area;
      
	DL1[3*i+k] += m1[k] * s2_d1 + m2[k] * s3_d1;
	DL2[3*i+k] += m1[k] * s2_d2 + m2[k] * s3_d2;
	DL3[3*i+k] += m1[k] * s2_d3 + m2[k] * s3_d3;
      }
      
    }

    double xi = 0.0; // target point is outside the boundary
    if ( abs(Target[i].b) < 1.0e-12 ) { // target point is on the boundary
      xi = 0.5; 
    }
    else if ( Target[i].b < 0.0 ) { // target point is inside the boundary
      xi = 1.0;
    }

    for (int k=0; k<3; k++) {
      DL1[3*i+k] = DL1[3*i+k] * EightPI6 + xi * g0[k];
      DL2[3*i+k] = DL2[3*i+k] * EightPI6 + xi * g0[k];
      DL3[3*i+k] = DL3[3*i+k] * EightPI6 + xi * g0[k];
    }
  }
}
