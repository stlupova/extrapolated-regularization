/* 
   this program reproduces results in Fig. 13
   in "Extrapolated regularization..."
   for the molecular surface
   values are in Table 16 of original version on arxiv
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <time.h>

#include "stokesIntegrals.h"

using namespace std;

void test_Stokes_Sum_5ord(double h,
			  int N_quad, vector<Surf_point>& Surfc,
			  int N_target, vector<Target_point>& Target);

void set_rho_del(double h,
		 vector<double>& rho,
		 vector<double>& DEL);

void solve_3by3(double h, const vector<double>& rho, const vector<double>& DEL,
		int N_target, const vector<Target_point>& Target,
		const vector<double>& S_d1,
		const vector<double>& S_d2,
		const vector<double>& S_d3,
		vector<double>& S_comp);

double I0(double lam);
double I2(double lam);

void cout_norms(int n,
		const vector<double>& S_ex,
		const vector<double>& S_comp);

//*********************************************************************

int main(int argc, char** argv) {
  int N = N_gridlines;
  double L = 1.1;   // bounding box: [-L,L]^3
  double h = 1.0 / N; // grid spacing

  clock_t tm = clock();
  
  vector<Surf_point> Surfc; // array of surface (quadrature) points                    
  Generate_Surface(N, h, &Surfc);
  int N_quad = Surfc.size();
  cout << "Number of quadrature points = " << N_quad << endl;

  
  vector<Target_point> Target; // array of target points
  Generate_Targets(h,N_quad,Surfc,&Target);
  int N_target = Target.size();
  cout << "Number of target points = " << N_target << endl;
  

  test_Stokes_Sum_5ord(h, N_quad, Surfc, N_target, Target);

  
  tm = clock() - tm;
  cout << "CPU time = " << ((float)tm)/CLOCKS_PER_SEC << " seconds" << endl;
  
  return 0;
}

//*********************************************************************
//*********************************************************************

void test_Stokes_Sum_5ord(double h,
			  int N_quad, vector<Surf_point>& Surfc,
			  int N_target, vector<Target_point>& Target) {
  
  cout << "Computing Stokes SUM of SL and DL..." << endl;
  
  int N_quad3 = 3 * N_quad;
  int N_target3 = 3 * N_target;
  
  vector<double> S_ex(N_target3,0), S_comp(N_target3,0);

  for (int i=0; i<N_quad; i++) {
    Surfc[i].f = stokeslet_density(Surfc[i].x, Surfc[i].Nrml);
    Surfc[i].g = stresslet_density(Surfc[i].x, Surfc[i].Nrml);
  }

  Exact_solution(N_target, Target);
  for (int i=0; i<N_target; i++) {
    S_ex[3*i  ] = Target[i].u_ex;
    S_ex[3*i+1] = Target[i].v_ex;
    S_ex[3*i+2] = Target[i].w_ex;
  }

  vector<double> rho(3,0), DEL(3,0);
  set_rho_del(h, rho, DEL);

  // Compute the regularized SL and DL for each delta  
  vector<double> SL_del1(N_target3,0), SL_del2(N_target3,0), SL_del3(N_target3,0);
  vector<double> DL_del1(N_target3,0), DL_del2(N_target3,0), DL_del3(N_target3,0);
  
  eval_Stokes_SL_offSurf_3del(h, DEL[0], DEL[1], DEL[2],
			      N_quad, Surfc,
			      N_target, Target,
			      SL_del1,
			      SL_del2,
			      SL_del3,
			      1);  // use subtraction
  
  eval_Stokes_DL_offSurf_3del(h, DEL[0], DEL[1], DEL[2],
			      N_quad, Surfc,
			      N_target, Target,
			      DL_del1,
			      DL_del2,
			      DL_del3); 
  
  // Set up and solve 3x3 system for S at each target point
  vector<double> SL_comp(N_target3,0), DL_comp(N_target3,0);
  solve_3by3(h, rho, DEL,
	     N_target, Target,
	     SL_del1,
	     SL_del2,
	     SL_del3,
	     SL_comp);

  solve_3by3(h, rho, DEL,
	     N_target, Target,
	     DL_del1,
	     DL_del2,
	     DL_del3,
	     DL_comp);

  for (int k=0; k<N_target3; k++)  S_comp[k] = SL_comp[k] + DL_comp[k];  
					
  // display norms of solution and error
  cout_norms(N_target, S_ex, S_comp);

}

//*********************************************************************

void set_rho_del(double h,
		 vector<double>& rho,
		 vector<double>& DEL) {
  
  if (del_h == 1) {
    rho[0] = 3.0; rho[1] = 4.0; rho[2] = 5.0;
    //rho[0] = 2.0; rho[1] = 3.0; rho[2] = 4.0;
  }
  else {
    rho[0] = 3.0; rho[1] = 4.0; rho[2] = 5.0;
    //rho[0] = 2.0; rho[1] = 3.0; rho[2] = 4.0;
  }
  
  for (int i=0; i<3; i++) {
    if (del_h == 1) {
      DEL[i] = rho[i] * h;
    }
    else {
      DEL[i] = rho[i] * pow(1.0/64.0,1.0/5.0) * pow(h,4.0/5.0);
    }
  }
}

//*********************************************************************

void solve_3by3(double h, const vector<double>& rho, const vector<double>& DEL,
		int N_target, const vector<Target_point>& Target,
		const vector<double>& S_d1,
		const vector<double>& S_d2,
		const vector<double>& S_d3,
		vector<double>& S_comp) {

  for (int k=0; k<N_target; k++) {
    double b = Target[k].b;
    
    vector<double> lam(3,0);    
    for (int i=0; i<3; i++)  lam[i] = b / DEL[i];
    
    vector<double> Line(3,0), RHS(3,0), SOL(3,0);
    vector< vector<double> > Mat(3,Line);

    for (int i=0; i<3; i++) {
      Mat[i][0] = 1.0;
      Mat[i][1] = rho[i] * I0(lam[i]);
      Mat[i][2] = rho[i] * rho[i] * rho[i] * I2(lam[i]);
    }

    vector<int> p(3,0);
    int result = LU_factorization(Mat, p);

    // first component
    RHS[0] = S_d1[3*k];
    RHS[1] = S_d2[3*k];
    RHS[2] = S_d3[3*k];
    SOL = solveByLU(Mat, p, RHS);

    S_comp[3*k] = SOL[0];

    // second component
    RHS[0] = S_d1[3*k+1];
    RHS[1] = S_d2[3*k+1];
    RHS[2] = S_d3[3*k+1];
    SOL = solveByLU(Mat, p, RHS);

    S_comp[3*k+1] = SOL[0];

    // third component
    RHS[0] = S_d1[3*k+2];
    RHS[1] = S_d2[3*k+2];
    RHS[2] = S_d3[3*k+2];
    SOL = solveByLU(Mat, p, RHS);

    S_comp[3*k+2] = SOL[0];
  }
}
  
//*********************************************************************

double I0(double lam) {
  return exp(-lam*lam)/rootPI - abs(lam) * erfc(abs(lam));
}

double I2(double lam) {
  return 2.0/3.0 * ((0.5-lam*lam) * exp(-lam*lam)/rootPI
		    + abs(lam)*lam*lam * erfc(abs(lam)));
}

//*********************************************************************

void cout_norms(int n,
		const vector<double>& S_ex,
		const vector<double>& S_comp) {
		
  double S_max, S_l2;
  int ind_S_max;  
  norms_3d_vector(n, S_ex, S_max, ind_S_max, S_l2);
  cout << "max exact value = " << S_max << "   l2 exact values = " << S_l2 << endl;
  
  double err_max, err_l2;
  int ind_max;  
  error_3d_vector(n, S_comp, S_ex, err_max, ind_max, err_l2);
  cout << "max error = " << err_max << "   l2 error = " << err_l2 << endl;

}
