#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

static const int N_gridlines = 32;
static const int del_h = 2; // 1:    del=rho*h,       rho = (3,4,5)
                            // else: del=rho*h^(4/5), rho=(2,3,4)*(1/64)^(1/5)


static const double PI = 3.14159265358979323846;
static const double rootPI = sqrt(PI);
static const double theta = 70.0*PI/180.0; //angle in Beale et al (2.2)
static const double tol = 1e-14; //tolerance in the search of quadrature points (Newton's/bisection method)

struct Surf_point
{
  Surf_point() : x(3,0), f(3,0), g(3,0), Nrml(3,0), Area(0) {}
  vector<double> x;
  vector<double> f;  // usually used for Stokeslet density
  vector<double> g;  // usually used for stresslet density
  vector<double> Nrml;
  double Area;
};

struct Target_point
{
  Target_point() : flag(0), x(0), y(0), z(0), nrst_x(0), nrst_y(0), nrst_z(0), b(0), nrst_Nrml(3,0), nrst_f(3,0), nrst_g(3,0), S1(0), S2(0), S3(0), T1(0), T2(0), T3(0), u(0), v(0), w(0), u_ex(0), v_ex(0), w_ex(0) {}
  int flag;
  double S1, S2, S3;
  double T1, T2, T3;
  double x, y, z;
  double nrst_x, nrst_y, nrst_z, b;
  vector<double> nrst_Nrml, nrst_f, nrst_g;
  double u, v, w;
  double u_ex, v_ex, w_ex;
};

double phi(const vector<double>& x);
double Dphi(int i, const vector<double>& x);
void D2phi(const vector<double>& x, double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23, double& phi31,
	   double& phi32, double& phi33);

void H_gauss_opt2(double r, double d, double& H1, double& H2);
void s_opt2(double r, double& s1, double& s2);

void Stresslet_orig(double r, double d, double& s2, double& s3);

vector<double> stokeslet_density(const vector<double>& pt,
				 const vector<double>& nl);
vector<double> stresslet_density(const vector<double>& pt,
				 const vector<double>& nl);

void Exact_solution(int N_target, vector<Target_point>& Target);

double distance(const vector<double>& x,
		const vector<double>& y,
		const vector<double>& normal);

double dot_product(const vector<double>& x,
		   const vector<double>& y);

vector<double> cross_product(const vector<double>& x,
			     const vector<double>& y);

void orthogonalize(vector<double>& x,
		   vector<double>& y);

int LU_factorization(vector< vector<double> >& A,
		     vector<int>& p);

vector<double> solveByLU(const vector< vector<double> >& A,
			 const vector<int>& p,
			 const vector<double>& b);

void initialize_vector(int n, vector<double>& vec);

void error_vector(int n,
		  const vector<double>& x_comp,
		  const vector<double>& x_ex,
		  double &err_max,
		  int &ind_max,
		  double &err_l2);

void error_3d_vector(int n,
		     const vector<double>& x_comp,
		     const vector<double>& x_ex,
		     double &err_max,
		     int &ind_max,
		     double &err_l2);

void norms_3d_vector(int n,
		     const vector<double>& x,
		     double &x_max,
		     int &ind_max,
		     double &x_l2);
