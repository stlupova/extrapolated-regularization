#include <iostream>
#include <vector>

#include "utilities.h"

using namespace std;

//*********************************************************************

double phi(const vector<double>& x) {

  vector<double> x1(3,0), x2(3,0), x3(3,0), x4(3,0);
  x1[0] = sqrt(3.0)/3.0;  x1[1] = 0.0;  x1[2] =-sqrt(6.0)/12.0;
  x2[0] =-sqrt(3.0)/6.0;  x2[1] = 0.5;  x2[2] =-sqrt(6.0)/12.0;
  x3[0] =-sqrt(3.0)/6.0;  x3[1] =-0.5;  x3[2] =-sqrt(6.0)/12.0;
  x4[0] = 0.0;            x4[1] = 0.0;  x4[2] = sqrt(6.0)/4.0;
  
  double r2 = 0.5*0.5;
  double c = 0.6;

  double a1, a2, a3, a4;
  a1 = (x[0]-x1[0])*(x[0]-x1[0]) + (x[1]-x1[1])*(x[1]-x1[1]) + (x[2]-x1[2])*(x[2]-x1[2]); 
  a2 = (x[0]-x2[0])*(x[0]-x2[0]) + (x[1]-x2[1])*(x[1]-x2[1]) + (x[2]-x2[2])*(x[2]-x2[2]); 
  a3 = (x[0]-x3[0])*(x[0]-x3[0]) + (x[1]-x3[1])*(x[1]-x3[1]) + (x[2]-x3[2])*(x[2]-x3[2]); 
  a4 = (x[0]-x4[0])*(x[0]-x4[0]) + (x[1]-x4[1])*(x[1]-x4[1]) + (x[2]-x4[2])*(x[2]-x4[2]); 

  double val = c - exp(-a1/r2) - exp(-a2/r2) - exp(-a3/r2) - exp(-a4/r2);
  return val;
}

//*********************************************************************

// D_phi/D_x(i): i-th derivative of phi

double Dphi(int i, const vector<double>& x) {

  vector<double> x1(3,0), x2(3,0), x3(3,0), x4(3,0);
  x1[0] = sqrt(3.0)/3.0;  x1[1] = 0.0;  x1[2] =-sqrt(6.0)/12.0;
  x2[0] =-sqrt(3.0)/6.0;  x2[1] = 0.5;  x2[2] =-sqrt(6.0)/12.0;
  x3[0] =-sqrt(3.0)/6.0;  x3[1] =-0.5;  x3[2] =-sqrt(6.0)/12.0;
  x4[0] = 0.0;            x4[1] = 0.0;  x4[2] = sqrt(6.0)/4.0;
  
  double r2 = 0.5*0.5;

  double a1, a2, a3, a4;
  a1 = (x[0]-x1[0])*(x[0]-x1[0]) + (x[1]-x1[1])*(x[1]-x1[1]) + (x[2]-x1[2])*(x[2]-x1[2]); 
  a2 = (x[0]-x2[0])*(x[0]-x2[0]) + (x[1]-x2[1])*(x[1]-x2[1]) + (x[2]-x2[2])*(x[2]-x2[2]); 
  a3 = (x[0]-x3[0])*(x[0]-x3[0]) + (x[1]-x3[1])*(x[1]-x3[1]) + (x[2]-x3[2])*(x[2]-x3[2]); 
  a4 = (x[0]-x4[0])*(x[0]-x4[0]) + (x[1]-x4[1])*(x[1]-x4[1]) + (x[2]-x4[2])*(x[2]-x4[2]); 

  double d1, d2, d3, d4;
  d1 = x[i] - x1[i];
  d2 = x[i] - x2[i];
  d3 = x[i] - x3[i];
  d4 = x[i] - x4[i];

  double val = 2.0/r2* (exp(-a1/r2)*d1 + exp(-a2/r2)*d2 + exp(-a3/r2)*d3 + exp(-a4/r2)*d4);
  return val;
}

//*********************************************************************

// Second derivatives of phi
void D2phi(const vector<double>& x, double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23, double& phi31,
	   double& phi32, double& phi33) {
  
  vector<double> x1(3,0), x2(3,0), x3(3,0), x4(3,0);
  x1[0] = sqrt(3.0)/3.0;  x1[1] = 0.0;  x1[2] =-sqrt(6.0)/12.0;
  x2[0] =-sqrt(3.0)/6.0;  x2[1] = 0.5;  x2[2] =-sqrt(6.0)/12.0;
  x3[0] =-sqrt(3.0)/6.0;  x3[1] =-0.5;  x3[2] =-sqrt(6.0)/12.0;
  x4[0] = 0.0;            x4[1] = 0.0;  x4[2] = sqrt(6.0)/4.0;
  
  double r2 = 0.5*0.5;

  vector<double> d1 = x; d1[0]-=x1[0]; d1[1]-=x1[1]; d1[2]-=x1[2];
  vector<double> d2 = x; d2[0]-=x2[0]; d2[1]-=x2[1]; d2[2]-=x2[2];
  vector<double> d3 = x; d3[0]-=x3[0]; d3[1]-=x3[1]; d3[2]-=x3[2];
  vector<double> d4 = x; d4[0]-=x4[0]; d4[1]-=x4[1]; d4[2]-=x4[2];
  double a1, a2, a3, a4;
  a1 = exp(-dot_product(d1,d1)/r2);
  a2 = exp(-dot_product(d2,d2)/r2);
  a3 = exp(-dot_product(d3,d3)/r2); 
  a4 = exp(-dot_product(d4,d4)/r2); 

  phi11 = 2.0/r2*(a1*(1.0-2.0/r2*d1[0]*d1[0]) + a2*(1.0-2.0/r2*d2[0]*d2[0]) +
		  a3*(1.0-2.0/r2*d3[0]*d3[0]) + a4*(1.0-2.0/r2*d4[0]*d4[0]));
  phi22 = 2.0/r2*(a1*(1.0-2.0/r2*d1[1]*d1[1]) + a2*(1.0-2.0/r2*d2[1]*d2[1]) +
		  a3*(1.0-2.0/r2*d3[1]*d3[1]) + a4*(1.0-2.0/r2*d4[1]*d4[1]));
  phi33 = 2.0/r2*(a1*(1.0-2.0/r2*d1[2]*d1[2]) + a2*(1.0-2.0/r2*d2[2]*d2[2]) +
		  a3*(1.0-2.0/r2*d3[2]*d3[2]) + a4*(1.0-2.0/r2*d4[2]*d4[2]));
  
  phi12 = -4.0/(r2*r2)*(a1*d1[0]*d1[1] + a2*d2[0]*d2[1] +
			a3*d3[0]*d3[1] + a4*d4[0]*d4[1]);
  phi13 = -4.0/(r2*r2)*(a1*d1[0]*d1[2] + a2*d2[0]*d2[2] +
			a3*d3[0]*d3[2] + a4*d4[0]*d4[2]);
  phi23 = -4.0/(r2*r2)*(a1*d1[1]*d1[2] + a2*d2[1]*d2[2] +
			a3*d3[1]*d3[2] + a4*d4[1]*d4[2]);

  phi21 = phi12;
  phi31 = phi13;
  phi32 = phi23;
}

//*********************************************************************
//*********************************************************************

// Gaussian regularization Option 2
void H_gauss_opt2(double r, double d, double& H1, double& H2) {
  if (r < 1e-14) {
    H1 = 2.0/rootPI/d;
    H2 = 4.0/3.0/rootPI/(d*d*d);
  }
  else {
    double s1, s2;
    s_opt2(r/d,s1,s2);
    H1 = s1/r;
    H2 = s2/(r*r*r);
  }
}
//*********************************************************************
void s_opt2(double r, double& s1, double& s2) {
  s1 = erf(r);
  s2 = s1 - 2.0*r*exp(-r*r)/rootPI; 
}
//*********************************************************************

void Stresslet_orig(double r, double d, double& s2, double& s3) {
  if (r < 1e-14) {
    s2 = 4.0/3.0/rootPI/(d*d*d);
    s3 = s2 * 2.0/5.0/(d*d);
  }
  else {
    double r2 = r * r;
    double rd = r/d;
    double rd2 = rd * rd;
    double c1 = erf(rd);
    double c2 = 2.0 * rd * exp(-rd2) / rootPI; 
    s2 = ( c1 - c2 ) / (r2*r);
    s3 = ( c1 - c2 * (2.0/3.0*rd2 + 1.0) ) / (r2*r2*r);
  }
}

//*********************************************************************
//*********************************************************************

vector<double> stokeslet_density(const vector<double>& pt,
				 const vector<double>& nl) {
  vector<double> f(3,0);

  // TEST: Combined Single and Double Layer
  // Stokeslet density is the boundary value of traction = sigma*n
  // Putting a Stokeslet strength 4pi*(1,0,0) at (2,0,0)
  
  vector<double> b(3,0), ctr(3,0), dx(3,0);
  b[0] = 4.0 * PI;
  ctr[0] = 2.0; 
  dx[0] = pt[0] - ctr[0];
  dx[1] = pt[1] - ctr[1];
  dx[2] = pt[2] - ctr[2];
  double dist2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
  double dx_dot_b = dx[0]*b[0] + dx[1]*b[1] + dx[2]*b[2];
  double dx_dot_n = dx[0]*nl[0] + dx[1]*nl[1] + dx[2]*nl[2]; 
  double temp = -6.0*dx_dot_b*dx_dot_n/(dist2*dist2*sqrt(dist2)*8.0*PI);
  
  f[0] = dx[0]*temp;
  f[1] = dx[1]*temp;
  f[2] = dx[2]*temp;
  
  return f;
}

//*********************************************************************

vector<double> stresslet_density(const vector<double>& pt,
				 const vector<double>& nl) {
  vector<double> g(3,0);

  // TEST: Combined Single and Double Layer
  // Stresslet density is the boundary value of velocity
  // Putting a Stokeslet strength 4pi*(1,0,0) at (2,0,0)
  
  vector<double> b(3,0), ctr(3,0), dx(3,0);
  b[0] = 4.0 * PI;
  ctr[0] = 2.0; 
  dx[0] = pt[0] - ctr[0];
  dx[1] = pt[1] - ctr[1];
  dx[2] = pt[2] - ctr[2];
  double dist = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
  double dx_dot_b = dx[0]*b[0]+dx[1]*b[1]+dx[2]*b[2];
  double temp1 = 1.0/(8.0*PI*dist);
  double temp2 = dx_dot_b/(dist*dist*dist*8.0*PI);
  
  g[0] = b[0]*temp1 + dx[0]*temp2;
  g[1] = b[1]*temp1 + dx[1]*temp2;
  g[2] = b[2]*temp1 + dx[2]*temp2;

  return g;
}

//*********************************************************************
//*********************************************************************

void Exact_solution(int N_target, vector<Target_point>& Target) {
  
  for (int i=0; i<N_target; i++) {
    
    // TEST: Combined Single and Double Layer
      
    vector<double> temp_x(3,0);
    temp_x[0] = Target[i].x;
    temp_x[1] = Target[i].y;
    temp_x[2] = Target[i].z;
    
    double loc = phi(temp_x);
    double xi = 0.0; // target point is outside the boundary
    if ( abs(loc) < 1.0e-12 ) { // target point is on the boundary
      xi = 0.5; 
    }
    else if ( loc < 0.0 ) { // target point is inside the boundary
      xi = 1.0;
    }
    
    vector<double> temp_u = stresslet_density(temp_x,temp_x);
    Target[i].u_ex = xi*temp_u[0];
    Target[i].v_ex = xi*temp_u[1];
    Target[i].w_ex = xi*temp_u[2];     
  }
}

//*********************************************************************
//*********************************************************************

// Signed distance between two points in 3D

double distance(const vector<double>& x,
		const vector<double>& y,
		const vector<double>& normal) {
  
  return (x[0]-y[0])*normal[0] + (x[1]-y[1])*normal[1] + (x[2]-y[2])*normal[2];
}

//*********************************************************************

double dot_product(const vector<double>& x,
		   const vector<double>& y) {
  
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

//*********************************************************************

vector<double> cross_product(const vector<double>& x,
			     const vector<double>& y) {
  
  vector<double> cross_prod(3,0);
  cross_prod[0] = x[1]*y[2] - x[2]*y[1];
  cross_prod[1] = x[2]*y[0] - x[0]*y[2];
  cross_prod[2] = x[0]*y[1] - x[1]*y[0];
  return cross_prod;
}

//*********************************************************************

// assumes neither vector is 0 length
void orthogonalize(vector<double>& x,
		   vector<double>& y) {

  double x_dot_x = dot_product(x,x);
  double x_dot_y = dot_product(x,y);

  double coeff = x_dot_y / x_dot_x;

  y[0] -= coeff * x[0];
  y[1] -= coeff * x[1];
  y[2] -= coeff * x[2];

  double y_dot_y = dot_product(y,y);

  double norm_x = sqrt(x_dot_x);
  double norm_y = sqrt(y_dot_y);

  for (int i=0; i<3; i++) {
    x[i] = x[i] / norm_x;
    y[i] = y[i] / norm_y;
  }      
}

//*****************************************************************************
// Solve linear system through LU decomposition with pivoting.                
//                                                                           
// Factorize PA = LU with pivoting:                                           
//   The lower and upper triangular matrices are still stored in the original 
//   matrix and the permutation matrix "P" is stored in the vector "int *p".  
//******************************************************************************
 
int LU_factorization(vector< vector<double> >& A,
		     vector<int>& p) {
  int n = A.size();

  for (int j=0; j<n; j++) p[j] = j;

  for (int j=0; j<n; j++) {
    int k = j;
    double m = A[p[j]][j];
    // Search for maximum in this column
    for (int i=j+1; i<n; i++) {
      if (abs(A[p[i]][j]) > abs(m)) {
        m = A[p[i]][j];
        k = i;
      }
    }

    // "Swap" maximum row with current row (using the permutation vector)
    if (k != j) {
      int temp = p[j];
      p[j] = p[k];
      p[k] = temp;
    }

    double ajj = A[p[j]][j];
    if (abs(ajj) < 1.0e-15) {
      return (-1);
    }

    for (int i=j+1; i<n; i++) {
      double lij = A[p[i]][j] / ajj;
      A[p[i]][j] = lij; // lower triangular elements
      for (int k=j+1; k<n; k++) {
        A[p[i]][k] -= lij * A[p[j]][k]; // upper triangular elements
      }
    }
  }
  return 0;
}

//*****************************************************************************

vector<double> solveByLU(const vector< vector<double> >& A,
			 const vector<int>& p,
			 const vector<double>& b) {
  int n = A.size();
  vector<double> x(n,0);

  // Solve Ly=b by forward substitution
  x[0] = b[p[0]];
  for (int i=1; i<n; i++) {
    x[i] = b[p[i]];
    double rowsum = 0.0;
    for (int j=0; j<i; j++) {
      rowsum += A[p[i]][j] * x[j];
    }
    x[i] -= rowsum;
  }

  // Solve Ux=y by back substitution
  x[n-1] = x[n-1] / A[p[n-1]][n-1];
  for (int i=n-2; i>=0; i--) {
    double rowsum = 0.0;
    for (int j = n - 1; j > i; j--) {
      rowsum += A[p[i]][j] * x[j];
    }    
    x[i] = (x[i] - rowsum) / A[p[i]][i];
  }

  return x;
}

//*********************************************************************
//*********************************************************************

void initialize_vector(int n, vector<double>& vec) {
  for (int i=0; i<n; i++)  vec[i] = 0.0;
}

//*********************************************************************

void error_vector(int n,
		  const vector<double>& x_comp,
		  const vector<double>& x_ex,
		  double &err_max,
		  int &ind_max,
		  double &err_l2) {

  err_max = abs(x_comp[0] - x_ex[0]);
  ind_max = 0;
  err_l2 = err_max * err_max;
  
  for (int i=1; i<n; i++) {
    double err_temp = abs(x_comp[i] - x_ex[i]);
    if ( err_temp > err_max) {
      err_max = err_temp;
      ind_max = i;
    }
    err_l2 += err_temp * err_temp;
  }
  err_l2 = sqrt(err_l2 / n);

}
//*********************************************************************

void error_3d_vector(int n,
		     const vector<double>& x_comp,
		     const vector<double>& x_ex,
		     double &err_max,
		     int &ind_max,
		     double &err_l2) {

  double val = (x_comp[0] - x_ex[0]) * (x_comp[0] - x_ex[0])
             + (x_comp[1] - x_ex[1]) * (x_comp[1] - x_ex[1])
             + (x_comp[2] - x_ex[2]) * (x_comp[2] - x_ex[2]);
  
  err_max = sqrt(val);
  ind_max = 0;
  err_l2 = val;
  
  for (int i=1; i<n; i++) {
    double err_temp = (x_comp[3*i] - x_ex[3*i]) * (x_comp[3*i] - x_ex[3*i])
      + (x_comp[3*i+1] - x_ex[3*i+1]) * (x_comp[3*i+1] - x_ex[3*i+1])
      + (x_comp[3*i+2] - x_ex[3*i+2]) * (x_comp[3*i+2] - x_ex[3*i+2]);
    if ( sqrt(err_temp) > err_max) {
      err_max = sqrt(err_temp);
      ind_max = i;
    }
    err_l2 += err_temp;
  }
  err_l2 = sqrt(err_l2 / n);

}

//*********************************************************************

void norms_3d_vector(int n,
		     const vector<double>& x,
		     double &x_max,
		     int &ind_max,
		     double &x_l2) {

  double val = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];  

  x_max = sqrt(val);
  ind_max = 0;
  x_l2 = val;
  
  for (int i=1; i<n; i++) {
    double x_temp = x[3*i] * x[3*i] + x[3*i+1] * x[3*i+1] + x[3*i+2] * x[3*i+2];
    if ( sqrt(x_temp) > x_max) {
      x_max = sqrt(x_temp);
      ind_max = i;
    }
    x_l2 += x_temp;
  }
  x_l2 = sqrt(x_l2 / n);

}
