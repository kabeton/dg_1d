#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <vector>
#include <valarray>
#include <string>
#include "quad.h"
#include "rk.h"

typedef std::vector<Eigen::VectorXd> vect;

// numerical flux function
double g(std::vector<Eigen::VectorXd> u, int i);

double sol(vect &u, int i, double &x, int k) {
  double res = 0;
  for(int j = 0; j < k; i++) {
    res += u[i](j)*P(j, x);
  }
  return res;
}

double reconstruct(double x, int k, Eigen::VectorXd coef) {
  double res = 0;
  for(int i = 0; i < k; i++) {
    res += coef[i]*P(i, x);
  }
  return res;
}

void write_vtk(std::string name, int N, Eigen::VectorXd x) {
  std::cout << "writing " << name << std::endl;
  std::ofstream vtk_file;
  vtk_file.open(name.c_str(), std::ios::out);
  vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
  vtk_file << "DATASET POLYDATA\nPOINTS " << N + 1 << " float\n";
  for(int i = 0; i < N + 1; i++) {
    vtk_file << x(i) << " " << 0.0 << " " << 0.0 << "\n";
  }

}

Eigen::VectorXd valarray__to_vector(std::valarray<double> a);

std::valarray<double> vector_to_valarray(Eigen::VectorXd a);

int main() {  
  // 1-d grid, equidistant points
  // we assume that problem is initially defined on [0, L] 
  double L = 5.;
  double T0 = 0., T = 1.;
  int N;
  std::cout << "input number of sections N" << std::endl;
  std::cin >> N;
  Eigen::VectorXd x(N + 1);
  double step = L / N;
  for(int i = 0; i <= N; i++) {
    x(i) = step*i;
  }

  int k;
  // std::cout << "input max polynomial degree k" << std::endl;
  k = 5;
  // std::cin >> k; 
  k++;

  // precompute matrices
  Eigen::MatrixXd Minv = Eigen::MatrixXd::Zero(k, k);
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(k, k);

  for(int i = 0; i < k; i++) {
    for(int j = 0; j < k; j++) {
      if(i == j) Minv(i, i) = 1./gl_quad0d(j, j, k);
      K(i, j) = gl_quad1d(j, i, k);
      //std::cout << "(" << i << "," << j << ") " << Minv(i, j) << " " << K(i, j) << std::endl;
      //std::cout << "M_" << j << j << "= " << gl_quad0d(j, j, k) << std::endl;
    }
  }

  std::vector<vect> U(N);
  // U[0].resize(N);
  vect fluxv(N + 1, Eigen::VectorXd::Zero(k));

  // initial condition
  for(int i = 0; i < N; i++) {
    std::cout << "(" << x(i) << "," << x(i + 1) << ")\n";
    U[i].push_back(Eigen::VectorXd::Zero(k));
    auto u0 = [L, i, &x](double y) {return sin(4*M_PI*( (x(i + 1) - x(i))*y/2 + (x(i + 1) - x(i))/2 )/L);};
    for(int j = 0; j < k; j++) {
      U[i][0](j) = gl_quadcoef(u0, k, j);
      std::cout << U[i][0](j) << std::endl;
    }
  }
 
  // RK4 on each section
  for(int i = 0; i < N; i++) {
    double ci = 2./(x(i + 1) - x(i));

    Eigen::VectorXd G = fluxv[i];

    /* auto method = [&Minv, &K, G, ci](Eigen::VectorXd u)->Eigen::VectorXd {
      return ci*Minv*K*u - ci*Minv*G;
    }; */

    auto rhs = [&Minv, &G, &K, ci](double x, std::valarray<double> arr)->std::valarray<double> {
      Eigen::VectorXd temp1, temp2;
      temp1 = valarray__to_vector(arr);
      temp2 = ci*Minv*K*temp1 - ci*Minv*G;
      return vector_to_valarray(temp2);
    };

    rk4 Solver(T0, T, vector_to_valarray(U[0][i]), rhs);
    Solver.solve_precision(1e-3);
    U[i] = Solver.write_vector();
  }
 
  // save results

  return 0;
}