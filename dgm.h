#ifndef DGM_H
#define DGM_H

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>

typedef std::vector<Eigen::VectorXd> vect;
typedef std::function<double(double)> func;

class dgm {
private:
  int k, N, v, j;
  double t0, t, dt, nt; 
  double rt, a;
  double errfin;
  Eigen::VectorXd x;
  vect G;
  Eigen::MatrixXd Minv, K;
  std::vector<vect> U;
  bool v_def, k_def;

public:
  dgm(int _k, int _N, double T0, double T, double dT, double x0, double xf);
  void initial_condition(func u0, double _a = 1);
  double g(int i);
  double phi(double x);
  double scalar_product_basis(func f1);
  void update_flux();
  double reconstruct(double x, Eigen::VectorXd coef);
  void write(std::string name);
  void time_step(double h);
  void integrate(double eps);
  void solve(double h);
  double sol_dist(std::vector<vect> U1, std::vector<vect> U2);
  Eigen::VectorXd method(Eigen::VectorXd u);
  ~dgm(){};
};

#endif