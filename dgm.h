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
  int k, N, v;
  double t0, t, dt, nt; 
  double a;
  Eigen::VectorXd x;
  vect G;
  Eigen::MatrixXd Minv, K;
  std::vector<vect> U;
  bool v_def;

public:
  dgm(int _k, int _N, double T0, double T, double dT, double x0, double xf);
  void initial_condition(func u0, double _a = 1);
  double g(int i);
  double phi(double x, int j);
  double scalar_product_basis(func f1, int j);
  void update_flux();
  double reconstruct(double x, Eigen::VectorXd &coef);
  void write(std::string name);
  void time_step(double h);
  void time_step_mp(double h);
  void integrate(double eps);
  void solve(double h);
  void debug_drawt(std::string name);
  double sol_dist(std::vector<vect> &U1, std::vector<vect> &U2);
  Eigen::VectorXd method(const Eigen::VectorXd &u);
  ~dgm(){};
};

#endif