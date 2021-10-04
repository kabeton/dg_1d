#ifndef DGM_H
#define DGM_H

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include "quad.h"

typedef std::vector<Eigen::VectorXd> vect;

class dgm {
private:
  int k, N, v;
  double t0, t, dt; 
  double rt, a;
  Eigen::VectorXd x, G;
  Eigen::MatrixXd Minv, K;
  std::vector<vect> U;

public:
  dgm(int _k, int _N, double T0, double T, double dT, double x0, double xf);
  void initial_condition(std::function<double(double)> u0, double _a = 1);
  double g(int i);
  double reconstruct(double x, Eigen::VectorXd coef);
  void write(std::string name);
  void time_step(double h);
  void integrate(double eps);
  Eigen::VectorXd method(Eigen::VectorXd u);
  ~dgm();
};

#endif