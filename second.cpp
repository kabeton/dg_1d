#include <cmath>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include "dgm.h"

int main() {
  int k = 5, N = 10;
  double T0 = 0., T = 5., dt = 1.;
  double x0 = 0., xf = 5.;
  auto u0 = [xf](double x) {return sin(4*M_PI*x/xf);};
  auto u1 = [](double x) {return x;};
  dgm solver(k, N, T0, T, dt, x0, xf);
  solver.initial_condition(u1);
  //solver.solve((T - T0)/10);
  return 0;
}