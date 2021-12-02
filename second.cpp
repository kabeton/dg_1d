#include <cmath>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include "dgm.h"

int main() {
  int k = 4, N = 15;
  double t0 = 0., t = 15., dt = 1.;
  double x0 = 0., xf = 5.;
  auto u0 = [xf](double x) {return sin(4*M_PI*x/xf);};
  auto u1 = [](double x) {return 1;};
  dgm solver(k, N, t0, t, dt, x0, xf);
  solver.initial_condition(u0);
  // solver.solve((t - t0)/720);
  solver.integrate(1e-1);
  return 0;
}