#include <cmath>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include "dgm.h"

int main() {
  int k = 6, N = 10;
  double t0 = 0., t = 30., dt = 1.;
  double x0 = 0., xf = 20.;
  auto u0 = [xf](double x) {return sin(4*M_PI*x/xf);};
  auto u1 = [](double x) {return x > 5 && x < 15 ? 3 : 0;};
  dgm solver(k, N, t0, t, dt, x0, xf);
  solver.initial_condition(u0);
  solver.intgr();
  return 0;
}