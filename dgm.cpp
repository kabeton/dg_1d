#include "dgm.h"

dgm::dgm(int _k, int _N, double T0, double T, double dT, double x0, double xf) {
  // grid parameters
  k = _k;
  N = _N;
  t0 = T0;
  t = T;
  dt = dT;
  x = Eigen::VectorXd::Zero(N + 1);
  double step = (xf - x0)/N;
  for(int i = 0; i <= N; i++) {
    x[i] = step*i;
  }

  // allocate space for solution matrix and flux vector
  U.resize(N);
  G = Eigen::VectorXd::Zero(N);

  // precompute matrices
  Minv = Eigen::MatrixXd::Zero(k, k);
  K = Eigen::MatrixXd::Zero(k, k);

  for(int i = 0; i < k; i++) {
    for(int j = 0; j < k; j++) {
      if(i == j) Minv(i, i) = 1./gl_quad0d(j, j, k);
      K(i, j) = gl_quad1d(j, i, k);
    }
  }
}

void dgm::initial_condition(std::function<double(double)> u0, double _a = 1) {
  a = _a;

  // get initial coefficients from scalar product with basis function
  for(v = 0; v < N; v++) {
    U[v].push_back(Eigen::VectorXd::Zero(k));
    auto u0_stdint = [this, &u0](double y) {return u0((x(v + 1) - x(v))*y/2 + (x(v + 1) - x(v))/2);};
    for(int j = 0; j < k; j++) {
      U[v][0][j] = gl_quadcoef(u0, k, j);
    }
  }

  // build flux vector
  for(int v = 0; v < N; v++) {
    if(v % 2 == 0) G[v] = g(v + 1) + g(v);
    if(v % 2 != 0) G[v] = g(v + 1) - g(v);
  }
}

double dgm::reconstruct(double x, Eigen::VectorXd coef) {
  double res = 0;
  for(int i = 0; i < k; i++) {
    res += coef[i]*P(i, x);
  }
  return res;
}

void dgm::time_step(double h) {
  for(v = 0; v < N; v++) {
    Eigen::VectorXd f1, f2, f3, f4;
  }
}

void dgm::integrate(double eps) {
  double h = (t - t0)/10;
  rt = t0;
  while(rt < t) {
    time_step(h);
    rt += h;
  }
}

Eigen::VectorXd dgm::method(Eigen::VectorXd u) {
  double ci = 2.*a/(x(v + 1) - x(v));
  return ci*Minv*K*u - ci*Minv*G[v];
}

double dgm::g(int i) {
  if(i == 0) {
    return 0;
  } else if(i == N - 1) {
    return 0;
  } else {
    return 0;
  }
}