#include "dgm.h"
#include <float.h>
#include <fstream>
#include "quad.h"

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

void dgm::initial_condition(std::function<double(double)> u0, double _a) {
  a = _a;

  // get initial coefficients from scalar product with basis function
  for(v = 0; v < N; v++) {
    U[v].push_back(Eigen::VectorXd::Zero(k));
    auto u0_stdint = [this, &u0](double y) {return u0((x(v + 1) - x(v))*y/2 + (x(v + 1) - x(v))/2);};
    for(int j = 0; j < k; j++) {
      U[v][0][j] = gl_quadcoef(u0, k, j);
    }
  }

  update_flux();
}

void dgm::update_flux() { //TODO: fix this
  for(int v = 0; v < k; v++) {
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
    Eigen::VectorXd tmp = U[v].back();
    f1 = method(tmp);
    f2 = method(tmp + h*f1/2);
    f3 = method(tmp + h*f2/2);
    f4 = method(tmp + h*f3);
    U[v].push_back(tmp + h/6*(f1 + 2*f2 + 2*f3 + f4));
  }
}

void dgm::solve(double h) {
  rt = t0;
  nt = t0 + dt;
  while(rt < t) {
    time_step(h);
    update_flux();
    if(rt <= nt && rt + h >= nt) {
      write("out.vtk");
    }
    rt += h;
  }
}

void dgm::integrate(double eps) {
  double h = (t - t0)/10;
  solve(h);
  std::vector<vect> Utmp = U;
  h /= 2;
  solve(h);
  errfin = sol_dist(U, Utmp);
  int iters = 1;
  while(errfin < eps) {
    std::cout << "solving with h = " << h << " error = " << errfin << std::endl;
    Utmp = U;
    h /= 2;
    solve(h);
    iters++;
  }
  errfin = sol_dist(U, Utmp);
  std::cout << "finished in " << iters << " iterations, final error = " << errfin << std::endl;
}

double dgm::sol_dist(std::vector<vect> U1, std::vector<vect> U2) {
  double t = 0, tm = -DBL_MAX;
  int total_points = U1[0].size();
  for(int v = 0; v < N; v++) {
    for(int i = 0; i < total_points; i+= total_points/10) {
      t = (U1[v][i] - U2[v][i*2]).array().abs().maxCoeff();
      if(t > tm) tm = t;
    }
  }
  return tm;
}

Eigen::VectorXd dgm::method(Eigen::VectorXd u) { //TODO: fix flux
  double ci = 2.*a/(x(v + 1) - x(v));
  std::cout << (Minv*K*u).rows() << " " << (Minv*K*u).cols() << std::endl;
  std::cout << (Minv*G[v]).rows() << " " << (Minv*G[v]).cols() << std::endl;
  return ci*((Minv*K)*u) - ci*(Minv*G);
}

double dgm::g(int i) {
  if(i == 0) {
    return reconstruct(x[N], U[N - 1].back());
  } else {
    return reconstruct(x(i), U[i - 1].back());
  }
}

void dgm::write(std::string name) {
  std::cout << "writing " << name << std::endl;
  std::string nname = "output/" + name;
  std::ofstream vtk_file;
  vtk_file.open(nname.c_str(), std::ios::out);
  vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
  vtk_file << "DATASET POLYDATA\nPOINTS " << N + 1 << " float\n";
  for(int i = 0; i < N + 1; i++) {
    vtk_file << x(i) << " " << 0.0 << " " << 0.0 << "\n";
  }
	vtk_file << "\nPOINT_DATA " << N + 1 << "\n";
	vtk_file << "SCALARS v FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
  for(int i = 0; i < N; i++) {
    vtk_file << reconstruct(x(i), U[i].back()) << "\n";
    vtk_file << reconstruct(x(i + 1), U[i].back()) << "\n";
  }
  vtk_file.close();
}