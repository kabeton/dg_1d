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
  G.assign(N, Eigen::VectorXd::Zero(k));

  // precompute matrices
  Minv = Eigen::MatrixXd::Zero(k, k);
  K = Eigen::MatrixXd::Zero(k, k);

  for(int i = 0; i < k; i++) {
    for(j = 0; j < k; j++) {
      if(i == j) Minv(i, i) = 1./gl_quad0d(j, j, k);
      K(i, j) = gl_quad1d(j, i, k);
    }
  }
  v_def = false;
  k_def = false;
}

void dgm::initial_condition(func u0, double _a) {
  a = _a;

  // get initial coefficients from scalar product with basis function
  for(v = 0; v < N; v++) {
    v_def = true;
    U[v].push_back(Eigen::VectorXd::Zero(k));
    for(j = 0; j < k; j++) {
      k_def = true;
      U[v][0][j] = 0.25*scalar_product_basis(u0); // idk wtf this 0.25 does but it seems right with it
    }
    std::cout << u0(x(v)) << " " << reconstruct(x(v), U[v][0]) << std::endl;
    std::cout << reconstruct(x(v), U[v][0]) - u0(x(v)) << std::endl;
    std::cout << "-----" << std::endl;
  }

  v_def = false;
  k_def = false;
  update_flux();
}

double dgm::scalar_product_basis(func f1) {
  if(!(v_def)) {
    std::cout << "warning: " << v_def << std::endl;
    return 0;
  }

  auto f1_stdint = [this, &f1](double y) {return f1((x(v + 1) - x(v))*y/2 + (x(v + 1) + x(v))/2);};
  auto f2_stdint = [this](double y) {return phi((x(v + 1) - x(v))*y/2 + (x(v + 1) + x(v))/2);};

  Eigen::ArrayXd roots(k);
  roots = get_roots(k);
  Eigen::ArrayXd weights(k);
  weights = get_weights(k);
  double res = 0;

  for(int i = 0; i < k; i++) {
    res += weights(i)*f1_stdint(roots(i))*f2_stdint(roots(i));
  }

  return res;
}

void dgm::update_flux() { 
  for(v = 0; v < N; v++) {
    v_def = true;
    for(int l = 0; l < k; l++) {
      k_def = true;
      if(l % 2 == 0) G[v][l] = g(v + 1) + g(v);
      if(l % 2 != 0) G[v][l] = g(v + 1) - g(v);
    }
  }
  v_def = false;
  k_def = false;
}

double dgm::reconstruct(double p, Eigen::VectorXd coef) {
  double res = 0;
  int nc = 0;
  if(k_def == false) {k_def = true; nc++;}
  for(j = 0; j < k; j++) {
    res += coef[j]*phi(p);
  }
  if(nc) k_def = false;
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
      std::string oname = std::to_string(nt) + ".vtk";
      write(oname);
      nt += dt;
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
  std::cout << errfin << std::endl;
  int iters = 1;
  while(errfin > eps) {
    std::cout << "solving with h = " << h << " error = " << errfin << std::endl;
    Utmp = U;
    h /= 2;
    solve(h);
    iters++;
    errfin = sol_dist(U, Utmp);
  }
  errfin = sol_dist(U, Utmp);
  std::cout << "finished in " << iters << " iterations, final error = " << errfin << std::endl;
}

double dgm::sol_dist(std::vector<vect> U1, std::vector<vect> U2) {
  double t = 0, tm = -DBL_MAX;
  int total_points = U2[0].size();
  for(int l = 0; l < N; l++) {
    for(int i = 0; i < total_points; i += total_points/10) {
      t = (U1[l][i*2] - U2[l][i]).array().abs().maxCoeff();
      if(t > tm) tm = t;
    }
  }
  return tm;
}

Eigen::VectorXd dgm::method(Eigen::VectorXd u) { //TODO: fix flux
  double ci = 2.*a/(x(v + 1) - x(v));
  return ci*((Minv*K)*u) - ci*(Minv*G[v]);
}

double dgm::g(int i) {
  if(i == 0) {
    return reconstruct(x(Eigen::last), U.back().back());
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
  for(v = 0; v < N + 1; v++) {
    vtk_file << x(v) << " " << 0.0 << " " << 0.0 << "\n";
  }
	vtk_file << "\nPOINT_DATA " << N + 1 << "\n";
	vtk_file << "SCALARS v FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
  v_def = true;
  for(v = 0; v < N; v++) {
    vtk_file << reconstruct(x(v), U[v].back()) << "\n";
  }
  vtk_file.close();
  v_def = false;
}

double dgm::phi(double p) {
  if(!(v_def && k_def)) {
    std::cout << "warning: " << v_def << " " << k_def << std::endl;
  }
  double hv = abs(x(v + 1) - x(v));
  return 1./sqrt(hv)*sqrt(2*j + 1)*P(j, 2*(p - x(v))/hv - 1);
}