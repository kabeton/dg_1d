#include "dgm.h"
#include <float.h>
#include <fstream>
#include "quad.h"

#define UPWIND

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
      //if(i == j) Minv(i, i) = 1./gl_quad0d(i, i, k);
      if(i == j) Minv(i, i) = 1.;
      //K(i, j) = gl_quad1d(j, i, k);
      K(i, j) = sqrt((2*i + 1)*(2*j + 1))*gl_quad1d(j, i, k);
    }
  }
  v_def = false;
  k_def = false;
}

void dgm::initial_condition(func u0, double _a) {
  a = _a;
  double ratio = (double)k/N/4; // idk wtf is this and why we need it

  // get initial coefficients from scalar product with basis function
  for(v = 0; v < N; v++) {
    v_def = true;
    U[v].push_back(Eigen::VectorXd::Zero(k));
    for(j = 0; j < k; j++) {
      k_def = true;
      U[v][0][j] = ratio*scalar_product_basis(u0); 
    }

  }

  v_def = false;
  k_def = false;
  update_flux();
  write("dg_" + std::to_string(0));
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
    double hv = sqrt(x(v + 1) - x(v));
    for(j = 0; j < k; j++) {
      k_def = true;
      G[v][j] = (g(v + 1)*phi(x(v + 1)) - g(v)*phi(x(v)));
      //if(l == 0) G[v][l] = g(v + 1);
      //if(l == k - 1) G[v][l] = g(v + 1);
    } 
  }
  v_def = false;
  k_def = false;
}

double dgm::g(int i) {
  // indexing for x(i) works because we have v fixed and reconstruct 
  // maps x(v + 1) to 1 for any i so x(v + 1) returns the rightmost point 
  // of corresponding section
  #ifdef UPWIND
  if(i == v) {
    if(v == 0) {
      return reconstruct(x(v + 1), U.back().back());
    } else {
      return reconstruct(x(v + 1), U[v - 1].back());
    }
  } else {
    return reconstruct(x(v + 1), U[v].back());
  }
  #endif

  #ifdef CENTERED
    if(i == 0) {
      return (reconstruct(x(i + 1), U[N - 1].back()) + reconstruct(x(i), U[0].back()))/2;
    } else if(i == N) {
      return (reconstruct(x(i), U[i - 1].back()) + reconstruct(x(i - 1), U[0].back()))/2;
    } else {
      return (reconstruct(x(i + 1), U[i - 1].back()) + reconstruct(x(i), U[i].back()))/2;
    }
  #endif
}

double dgm::reconstruct(double p, Eigen::VectorXd coef) {
  double res = 0;
  int j_old = j;
  int nc = 0;
  if(k_def == false) {k_def = true; nc++;}
  for(j = 0; j < k; j++) {
    res += coef[j]*phi(p);
  }
  if(nc) k_def = false;
  j = j_old;
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
      std::string oname = "dg_" +  std::to_string((int)(nt*1e4));
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

Eigen::VectorXd dgm::method(Eigen::VectorXd u) { 
  double ci = a/(x(v + 1) - x(v));
  return ci*K*u + a*G[v];
}

void dgm::write(std::string name) {
  int r = 5;
  int N_out = N*r;
  std::cout << "writing " << name << std::endl;
  std::string nname = "output/" + name + ".vtk";
  std::string gpname = "gp/" + name + ".png";
  debug_drawt(gpname);
  std::ofstream vtk_file;
  vtk_file.open(nname.c_str(), std::ios::out);
  vtk_file << "# vtk DataFile Version 3.0\nVx data\nASCII\n\n";
  vtk_file << "DATASET POLYDATA\nPOINTS " << N_out + 1 << " float\n";
	for (int i = 0; i < N; i++)
	{
    for(int j = 0; j < r; j++) {
      double s = (x(i + 1) - x(i))/r;
		  vtk_file << x(i) + j*s << " " << 0.0 << " "  << 0.0 << "\n";
    }
	}
  vtk_file << x(N) << " " << 0. << " " << 0. <<  "\n";
	vtk_file << "\nLINES " << N_out << " " << N_out*3 << "\n";
	for (int i = 0; i < N_out; i++)
	{
		vtk_file << 2 << " " << i << " " << i+1 << "\n";
	}
	vtk_file << "\nPOINT_DATA " << N_out + 1 << "\n";
	vtk_file << "SCALARS v FLOAT\n";
	vtk_file << "LOOKUP_TABLE default\n";
  v_def = true;
  for(v = 0; v < N; v++) {
    for(int j = 0; j < r; j++) {
      double s = (x(v + 1) - x(v))/r;
      vtk_file << reconstruct(x(v) + j*s, U[v].back()) << "\n";
    }
  }
  v--;
  vtk_file << reconstruct(x(v + 1), U.back().back()) << "\n";
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

void dgm::debug_drawt(std::string name) {
  std::ofstream ofi("temp.dat");
  std::ofstream gfi("temp2.dat");
  std::ofstream ggfi("temp3.dat");
  v_def = true;
  int pp = 5;
  for(v = 0; v < N; v++) {
    double s = (x(v + 1) - x(v))/pp;
    for(int i = 0; i <= pp; i++) {
      double po = x(v) + i*s;
      ofi << po << " " << reconstruct(po, U[v].back()) << std::endl;
    }
    gfi << x(v) << " " << g(v) << std::endl;
    ggfi << x(v + 1) << " " << g(v + 1) << std::endl;
  }
  v_def = false;

  FILE *gp = popen("gnuplot", "w");
  fprintf(gp, "set terminal png size 1300, 800\n");
  fprintf(gp, "set output \"%s\"\n", name.c_str());
  fprintf(gp, "set grid x y\n");
  fprintf(gp, "show grid \n");
  fprintf(gp, "plot 'temp.dat' w l, 'temp2.dat', 'temp3.dat'\n");
  pclose(gp);
}