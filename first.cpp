#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include "quad.h"

int main() {  
  //1-d grid (equidistant points for now)
  //we assume that problem is initially defined on [0, L] 
  double L = 5.;
  int N;
  std::cout << "input number of points N" << std::endl;
  std::cin >> N;
  Eigen::ArrayXd x(N);
  double step = L / N;
  for(int i = 0; i < N; i++) {
    x(i) = step*i;
  }

  int k;
  std::cout << "input max polynomial degree k" << std::endl;
  k = 5;
  //std::cin >> k; 
  k++;

  //precompute matrices
  Eigen::MatrixXd M(k, k);
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(k);

  M << gl_quad0d(0, 0, k), 0, 0,
       0, gl_quad0d(1, 1, k), 0, 
       0, 0, gl_quad0d(2, 2, k);                  

  for(int i = 0; i < k; i++) {
    for(int j = 0; j < k; j++) {
      K(i, j) = gl_quad1d(j, i, k);
    }
  }

  //cycle over each section
  for(int i = 0; i < N - 1; i++) {

  }

  //save results

  return 0;
}