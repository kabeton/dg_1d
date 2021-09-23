#ifndef QUAD_H
#define QUAD_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <functional>

double P(int n, double x) {
  if(n == 0) return 1;
  if(n == 1) return x;
  return ((2*n - 1)*x*P(n - 1, x) - (n - 1)*P(n - 2, x))/n;
}

double dP(int n, double x) {
  if(n == 0) return 0;
  if(n == 1) return 1;
  return P(n - 1, x)*n + x*dP(n - 1, x);
}

Eigen::ArrayXd get_roots(int k) { 
  Eigen::ArrayXd roots(k);
  double eps = 1e-12;
  int i = 0;

  for(double x = -1.; x <= 1.; x += 0.01) {
    if(P(k, x)*P(k, x + 0.01) < 0) {
      double left = x;
      double right = x + 0.01;
      double c = (left + right)/2;
      while(fabs(left - right) > eps) {
        if(P(k, c)*P(k, right) < 0) {
          left = c;
        } else {
          right = c;
        }
        c = (left + right)/2;
      }
      c = (left + right)/2;
      roots(i) = c;
      i++;
    } 
  }

  return roots;
}

Eigen::ArrayXd get_weights(int k) {
  Eigen::ArrayXd roots(k); 
  Eigen::ArrayXd weights(k);
  roots = get_roots(k);

  for(int i = 0; i < k; i++) {
    double pp = dP(k, roots(i));
    weights(i) = 2./(1 - roots(i)*roots(i))/pp/pp;
  }

  return weights;
}

double gl_quad0d(int m, int n, int k) {
  Eigen::ArrayXd roots(k); 
  Eigen::ArrayXd weights(k);
  double res = 0;
  
  for(int i = 0; i < k; i++) {
    res += weights(i)*P(i, roots(m))*P(n, roots(i));
  }
  return res;
}

double gl_quad1d(int m, int n, int k) {
  Eigen::ArrayXd roots(k); 
  Eigen::ArrayXd weights(k);
  double res = 0;
  
  for(int i = 0; i < k; i++) {
    res += weights(i)*dP(i, roots(m))*P(n, roots(i));
  }
  return res;
}
#endif