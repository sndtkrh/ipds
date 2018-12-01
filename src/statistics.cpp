#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include "linear_algebra.hpp"
#include "statistics.hpp"

namespace ipds {
  double covar(const std::vector<double> & x, const std::vector<double> & y) {
    if( x.size() != y.size() ) {
      throw "size is different";
    }
    int n = x.size();
    double ex = 0, ey = 0;
    for(int i = 0; i < n; i++) {
      ex += x[i];
      ey += y[i];
    }
    ex = ex / n;
    ey = ey / n;
    double cov = 0;
    for(int i = 0; i < n; i++) {
      cov += (x[i] - ex) * (y[i] - ey);
    }
    cov = cov / n;
    return cov;
  }

  double average(const std::vector<double> & x) {
    double sum = 0;
    for(double v : x) {
      sum += v;
    }
    return sum / x.size();
  }

  double standard_deviation(const std::vector<double> & x) {
    return std::sqrt(covar(x, x));
  }

  std::vector<double> standardize(const std::vector<double> & x) {
    double avg = average(x);
    double sigma = standard_deviation(x);
    std::vector<double> y;
    for(double v : x) {
      y.push_back( (v - avg) / sigma );
    }
    return y;
  }

  SymmetricMatrix var_covar_matrix(const std::vector<std::vector<double>> & data) {
    int n = data.size();
    std::vector<std::vector<double>> varcovar(n, std::vector<double>(n));
    for(int i = 0; i < n; i++) {
      for(int j = i; j < n; j++) {
        varcovar[i][j] = varcovar[j][i] = covar(data[i], data[j]);
      }
    }
    return SymmetricMatrix(varcovar);
  }

  SymmetricMatrix corel_matrix(const std::vector<std::vector<double>> & data) {
    int n = data.size();
    std::vector<std::vector<double>> corel(n, std::vector<double>(n));
    std::vector<double> sigma(n);
    for(int i = 0; i < n; i++) {
      sigma[i] = std::sqrt(covar(data[i], data[i]));
    }
    for(int i = 0; i < n; i++) {
      for(int j = i; j < n; j++) {
        corel[i][j] = corel[j][i] = covar(data[i], data[j]) / sigma[i] / sigma[j];
      }
    }
    return SymmetricMatrix(corel);
  }
}
