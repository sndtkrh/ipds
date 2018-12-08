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

  std::vector<std::vector<double>> standardize(const std::vector<std::vector<double>> & data) {
    std::vector<std::vector<double>> standardized_data;
    for(const auto & x : data) {
      standardized_data.push_back(standardize(x));
    }
    return standardized_data;
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
    std::cout << "n=" << n << std::endl;
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

  std::vector<std::tuple<double, std::vector<double>>> pca(const std::vector<std::vector<double>> & data) {
    ipds::SymmetricMatrix seiseki_corel = ipds::corel_matrix(data);
    return seiseki_corel.get_eigen_householder_binarysearch();
  }

  void print_pca_result(const std::vector<std::tuple<double, std::vector<double>>> & pca_result) {
    double sum_of_lambda = 0;
    for(const auto & eigen : pca_result) {
      sum_of_lambda += std::get<0>(eigen);
    }
    double cumulative_proportion = 0;
    for(const auto & eigen : pca_result) {
      const auto & e_lambda = std::get<0>(eigen);
      const auto & e_vector = std::get<1>(eigen);
      cumulative_proportion += e_lambda / sum_of_lambda;
      std::cout
      << "\x1b[36mStandard deviation\x1b[39m: " << std::sqrt(e_lambda)
      << ", \x1b[36mProportion of Variance\x1b[39m: " << e_lambda / sum_of_lambda
      << ", \x1b[36mCumulative Proportion\x1b[39m: " << cumulative_proportion
      << std::endl;
      std::cout << "\x1b[36mEigen vector\x1b[39m: ";
      for(std::size_t i = 0; i < e_vector.size(); i++) {
        std::cout << e_vector[i] << " ";
      }
      std::cout << std::endl;
      std::cout << std::endl;
    }
  }
}
