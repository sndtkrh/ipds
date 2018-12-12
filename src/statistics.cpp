#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>

#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#define BOOST_NUMERIC_BINDINGS_USE_CLAPACK
#include <boost/numeric/bindings/lapack/syevd.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#undef  BOOST_NUMERIC_BINDINGS_USE_CLAPACK
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

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

  boost::numeric::ublas::symmetric_matrix<double> var_covar_matrix(const std::vector<std::vector<double>> & data) {
    int n = data.size();
    boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> A(n, n);
    for(int i = 0; i < n; i++) {
      for(int j = i; j < n; j++) {
        A(i, j) = covar(data[i], data[j]);
      }
    }
    return A;
  }

  boost::numeric::ublas::symmetric_matrix<double> corel_matrix(const std::vector<std::vector<double>> & data) {
    int n = data.size();
    boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> A(n, n);
    std::vector<double> sigma(n);
    for(int i = 0; i < n; i++) {
      sigma[i] = std::sqrt(covar(data[i], data[i]));
    }
    for(int i = 0; i < n; i++) {
      for(int j = i; j < n; j++) {
        A(i, j) = covar(data[i], data[j]) / sigma[i] / sigma[j];
      }
    }
    return A;
  }

  std::vector<std::tuple<double, std::vector<double>>> pca(const std::vector<std::vector<double>> & data) {
    auto [eigenvalues, eigenvectors] = eigen(ipds::corel_matrix(data));
    std::vector<std::tuple<double, std::vector<double>>> eigen_;
    for(std::size_t i = 0; i < eigenvectors.size1(); i++) {
      auto u = boost::numeric::ublas::column(eigenvectors, i);
      std::vector<double> v(u.size());
      for(std::size_t i = 0; i < u.size(); i++) {
        v[i] = u(i);
      }
      eigen_.emplace_back(eigenvalues[i], v);
    }
    std::sort(eigen_.begin(), eigen_.end(),
      [](const std::tuple<double, std::vector<double>> & a, const std::tuple<double, std::vector<double>> & b) {
        return std::abs(std::get<0>(a)) > std::abs(std::get<0>(b));
      } );
    return eigen_;
  }

  void print_pca_result(const std::vector<std::tuple<double, std::vector<double>>> & pca_result, int limit) {
    double sum_of_eigenvalues = 0;
    for(const auto &  eigen : pca_result) {
      sum_of_eigenvalues += std::get<0>(eigen);
    }
    double cumulative_proportion = 0;
    for(std::size_t i = 0; i < ((limit == -1) ? pca_result.size() : static_cast<std::size_t>(limit)); i++) {
      double eigenvalue = std::get<0>(pca_result[i]);
      const auto & eigenvector = std::get<1>(pca_result[i]);
      cumulative_proportion += eigenvalue / sum_of_eigenvalues;
      std::cout
      << "\x1b[36mStandard deviation\x1b[39m: " << std::sqrt(eigenvalue)
      << ", \x1b[36mProportion of Variance\x1b[39m: " << eigenvalue / sum_of_eigenvalues
      << ", \x1b[36mCumulative Proportion\x1b[39m: " << cumulative_proportion
      << std::endl;
      std::cout << "\x1b[36mEigen vector\x1b[39m: ";
      for(auto v : eigenvector) {
        std::cout << v << " ";
      }
      std::cout << std::endl;
      std::cout << std::endl;
    }
  }

  void print_pca_result(const std::vector<double> & pca_result, int limit) {
    double sum_of_eigenvalues = 0;
    const auto & eigenvalues = pca_result;
    for(double eigenvalue : eigenvalues) {
      sum_of_eigenvalues += eigenvalue;
    }
    double cumulative_proportion = 0;
    for(std::size_t i = 0; i < ((limit == -1) ? eigenvalues.size() : static_cast<std::size_t>(limit)); i++) {
      double eigenvalue = eigenvalues[i];
      cumulative_proportion += eigenvalue / sum_of_eigenvalues;
      std::cout
      << "\x1b[36mStandard deviation\x1b[39m: " << std::sqrt(eigenvalue)
      << ", \x1b[36mProportion of Variance\x1b[39m: " << eigenvalue / sum_of_eigenvalues
      << ", \x1b[36mCumulative Proportion\x1b[39m: " << cumulative_proportion
      << std::endl;
    }
  }
}
