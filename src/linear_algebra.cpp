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

namespace ipds {

  const double eps = 1e-13;

  double inner_product(const std::vector<double> & v, const std::vector<double> & u) {
    if( v.size() != u.size() ) throw "inner_product: different dimmension!";
    double a = 0;
    for(std::size_t i = 0; i < v.size(); i++) {
      a += v[i] * u[i];
    }
    return a;
  }

  double norm(const std::vector<double> & v) {
    return std::sqrt(inner_product(v,v));
  }


  std::vector<double> operator + (const std::vector<double> & a, const std::vector<double> & b) {
    if( a.size() != b.size() ) throw "operator +: size is different";
    std::vector<double> c;
    for(std::size_t i = 0; i < a.size(); i++) {
      c.push_back(a[i] + b[i]);
    }
    return c;
  }

  std::vector<double> operator * (const double lambda, const std::vector<double> & a) {
    std::vector<double> b;
    for(std::size_t i = 0; i < a.size(); i++) {
      b.push_back(lambda * a[i]);
    }
    return b;
  }

  std::vector<std::vector<double>> orthonormalize(const std::vector<std::vector<double>> & v) {
    // Gram–Schmidt orthonormalization
    // orthogonalize
    std::vector<std::vector<double>> u(v.size());
    std::vector<double> norms(v.size());
    u[0] = v[0];
    norms[0] = norm(u[0]);
    for(std::size_t i = 1; i < v.size(); i++) {
      u[i] = v[i];
      for(std::size_t j = 0; j < i; j++) {
        double coeff = inner_product(u[j], v[i]) / norms[j];
        for(std::size_t k = 0; k < u[i].size(); k++) {
          u[i][k] -= coeff * u[j][k];
        }
      }
      norms[i] = norm(u[i]);
    }
    // normalize
    for(std::size_t i = 0; i < u.size(); i++) {
      u[i] = (1/norms[i]) * u[i];
    }
    return u;
  }

  std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> & a) {
    std::vector<std::vector<double>> aT(a[0].size(), std::vector<double>(a.size()));
    for(std::size_t i = 0; i < a.size(); i++) {
      if( a[0].size() != a[i].size() ) throw "transpose: size is different";
      for(std::size_t j = 0; j < a[i].size(); j++) {
        aT[j][i] = a[i][j];
      }
    }
    return aT;
  }

  std::vector<double> project(const std::vector<double> & point, const std::vector<std::vector<double>> & basis_of_subspace) {
    auto e = orthonormalize(basis_of_subspace);
    std::vector<double> projected_point(e.size());
    for(std::size_t k = 0; k < e.size(); k++) {
      projected_point[k] = inner_product(point, e[k]);
    }
    return projected_point;
  }

  std::vector<std::vector<double>> project(const std::vector<std::vector<double>> & points, const std::vector<std::vector<double>> & basis_of_subspace) {
    std::vector<std::vector<double>> projected_points;
    for(const auto & p : points) {
      projected_points.push_back(project(p, basis_of_subspace));
    }
    return projected_points;
  }

  std::tuple<std::vector<double>, boost::numeric::ublas::matrix<double>> eigen(const boost::numeric::ublas::symmetric_matrix<double> & A) {
    namespace ublas = boost::numeric::ublas;
    std::vector<double> v(A.size1());
    ublas::matrix<double, ublas::column_major> Q(A.size1(), A.size2());
    int info;
    for (std::size_t i = 0; i < A.size1(); ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        Q(j,i) = Q(i, j) = A(i, j);
      }
    }
    info = boost::numeric::bindings::lapack::syevd('V', 'L', Q, v, boost::numeric::bindings::lapack::optimal_workspace());
    BOOST_UBLAS_CHECK(info == 0, ublas::internal_logic());
    return {v, Q};
  }
}
