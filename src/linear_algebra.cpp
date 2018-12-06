#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include "linear_algebra.hpp"

namespace ipds {

  const double eps = 1e-13;

  std::ostream & operator<<(std::ostream & os, const Matrix & mat) {
    for(std::size_t i = 0; i < mat.rows(); i++) {
      for(std::size_t j = 0; j < mat.cols(); j++) {
        os << mat[i][j] << " ";
      }
      if(i != mat.rows() - 1) {
        os << "\n";
      }
    }
    return os;
  }

  Matrix::Matrix(const std::vector<std::vector<double>> & mat_) {
    mat = mat_;
    row = mat_.size();
    col = mat_[0].size();
    // check whether mat_ is well-formed
    for(std::size_t i = 0; i < row; i++) {
      if(row != mat_[i].size()) {
        throw "matrix error: ill-formed";
      }
    }
  }
  const std::vector<double> & Matrix::operator[] (int i) const { return mat[i]; }
  std::size_t Matrix::rows() const { return row; }
  std::size_t Matrix::cols() const { return col; }

  std::vector<double> Matrix::operator() (const std::vector<double> & v) const {
    if( v.size() != col ) throw "matrix error: size is different";
    std::vector<double> u(row, 0);
    for(std::size_t i = 0; i < row; i++) {
      for(std::size_t j = 0; j < col; j++) {
        u[i] += v[i] * mat[i][j];
      }
    }
    return u;
  }

  SquareMatrix::SquareMatrix(const std::vector<std::vector<double>> & mat_) : Matrix(mat_) {
    // check whether mat_ is square matrix
    if(row != col) {
      throw "matrix error: not square";
    }
  }

  SymmetricMatrix::SymmetricMatrix(const std::vector<std::vector<double>> & mat_) : SquareMatrix(mat_) {
    // check whether mat_ is symmetric
    for(std::size_t i = 0; i < row; i++) {
      for(std::size_t j = i + 1; j < col; j++) {
        if(std::abs(mat_[i][j] - mat_[j][i]) > eps) {
          throw "matrix error: not symmetric";
        }
      }
    }
  }

  std::vector<std::tuple<double, std::vector<double>>> SymmetricMatrix::get_eigen() const {
    int n = row;
    if( n == 1 ) { return { {mat[0][0], {1}} }; }
    std::vector<std::vector<double>> a = mat;
    std::vector<std::vector<double>> eigen_vectors(n, std::vector<double>(n, 0));
    for(int i = 0; i < n; i++){ eigen_vectors[i][i] = 1; }

    auto givens_rotate = [&]() {
      int p, q;
      double max_value = -1;
      // find argmax_{(p,q), p<q} {|a[p][q]|}
      for(int i = 0; i < n; i++) {
        for(int j = i + 1; j < n; j++) {
          if( max_value < std::abs(a[i][j]) ) {
            p = i;
            q = j;
            max_value = std::abs(a[i][j]);
          }
        }
      }
      std::vector<double> ap = a[p], aq = a[q];
      double app = a[p][p], aqq = a[q][q], apq = a[p][q];
      double theta = 0.5 * std::atan(-2 * apq / (app - aqq)); // angle of rotation
      double sin_theta = std::sin(theta), sin_2theta = std::sin(2 * theta);
      double cos_theta = std::cos(theta), cos_2theta = std::cos(2 * theta);
      // rotate
      for(int i = 0; i < n; i++) {
        a[i][p] = a[p][i] = ap[i] * cos_theta - aq[i] * sin_theta;
        a[i][q] = a[q][i] = ap[i] * sin_theta + aq[i] * cos_theta;
      }
      a[p][p] = (app + aqq) / 2 + ((app - aqq) / 2) * cos_2theta - apq * sin_2theta;
      a[q][q] = (app + aqq) / 2 - ((app - aqq) / 2) * cos_2theta + apq * sin_2theta;
      a[p][q] = a[q][p] = 0;
      // update the approximation of the eigen vectors
      std::vector<double> ev_p(n), ev_q(n);
      for(int i = 0; i < n; i++) {
        ev_p[i] = eigen_vectors[i][p];
        ev_q[i] = eigen_vectors[i][q];
      }
      for(int i = 0; i < n; i++) {
        eigen_vectors[i][p] = ev_p[i] * cos_theta - ev_q[i] * sin_theta;
        eigen_vectors[i][q] = ev_p[i] * sin_theta + ev_q[i] * cos_theta;
      }
    };

    while(true){
      givens_rotate();
      double max_non_diag = 0;
      for(int i = 0; i < n; i++) {
        for(int j = i + 1; j < n; j++) {
          max_non_diag = std::max(max_non_diag, a[i][j]);
        }
      }
      if( max_non_diag < eps ) break;
    }
    std::vector<std::tuple<double, std::vector<double>>> eigen;
    for(int i = 0; i < n; i++) {
      std::vector<double> ev(n);
      for(int k = 0; k < n; k++) {
        ev[k] = eigen_vectors[k][i];
      }
      eigen.emplace_back(a[i][i], ev);
    }

    auto comp = [](const std::tuple<double, std::vector<double>> & a, const std::tuple<double, std::vector<double>> & b) -> bool {
      return std::abs(std::get<0>(a)) < std::abs(std::get<0>(b));
    };
    std::sort(eigen.begin(), eigen.end(), comp);
    std::reverse(eigen.begin(), eigen.end());
    return eigen;
  }

  std::tuple<SymmetricMatrix, std::vector<std::vector<double>>> SymmetricMatrix::tridiagonalize() const {
    std::size_t n = mat.size();
    auto T = mat;
    std::vector<std::vector<double>> vs;
    for(std::size_t i = 0; i < n - 2; i++) {
      auto [s, v] = householder(T[i], i + 1);
      vs.push_back(v);
      // transform T
      T[i][i + 1] = T[i + 1][i] = s;
      for(std::size_t j = i + 2; j < T.size(); j++) {
        T[i][j] = T[j][i] = 0;
      }
      std::vector<double> d(T.size() - i - 1, 0);
      std::vector<double> g(T.size() - i - 1, 0);
      for(std::size_t k = 0; k < d.size(); k++) {
        for(std::size_t l = 0; l < d.size(); l++) {
          d[k] += T[i + 1 + k][i + 1 + l] * v[i + 1 + l];
        }
      }
      double tmp = 0;
      for(std::size_t k = 0; k < d.size(); k++) {
        tmp += v[i + 1 + k] * d[k];
      }
      for(std::size_t k = 0; k < g.size(); k++) {
        g[k] = 2 * (d[k] - tmp * v[i + 1 + k]);
      }
      for(std::size_t k = 0; k < g.size(); k++) {
        for(std::size_t l = 0; l < g.size();l++) {
          T[i + 1 + k][i + 1 + l] -= g[k] * v[i + 1 + l] + v[i + 1 + k] * g[l];
        }
      }
    }
    return {SymmetricMatrix(T), vs};
  }

  std::tuple<double, std::vector<double>> householder(const std::vector<double> & u, const std::size_t idx) {
    std::vector<double> v = u;
    for(std::size_t i = 0; i < idx; i++) v[i] = 0;
    double s = ((v[idx] >= 0) ? 1 : -1) * norm(v);
    v[idx] -= s;
    double n = norm(v);
    for(auto & a : v) { a /= n; }
    return {s, v};
  }

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
      /*
      for(std::size_t k = 0; k < u[i].size(); k++) {
        u[i][k] /= norms[i];
      }
      */
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
}
