#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <cmath>
#include "matrix.hpp"

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
