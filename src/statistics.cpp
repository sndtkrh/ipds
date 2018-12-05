#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include "linear_algebra.hpp"
#include "statistics.hpp"
#include "svg/svg.hpp"

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

  std::vector<std::tuple<double, std::vector<double>>> pca(const std::vector<std::vector<double>> & data) {
    ipds::SymmetricMatrix seiseki_corel = ipds::corel_matrix(data);
    return seiseki_corel.get_eigen();
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

  void plot_point_2D (std::tuple<double, double> point, SVGcanvas & svg, double scale, std::string color, double point_radius) {
    auto [w, h] = svg.get_canvas_size();
    auto [x, y] = point;
    svg.circles.emplace_back(w / 2 + scale * x, h / 2 + scale * y, point_radius, color);
  }

  void plot_line_2D (std::tuple<double, double> point_begin, std::tuple<double, double> point_end, SVGcanvas & svg, double scale, std::string color, double stroke_width) {
    auto [w, h] = svg.get_canvas_size();
    auto [x1, y1] = point_begin;
    auto [x2, y2] = point_end;
    svg.lines.emplace_back(w / 2 + scale * x1, h / 2 + scale * y1, w / 2 + scale * x2, h / 2 + scale * y2, color, stroke_width);
  }
}
