#ifndef IPDS_STATISTICS
#define IPDS_STATISTICS
#include "svg/svg.hpp"

namespace ipds {
  double covar(const std::vector<double> & x, const std::vector<double> & y);
  double average(const std::vector<double> & x);
  double standard_deviation(const std::vector<double> & x);
  std::vector<double> standardize(const std::vector<double> & x);
  SymmetricMatrix var_covar_matrix(const std::vector<std::vector<double>> & data);
  SymmetricMatrix corel_matrix(const std::vector<std::vector<double>> & data);
  /* pca's return value: the sorted vector of pairs of eigenvalue and eigenvector */
  std::vector<std::tuple<double, std::vector<double>>> pca(const std::vector<std::vector<double>> & data);
  void print_pca_result(const std::vector<std::tuple<double, std::vector<double>>> & pca_result);
  void plot_point_2D (std::tuple<double, double> point, SVGcanvas & svg, double scale = 1, std::string color = "black", double point_radius = 1);
  void plot_line_2D (std::tuple<double, double> point_begin, std::tuple<double, double> point_end, SVGcanvas & svg, double scale = 1, std::string color = "black", double stroke_width = 1);
}

#endif
