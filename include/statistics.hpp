#ifndef IPDS_STATISTICS
#define IPDS_STATISTICS

namespace ipds {
  double covar(const std::vector<double> & x, const std::vector<double> & y);
  double average(const std::vector<double> & x);
  double standard_deviation(const std::vector<double> & x);
  std::vector<double> standardize(const std::vector<double> & x);
  SymmetricMatrix var_covar_matrix(const std::vector<std::vector<double>> & data);
  SymmetricMatrix corel_matrix(const std::vector<std::vector<double>> & data);
  std::vector<std::tuple<double, std::vector<double>>> pca(const std::vector<std::vector<double>> & data);
}

#endif
