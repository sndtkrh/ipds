#ifndef IPDS_STATISTICS
#define IPDS_STATISTICS

namespace ipds {
  double covar(const std::vector<double> & x, const std::vector<double> & y);
  double average(const std::vector<double> & x);
  double standard_deviation(const std::vector<double> & x);
  std::vector<double> standardize(const std::vector<double> & x);
  std::vector<std::vector<double>> standardize(const std::vector<std::vector<double>> & data);
  boost::numeric::ublas::symmetric_matrix<double> var_covar_matrix(const std::vector<std::vector<double>> & data);
  boost::numeric::ublas::symmetric_matrix<double> corel_matrix(const std::vector<std::vector<double>> & data);
  /* pca's return value: the sorted vector of pairs of eigenvalue and eigenvector */
  std::vector<std::tuple<double, std::vector<double>>> pca(const std::vector<std::vector<double>> & data);
  void print_pca_result(const std::vector<std::tuple<double, std::vector<double>>> & pca_result, int limit = -1);
  void print_pca_result(const std::vector<double> & pca_result, int limit = -1);
}

#endif
