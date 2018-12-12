#ifndef IPDS_MATRIX
#define IPDS_MATRIX

namespace ipds {
  std::vector<double> operator + (const std::vector<double> & a, const std::vector<double> & b);
  std::vector<double> operator * (const double lambda, const std::vector<double> & a);

  double inner_product(const std::vector<double> & v, const std::vector<double> & u);
  double norm(const std::vector<double> & v);
  std::vector<std::vector<double>> orthonormalize(const std::vector<std::vector<double>> & v);
  std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> & a);
  std::vector<double> project(const std::vector<double> & point, const std::vector<std::vector<double>> & basis_of_subspace);
  std::vector<std::vector<double>> project(const std::vector<std::vector<double>> & points, const std::vector<std::vector<double>> & basis_of_subspace);

  std::tuple<std::vector<double>, boost::numeric::ublas::matrix<double>> eigen(const boost::numeric::ublas::symmetric_matrix<double> & A);
}

#endif
