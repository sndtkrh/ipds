#ifndef IPDS_DATAMANIPULATION
#define IPDS_DATAMANIPULATION

namespace ipds {
  std::vector<double> pointwise(const std::function<double(double, double)> & op, const std::vector<double> v, const std::vector<double> & u);
  std::vector<double> pointwise_mult(const std::vector<double> v, const std::vector<double> & u);
  std::vector<std::vector<double>> polynomial_features(const std::vector<std::vector<double>> & data, unsigned int deg = 2);
}

#endif
