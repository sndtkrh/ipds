#include <vector>
#include <functional>
#include "data_manipulation.hpp"

namespace ipds {
  std::vector<double> pointwise(const std::function<double(double, double)> & op, const std::vector<double> v, const std::vector<double> & u) {
    if( v.size() != u.size() ) throw "pointwise: size is different";
    std::vector<double> a(v.size());
    for(std::size_t i = 0; i < v.size(); i++) {
      a[i] = op(v[i], u[i]);
    }
    return a;
  }

  std::vector<double> pointwise_mult(const std::vector<double> v, const std::vector<double> & u) {
    return pointwise([](double x, double y){return x * y;}, v, u);
  }

  std::vector<std::vector<double>> polynomial_features(const std::vector<std::vector<double>> & data, unsigned int deg) {
    std::vector<std::vector<double>> newdata = data;
    std::size_t n = data.size();
    std::size_t pre_poly_begin = 0, pre_poly_end = data.size();
    std::vector<std::size_t> pre_index(data.size(), 1);
    for(unsigned int d = 2; d <= deg; d++) {
      std::size_t idx = 0;
      std::vector<std::size_t> new_pre_index;
      for(std::size_t i = 0; i < n; i++) {
        for(std::size_t j = pre_poly_begin + idx; j < pre_poly_end; j++) {
          newdata.push_back(pointwise_mult(newdata[i], newdata[j]));
        }
        new_pre_index.push_back(pre_poly_end - pre_poly_begin - idx);
        idx += pre_index[i];
      }
      pre_index = new_pre_index;
      pre_poly_begin = pre_poly_end;
      pre_poly_end = newdata.size();
    }
    return newdata;
  }
}
